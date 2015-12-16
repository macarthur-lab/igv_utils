import logging
import os
import re
import socket
import subprocess
import tempfile


class _IGVRobot(object):
    """Parent class that defines the IGV command interface (http://www.broadinstitute.org/software/igv/PortCommands)
     and expects subclasses to implement how these commands are actually sent to IGV.
     """

    def __init__(self, verbose=False):
        """
        Args:
          verbose: whether to print each command to stdout
        """
        self._verbose = verbose
        if verbose:
            logging.basicConfig(level=logging.INFO)

        self.__command_queue = []

    def new_session(self):
        """Creates a new session. Unloads all tracks except the default genome annotations."""
        self.command("new")

    def load(self, filenames, collapse=False, squish=False, expand=False, viewaspairs=False):
        """Loads data or session files.

        Args:
          filenames: list or comma-separated string of one or more paths or URLs to load.
        """
        if collapse + squish + expand > 1:
            raise ValueError("Contradictory options set to True: collapse (%(collapse)s), squish (%(squish)s), "
                             "expand (%(expand)s)" % locals())

        if type(filenames) is str:
            filename_list = filenames.split(",")
        elif type(filenames) is list:
            filename_list = filenames
        else:
            raise ValueError("Unexpected value: %(filenames)s" % locals())

        filenames = ",".join([f.strip() for f in filename_list])
        self.command("load %(filenames)s" % locals())

        for filename in filename_list:
            track_name = self._get_track_name(filename)
            if collapse:
                self.command("collapse %(track_name)s" % locals())
            elif squish:
                self.command("squish %(track_name)s" % locals())
            elif expand:
                self.command("expand %(track_name)s" % locals())

            if viewaspairs:
                self.command("viewaspairs %(track_name)s" % locals())

    def collapse_all_tracks(self):
        """Switch all tracks to collapsed."""
        self.command("collapse")

    def squish_all_tracks(self):
        """Switch all tracks to squished."""
        self.command("squish")

    def expand_all_tracks(self):
        """Switch all tracks to expanded."""
        self.command("expand")

    def view_all_tracks_as_pairs(self):
        """Switches all BAM tracks to view-as-pairs mode"""
        self.command("viewaspairs")

    def genome(self, genome_id="hg19"):
        """Switch to the given genome.

        Args:
            genome_id: Either the local path of a .genome file, or a genome id chosen from the right-most column of
            http://igv.broadinstitute.org/genomes/genomes.txt
        """
        self.command("genome %(genome_id)s" % locals())

    def max_panel_height(self, height):
        """Sets the number of vertical pixels (height) of each panel. Images are not limited to the
        data visible on the screen. Stated another way, images can include the entire panel not just the portion visible
        in the scrollable screen area. The default value for this setting is 1000, increase it to see more data,
        decrease it to create smaller images.

        Args:
            height: in pixels
        """
        self.command("maxPanelHeight %(height)s" % locals())

    def goto(self, locus):
        """Jump to a locus (eg. "chr3:12345") or region (eg. "chr3:12345-12543").
        Any syntax that is valid in the IGV search box can be used.

        Args:
            locus: genomic position string
        """
        if self._match_locus_string(locus) or self._match_region_string(locus):
            self.locus(locus)
        else:
            raise ValueError("Invalid locus or region: %(locus)s" % locals())

    def locus(self, locus):
        """Scrolls to the given locus or region. Use any syntax that is valid in the IGV search box.

        Args:
            locus: genomic position string (eg. "1:12345" or "chr3:12345-54321")
        """
        if not self._match_locus_string(locus) and not self._match_region_string(locus):
            raise ValueError("Invalid locus: %(locus)s" % locals())

        self.command("goto %(locus)s" % locals())

    def region(self, region_string_or_chrom, start=None, end=None):
        """Defines a region of interest bounded by the two loci (e.g., region chr1 100 200).

        Args:
          region_string_or_chrom: Entire region string (eg. "chr1:100-1000") or just the chrom (can be "chr1" or "1")
          start: if the 1st argument is a chromosome, this must be the start position. It can be a string or integer.
          end: if the 1st argument is a chromosome, this must be the end position. It can be a string or integer.
        """
        region_match = self._match_region_string(region_string_or_chrom)
        if region_match:
            chrom, start, end = region_match.group(1), region_match.group(2), region_match.group(3)
        else:
            chrom = region_string_or_chrom
            try:
                start = int(start)
                end = int(end)
            except Exception as e:
                raise ValueError(
                    ("Since %(region_string_or_chrom)s is not a valid region string, it should be a chromosome, and "
                    "start, stop args should be valid genomic positions. %(e)s") % locals())

        self.command("region %(chrom)s %(start)s %(end)s" % locals())

    def loci_split_screen(self, loci):
        """Displays loci in split screen. Loci can be specified using any syntax that is valid in the IGV search box.

        Args:
            loci: list of strings or a single string of comma-separated genomic position(s) (eg. "1:12345, chr3:12345")
        """
        if type(loci) is str:
            loci_list = loci.split(",")
        elif type(loci) is list:
            loci_list = loci
        else:
            raise ValueError("Unexpected arg type: %(loci)s" % locals())

        loci_list = [l.strip() for l in loci_list]
        for locus in loci_list:
            locus_match = self._match_locus_string(locus)
            if not locus_match:
                raise ValueError("Invalid locus: %(locus)s" % locals())

        loci = ",".join(loci_list)
        self.command("goto %(loci)s" % locals())


    def snapshot(self, filename=None):
        """Saves a snapshot of the IGV window to an image file.

        Args:
            filename: The image filename where the extension determines the image file format (must
            be .png, .jpg, or .svg). If filename is omitted, this writes a PNG file with a filename generated based
            on the locus.
        """
        if filename:
            self.command("snapshot %(filename)s" % locals())
        else:
            self.command("snapshot")

    def screenshot(self, filename=None):
        """Alias for snapshot"""
        self.snapshot(filename)

    def preference(self, key, value):
        """
        Sets an IGV preference. Many of these are from the IGV Settings panel.

        Some interesting preference keys and their default IGV values are:

            "DEFAULT_FONT_SIZE", "10"
            "DEFAULT_FONT_FAMILY", "Arial"
            "NAME_PANEL_WIDTH", "160"
            "BACKGROUND_COLOR", "250,250,250"
            "SASHIMI_SHOW_COVERAGE", "true"
            "NORMALIZE_COVERAGE", "false"
            "TRACK_HEIGHT_KEY", "15"
            "COLOR_MUTATIONS", "false"
            "SHOW_SINGLE_TRACK_PANE_KEY", "false"
            "SHOW_ATTRIBUTE_VIEWS_KEY", "true"

            "SAM_SHOW_COV_TRACK", "true"
            "SAM_SHOW_REF_SEQ", "false"
            "SAM_SHOW_GROUP_SEPARATOR", "true"
            "SAM_SHOW_CENTER_LINE", "true"
            "SAM_SHOW_DUPLICATES", "false"
            "SAM_SHOW_SOFT_CLIPPED", "false"
            "SAM_FLAG_UNMAPPED_PAIR", "false"
            "SAM_SHADE_CENTER", "true"
            "SAM_DOWNSAMPLE_READS", "true"
            "SAM_SAMPLING_COUNT", "100"
            "SAM_BASE_QUALITY_MIN", "5"
            "SAM_BASE_QUALITY_MAX", "20"
            "SAM_QUALITY_THRESHOLD", "0"
            "SAM_ALLELE_THRESHOLD", "0.2f"
            "SAM_ALLELE_USE_QUALITY", "true"
            "SAM_MIN_INSERT_SIZE_THRESHOLD", "50"
            "SAM_MAX_INSERT_SIZE_THRESHOLD", "1000"
            "SAM_MIN_INSERT_SIZE_PERCENTILE", "0.5"
            "SAM_MAX_INSERT_SIZE_PERCENTILE", "99.5"
            "SAM_COLOR_BY", "UNEXPECTED_PAIR"
            "SAM_COMPUTE_ISIZES", "true"
            "SAM_SHOW_JUNCTION_TRACK", "false"
            "SAM_JUNCTION_MIN_FLANKING_WIDTH", "0"
            "SAM_JUNCTION_MIN_COVERAGE", "1"
            "SAM_SHOW_JUNCTION_FLANKINGREGIONS", "true"


        For the full list of keys and default values see initDefaultValues() in org.broad.igv.PreferenceManager:
            https://github.com/broadinstitute/IGV/blob/master/src/org/broad/igv/PreferenceManager.java
        Also, there are relevant discussions on the IGV mailing list:
            https://groups.google.com/forum/#!topic/igv-help/p50O414MiTg
        """
        self.command("preference %(key)s %(value)s" % locals())

    def exit_igv(self):
        """Exit IGV"""
        self.command("exit")

    def command(self, command_string):
        """Queue an arbitrary IGV command to be executed when execute() is called.

        Args:
          command_string: the full command to send to IGV (see http://www.broadinstitute.org/software/igv/PortCommands).
        """
        self.__command_queue.append(command_string)

    def execute(self):
        """Execute queued-up IGV commands."""

        if self.__command_queue:
            logging.info("Executing commands:")
            for c in self.__command_queue:
                logging.info("  " + c)
        else:
            logging.info("0 commands to execute.")

        self._execute_impl(self.__command_queue)

        self.__command_queue = []

    def _execute_impl(self, commands):
        """Abstract method that must be overridden by subclasses.

        Executes a list of IGV commands."""
        raise NotImplementedError("abstract method")

    def _get_track_name(self, filename):
        """Returns the track name IGV would assign to a track with the given filename."""
        return os.path.basename(filename)

    def _match_locus_string(self, locus_str):
        """If this is a valid locus_str (eg. 'chrX:1000'), return a reg-exp Match object where match.group(1) is the
        chromosome and match.group(2) is the position. Otherwise, return None"""
        return re.match("([a-zA-Z0-9]{1,5}):([0-9]{1,10})", locus_str)

    def _match_region_string(self, region_str):
        """If region_str is valid (eg. 'chrX:10-100'), return a reg-exp Match object where match.group(1) is the
        chromosome and match.group(2) and match.group(3) are the start and stop positions. Otherwise, return None."""
        return re.match("([a-zA-Z0-9]{1,5}):([0-9]{1,10})-([0-9]{1,10})", region_str)

    def _get_command_queue(self):
        """Used for testing"""
        return self.__command_queue


class IGVSocketRobot(_IGVRobot):
    """Sends commands to IGV through a socket. This is useful for controlling an IGV instance that's already running,
    either locally or on a remote machine."""

    def __init__(self, host="localhost", port=60151, verbose=False):
        """Creates an IGVSocketRobot that talks to a local or a remote IGV instance via the given port.

        Args:
            host: name or ip address of computer where IGV is running. Use 'localhost' or '127.0.0.1' to connect to IGV
              running on the same computer as this script.
            port: the port IGV is listening to. Shortly after it starts up, IGV prints this to the console in a message
              like 'Listening on port X'.
            verbose: whether to print each command.
        """
        super(IGVSocketRobot, self).__init__(verbose=verbose)

        self.host = host
        self.port = port

    def _execute_impl(self, commands):
        """Executes the list of IGV commands."""
        conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        conn.connect((self.host, self.port))
        for c in commands:
            conn.sendall(c)
            conn.recv(4096)
        conn.close()


class IGVCommandLineRobot(_IGVRobot):
    """Sends commands to IGV through a batch file. This class should be used unless IGV is already running."""

    def __init__(self,
                 igv_jar_path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin/igv.jar"),
                 hide_igv_window=True,
                 igv_window_width=800,
                 igv_window_height=600,
                 max_memory="1024M",
                 verbose=False):
        """
        Args:
            igv_jar_path: Either an absolute or relative path to the igv.jar (relative to the location of this script).
            hide_igv_window: runs IGV without its window being displayed. This is necessary when running IGV
                on a server that doesn't have a display.
            igv_window_width: IGV window width in pixels
            igv_window_height: IGV window height in pixels
            max_memory: java maximum memory limit when running IGV (eg. 512M or 2G)
            verbose: whether to print out every command
        """
        super(IGVCommandLineRobot, self).__init__(verbose=verbose)

        self.igv_jar_path = igv_jar_path
        self.hide_igv_window = hide_igv_window
        self.igv_window_width = igv_window_width
        self.igv_window_height = igv_window_height

        if max_memory and not re.match("^[0-9]+[kmg]?$", max_memory, re.IGNORECASE):
            raise ValueError("Invalid max_memory arg: %s" % max_memory)
        self.max_memory = max_memory



    def set_igv_window_size(self, width=800, height=600):
        """Set IGV window size.
          igv_window_width: IGV display width in pixels
          igv_window_height: IGV window height in pixels
        """
        self.set_igv_window_width(width)
        self.set_igv_window_height(height)

    def set_igv_window_width(self, width=800):
        """Set IGV window width.
        Args:
          igv_window_width: IGV window width in pixels.
        """
        self.igv_window_width = int(width)

    def set_igv_window_height(self, height=600):
        """Set IGV window height.
        Args:
          igv_window_height: IGV window height in pixels.
        """
        self.igv_window_height = int(height)

    def _execute_impl(self, commands):
        """Executes the list of IGV commands"""

        # save commands to batch file
        if commands:
            with self._create_temp_batch_file() as f:
                batch_filename = f.name
                for c in commands:
                    f.write(("%s\n" % c).encode('UTF-8'))
        else:
            batch_filename = None

        # start IGV
        if self.hide_igv_window:
            logging.info("Hiding IGV window.")
            import xvfbwrapper
            fake_display = xvfbwrapper.Xvfb(
                width=self.igv_window_width, height=self.igv_window_height, colordepth=24)
            fake_display.start()

        self.launch_igv(batch_filename=batch_filename)

        if self.hide_igv_window:
            fake_display.stop()

    def _create_temp_batch_file(self):
        """Returns a new temp file, open for writing."""
        return tempfile.NamedTemporaryFile(delete=False)

    def launch_igv(self, batch_filename=None):
        """Launches IGV, optionally passing it the given batch file."""

        # http://www.broadinstitute.org/software/igv/startingIGV
        igv_jar = self.igv_jar_path
        igv_command = "java "
        if self.max_memory:
            igv_command += " -Xmx%s" % self.max_memory
        igv_command += " -jar %(igv_jar)s " % locals()
        if batch_filename is not None:
            igv_command += " -b " + batch_filename

        logging.info("Launching IGV: " + igv_command)
        s = subprocess.Popen(igv_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while True:
            line = s.stdout.readline().decode('UTF-8').strip('\n')
            if line:
                logging.info(line)
            elif s.poll() is not None:
                break

        if s.returncode != 0:
            raise IGVException("IGV exited with non-zero exit code: %s" % s.returncode)

        logging.info("Finished.")

class IGVException(Exception):
    pass
