from igv_api import _IGVRobot, IGVCommandLineRobot
import logging
import sys
import unittest


class IGVRobotStub(_IGVRobot):
    def __init__(self):
        super(IGVRobotStub, self).__init__()

        self.execute_impl_called = False

    def _execute_impl(self, commands):
        self.execute_impl_called = True


class TestBasicUseCases(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # enable logging to simplify debugging
        logger = logging.getLogger()
        logger.level = logging.DEBUG
        stream_handler = logging.StreamHandler(sys.stdout)
        logger.addHandler(stream_handler)

    def setUp(self):
        if sys.version_info.major < 3:
            self.assertRaisesRegex = self.assertRaisesRegexp

    def test_igv_robot_parent_class(self):
        expected_command_queue = [
            "new",
            "load file1.bam",
            "load file1.bam,file2.bam",
            "collapse",
            "squish",
            "expand",
            "viewaspairs",
            "genome hg18",
            "maxPanelHeight 30",
            "goto Y:1000",
            "goto 3:10-1000",
            "goto M:1000",
            "goto chr14:1000,chr3:12345",
            "goto X:1000,Y:12345",
            "region X 10 1000",
            "region chrMT 20 2000",
            "snapshot",
            "snapshot file1.png",
            "snapshot file2.png",
            "preference DEFAULT_FONT_SIZE 20",
            "exit",
        ]

        robot = IGVRobotStub()
        robot.new_session()
        robot.load("file1.bam")
        robot.load(["file1.bam", "file2.bam"])
        robot.collapse_all_tracks()
        robot.squish_all_tracks()
        robot.expand_all_tracks()
        robot.view_all_tracks_as_pairs()
        robot.genome("hg18")
        robot.max_panel_height(30)
        robot.goto("Y:1000")
        robot.goto("3:10-1000")
        robot.locus("M:1000")
        robot.loci_split_screen("chr14:1000,chr3:12345")
        robot.loci_split_screen(["X:1000", "Y:12345"])
        robot.region("X", "10", "1000")
        robot.region("chrMT:20-2000")
        robot.snapshot()
        robot.snapshot("file1.png")
        robot.screenshot("file2.png")
        robot.preference("DEFAULT_FONT_SIZE", "20")
        robot.exit_igv()

        self.assertListEqual(robot._get_command_queue(), expected_command_queue)

        self.assertFalse(robot.execute_impl_called)
        robot.execute()
        self.assertTrue(robot.execute_impl_called)

        self.assertEqual(len(robot._get_command_queue()), 0)  # check that execute() resets the command queue

        self.assertRaisesRegex(ValueError, "Unexpected value", robot.load, None)
        self.assertRaisesRegex(ValueError, "Invalid locus", robot.locus, "12345")
        self.assertRaisesRegex(ValueError, "Invalid locus", robot.loci_split_screen, "12345")
        self.assertRaisesRegex(ValueError, "Unexpected arg type", robot.loci_split_screen, None)
        self.assertRaisesRegex(ValueError, "start, stop args should be valid genomic positions", robot.region, "12345")

    def test_execute_method(self):
        robot = IGVRobotStub()

        self.assertFalse(robot.execute_impl_called)
        robot.execute()
        self.assertTrue(robot.execute_impl_called)

        robot.new_session()
        self.assertListEqual(robot._get_command_queue(), ["new"])
        robot.exit_igv()
        self.assertListEqual(robot._get_command_queue(), ["new", "exit"])
        robot.execute()
        self.assertListEqual(robot._get_command_queue(), [])
        robot.new_session()
        self.assertListEqual(robot._get_command_queue(), ["new"])

    def test_load_method(self):
        robot = IGVRobotStub()
        self.assertRaisesRegex(ValueError, "collapse.*squish.*expand", robot.load, "foo.bam", expand=True, squish=True)

        robot.load("f.bam", collapse=True)
        self.assertEqual(robot._get_command_queue(), ["load f.bam", "collapse f.bam"])
        robot.execute()

        robot.load("f.bam", squish=True)
        self.assertEqual(robot._get_command_queue(), ["load f.bam", "squish f.bam"])
        robot.execute()

        robot.load("f.bam", expand=True)
        self.assertEqual(robot._get_command_queue(), ["load f.bam", "expand f.bam"])
        robot.execute()

        robot.load("f.bam", expand=True, viewaspairs=True)
        self.assertEqual(robot._get_command_queue(), ["load f.bam", "expand f.bam", "viewaspairs f.bam"])
        robot.execute()

    def test_goto_locus_or_region(self):
        robot = IGVRobotStub()
        self.assertRaisesRegex(ValueError, "Invalid locus", robot.goto, "X")
        self.assertRaisesRegex(ValueError, "Invalid locus", robot.goto, "chrX: 123")
        self.assertRaisesRegex(ValueError, "Invalid locus", robot.goto, "chrX: 123")
        self.assertRaisesRegex(ValueError, "Invalid locus", robot.goto, "chrX: 123")

    def test_igv_command_line_robot_constructor(self):
        self.assertRaisesRegex(ValueError, "Invalid max_memory arg", IGVCommandLineRobot, max_memory="abcde")
        self.assertRaisesRegex(ValueError, "Invalid max_memory arg", IGVCommandLineRobot, max_memory="1kb")
        self.assertRaisesRegex(ValueError, "Invalid max_memory arg", IGVCommandLineRobot, max_memory="0.1k")
        for value in ["10k", "10m", "10G"]:
            r = IGVCommandLineRobot(max_memory=value)
            self.assertEquals(r.max_memory, value)

if __name__ == '__main__':
    unittest.main()
