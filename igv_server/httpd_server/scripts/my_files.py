#!/usr/bin/python

from __future__ import print_function

import hashlib
import os
import sys
import urlparse

WWWROOT_PATH = "/local/httpd/wwwroot/"
BASE_URL = "http://maclab-utils:8000"

print("Content-type: text/html\n\n")

def log(msg):
    """Writes a log message to the Apache httpd error log"""
    print(msg, file=sys.stderr)  # will go into logs/http-error.log

def compute_md5(path):
    """Returns a fixed-width md5 hash string for the given directory path"""
    h = hashlib.md5(path)
    return h.digest().encode('base64')[:10]   # even with 10 chars, collisions are highly improbable

def create_error_xml(directory, error_message):
    """Generates an xml message that can't be sent to IGV to report an error"""
    log("%s: %s" % (directory, error_message))

    top_level_label = os.path.basename(directory) or "my_files"
    result = ""
    #result += "<?xml version='1.0' encoding='UTF-8'?>\n"
    result += "<Global name='%(top_level_label)s'>\n" % locals()
    result += "<Resource name='%(error_message)s' />\n" % locals()
    result += "</Global>"
    return result

request_url = os.environ.get('REQUEST_URI')
url_args = os.environ.get('QUERY_STRING')
if not url_args:
    #print("""<b>ERROR:</b> 'genome' and 'directory' URL parameters not provided.<br><br>
    #URL must be of the form: %(BASE_URL)s/scripts/my_files.py?genome=hg19&directory=/some/directory/with/files""" % locals())
    print(create_error_xml('', "ERROR: invalid registry url: %(BASE_URL)s%(request_url)s. genome and directory args not provided.\nPlease check your IGV preferences and make sure that the Data Registry URL under the Advanced tab looks like http://maclab-utils:8000/scripts/dataServerRegistry.py?genome=$$&directory=/humgen/atgu1/fs03/my/directory/with/files/to/show/in/igv) " % locals()))
    sys.exit(0)

url_args_dict = urlparse.parse_qs(url_args)
genome = url_args_dict.get('genome', [''])[0]   # IGV genome: eg. 'hg19' or 'b37'
directory = url_args_dict.get('directory', [''])[0]
if not directory:
    print(create_error_xml(directory, "ERROR: directory not specified in registry url: %(BASE_URL)s%(request_url)s.\nThis can be fixed by going into your IGV preferences and editing the Data Registry URL under the Advanced tab so that it looks like http://maclab-utils:8000/scripts/dataServerRegistry.py?genome=$$&directory=/humgen/atgu1/fs03/my/directory/with/files/to/show/in/igv" % locals()))
    sys.exit(0)
if not os.path.isdir(directory) or directory.startswith("."):
    print(create_error_xml(directory, "ERROR: invalid directory: %(directory)s\nThis can be fixed either by making sure %(directory)s exists on the cluster (eg. on gold) and is world-readable (eg. has chmod permissions = 777), or alternatively by going into your IGV preferences and editing the Data Registry URL under the Advanced tab so that the directory at the end of that url is correct (eg. it must be the directory that contains the files you want to view in IGV, and must be world-readable (eg. chmod permissions = 777)). " % locals()))
    #print("""<b>ERROR:</b> Invalid directory: %s""" % directory)
    sys.exit(0)


# create a symlink to the given directory if it doesn't exist already
WWWROOT_MD5_SIMLINK = compute_md5(directory)
simlink_path = os.path.join(WWWROOT_PATH, WWWROOT_MD5_SIMLINK)
if not os.path.islink(simlink_path):
    os.symlink(directory, simlink_path) 


# file formats viewable in IGV, mapped to validators
IGV_SUPPORTED_EXTENSIONS = set([".sam", ".bam", ".bed", ".bedgraph", ".bw", ".bb", 
    ".birdseye_canary_calls", ".broadPeak", ".seg", ".cbs", ".cn", ".expr", 
    ".igv", ".fasta", ".gct", ".gff", ".vcf"])

EXTENSIONS_TO_IGNORE = set([".idx", ".tbi", ".bai"])


def create_resources_xml(path, subpath="", level=0):
    """Gives IGV access to all files underneath this path by recursively 
    generating the contents of an IGV resources.xml file.

    Args:
        path: A directory to crawl for files that can be openend in IGV
        subpath: A subdirectory relative to path
    Returns:
        The contents of the xml file as a string.
    """
    result = ""
    if level == 0:
        top_level_label = os.path.basename(path) or "my_files"
        result = ""
        #result += "<?xml version='1.0' encoding='UTF-8'?>\n"
        result += "<Global name='%(top_level_label)s'>\n" % locals()

    has_igv_file_or_dir = has_any_file_or_dir = False
    for file_or_dir in os.listdir(os.path.join(path, subpath)):
        new_subpath = os.path.join(path, subpath, file_or_dir)
        if os.path.isdir(new_subpath):
            has_igv_file_or_dir = has_any_file_or_dir = True
            # process the directory
            category_label = file_or_dir
            result += "\t"*level + "<Category name='%(category_label)s'>\n" % locals()
            result += create_resources_xml(path, os.path.join(subpath, file_or_dir), level + 1)
            result += "\t"*level + "</Category>\n"
        else:
            # process the file
            has_any_file_or_dir = True
            if not any(e for e in IGV_SUPPORTED_EXTENSIONS if e in file_or_dir) or any(e for e in EXTENSIONS_TO_IGNORE if e in file_or_dir):
                continue
            has_igv_file_or_dir = True
            has_file_or_dir = True
            resource_label = file_or_dir
            resource_url = os.path.join(BASE_URL, WWWROOT_MD5_SIMLINK, subpath, file_or_dir)
            result += "\t"*level + "<Resource name='%(resource_label)s' path='%(resource_url)s' />\n" % locals()

    if not has_any_file_or_dir or not has_igv_file_or_dir:
        result += "\t"*level + "<Resource name='ERROR: " + ("no IGV-compatible" if has_any_file_or_dir else "no") + " files found in this directory' />\n" % locals()

    if level == 0:
        result += "</Global>"

    return result



# generate the xml
xml = create_resources_xml(directory)
print(xml)

#with open("out.txt", "w") as f:
#    f.write(xml)


