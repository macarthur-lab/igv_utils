#!/usr/bin/python

from __future__ import print_function

import os
import sys
import urlparse

print("Content-type: text/html\n\n")

def log(msg):
    print(msg, file=sys.stderr)  # will go into logs/http-error.log

url_args = os.environ.get('QUERY_STRING')
url_args_dict = urlparse.parse_qs(url_args)
genome = url_args_dict.get('genome', [''])[0]   # IGV genome: eg. 'hg19' or 'b37'
directory = url_args_dict.get('directory', [''])[0]
server_name = os.environ.get('SERVER_NAME', 'maclab-utils')
server_port = os.environ.get('SERVER_PORT', '8000')
log(url_args)

print("""http://%(server_name)s:%(server_port)s/scripts/dataFiles.py?%(url_args)s
http://igv.broadinstitute.org/annotations/hg19/hg19_annotations.xml
http://www.broadinstitute.org/igvdata/1KG/b37/1KG.s3.bams.xml
#http://igv.broadinstitute.org/data/hg19/BodyMap/BodyMap_hg19.xml
#http://www.broadinstitute.org/igvdata/encode/hg19/hg19_encode.xml
#http://gdac.broadinstitute.org/tap/igv/load/tcga_hg19.xml
#http://www.broadinstitute.org/igvdata/crome/crome.xml
#http://www.broadinstitute.org/igvdata/1KG/b37/1KG_b37.php
#http://www.broadinstitute.org/igvdata/tutorials/tutorials.hg19.xml
""" % locals())
