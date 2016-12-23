import glob
import igv_api
import logging
import os
import sys


try:
    from setuptools import setup
except ImportError:
    print("WARNING: setuptools not installed. Will try using distutils instead..")
    from distutils.core import setup


def launch_http_server(directory):
    assert os.path.isdir(directory)

    try:
        try:
            from SimpleHTTPServer import SimpleHTTPRequestHandler
        except ImportError:
            from http.server import SimpleHTTPRequestHandler

        try:
            import SocketServer
        except ImportError:
            import socketserver as SocketServer

        import socket

        for port in [80] + list(range(8000, 8100)):
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.bind(('localhost', port))
                s.close()
            except socket.error as e:
                logging.debug("Can't use port %d: %s" % (port, e.strerror))
                continue

            print("HTML coverage report now available at http://%s%s" % (
                socket.gethostname(), (":%s" % port) if port != 80 else ""))

            os.chdir(directory)
            SocketServer.TCPServer(("", port),
                SimpleHTTPRequestHandler).serve_forever()
        else:
            logging.debug("All network port. ")
    except Exception as e:
        logging.error("ERROR: while starting an HTTP server to serve "
                      "the coverage report: %s" % e)


command = sys.argv[-1]
if command == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()
elif command == "coverage":
    try:
        import coverage
    except:
        sys.exit("coverage.py not installed (pip install --user coverage)")
    setup_py_path = os.path.abspath(__file__)
    os.system('coverage run --source=igv_api ' + setup_py_path +' test')
    os.system('coverage report')
    os.system('coverage html')
    print("Done computing coverage")
    launch_http_server(directory="htmlcov")
    sys.exit()

long_description = ''
if command not in ['test', 'coverage']:
    long_description = open('README.rst').read()

setup(
    name='igv_plotter',
    version="0.9.8.7",
    description='python interface to IGV that simplifies creating screenshots of BAMs, VCFs, BEDs, etc for one-off '
                'spot checking or automated / scripted image collection',
    long_description=long_description,
    author='Ben Weisburd',
    author_email='weisburd@broadinstitute.org',
    url='https://github.com/macarthur-lab/igv_utils',
    py_modules=['igv_api'],
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'xvfbwrapper',  # creates a fake X server, allowing IGV to run without crashing with java.awt.HeadlessException
        'configargparse',  # argparse replacement that adds support for config files
        'flask',
        'requests',
    ],
    license="MIT",
    keywords='IGV, snapshots, screenshots, images, BAM, VCF, bioinformatics, sequencing, NGS, variants',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],

    test_suite='tests',
    scripts=['bin/igvweb_viewer', 'bin/igv_plotter', 'bin/igv'],
    data_files=[
        #('EGG-INFO/scripts/', ['lib/igv.jar']),
        ('static/js', glob.glob('static/js/*.*')), 
        ('static/css', glob.glob('static/css/*.*')), 
        ('static/css/img', glob.glob('static/css/img/*.*')), 
        ('static/fonts', glob.glob('static/fonts/*.*')), 
        ('bin/', ['lib/igv.jar']), # + glob.glob('bin/*')), 
        #('bin/js', glob.glob('static/js/*.*')), 
        #('bin/css', glob.glob('static/css/*.*')), 
        #('bin/css/img', glob.glob('static/css/img/*.*')), 
        #('bin/fonts', glob.glob('static/fonts/*.*')), 
    ],    
)
