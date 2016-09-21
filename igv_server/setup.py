import logging
import os
import sys


try:
    from setuptools import setup
except ImportError:
    print("WARNING: setuptools not installed. Will try using distutils instead..")
    from distutils.core import setup


command = sys.argv[-1]
#if command == 'publish':
#    os.system('python setup.py sdist upload')
#    sys.exit()

long_description = open('README.md').read()

setup(
    name='add_to_igv_server',
    version="1.0",
    description='convenient way to add files to your $IGV_SERVER_DIRECTORY',
    long_description=long_description,
    author='MacArthur Lab',
    author_email='macarthurlab@gmail.com',
    url='https://github.com/macarthur-lab',
    #py_modules=['add_to_igv_server'],
    include_package_data=True,
    install_requires=[
        'configargparse',  # argparse replacement that adds support for config files
    ],
    license="MIT",
    keywords='IGV',
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
    scripts=['add_to_igv_server'],
)
