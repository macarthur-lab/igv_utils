Overview
~~~~~~~~

Python interface to IGV - inspired by an IGV_plotter tool originally written by `@tmfromer <https://github.com/tmfromer>`_.

It includes the following executables:

* :code:`igv`    Launches IGV from the command line and optionally makes it load some file(s) and jump to some locus. 
* :code:`igv_plotter`   Automates taking IGV screenshots of one or more data files at one or more loci.
* :code:`igvweb_viewer`  Allows bam, vcf, and/or bed file tracks to be viewed in a web browser using `igv.js <https://github.com/jrobinso>`_.


Install
~~~~~~~~

* To install for all users, run:   
  :code:`pip install igv_plotter`   

* To install in your home directory (eg. if you don't have permissions to install for all users), add ~/.local/bin to your PATH and run:
  :code:`pip install igv_plotter --user` 


Run
~~~

To see all command line options, you can do:

:code:`igvweb_viewer -h`

:code:`igv_plotter -h`

:code:`igv -h`

**igvweb_viewer script** - usage example:

This starts a webserver for viewing 3 files at 2 loci:

:code:`igvweb_viewer my_file1.vcf  my_file2.bam  my_file3.bed 1:12345 chrX:12345`

After starting this script, open your web browser to 127.0.0.1:8000 for an interactive
browser-based IGV view of these files.

**igv_plotter script** - usage example:

This loads 3 files, and takes 2 snapshots:

:code:`igv_plotter  my_file1.vcf  my_file2.bam  my_file3.bed 1:12345 chrX:12345`

**igv script** - usage example:

This launches IGV with 2 files loaded at locus 1:12345:

:code:`igv  my_file1.vcf  my_file3.bed 1:12345`


Configure
~~~~~~~~~

An :code:`~/.igv_plotter` config file can be created and used to persitantly set any :code:`igv_plotter` or :code:`igv` command line options, by using lines like:
:code:`<command line option (without --)> = <value>`


Contribute
~~~~~~~~~~

.. image:: https://travis-ci.org/macarthur-lab/igv_utils.svg?branch=master
   :target: https://travis-ci.org/macarthur-lab/igv_utils
    

Issue reporting and pull requests are appreciated.

Unit tests can be run with:

:code:`python setup.py test`

    
This code is open sourced under the MIT license. 



