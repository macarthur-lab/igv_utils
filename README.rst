Overview
~~~~~~~~

Python interface to IGV - inspired by a tool written by `@monkollek
<https://github.com/monkollek>`_.

It includes the following executables:

* :code:`igv`    Launches IGV from the command line and optionally makes it load some file(s) and jump to some locus. 
* :code:`igv_plotter`   Automates taking IGV screenshots of one or more data files at one or more loci.

Install
~~~~~~~~

* Install :code:`igv_plotter`:

.. code:: py
    
    git clone https://github.com/macarthur-lab/igv_plotter.git
    cd igv_plotter
    python setup.py install   # add --user to install in your home directory
    
* If you don't have IGV installed yet, `download the IGV binary distribution <https://www.broadinstitute.org/software/igv/download>`_ and unzip it in some directory (for example :code:`~/bin`). If you do have it already, find the install directory. 


Configure
~~~~~~~~~

For the :code:`igv_plotter` or :code:`igv` scripts to work, they need to know the path of the :code:`igv.jar` file (see Install section). 
This can be set using either the :code:`--igv-jar-path` command line option or by setting it in the 
:code:`~/.igv_plotter` config file using the following line:

:code:`igv-jar-path = <full path of igv.jar>`

Similarly, any of the other :code:`igv_plotter` or :code:`igv` command line options can be set in this config
file by putting:

:code:`<command line option (without --)> = <value>`

Run
~~~

**igv_plotter usage example:**

To run :code:`igv_plotter` so that it loads 3 files, and takes 2 snapshots, do:

:code:`igv_plotter  my_file1.vcf  my_file2.bam  my_file3.bed 1:12345 chrX:12345`

For all options, see:

:code:`igv_plotter -h`


**igv usage example:**

To launch IGV from the command line and make it load 2 files and go to 1:12345, do:

:code:`igv  my_file1.vcf  my_file3.bed 1:12345`

For all options, see:

:code:`igv -h`

Contribute
~~~~~~~~~~

.. image:: https://api.travis-ci.org/macarthur-lab/igv_plotter.svg?branch=master
   :target: https://travis-ci.org/macarthur-lab/igv_plotter
    

Issue reporting and pull requests are appreciated.

Unit tests can be run with:

:code:`python setup.py test`

    
This code is open sourced under the MIT license. 



