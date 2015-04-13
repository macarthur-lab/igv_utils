Overview
~~~~~~~~

Python interface to IGV inspired by a tool written by `@monkollek
<https://github.com/monkollek>`_.


Install
~~~~~~~~

.. code:: py
    
    git clone https://github.com/macarthur-lab/igv_plotter.git
    cd igv_plotter
    python setup.py install --user

Also, if you haven't already, you should `download IGV
<https://github.com/broadinstitute/IGV/releases/>`_ and unzip it in some directory like :code:`~/bin`.

Configure
~~~~~~~~~

Create a config file in your home directory: 

:code:`~/.igv_plotter`

and add a line specifying the full path of :code:`igv.jar` which is inside your IGV install dir.

:code:`igv-jar-path=<path of igv.jar>`


Run IGV plotter
~~~~~~~~~~~~~~~

**Example:**

To load 3 files, and take 2 snapshots, do:

:code:`igv_plotter  my_file1.vcf  my_file2.bam  my_file3.bed 1:12345 chrX:12345-54321`

For all options, see:

:code:`igv_plotter -h`

Development and testing
~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://api.travis-ci.org/macarthur-lab/igv_plotter.svg?branch=master
   :target: https://travis-ci.org/macarthur-lab/igv_plotter
    
    

To run tests:

:code:`python setup.py test`




