Overview
~~~~~~~~

Python interface to IGV. Based on a tool by @monkollek.


Install
~~~~~~~~

.. code:: py
    
    git pull https://github.com/macarthur-lab/igv_plotter.git
    cd igv_plotter
    python setup.py install --user

Also, if you haven't already, please `download a version of IGV
<https://github.com/broadinstitute/IGV/releases/>`_ and unzip it.

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

