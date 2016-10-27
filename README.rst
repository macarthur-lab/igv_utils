**igv_plotter**

python package with scripts for opening files remotely in either desktop IGV or a web viewer (igv.js), creating screenshots, etc.
It includes the following scripts:

* :code:`igv`    Launches IGV from the command line and optionally makes it load some file(s) and jump to some locus. 
* :code:`igv_plotter`   Automates taking IGV screenshots of one or more data files at one or more loci.
* :code:`igvweb_viewer`  Allows bam, vcf, and/or bed file tracks to be viewed in a web browser using `igv.js <https://github.com/jrobinso>`_.

To install it, run 
pip install igv_plotter

**igv_server**

A tool that makes it easy to open remote files on the cluster (bams, vcfs, etc.) in an IGV instance that's running on your laptop.

Setting it up involves running an Apache HTTP server on the cluster with mod-wsgi and the scripts here. 

After you set it up:

* on the cluster, go to the files you want to view and run *add_to_igv_server* *my_file1.bam* *my_file2.vcf.gz* ..  
* on your laptop, open IGV and click *File > Load from Server...*, then select *my_file1.bam* and/or *my_file2.vcf.gz* in the dialog that comes up.  


