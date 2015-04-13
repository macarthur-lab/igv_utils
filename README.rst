Overview
~~~~~~~~

Python interface to IGV. Based on a tool by @monkollek.


Install
~~~~~~~~

`git pull https://github.com/macarthur-lab/igv_plotter.git`  
`cd igv_plotter`  
`python setup.py install --user`  

Also, you must have IGV installed on your machine. 
If you don't, please download the latest version from: 

https://github.com/broadinstitute/IGV/releases/


Configure
~~~~~~~~~

Create a config file in your home directory: `~/.igv_plotter` 
Add to it the location of your IGV jar:

```
igv-jar-path=<path to jar>
```


Run IGV plotter
~~~~~~~~~~~~~~~

*Example:* To load 3 files, and take 2 snapshots, do:

`igv_plotter  my_file1.vcf  my_file2.bam  my_file3.bed 1:12345 chrX:12345-54321`

For all options, see:

`igv_plotter -h`

