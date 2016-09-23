This tool makes it easy to open remote files on the cluster (bams, vcfs, etc.) in an IGV instance that's running on your laptop.

After you set it up, all you have to do is:

1. on the cluster, go to the files you want to view and run
   `add_to_igv_server my_file1.bam  my_file2.vcf.gz ..` 

2. on your laptop, open IGV and click File > Load from Server..., then
   select `my_file1.bam` and/or `my_file2.vcf.gz` in the dialog that comes up.


HOW TO INSTALL 
--------------

There are 3 required and 3 optional install steps:  

**On the cluster (eg. after you ssh gold):**  
1. create a directory somewhere on the cluster to store the files and directories you
   want to view in IGV (or actually symlinks to them).
   It can be anywhere as long as its file permissions aren't restricted to only you or only your group (eg. they have to be chmod 755).

_Steps 2 and 3 are optional:_  
2. edit your ~/.my.bashrc and add this line:  `export IGV_SERVER_DIRECTORY=<full path of directory created in step 1>`  
3. run `python setup.py install --user` in the same directory as this README file.  
   This installs the `add_to_igv_server` script which makes it easier to add new files to your $IGV_SERVER_DIRECTORY.  

**On your laptop:**  
4. open IGV and go to View > Preferences...  the "Advanced" tab.  
5. click the "Edit server properties" checkbox, and make Data Registry URL = the url below (after modifying the
   directory path on the end to be the same as that created in step 1):  
```
http://maclab-utils:8000/scripts/dataServerRegistry.py?genome=$$&directory=/path/to/directory/created/in/step1
```
   click Ok. 

_Optional:_  
6. Add &sort=alphabetical (including the &) to the end of the URL in
   step 5. This will make files appear in alphabetical order in the IGV dialog.


_Troubleshooting:_  
If you're not seeing the files you expect when you go to File > Load from Server...:  
1. Open the URL from step 5 in Chrome, putting `view-source:` in front of the URL.  
For example:    
```
view-source:http://maclab-utils:8000/scripts/dataServerRegistry.py?genome=$$&directory=/path/to/directory/created/in/step1&sort=alphabetical
```  
You should see a list of URLs, with the 1st one being like: `http://maclab-utils:8000/scripts/dataFiles.py?...`  
2. open this 1st url, again putting `view-source:` in front of it.  
You should see XML tags corresponding to the files and directories in your cluster directory.  


HOW TO ADD AND VIEW FILES
-------------------------

Any time you want to view a new file in IGV  

**On the cluster:**  
you can either

   a. manually add (a symlink to) the file(s) or directory into your
      $IGV_SERVER_DIRECTORY . If you add .bam or .vcf files, remember
      to also add the .bai or .tbi index files.  

   or even simpler:  

   b. run `add_to_igv_server  path/to/my_file1.bam`      

   where the arg(s) can be a single file, multiple files (eg. *.vcf.gz) and/or directories (eg. my_bams/).   
   
   *NOTE:* If you pass it .bam or .vcf files, add_to_igv_server will automatically also add their .bai or .tbi index files (assuming they exist).   
   If passing director(ies), make sure they don't contain too many deeply nested subdirectories - otherwise 
   File > Load from Server... may be slow to come up because the script will be going through all the subdirectories
   and checking for new files to show.   

**On your laptop:**  
2. in IGV, go to File > Load from Server... and select some_bam_file.bam 
   (it will be under a category named like $IGV_SERVER_DIRECTORY)  

   click Ok.   


HOW TO REMOVE OLD FILES
-----------------------

**On the cluster:**  
1. cd $IGV_SERVER_DIRECTORY  
2. delete files or directories you no longer want to see.  (if
   deleting symlinks to directories, be careful to delete the symlink and not the directory itself - eg. do `rm symlink`, but never `rm -rf <symlink>/`)  



HOW IT WORKS
------------

On your laptop, setting the Data Registry URL in IGV Preferences to
```
http://maclab-utils:8000/scripts/dataServerRegistry.py?genome=$$&directory=/humgen/atgu1/fs03/your_data_dir/directory/with/files/you/want/to/view/in/igv
```   
tells IGV that, whenever you click on File > Load from Server... it should send a request to 
```
http://maclab-utils:8000/scripts/dataServerRegistry.py 
```
to get the list of files it should show in that dialog. 

That request goes to an Apache httpd server on maclab-utils and causes a python script to run. The python script takes the directory out of the last part of the url (eg. &directory=/humgen/atgu1/fs03/your_data_dir/igv_files/), 
figures out all the files in that directory, and then sends that list back to IGV.

Then, after you select some files and click ok, IGV sends requests to http://maclab-utils:8000/ , this time requesting the
subsections of selected files which are visible in the current window. These
requests are handled by the same Apache server on maclab-utils, where it can retrieve the data directly from the directories on the cluster.
