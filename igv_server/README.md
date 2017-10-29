Instructions for how to view .bam, .vcf and other files on the Broad cluster in an IGV instance running on your laptop.

SIMPLE SETUP
------------

**On the cluster:**  
1. ssh gold.broadinstitute.org
2. create a directory like: ```mkdir -p /humgen/atgu1/fs03/${USER}/igv_server_dir```
3. create symlinks to the files you want to view in IGV: 
      ```
      cd /humgen/atgu1/fs03/${USER}/igv_server_dir
      ln -s /path/to/my_file.bam
      ln -s /path/to/my_file.bam.bai
      ```
   *NOTE:* both the data and .bai or .tbi index files should be added to igv_server_dir   
   
   *NOTE:* symlinking to directories rather than individual files also works, but if be careful with adding deeply-nested directories, 
   because then **File > Load from Server...** will be slow to come up in IGV as IGV server will need to go through all the deeply-nested subdirectories to check for files to show.    


**On your laptop:**  
4. download/run IGV from https://software.broadinstitute.org/software/igv/download  
5. In IGV, go to **View > Preferences... | Advanced** , click the "Edit server properties" checkbox, and set Data Registry URL to the url below - but replacing **<<< USER >>>** with your actual Broad username:   
   
   ```
   http://maclab-utils:8000/scripts/dataServerRegistry.py?genome=$$&directory=/humgen/atgu1/fs03/<<< USER >>>/igv_server_dir&sort=alphabetical
   ```
6. In IGV, go to **File > Load from Server..** and select the files you want to view.

*NOTE:* Only steps 3 and 6 need to be repeated after the initial setup. 

ADVANCED SETUP
--------------

In addition to the above steps, you can setup the *add_to_igv_server* script to simplify adding new files to your igv_server_dir.


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
