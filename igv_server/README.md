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
5. in IGV, go to **View > Preferences... | Advanced** , click the *"Edit server properties"* checkbox, and set the *"Data Registry URL"* to the url below after replacing **<<< USER >>>** with your actual Broad username:   
   
   ```
   http://xbrowse-bams:8000/scripts/dataServerRegistry.py?genome=$$&directory=/humgen/atgu1/fs03/<<< USER >>>/igv_server_dir&sort=alphabetical
   ```
6. in IGV, go to **File > Load from Server..** and select the files you want to view.

*NOTE:* Only steps 3 and 6 need to be repeated after the initial setup. 

ADVANCED SETUP
--------------

In addition to the above steps, you can install the *add_to_igv_server* script which makes it easier to add new files to your igv_server_dir on the cluster. 

1. edit your `~/.my.bashrc` and add: 
```export IGV_SERVER_DIRECTORY=/humgen/atgu1/fs03/${USER}/igv_server_dir```
2. run 
```
git clone git@github.com:macarthur-lab/igv_utils.git
cd igv_utils/igv_server/
python setup.py install --user
```
3. to add new files or directories, run:
```
add_to_igv_server my_file.bam
```

TROUBLESHOOTING
---------------

If you're not seeing the files you expect when you go to File > Load from Server...:  
1. Open the URL from step 5 in Chrome, putting `view-source:` in front of the URL.  
For example:    
```
view-source:http://xbrowse-bams:8000/scripts/dataServerRegistry.py?genome=$$&directory=/path/to/directory/created/in/step1&sort=alphabetical
```  
You should see a list of URLs, with the 1st one being similar to: `http://xbrowse-bams:8000/scripts/dataFiles.py?...`  
2. open this url, also with `view-source:` in front of it.  
You should see XML tags corresponding to the files and directories in your *igv_server_dir* directory.  


REMOVING OLD FILES
-----------------------

**On the cluster:**  
1. cd $IGV_SERVER_DIRECTORY  
2. delete files or directories you no longer want to see.  (if
   deleting symlinks to directories, be careful to delete the symlink and not the directory itself - eg. do `rm symlink`, and not `rm -rf <symlink>/`)  



HOW IT WORKS
------------

On your laptop, setting the Data Registry URL in IGV Preferences to
```
http://xbrowse-bams:8000/scripts/dataServerRegistry.py?genome=$$&directory=/humgen/atgu1/fs03/${USER}/igv_server_dir
```   
tells IGV that, whenever you click on File > Load from Server... it should send an HTTP request to 
```
http://xbrowse-bams:8000/scripts/dataServerRegistry.py 
```
to get the list of files it should show in that dialog. 

That request goes to an Apache httpd server that's running on the `xbrowse-bams` VM and causes it to run the `dataServerRegistry.py` python script. This script takes the directory in the last part of the url (eg. &directory=/humgen/atgu1/fs03/your_data_dir/igv_files/), 
walks through all the files in that directory, and sends that list back to IGV.

Then, after you select some files and click ok, IGV sends new requests to http://xbrowse-bams:8000/ , this time requesting the
regions of these files that are visible in the current window. The httpd server on xbrowse-bams handles these as static file requests which it handles by serving the files from your igv_server_dir directory.
