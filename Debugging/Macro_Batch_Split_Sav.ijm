macro split_stack{ 
dir1 = getDirectory("Choose source directory"); 
list = getFileList(dir1); 
dir2 = getDirectory("Choose destination directory"); 
a = 1; 
Dialog.create("Input"); 
Dialog.addNumber("Number of images per substack",a); 
Dialog.show(); 
a = Dialog.getNumber(); 
setBatchMode(true); 
for (i=0; i<list.length; i++) { 
 open(dir1+list[i]); 
 getDimensions(width, height, channels, slices, frames); 
 counter=1; 
 title = File.nameWithoutExtension; 
 for (j=1; j< slices;j+=a){ 
  k = j+a-1; 
  if (k>slices)  k = slices;
	
  run("Make Substack...", "channels=1-"+channels+" slices="+j+"-"+k); 
  saveAs("tif", dir2+title+"-stack"+counter+".tiff"); 
  close(); 
  counter+=1; 
 } 
 close(); 
} 
showMessage("Created and saved substacks for all files in directory. Process is finished."); 
} 
