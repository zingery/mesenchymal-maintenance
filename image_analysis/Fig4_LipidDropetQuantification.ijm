//-------------------------------------------------------------------------------------
// Zinger Yang Loureiro <zinger.yang@umassmed.edu>
// 
// Last updated: December 5, 2020
//
// About : this ImageJ macro file analyzes the lipid droplets and outputs quantitation 
// in CSV
//-------------------------------------------------------------------------------------

//----- Settings
input_dir = getDirectory("Directory containing the images to be processed"); 

//----- Main
output_dir = input_dir + "/DROPLET_QUANT/";  //create a subdirectory for output
File.makeDirectory(output_dir); 
processFiles(input_dir); 


//----- Browse through a directory and call dropletQuant for each TIFF file
function processFiles(input_dir) {
  file_count = 1;  //this is a counter for # of TIFF files.
  list = getFileList(input_dir);
  for (i=0; i<list.length; i++) {
      if (endsWith(list[i], "tif") && !startsWith(list[i], "Fiducial")){
        print((file_count++) + ": " + input_dir+list[i]);
        dropletQuant(input_dir,list[i]);
      }
   }
}

//----- Lipid droplet quantitation
function dropletQuant(input_dir,path){
  output_path = output_dir + path;
  open(input_dir+path);
  run("RGB Color");
  run("8-bit");
  close("\\Others"); 
  run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
  run("Subtract Background...", "rolling=50 light");
  run("Enhance Local Contrast (CLAHE)", "blocksize=125 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
  setAutoThreshold("Default dark");
  setThreshold(132, 255); 
  run("Convert to Mask");
  run("Analyze Particles...", "  circularity=0.50-1.00 show=[Bare Outlines] display clear");
  saveAs("Results", output_path+".txt");
  close();
  run("Close");
}

