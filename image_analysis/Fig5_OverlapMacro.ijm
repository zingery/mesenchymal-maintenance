setOption("BlackBackground", false);
run("Convert to Mask", "method=Huang background=Dark");
run("Set Scale...", "distance=0 known=0 unit=pixel");
run("Set Measurements...", "area mean integrated stack display redirect=Dec9_MGP_bckadjusted.tif decimal=3");
run("Analyze Particles...", "size=500-Infinity show=Nothing display stack");

 
