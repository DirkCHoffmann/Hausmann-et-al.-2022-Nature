//opens folder
#@ File (label = "Input directory", style = "directory") input
list = getFileList(input);
list = Array.sort(list);
output = input + "_results_fiji";
//output = "/Users/davidhausmann/Desktop/results_fiji";
File.makeDirectory(output);

for (i = 0; i < list.length; i++) { 

	xx = i +1;
	print("Time Series " + xx + " of " + list.length);

	// empties Roi manager
	nROIs = roiManager("count");
	if (nROIs > 0)
	{
		roiManager("Deselect");
		roiManager("Delete");
	}

	//opens file and creates z-projection
	open(list[i]);
	filename = list[i];
	//run("Slice Remover", "first=2 last=10000 increment=2");
	run("8-bit");
	run("Z Project...", "projection=[Max Intensity]");

	//Converts to grayscale
	run("Conversions...", "scale");

	run("Duplicate...", " ");
	run("Duplicate...", " ");
	index = lastIndexOf(filename, "."); 
  	if (index!=-1) filename2 = substring(filename, 0, index);
	saveAs("png", output + "/" + filename2 + "_max projection.png");	//save zstack
	close();
	
	//Threshhold (local threshold is wahrscheinlich der Beste von den Dreien,
	//dauert aber auch deutlich länger):
	//setAutoThreshold("Default dark");
	//setThreshold(50, 255);
	//run("Convert to Mask", "method=Default background=Dark calculate black");
	run("Auto Local Threshold", "method=Otsu radius=35 parameter_1=0 parameter_2=0 white");
	setOption("BlackBackground", true);
	
	
	//Watershed
	run("Watershed", "stack");
	
	//Sets ROIs only in first slice
	run("Analyze Particles...", "size=60-1200 add");
	

	// sets ROIs on stack
	close();
	roiManager("Show All");

	//scales all ROIS with factor 0.6
	for (j = 0; j < roiManager("count"); j++) {
		roiManager("select", j);
		run("Scale... ", "x=0.6 y=0.6 centered");
		roiManager("update")
		}

	//waitForUser("User input required", "Were ROIs set correctly? If not, move, add or delete ROIs.");

	//Measures and saves location of each ROIs
	count=roiManager("count"); 
	array=newArray(count); 
	for(k=0; k<count;k++) { 
        array[k] = k; 
	} 
	roiManager("Select", array);
	run("Select All");
	run("Set Measurements...", "center redirect=None decimal=0");
	roiManager("Multi Measure");
  	selectWindow("Results");
	saveAs("text", output + "/" + filename2 + "_center.txt");
	close("Results");

   
	//Measures mean intensities in each ROIs over time
	close();
	roiManager("Show None");
	roiManager("Show All");
	run("Set Measurements...", "mean redirect=None decimal=0");
	roiManager("Multi Measure");
 
	//saves Results as .txt file in output folder and closes results
  	selectWindow("Results");
	saveAs("text", output + "/" + filename2 + "_mean intensities.txt");
	close("Results");

	close("*");


}