
/*
Advanced Optical Microscopy Unit
Scientific and Technological Centers. Clinic Medicine Campus
UNIVERSITY OF BARCELONA
C/ Casanova 143
Barcelona 08036 
Tel: 34 934037159

------------------------------------------------
Gemma Martin (gemmamartin@ub.edu) , Maria Calvo (mariacalvo@ub.edu)
------------------------------------------------

Name of Macro: Colocalization_LipidDroplets.ijm (144_Colocalization_Lipid_droplets_RedandGreen_v5.ijm)


Date: 1 April 2020														

Objective: Quantify localization of two fluorescent labelled proteins in fluorescent labelled lipid droplets from cultured cells in confocal microscopy stack z section images.

Input: Z section confocal images of cells labelled, 4 fluorescent channels (X protein labelling, Y protein labelling, lipid droplet labelling and DNA DAPI labelling)

Output: Sum of intensities of protein X in lipid droplets mask, sum of intensities of protein X in whole cell
        Sum of intensities of protein Y in lipid droplets mask, sum of intensities of protein Y in whole cell
        (Sum of intensities of protein X in lipid droplets mask/sum of intensities of protein X in whole cell)*100
        (Sum of intensities of protein Y in lipid droplets mask, sum of intensities of protein Y in whole cell)*100

*/


if(isOpen("Results")){
    IJ.deleteRows(0, nResults);
}
run("ROI Manager...");
roiManager("reset"); //to delete previous ROIs
IJ.deleteRows(0, nResults);


dir = getDirectory("Choose images folder");
list=getFileList(dir);
dirRes=dir+"Results"+File.separator;
File.makeDirectory(dirRes);
run("ROI Manager...");

run("Set Measurements...", "area mean min integrated limit display redirect=None decimal=5");
run("Options...", "iterations=1 count=1 do=Nothing");

path = dirRes+"Results.txt";

File.append( "Image \t Serie \t Cell \t Area Cell \t Area LD (um2) \t # LD total \t # LD colocalized Red \t % LD colocalized Red \t # LD colocalized Green \t % LD colocalized Green \t Int. Red total \t Int. Red colocalized \t % Coloc Red \t Int. Green total \t Int. Green colocalized \t % Coloc Green ", path);

for(i=0;i<list.length;i++){
	
	if(endsWith(list[i],".lif")){

		run("Bio-Formats Importer", "open=[" + dir + list[i] + "] autoscale color_mode=Default open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

		numseries=nImages;
		run("Close All");

		
		for (k=1; k<=numseries; k++) {

			run("Bio-Formats Importer", "open=[" + dir + list[i] + "] autoscale color_mode=Default rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT series_"+k);
			
			
			a=getTitle();
			
			t=substring(a,0,lengthOf(a));

			
		    dirResimages=dirRes+"Images_"+t+File.separator;
		        
			File.makeDirectory(dirResimages);

			nimages=nImages;

			for (j=1; j<=nimages; j++) {
	        
		        selectImage(j);
		        title=getTitle();
		        
		        t=substring(title,0,lengthOf(title)-6);

		        
		        if (matches(title,".*C=0.*")==1){  
					red=getImageID();
					run("Red");
					rename("Red");

					run("Duplicate...", "duplicate");

					green=getImageID();
					run("Green");
					rename("Green");
					

					
					//run("Subtract Background...", "rolling=50 stack");
				}
		        if (matches(title,".*C=1.*")==1){   
		        	blue=getImageID();
					run("Blue");
					rename("Blue");
					//run("Subtract Background...", "rolling=50 stack");

				}
				
			}

			selectImage(red);

			slices=nSlices;
			
			Intensity =newArray(slices+1);
			Intensity[0]=0;

			for (j=0;j<slices;j++) {
				setSlice(j+1);
				run("Measure");
				Intensity[j+1]=getResult("RawIntDen",j);
				if (Intensity[j+1]>Intensity[j]){
					slicemax=j+1;
					
				}
				
			}

			IJ.deleteRows(0, nResults);

			setSlice(slicemax);

			selectImage(blue);
			setSlice(slicemax);

			selectImage(green);
			setSlice(slicemax);

			selectImage(red);
			
			setTool("freehand");
			run("Brightness/Contrast...");
			
			waitForUser("Select ROI from each cell and press T. Once finished press OK");

			slicemax=getSliceNumber();

			numcells=roiManager("count");


			if(numcells>0){

				roiManager("Save", dirResimages+"ROIs_cells.zip");

				Areacell =newArray(numcells);
				
				for (j = 0; j < numcells; j++) {
					
					selectImage(red);
					setSlice(slicemax);
					roiManager("select", j);
					
					run("Measure");

					Areacell[j]=getResult("Area",0);
					IJ.deleteRows(0, nResults);

					
					run("Duplicate...", "title=red_"+j);
					setBackgroundColor(0, 0, 0);
					run("Clear Outside");
					run("Select None");
					
					saveAs("Tiff", dirResimages + "Red_"+j+".tif");

					selectImage(blue);
					setSlice(slicemax);
					roiManager("select", j);
					run("Duplicate...", "title=blue_"+j);
					setBackgroundColor(0, 0, 0);
					run("Clear Outside");
					run("Select None");
					saveAs("Tiff", dirResimages+ "Blue_"+j+".tif");

					
					selectImage(green);
					setSlice(slicemax);
					roiManager("select", j);
					run("Duplicate...", "title=green_"+j);
					setBackgroundColor(0, 0, 0);
					run("Clear Outside");
					run("Select None");
					saveAs("Tiff", dirResimages+ "Green_"+j+".tif");


				}

				roiManager("reset");

				for (m = 0; m < numcells; m++) {

					selectWindow("Blue_"+m+".tif");

					run("Duplicate...", "title=blue2_"+m);
	
					run("Threshold...");
					setAutoThreshold("Li dark");
		
					setOption("BlackBackground", false);
					run("Convert to Mask");

					run("Median...", "radius=2");
					run("Watershed");

					run("Duplicate...", "title=cells_"+m);
					selectWindow("blue2_"+m);


					run("Dilate");
					run("Dilate");

					run("Create Selection");
					roiManager("add");

					run("Measure");

					AreaLD=getResult("Area",0);


					//********************* Area Red calculation *************************************

					
					selectWindow("Red_"+m+".tif");

					run("Select None");

					run("Threshold...");
					setAutoThreshold("Moments dark");

					roiManager("Select", 0);

					run("Measure");

					Areavermell=getResult("Area",1);
					Intvermell=getResult("RawIntDen",1);

					run("Select None");

					setAutoThreshold("Moments dark");

					run("Measure");

					Intvermelltotal=getResult("RawIntDen",2);
					

					//********************* Area Green calculation ***********************************

					selectWindow("Green_"+m+".tif");

					run("Select None");

					run("Threshold...");
					setAutoThreshold("Otsu dark");

					roiManager("Select", 0);

					run("Measure");

					Areagreen=getResult("Area",3);
					Intgreen=getResult("RawIntDen",3);

					run("Select None");

					setAutoThreshold("Otsu dark");

					run("Measure");

					Intgreentotal=getResult("RawIntDen",4);

					IJ.deleteRows(0, nResults);

					//******************** LP segmentation *******************************************
					

					selectWindow("cells_"+m);

					roiManager("reset");

					run("Analyze Particles...", "size=10-Infinity pixel add");

					roiManager("Save", dirResimages+"ROIs_LP_"+m+".zip");

					count=roiManager("count");
					
					countred=0;

					for (l = 0; l < count; l++) {
						
						selectWindow("Red_"+m+".tif");
						run("Select None");
						setAutoThreshold("Moments dark");
						
						roiManager("Select", l);
						
						
						run("Measure");
						
						if(getResult("Area",l)>0){
							countred=countred+1;
						}
						
						
					}

					countgreen=0;

					for (l = 0; l < count; l++) {
						
						selectWindow("Green_"+m+".tif");
						run("Select None");
						setAutoThreshold("Otsu dark");
						
						roiManager("Select", l);
						
						
						run("Measure");
						
						if(getResult("Area",count+l)>0){
							countgreen=countgreen+1;
						}
						
						
					}



					File.append( t + "\t" + k + "\t" + m+1 + "\t" + Areacell[m] + "\t" +  AreaLD + "\t" +  count + "\t" + countred + "\t" + countred/count*100 + "\t" + countgreen + "\t" + countgreen/count*100 +  "\t" + Intvermelltotal + "\t" + Intvermell + "\t" + Intvermell/Intvermelltotal*100 + "\t" + Intgreentotal + "\t" + Intgreen + "\t" + Intgreen/Intgreentotal*100 , path);

					IJ.deleteRows(0, nResults);
					roiManager("reset");
					
					
				}

				
				
			}

			closeImagesWindows();
			
	
		}


	}


}

waitForUser("Macro has finished");

function closeImagesWindows(){
	run("Close All");
	if(isOpen("Results")){
		selectWindow("Results");
		run("Close");
	}
	if(isOpen("ROI Manager")){
		selectWindow("ROI Manager");
		run("Close");
	}
	if(isOpen("Threshold")){
		selectWindow("Threshold");
		run("Close");
	}
	if(isOpen("Summary")){
		selectWindow("Summary");
		run("Close");
	}
	if(isOpen("B&C")){
		selectWindow("B&C");
		run("Close");
	}

}


	
