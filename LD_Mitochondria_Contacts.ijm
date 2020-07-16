
/*
Advanced Optical Microscopy Unit
Scientific and Technological Centers. Clinic Medicine Campus
UNIVERSITY OF BARCELONA
C/ Casanova 143
Barcelona 08036 
Tel: 34 934037159

------------------------------------------------
Gemma Martin (gemmamartin@ub.edu) , Maria Calvo (mariacalvo@ub.edu).
------------------------------------------------

Name of Macro: LD_Mitochondria_Contacts (147_LD_and_Mithocondria_V11.ijm)


Date: 26 March 2020

Objective: Analyse contacts between fluorescently labelled mitochondria and lipid droplets in X protein expressing cells or in control cells from cultured cells in confocal microscopy images. 

Input: Confocal microscopy images of cells labelled, 4 fluorescent channels (X protein labelling, lipid droplet labelling, mitochondria labelling and DNA DAPI labelling)

Output: Mitochondria Area and Intensity parameters from each cell are quantified 
	    Contact perimeter and contact counts (a contact is defined as a continuous contact line) between mitochondria and lipid droplets are quantified from each cell and stored in the results table.
	    Mean X protein intensity is quantified from each cell.
Requirements: 
- Colocalization Highlighter Plugin (Pierre Bourdoncle, Institut Jacques Monod, Service Imagerie, Paris) 
- Advanced Weka Segmentation PMID 28369169, doi:10.1093/bioinformatics/btx180

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

run("Set Measurements...", "area mean min shape integrated display redirect=None decimal=5");

run("Options...", "iterations=1 count=1 do=Nothing");



path = dirRes+"Results.txt";

File.append( "Image \t # Cells \t Area Cell \t Area Mito (um2) \t Mean Grey Value \t Std desv. Grey Value \t RawIntDen \t Median \t # LD \t Mean Area LD (um2) \t Std desv. Area LD \t Mean Circ. LD \t Mean Perim. LD (um) \t Perim. total LD (um) \t # Contact \t Mean Length contact (um) \t Total length contact (um) \t Raw Int Den. Transfection", path);



for(i=0;i<list.length;i++){
	
	if(endsWith(list[i],".czi")){

		run("Bio-Formats Importer", "open=[" + dir + list[i] + "] autoscale color_mode=Composite open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

		numseries=nImages;
		
		run("Close All");

		
		for (m=1; m<=numseries; m++) {

			run("Bio-Formats Importer", "open=[" + dir + list[i] + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+m);

			title=getTitle();
			t=title;

			ImageID=getImageID();
			run("Duplicate...", "duplicate");

			run("Split Channels");

	        for (j=1; j<=nImages; j++) {

	        	selectImage(j);
	        	title=getTitle();
	        	

		        if (matches(title,".*C1.*")==1){
		        	green=getImageID();
		        	rename("Green");
		        	run("Green");
		        	run("Subtract Background...", "rolling=30");
		   	 		getPixelSize(unit, pixelWidth, pixelHeight);

	
				}	

				if (matches(title,".*C3.*")==1){
		        	magenta=getImageID();
		        	rename("Magenta");
		        	run("Magenta");

	
				}	


				if (matches(title,".*C4.*")==1){
		        	blue=getImageID();
		        	rename("Blue");
		        	run("Blue");
		        	

	
				}	
				
		        if (matches(title,".*C2.*")==1){
					red=getImageID();
					run("Red");
					rename("Red");
					run("Red");
					run("Subtract Background...", "rolling=30");

				}

	        }


			// ************************* Cells Count ***************************************

			selectImage(blue);
			run("Duplicate...", "title=DAPI");



			//Segment nuclei
			run("Mean...", "radius=7");
			setAutoThreshold("Huang dark");
			setOption("BlackBackground", true);
			run("Convert to Mask");
			run("Fill Holes");
			run("Duplicate...", "title=NucleiMaskUnsegmented");
			run("Distance Map");
			run("Find Maxima...", "prominence=5 output=[Segmented Particles]");
			imageCalculator("Min create", "DAPI","NucleiMaskUnsegmented Segmented");
			run("Analyze Particles...", "size=2000-Infinity pixel show=Masks");
			run("Invert LUT");
			rename("Nuclei");


			//Calculate addition 0.2LD+ 0.2Mitochondria
			selectImage(green);
			run("Duplicate...", " ");
			rename("LipidDroplet");
			run("Multiply...", "value=0.20000");



			selectImage(red);
			run("Duplicate...", " ");
			rename("Mitochondria");
			run("Multiply...", "value=0.20000");

			
			
			imageCalculator("Add create", "LipidDroplet","Mitochondria");
			rename("LD_Mitochondria");


			selectImage(magenta);
			run("Duplicate...", " ");
			rename("HCSCellMask");
			run("Subtract Background...", "rolling=50");
			run("Multiply...", "value=0.60000");



			//Addition: 0.6HCSCellMask+0.2LD+0.2Mitochondria Gamma and Filter= Cell Mask
			imageCalculator("Add create", "LD_Mitochondria","HCSCellMask");
			run("Gamma...", "value=0.50");
			run("Mean...", "radius=7");
			run("Subtract...", "value=50");
			rename("CellMask");
	

			//Addition: Binary Nuclei and Cell Mask; Cell Territories Segmented Particles
			imageCalculator("Max create", "Nuclei","CellMask");
			run("Find Maxima...", "prominence=100 output=[Segmented Particles]");
			run("Analyze Particles...", "size=1000-Infinity pixel show=Masks add");

		
			// Cells Limits thresholded from the original 4 channelsaddition 
			selectImage(ImageID);
			run("RGB Color");
			run("8-bit");
			run("Median...", "radius=2");
			setAutoThreshold("Huang dark");
			setOption("BlackBackground", true);
			run("Convert to Mask");
			rename("4channels");

		
			//Segmented Cells:  Cell Territories MIN Cell limits thresholded
			imageCalculator("Min create", "Result of Nuclei Segmented","4channels");
			run("Analyze Particles...", "size=7000-Infinity pixel show=Masks clear add");



			selectImage(ImageID);
			
			run("Channels Tool...");

			roiManager("Show All");

			waitForUser("Check Cells ROIs. Erase the excluded for measurement. ClicK OK when finished");

			
			numcells=roiManager("count");

			AreaCell =newArray(numcells);

			for(k=0;k<numcells;k++){
				AreaCell[k]=getResult("Area",k);
			}

			roiManager("save", dirRes+"_ROI_Cells_"+t +".zip");

			
			print("numcells");
			print(numcells);


			IJ.deleteRows(0, nResults);

				
			// *************************  Lipid Droplet segmentation  (WEKA)  *******************************

			selectImage(green);

			run("Duplicate...", " ");

			run("Advanced Weka Segmentation");
			wait(3000);
			call("trainableSegmentation.Weka_Segmentation.loadClassifier", dir+"model\\classifierLD.model");
			wait(5000);
			call("trainableSegmentation.Weka_Segmentation.loadData", dir+"model\\dataLD.arff");
			wait(5000);
			call("trainableSegmentation.Weka_Segmentation.trainClassifier");
			wait(10000);
			call("trainableSegmentation.Weka_Segmentation.getResult");
			wait(3000);
			
			selectWindow("Classified image");
	
			run("8-bit");
			run("Threshold...");
			setThreshold(0, 135);
			setOption("BlackBackground", false);
			run("Convert to Mask");
	

			run("Set Measurements...", "area mean min perimeter shape integrated display redirect=None decimal=5");


			
			CircLDCell =newArray(numcells);
			AreaLDCell =newArray(numcells);
			PerimLDCell =newArray(numcells);
			CircdvLDCell =newArray(numcells);
			AreadvLDCell =newArray(numcells);
			PerimdvLDCell =newArray(numcells);
			numLDs =newArray(numcells);

			
			
			for(k=0;k<numcells;k++){

				roiManager("select", k);

				run("Analyze Particles...", "size=5-Infinity pixel display exclude clear");
				numLD=nResults;

								
				CircLD =newArray(numLD);
				AreaLD =newArray(numLD);
				PerimLD =newArray(numLD);

				for (j=0;j<numLD;j++){

					CircLD[j]=getResult("Circ.",j);
					AreaLD[j]=getResult("Area",j);
					PerimLD[j]=getResult("Perim.", j);

				
				}

				Array.getStatistics(AreaLD, min, max, meanareaLD, stddesvLD);
				Array.getStatistics(CircLD, min, max, meancircLD, stddesvcircLD);
				Array.getStatistics(PerimLD, min, max, meanperimLD, stddesvperimLD);



				CircLDCell[k] =meancircLD;
				AreaLDCell[k] = meanareaLD;
				PerimLDCell[k] =meanperimLD;
				CircdvLDCell[k] =stddesvcircLD;
				AreadvLDCell[k] =stddesvLD;
				PerimdvLDCell[k] =stddesvperimLD;
				numLDs[k] =numLD;

				
			}
			
			roiManager("Deselect");
			//roiManager("Show None");
			//roiManager("Show All");
				
			IJ.deleteRows(0, nResults);
			
			
			// *************************  Mitochondria segmentation    *******************************
			
			
			selectImage(red);
			
			run("Duplicate...", " ");


			AreaMitoCell =newArray(numcells);
			MeangreyMitoCell =newArray(numcells);
			stddesvMitoCell =newArray(numcells);
			RawIntDenMitoCell =newArray(numcells);
			medianMitoCell =newArray(numcells);
			

			for(k=0;k<numcells;k++){

				roiManager("select", k);

				run("Set Measurements...", "area limit mean standard modal min perimeter shape integrated median display redirect=None decimal=5");
				setAutoThreshold("Otsu dark");

				roiManager("Measure");

				
				AreaMitoCell[k]=getResult("Area",k);
				MeangreyMitoCell[k]=getResult("Mean",k);
				stddesvMitoCell[k]=getResult("StdDev",k);
				RawIntDenMitoCell[k]=getResult("RawIntDen",k);
				medianMitoCell[k]=getResult("Median",k);

			}

			
			IJ.deleteRows(0, nResults);


			// ******************************** Magenta intensity ***************************************

			selectImage(magenta);

			run("Subtract Background...", "rolling=50");

			Meantransf =newArray(numcells);

			for(k=0;k<numcells;k++){



				run("Set Measurements...", "mean redirect=None decimal=5");

				roiManager("select", k);
				roiManager("Measure");

				Meantransf[k]=getResult("Mean",k);

			}

			
			IJ.deleteRows(0, nResults);

			// ******************************** Colocalization Mitochondria Lipid Droplets ***************************************


			run("Set Measurements...", "area limit mean standard modal min perimeter shape integrated median display redirect=None decimal=5");

			
			selectImage(green);

			setAutoThreshold("Yen dark");

			getThreshold(thresholdgreen,thresholdgreen2);


			selectImage(red);


			setAutoThreshold("Otsu dark");

			getThreshold(thresholdred,thresholdred2);
			

			run("Colocalization Highligter", "channel_1=[Red] channel_2=[Green] ratio=50 threshold_channel_1="+thresholdred+" threshold_channel_2="+thresholdgreen+" display=255");


			run("Split Channels");
			
			selectWindow("Colocalized points (RGB)  (blue)");

			run("Invert LUT");

			run("Skeletonize");


			colocpixels =newArray(numcells);
			colocmicron =newArray(numcells);
			numcontact =newArray(numcells);

			for(k=0;k<numcells;k++){

				roiManager("select", k);

				run("Analyze Particles...", "size=0-Infinity display exclude clear");

				numcontacts=nResults;
				numcontact[k]=numcontacts;

				Areacontact =newArray(numcontacts);


				for(n=0;n<numcontacts;n++){

					Areacontact[n]=getResult("Area",n);
					
				}

				Array.getStatistics(Areacontact, min, max, meanareacontact, stddesvcontact);
				
				colocpixels[k]=meanareacontact;
				colocmicron[k]=colocpixels[k]*pixelWidth;
			}

			roiManager("reset"); //to delete previous ROIs

			run("Create Selection");
			roiManager("Add");	

			roiManager("save", dirRes+"_ROI_contact_"+t +".zip");
			IJ.deleteRows(0, nResults);

			
			// ********************************  measurements ****************************************

			for(k=0;k<numcells;k++){
				
				File.append( t + "\t" + k+1 + "\t" + AreaCell[k] + "\t" + AreaMitoCell[k] + "\t" + MeangreyMitoCell[k] + "\t" + stddesvMitoCell[k] + "\t" + RawIntDenMitoCell[k] + "\t" + medianMitoCell[k] + "\t" + numLDs[k] + "\t" + AreaLDCell[k] + "\t" + AreadvLDCell[k] + "\t" + CircLDCell[k] + "\t" + PerimLDCell[k] + "\t" + PerimLDCell[k]*numLDs[k] + "\t" +  numcontact[k] + "\t" + colocmicron[k] + "\t" + colocmicron[k]*numcontact[k] + "\t" + Meantransf[k], path);
				//File.append( "Image \t # Cells \t Area Cell \t Area Mito (um2) \t Mean Grey Value \t Std desv. Grey Value \t RawIntDen \t Median \t # LD \t Mean Area LD (um2) \t Std desv. Area LD \t Mean Circ. LD \t Mean Perim. LD (um) \t Perim. total LD (um) \t # Contact \t Mean Length contact (um) \t Total length contact (um) \t Mean. Transfection", path);


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


	
