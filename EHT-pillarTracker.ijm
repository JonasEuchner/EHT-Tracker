// Ask the user to select a directory
dir = getDirectory("Choose a Directory ");
// Start processing files in the selected directory
processFiles(dir);

// Recursive function to process files and subdirectories
function processFiles(dir) {
	list = getFileList(dir);
	for (i=0; i<list.length; i++) {
		// If it's a subdirectory, recursively process it
		if (endsWith(list[i], "/"))
			processFiles(""+dir+list[i]);
		else {
			// If it's a file, process it
			path = dir+list[i];
			processFile(dir,path);
		}
	}
}

// Function to process a single file
function processFile(dir,path) {
	if (endsWith(path, ".tif")) {
		open(path);
		// Call the PillarTracker function to analyze it
		PillarTracker(dir);
	}
}

// Main function for tracking pillars in an image
function PillarTracker(dir){	
	ImageID=getImageID();
	title=getTitle();
	slices=nSlices;
	run("8-bit");
	waitForUser("Select Pillars with Multi-Point Tool");   
	RegionMatching=70;
	SearchRadiusX=15;
	SearchRadiusY=5;
	getPixelSize(unit, pixelWidth, pixelHeight);
	
	// Create a dialog to get analysis parameters from the user
	Dialog.create("Analysis parameter wizard");
		Dialog.addMessage("Please make sure both pillars have a tracking pointer before continuing");
		Dialog.addString("Maximal X-axis frame to frame movement in px", SearchRadiusX);
		Dialog.addString("Maximal Y-axis frame to frame movement in px", SearchRadiusY);
		Dialog.addString("Pattern matching size in px", RegionMatching);
		Dialog.addString("Pixelsize [mm/px]", pixelWidth);	
		Dialog.show();
	
	// Retrieve analysis parameters from the user
	SearchRadiusX = Dialog.getString();
	SearchRadiusY = Dialog.getString();
	RegionMatching = Dialog.getString();
	pixelWidth = Dialog.getString();
	RegionMatching=2*RegionMatching+1;
	
	getSelectionCoordinates(x, y);
	setBatchMode(true);
	for (p = 0;p < x.length; p++) {
		yRef=y[p];
		xRef=x[p];
		xPrev=xRef;
		yPrev=yRef;
		xTotal=0;
		yTotal=0;
		selectImage(ImageID);
		setSlice(1);
		run("Select None");
		// Create a rectangle around the tracking point
		makeRectangle(x[p]-floor(RegionMatching/2),y[p]-floor(RegionMatching/2), RegionMatching, RegionMatching);
		run("Duplicate...", "use");
		rename("Reference");
		RefID=getImageID();
		print("Distance;xPos;yPos");
		// Loop through slices of the image
		for (i = 1; i <= slices; i++) {
			selectImage(ImageID);
			setSlice(i);
			// Create a search area around the tracking point
			makeRectangle(x[p]-floor(RegionMatching/2)-SearchRadiusX, y[p]-floor(RegionMatching/2)-SearchRadiusY, (RegionMatching+2*SearchRadiusX), (RegionMatching+2*SearchRadiusY));
			run("Duplicate...", "use");
			rename("tempSearch");
			SearchID=getImageID();
			newImage("tempStack", "8-bit grayscale-mode", RegionMatching, RegionMatching, 1, 1, ((2*SearchRadiusX)+1)*((2*SearchRadiusY)+1));
			StackID=getImageID();
			// Create a stack from the search area of all possibly ref positions
			for (j = 0; j < (2*SearchRadiusX+1); j++) {
				for (k = 0; k < (2*SearchRadiusY+1); k++) {
					selectImage(SearchID);
					makeRectangle(j, k, RegionMatching, RegionMatching);
					run("Copy");
					selectWindow("tempStack");
					setSlice((k+1)+(j*((2*SearchRadiusY)+1)));
					run("Paste");
				}
			}
			// Calculate the difference between the search area stack and the reference
			imageCalculator("Difference create stack", "tempStack","Reference");
			selectWindow("Result of tempStack");
			rename("tempResults");
			ResultsID=getImageID();
			meanMin=255;
			// Find the frame that best matches the reference (best match)
			for (j = 1; j <= ((2*SearchRadiusX)+1)*((2*SearchRadiusY)+1); j++) {
				setSlice(j);
				getStatistics(area, mean);
				if (mean<meanMin){
					minJ=j-1;
					meanMin=mean;
				}
			}
			selectImage(StackID);
			// Calculate displacement and position of marker
			yCur=minJ%(2*SearchRadiusY+1);
			xCur=((minJ-yCur)/(2*SearchRadiusY+1))+x[p]-SearchRadiusX;
			yCur=yCur+y[p]-SearchRadiusY;
			xDist=(xCur-xPrev);
			yDist=(yCur-yPrev);
	
			xTotal=xDist+xTotal;
			yTotal=yDist+yTotal;
			distance=sqrt((xTotal*xTotal*pixelWidth*pixelWidth)+(yTotal*yTotal*pixelWidth*pixelWidth));
			xPrev=xCur;
			yPrev=yCur;
			// Print the distance and current position
			print(distance+";"+xCur+";"+yCur);
			selectImage(ImageID);
			close("temp*");
			x[p]=x[p]+xDist;
			y[p]=y[p]+yDist;
		}
		close("Ref*");
		// Save the results as a text file
		selectWindow("Log");
		name=split(title, ".");
		root=File.getParent(dir);
		saveAs("Text", ""+root+"/"+name[0]+"_Pillar"+(p+1)+".txt");
		close("Log");
	}
	run("Close All");
	setBatchMode(false);
}