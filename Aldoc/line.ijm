nR = roiManager("Count"); 
save_path = "D:/mouse_T357.txt"
roiManager("Select", 0); 
rName = Roi.getName(); 
//rName = split(rName, "-");
//rName = rName[0];
line_id = rName + "+" + 0;
//run("Interpolate", "interval=0.1");
getSelectionCoordinates(xpoints1, ypoints1);
a11 = newArray(xpoints1.length);
for (i=0; i<a11.length; i++)
  a11[i] = line_id;

if(nR>1){
	for (j=1; j<nR; j++) { 
		roiManager("Select", j); 
		rName = Roi.getName(); 
		//rName = split(rName, "-");
		//rName = rName[0];
		line_id = rName+ "+" + j;
		getSelectionCoordinates(xpoints, ypoints);
		a1 = newArray(xpoints.length);
		for (i=0; i<a1.length; i++)
		  a1[i] = line_id;
		xpoints1 = Array.concat(xpoints1,xpoints);
		ypoints1 = Array.concat(ypoints1,ypoints);
		a11 = Array.concat(a11,a1);
	}
}
Table.create("Points");
Table.setColumn("x", xpoints1);
Table.setColumn("y", ypoints1);
Table.setColumn("line_id", a11)
Table.save(save_path);
//roiManager("deselect");
//roiManager("delete");