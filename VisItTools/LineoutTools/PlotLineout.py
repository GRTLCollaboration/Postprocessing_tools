#
## Run this using 'visit -nowin -cli -s PlotLineout.py'
# --------------------------------------------------------------------------
## Inputs - you need to update these!
# --------------------------------------------------------------------------

# file details
begin_file = 0
end_file = 100
file_step = 20
plt_prefix = "name_of_your_hdf5_files"
path_to_hdf5_files = "path_to_hdf5s"

# plot details
# select variable 
plot_variable = "phi"
outputfile_prefix = plot_variable + '_'
output_directory = "."
# slice origin and normal direction
origin_point_x = 0 # middle of the box is [L/2, L/2, L/2]
origin_point_y = 0
origin_point_z = 0
normal_in_x = 0
normal_in_y = 0
normal_in_z = 1
# specify a begin and end point for the lineout, e.g. two (x,y) points if you sliced normal to z
[u_min, v_min] = [0, 0]
[u_max, v_max] = [40, 40]

# set the range of your plots manually, if you want
manual_range = 0 # 0 for automatic range, 1 for manual range
range_min = -.003
range_max = 0.07

# --------------------------------------------------------------------------
# From here you only need to amend if you want to do something non standard
# --------------------------------------------------------------------------

# Loop through the files, and plot each slice
for i in range ( begin_file , end_file + file_step , file_step) :
	hdf5filename = plt_prefix + "%06i.3d.hdf5"%i
	print("Plotting file " + hdf5filename)
	OpenDatabase(path_to_hdf5_files + hdf5filename)
	
	# annotation settings
	AnnotationAtts = AnnotationAttributes()
	AnnotationAtts.axes3D.visible = 1
	AnnotationAtts.userInfoFlag = 0
	AnnotationAtts.databaseInfoFlag = 1
	AnnotationAtts.timeInfoFlag = 1
	AnnotationAtts.legendInfoFlag = 1
	AnnotationAtts.axesArray.visible = 0
	AnnotationAtts.backgroundColor = (0, 0, 0, 255)
	AnnotationAtts.foregroundColor = (255, 255, 255, 255)
	AnnotationAtts.axes3D.bboxFlag = 1
	AnnotationAtts.axes3D.triadFlag = 0
	SetAnnotationAttributes(AnnotationAtts)
	
	# save settings
	InvertBackgroundColor()
	saveAtts = SaveWindowAttributes()
	saveAtts.format = saveAtts.PNG
	saveAtts.resConstraint = saveAtts.NoConstraint
	saveAtts.width = 1280
	saveAtts.height = 620
	SetSaveWindowAttributes(saveAtts)
	
	# add pseudocolour plot
	AddPlot("Pseudocolor", plot_variable, 1, 1)
	PseudocolorAtts = PseudocolorAttributes()
	PseudocolorAtts.scaling = PseudocolorAtts.Linear # Linear, Log, Skew
	PseudocolorAtts.skewFactor = 1
     # limitsMode =  OriginalData, CurrentPlot
	PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData
	#PseudocolorAtts.colorTableName = "BuGn"
	PseudocolorAtts.invertColorTable = 0
     # opacityType = ColorTable, FullyOpaque, Constant, Ramp, VariableRange
	PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque 
	PseudocolorAtts.smoothingLevel = 0
     # centering = Natural, Nodal, Zonal
	PseudocolorAtts.centering = PseudocolorAtts.Nodal
	SetPlotOptions(PseudocolorAtts)
	
	# slice the pseudocolour plot
	AddOperator("Slice", 1)
	SliceAtts = SliceAttributes()
     # originType = Point, Intercept, Percent, Zone, Node
	SliceAtts.originType = SliceAtts.Point  
	SliceAtts.originPoint = (origin_point_x, origin_point_y, origin_point_z)
     # axisType = XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.axisType = SliceAtts.Arbitrary  
	SliceAtts.normal = (normal_in_x, normal_in_y, normal_in_z)
	SliceAtts.project2d = 1
	SliceAtts.flip = 0
	SetOperatorOptions(SliceAtts, 1)

	# plot all levels
	silr = SILRestriction()
	silr.TurnOnAll()
	SetPlotSILRestriction(silr ,1)
	
	# engage!
	DrawPlots()

	Lineout((u_min, v_min), (u_max, v_max))
	SetActiveWindow(2)

	# make the line black
	CurveAtts = CurveAttributes()
	CurveAtts.curveColorSource = CurveAtts.Custom  # Cycle, Custom
	CurveAtts.showLabels = 0
	CurveAtts.curveColor = (0, 0, 0, 255)
	SetPlotOptions(CurveAtts)

	if manual_range == 1:
		ViewCurveAtts = ViewCurveAttributes()
		ViewCurveAtts.domainCoords = (0, ((u_max - u_min)**2. + (v_max - v_min)**2.)**(1/2.))
		ViewCurveAtts.rangeCoords = (range_min, range_max)
		SetViewCurve(ViewCurveAtts)
	
	# Save Me
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputDirectory = output_directory
	SaveWindowAtts.fileName = outputfile_prefix
	SaveWindowAtts.family = 1
     # format = BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, 
     # RGB, STL, TIFF, ULTRA, VTK, PLY
	SaveWindowAtts.format = SaveWindowAtts.PNG
	SaveWindowAtts.width = 1024
	SaveWindowAtts.height = 1024
	SaveWindowAtts.quality = 80
	# resConstraint = NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.resConstraint = SaveWindowAtts.EqualWidthHeight 
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()

	# clear plots ready for next timestep in loop
	DeleteAllPlots()
	SetActiveWindow(1)
	DeleteAllPlots()
	CloseDatabase(path_to_hdf5_files + hdf5filename)

# loop ends
print ("Done and dusted!")