# JCAurre 10/12/2019
begin = 0
end = 350

DefineVectorExpression("Coords", "coord(Mesh)")
DefineScalarExpression("X", "Coords[0]")
DefineScalarExpression("Y", "Coords[1]")
DefineScalarExpression("Dist", "sqrt(X*X+Y*Y)")
DefineScalarExpression("rPsi4", "Dist*Weyl4_Re")

def annot_ini():

	AnnotationAtts = AnnotationAttributes()
	AnnotationAtts.axes3D.visible = 0
	AnnotationAtts.userInfoFlag = 0
	AnnotationAtts.databaseInfoFlag = 0
	AnnotationAtts.timeInfoFlag = 1
	AnnotationAtts.legendInfoFlag = 0
	AnnotationAtts.axesArray.visible = 0
	AnnotationAtts.backgroundColor = (0, 0, 0, 255)
	AnnotationAtts.foregroundColor = (255, 255, 255, 255)
	AnnotationAtts.axes3D.bboxFlag = 1
	AnnotationAtts.axes3D.triadFlag = 0
	SetAnnotationAttributes(AnnotationAtts)

def annot_evo():

	AnnotationAtts = AnnotationAttributes()
	AnnotationAtts.axes3D.visible = 0
	AnnotationAtts.userInfoFlag = 0
	AnnotationAtts.databaseInfoFlag = 0
	AnnotationAtts.timeInfoFlag = 1
	AnnotationAtts.legendInfoFlag = 0
	AnnotationAtts.axesArray.visible = 0
	AnnotationAtts.backgroundColor = (0, 0, 0, 255)
	AnnotationAtts.foregroundColor = (255, 255, 255, 255)
	AnnotationAtts.axes3D.bboxFlag = 0
	AnnotationAtts.axes3D.triadFlag = 0
	SetAnnotationAttributes(AnnotationAtts)

def window_options():
	InvertBackgroundColor()
	saveAtts = SaveWindowAttributes()
	saveAtts.format = saveAtts.PNG
	saveAtts.resConstraint = saveAtts.NoConstraint
	saveAtts.width = 1280
	saveAtts.height = 620
	SetSaveWindowAttributes(saveAtts)

	
def rho():
    AddPlot("Volume", "rho", 1, 0)
    VolumeAtts = VolumeAttributes()
    VolumeAtts.resampleTarget = 1000000
    VolumeAtts.useColorVarMin = 1
    VolumeAtts.colorVarMin = -0
    VolumeAtts.useColorVarMax = 1
    VolumeAtts.colorVarMax = 1e-06
    SetPlotOptions(VolumeAtts)

    #box
    AddOperator("Box", 0)
    BoxAtts = BoxAttributes()
    BoxAtts.minx = 0
    BoxAtts.maxx = 700
    BoxAtts.miny = 0
    BoxAtts.maxy = 700
    BoxAtts.minz = 0
    BoxAtts.maxz = 100
    SetOperatorOptions(BoxAtts, 0)

    #reflect
    AddOperator("Reflect", 0)
    ReflectAtts = ReflectAttributes()
    ReflectAtts.octant = ReflectAtts.PXPYPZ
    ReflectAtts.useXBoundary = 1
    ReflectAtts.specifiedX = 0 
    ReflectAtts.useYBoundary = 1
    ReflectAtts.specifiedY = 0 
    ReflectAtts.useZBoundary = 1
    ReflectAtts.specifiedZ = 0 
    ReflectAtts.reflections = (1, 1, 1, 1, 1, 1, 1, 1)
    SetOperatorOptions(ReflectAtts, 0)

def gws(time):

	AddPlot("Pseudocolor", "rPsi4", 1, 0)
	PseudocolorAtts = PseudocolorAttributes()
	PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
	PseudocolorAtts.skewFactor = 1
	PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
	PseudocolorAtts.minFlag = 1
	PseudocolorAtts.min = 0
	PseudocolorAtts.maxFlag = 1
	PseudocolorAtts.max = 1e-3
	PseudocolorAtts.colorTableName = "gws"
	PseudocolorAtts.invertColorTable = 0
	PseudocolorAtts.opacityType = PseudocolorAtts.ColorTable  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
	PseudocolorAtts.smoothingLevel = 0
	PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
	SetPlotOptions(PseudocolorAtts)


	AddOperator("Box", 0)
	BoxAtts = BoxAttributes()
	BoxAtts.minx = 0
	BoxAtts.maxx = 3000
	BoxAtts.miny = 0
	BoxAtts.maxy = 3000
	BoxAtts.minz = 0
	BoxAtts.maxz = 100
	SetOperatorOptions(BoxAtts, 0) 

	AddOperator("Cylinder", 0)
	CylinderAtts = CylinderAttributes()
	CylinderAtts.point1 = (0, 0, 0)
	CylinderAtts.point2 = (0, 0, 100)
	CylinderAtts.radius = time - 200
	CylinderAtts.inverse = 0
	SetOperatorOptions(CylinderAtts, 0)

	AddOperator("Clip", 0)
	ClipAtts = ClipAttributes()
	ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
	ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
	ClipAtts.center = (0, 0, 0)
	ClipAtts.radius = 100
	ClipAtts.sphereInverse = 0
	SetOperatorOptions(ClipAtts, 0)


	AddOperator("Slice", 0)
	SliceAtts = SliceAttributes()
	SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
	SliceAtts.originPoint = (0, 0, 0)
	SliceAtts.normal = (0, 0, 1)
	SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.upAxis = (0, 1, 0)
	SliceAtts.project2d = 1
	SetOperatorOptions(SliceAtts, 0)

 	AddOperator("Elevate", 0)
	ElevateAtts = ElevateAttributes()
	ElevateAtts.useXYLimits = 1
	ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
	ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
	ElevateAtts.minFlag = 1
	ElevateAtts.min = 0
	ElevateAtts.maxFlag = 1
	ElevateAtts.max = 9e-3
	ElevateAtts.zeroFlag = 0
	ElevateAtts.variable = "default"
	SetOperatorOptions(ElevateAtts, 0)

	AddOperator("Reflect", 0)
	ReflectAtts = ReflectAttributes()
	ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
	ReflectAtts.useXBoundary = 1
	ReflectAtts.specifiedX = 0
	ReflectAtts.useYBoundary = 1
	ReflectAtts.specifiedY = 0
	ReflectAtts.useZBoundary = 1
	ReflectAtts.specifiedZ = 0
	ReflectAtts.reflections = (1, 1, 1, 1, 0, 0, 0, 0)
	SetOperatorOptions(ReflectAtts, 0)

	AddOperator("Smooth", 0)
	SmoothOperatorAtts = SmoothOperatorAttributes()
	SmoothOperatorAtts.numIterations = 100
	SmoothOperatorAtts.relaxationFactor = 0
	SmoothOperatorAtts.convergence = 0
	SmoothOperatorAtts.maintainFeatures = 0
	SmoothOperatorAtts.featureAngle = 45
	SmoothOperatorAtts.edgeAngle = 15
	SmoothOperatorAtts.smoothBoundaries = 1
	SetOperatorOptions(SmoothOperatorAtts, 0)

	AddOperator("Transform", 0)
	TransformAtts = TransformAttributes()   																																					  
	TransformAtts.doTranslate = 1
	TransformAtts.translateZ = -250
	SetOperatorOptions(TransformAtts, 0)


def bh():

	AddPlot("Pseudocolor", "chi", 1, 0)
	PseudocolorAtts = PseudocolorAttributes()
	PseudocolorAtts.minFlag = 1
	PseudocolorAtts.min = 0.86
	PseudocolorAtts.maxFlag = 0
	PseudocolorAtts.max = 1
	PseudocolorAtts.colorTableName = "xray"
	PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
	SetPlotOptions(PseudocolorAtts)

	#box
	AddOperator("Box", 0)
	BoxAtts = BoxAttributes()
	BoxAtts.minx = 0
	BoxAtts.maxx = 700
	BoxAtts.miny = 0
	BoxAtts.maxy = 700
	BoxAtts.minz = 0
	BoxAtts.maxz = 100
	SetOperatorOptions(BoxAtts, 0)

	AddOperator("Isosurface", 0)
	IsosurfaceAtts = IsosurfaceAttributes()
	IsosurfaceAtts.contourValue = (0.9)
	IsosurfaceAtts.contourMethod = IsosurfaceAtts.Value  # Level, Value, Percent
	IsosurfaceAtts.minFlag = 0
	IsosurfaceAtts.min = 0
	IsosurfaceAtts.maxFlag = 0
	IsosurfaceAtts.max = 1
	SetOperatorOptions(IsosurfaceAtts, 0)

	#reflect
	AddOperator("Reflect", 0)
	ReflectAtts = ReflectAttributes()
	ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
	ReflectAtts.useXBoundary = 1
	ReflectAtts.specifiedX = 0 
	ReflectAtts.useYBoundary = 1
	ReflectAtts.specifiedY = 0 
	ReflectAtts.useZBoundary = 1
	ReflectAtts.specifiedZ = 0 
	ReflectAtts.reflections = (1, 1, 1, 1, 1, 1, 1, 1)
	SetOperatorOptions(ReflectAtts, 0)

	AddOperator("Transform", 0)
	TransformAtts = TransformAttributes()   																																					  
	TransformAtts.doTranslate = 1
	TransformAtts.translateZ = 0
	SetOperatorOptions(TransformAtts, 0)

def fly_evo(k,hdf5file):

	
	c0 = View3DAttributes()
	c0.viewNormal = (-0.15,-0.3, 0.2)
	c0.focus = (0,0,0)
	c0.viewUp = (0, 0, 1)
	c0.viewAngle = 10
	c0.parallelScale = 2660.43
	c0.nearPlane = -3000 #-5320.86
	c0.farPlane = 3000 #5320.86
	c0.imageZoom = 4
	c0.perspective = 1

	c1 = View3DAttributes()
	c1.viewNormal = (-0.15,-0.3, 0.2)
	c1.focus = (0,0,0)
	c1.viewUp = (0, 0, 1)
	c1.viewAngle = 10
	c1.parallelScale = 2660.43
	c1.nearPlane = -3000 #-5320.86
	c1.farPlane = 3000 #5320.86
	c1.imageZoom = 3
	c1.perspective = 1


	c2 = View3DAttributes()
	c2.viewNormal = (-0.,0, 1)
	c2.focus = (0,0,0)
	c2.viewUp = (0, 1, 0)
	c2.viewAngle = 10
	c2.parallelScale = 2660.43
	c2.nearPlane = -3000 #-5320.86
	c2.farPlane = 3000 #5320.86
	c2.imageZoom = 2
	c2.perspective = 1


	cpts = (c0,c2)
	x=[]
	for i in range(2):
		x = x + [float(i) / float(1.)]

	SetView3D(c0)
	if k > 270 and k < 341:
		t = float(k-271) / float(70 - 1)
		c = EvalCubicSpline(t, x, cpts)
		c.nearPlane = -3000
		c.farPlane = 3000
		SetView3D(c)
	if k>341:
		SetView3D(c2)



def meshgw(time,i):
	
	AddPlot("Mesh", "Mesh", 1, 0)
	MeshAtts = MeshAttributes()
	MeshAtts.showInternal = 0
	MeshAtts.meshColorSource = MeshAtts.Foreground  # Foreground, MeshCustom
	MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
	MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
	MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
	MeshAtts.opacity = 0.15*(i-141)/50
	SetPlotOptions(MeshAtts)

	AddOperator("MultiresControl", 0)

	AddOperator("Box", 0)
	BoxAtts = BoxAttributes()
	BoxAtts.minx = 0
	BoxAtts.maxx = 3000
	BoxAtts.miny = 0
	BoxAtts.maxy = 3000
	BoxAtts.minz = 0
	BoxAtts.maxz = 100
	SetOperatorOptions(BoxAtts, 1) 

	AddOperator("Cylinder", 0)
	CylinderAtts = CylinderAttributes()
	CylinderAtts.point1 = (0, 0, 0)
	CylinderAtts.point2 = (0, 0, 100)
	CylinderAtts.radius = time
	CylinderAtts.inverse = 0
	SetOperatorOptions(CylinderAtts, 0)

	AddOperator("Clip", 0)
	ClipAtts = ClipAttributes()
	ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
	ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
	ClipAtts.center = (0, 0, 0)
	ClipAtts.radius = 100
	ClipAtts.sphereInverse = 0
	SetOperatorOptions(ClipAtts, 0)

	AddOperator("Slice", 0)
	SliceAtts = SliceAttributes()
	SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
	SliceAtts.originPoint = (0, 0, 0)
	SliceAtts.normal = (0, 0, 1)
	SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.upAxis = (0, 1, 0)
	SliceAtts.project2d = 1
	SetOperatorOptions(SliceAtts, 0)

 	AddOperator("Elevate", 0)
	ElevateAtts = ElevateAttributes()
	ElevateAtts.useXYLimits = 1
	ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
	ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
	ElevateAtts.minFlag = 1
	ElevateAtts.min = 0
	ElevateAtts.maxFlag = 1
	ElevateAtts.max = 9e-3
	ElevateAtts.zeroFlag = 0
	ElevateAtts.variable = "rPsi4"
	SetOperatorOptions(ElevateAtts, 0)

	AddOperator("Reflect", 0)
	ReflectAtts = ReflectAttributes()
	ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
	ReflectAtts.useXBoundary = 1
	ReflectAtts.specifiedX = 0
	ReflectAtts.useYBoundary = 1
	ReflectAtts.specifiedY = 0
	ReflectAtts.useZBoundary = 1
	ReflectAtts.specifiedZ = 0
	ReflectAtts.reflections = (1, 1, 1, 1, 0, 0, 0, 0)
	SetOperatorOptions(ReflectAtts, 0)

	AddOperator("Smooth", 0)
	SmoothOperatorAtts = SmoothOperatorAttributes()
	SmoothOperatorAtts.numIterations = 100
	SmoothOperatorAtts.relaxationFactor = 0
	SmoothOperatorAtts.convergence = 0
	SmoothOperatorAtts.maintainFeatures = 0
	SmoothOperatorAtts.featureAngle = 45
	SmoothOperatorAtts.edgeAngle = 15
	SmoothOperatorAtts.smoothBoundaries = 1
	SetOperatorOptions(SmoothOperatorAtts, 0)

	AddOperator("Transform", 0)
	TransformAtts = TransformAttributes()   																																					  
	TransformAtts.doTranslate = 1
	TransformAtts.translateZ = -250
	SetOperatorOptions(TransformAtts, 0)

def meshgwflat(time,i):
	
	AddPlot("Mesh", "Mesh", 1, 0)
	MeshAtts = MeshAttributes()
	MeshAtts.showInternal = 0
	MeshAtts.meshColorSource = MeshAtts.Foreground  # Foreground, MeshCustom
	MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
	MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
	MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
	MeshAtts.opacity = 0.15*(i-141)/50
	SetPlotOptions(MeshAtts)

	AddOperator("MultiresControl", 0)

	AddOperator("Box", 0)
	BoxAtts = BoxAttributes()
	BoxAtts.minx = 0
	BoxAtts.maxx = 3000
	BoxAtts.miny = 0
	BoxAtts.maxy = 3000
	BoxAtts.minz = 0
	BoxAtts.maxz = 100
	SetOperatorOptions(BoxAtts, 1) 

	AddOperator("Clip", 0)
	ClipAtts = ClipAttributes()
	ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
	ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
	ClipAtts.center = (0, 0, 0)
	ClipAtts.radius = time
	ClipAtts.sphereInverse = 0
	SetOperatorOptions(ClipAtts, 0)

	AddOperator("Slice", 0)
	SliceAtts = SliceAttributes()
	SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
	SliceAtts.originPoint = (0, 0, 0)
	SliceAtts.normal = (0, 0, 1)
	SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.upAxis = (0, 1, 0)
	SliceAtts.project2d = 1
	SetOperatorOptions(SliceAtts, 0)

	AddOperator("Elevate", 0)
	ElevateAtts = ElevateAttributes()
	ElevateAtts.useXYLimits = 1
	ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
	ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
	ElevateAtts.minFlag = 1
	ElevateAtts.min = 0
	ElevateAtts.maxFlag = 1
	ElevateAtts.max = 9e-3
	ElevateAtts.zeroFlag = 1
	ElevateAtts.variable = "rPsi4"
	SetOperatorOptions(ElevateAtts, 0)

	AddOperator("Reflect", 0)
	ReflectAtts = ReflectAttributes()
	ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
	ReflectAtts.useXBoundary = 1
	ReflectAtts.specifiedX = 0
	ReflectAtts.useYBoundary = 1
	ReflectAtts.specifiedY = 0
	ReflectAtts.useZBoundary = 1
	ReflectAtts.specifiedZ = 0
	ReflectAtts.reflections = (1, 1, 1, 1, 0, 0, 0, 0)
	SetOperatorOptions(ReflectAtts, 0)

	AddOperator("Smooth", 0)
	SmoothOperatorAtts = SmoothOperatorAttributes()
	SmoothOperatorAtts.numIterations = 100
	SmoothOperatorAtts.relaxationFactor = 0
	SmoothOperatorAtts.convergence = 0
	SmoothOperatorAtts.maintainFeatures = 0
	SmoothOperatorAtts.featureAngle = 45
	SmoothOperatorAtts.edgeAngle = 15
	SmoothOperatorAtts.smoothBoundaries = 1
	SetOperatorOptions(SmoothOperatorAtts, 0)

	AddOperator("Transform", 0)
	TransformAtts = TransformAttributes()   																																					  
	TransformAtts.doTranslate = 1
	TransformAtts.translateZ = -250
	SetOperatorOptions(TransformAtts, 0)


def mesh():
	AddPlot("Mesh", "Mesh", 1, 0)
	AddOperator("Isosurface", 0)
	IsosurfaceAtts = IsosurfaceAttributes()
	IsosurfaceAtts.contourPercent = (20)
	IsosurfaceAtts.contourMethod = IsosurfaceAtts.Percent  # Level, Value, Percent
	IsosurfaceAtts.variable = "rho"
	SetOperatorOptions(IsosurfaceAtts, 0)																					   
	MeshAtts = MeshAttributes()
	MeshAtts.opacity = 0.05
	SetPlotOptions(MeshAtts)

	#box
	AddOperator("Box", 0)
	BoxAtts = BoxAttributes()
	BoxAtts.minx = 0
	BoxAtts.maxx = 700
	BoxAtts.miny = 0
	BoxAtts.maxy = 700
	BoxAtts.minz = 0
	BoxAtts.maxz = 100
	SetOperatorOptions(BoxAtts, 0)
	
	#reflect
	AddOperator("Reflect", 0)
	ReflectAtts = ReflectAttributes()
	ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
	ReflectAtts.useXBoundary = 1
	ReflectAtts.specifiedX = 0
	ReflectAtts.useYBoundary = 1
	ReflectAtts.specifiedY = 0
	ReflectAtts.useZBoundary = 1
	ReflectAtts.specifiedZ = 0
	ReflectAtts.reflections = (1, 1, 1, 1, 1, 1, 1, 1)
	SetOperatorOptions(ReflectAtts, 0)
	

def fly_ini():

	# Create the control points for the views.
	c0 = View3DAttributes()
	c0.viewNormal = (0, 0, 1)
	c0.focus = (0, 0, 0)
	c0.viewUp = (0, 0, 1)
	c0.viewAngle = 10
	c0.parallelScale = 2660.43
	c0.nearPlane = -3000 #-5320.86
	c0.farPlane = 3000 #5320.86
	c0.imageZoom = 2
	c0.perspective = 1

	c1 = View3DAttributes()
	c1.viewNormal = (0., -1, 0.3)
	c1.focus = (0, 0, 0)
	c1.viewUp = (0, 0, 1)
	c1.viewAngle = 1
	c1.parallelScale = 2660.43
	c1.nearPlane = -3000 #-5320.86
	c1.farPlane = 3000 #5320.86
	c1.imageZoom = 4
	c1.perspective = 1
 
	c2 = View3DAttributes()
	c2.viewNormal = (0.3, -0.5, 0.3)
	c2.focus = (0,0,0)
	c2.viewUp = (0, 0, 1)
	c2.viewAngle = 10
	c2.parallelScale = 2660.43
	c2.nearPlane = -3000 #-5320.86
	c2.farPlane = 3000 #5320.86
	c2.imageZoom = 3
	c2.perspective = 1
 
	c3 = View3DAttributes()
	c3.viewNormal = (0.5, -0.2, 0.3)
	c3.focus = (0,0,0)
	c3.viewUp = (0, 0, 1)
	c3.viewAngle = 10
	c3.parallelScale = 2660.43
	c3.nearPlane = -3000 #-5320.86
	c3.farPlane = 3000 #5320.86
	c3.imageZoom = 4
	c3.perspective = 1

	c4 = View3DAttributes()
	c4.viewNormal = (0.4, 1.0, 0.1)
	c4.focus = (0,400,0)
	c4.viewUp = (0, 0, 1)
	c4.viewAngle = 10
	c4.parallelScale = 2660.43
	c4.nearPlane = -3000 #-5320.86
	c4.farPlane = 3000 #5320.86
	c4.imageZoom = 10
	c4.perspective = 1

	c5 = View3DAttributes()
	c5.viewNormal = (0.2, 1.0, 0.0)
	c5.focus = (-200,400,0)
	c5.viewUp = (0, 0, 1)
	c5.viewAngle = 10
	c5.parallelScale = 2660.43
	c5.nearPlane = -3000 #-5320.86
	c5.farPlane = 3000 #5320.86
	c5.imageZoom = 12
	c5.perspective = 1

	c6 = View3DAttributes()
	c6.viewNormal = (-0.2,1.0, 0.2)
	c6.focus = (-200,400,0)
	c6.viewUp = (0, 0, 1)
	c6.viewAngle = 10
	c6.parallelScale = 2660.43
	c6.nearPlane = -3000 #-5320.86
	c6.farPlane = 3000 #5320.86
	c6.imageZoom = 40
	c6.perspective = 1

	c7 = View3DAttributes()
	c7.viewNormal = (-0.5,0.5, 0.2)
	c7.focus = (-400,400,0)
	c7.viewUp = (0, 0, 1)
	c7.viewAngle = 10
	c7.parallelScale = 2660.43
	c7.nearPlane = -3000 #-5320.86
	c7.farPlane = 3000 #5320.86
	c7.imageZoom = 20
	c7.perspective = 1

	c8 = View3DAttributes()
	c8.viewNormal = (-0.8,0.2, 0.1)
	c8.focus = (-400,200,0)
	c8.viewUp = (0, 0, 1)
	c8.viewAngle = 10
	c8.parallelScale = 2660.43
	c8.nearPlane = -3000 #-5320.86
	c8.farPlane = 3000 #5320.86
	c8.imageZoom = 12
	c8.perspective = 1

	c9 = View3DAttributes()
	c9.viewNormal = (-0.15,-0.3, 0.2)
	c9.focus = (0,0,0)
	c9.viewUp = (0, 0, 1)
	c9.viewAngle = 10
	c9.parallelScale = 2660.43
	c9.nearPlane = -3000 #-5320.86
	c9.farPlane = 3000 #5320.86
	c9.imageZoom = 4
	c9.perspective = 1

	c10 = View3DAttributes()
	c10.viewNormal = (-0.,-0, 1)
	c10.focus = (0,0,0)
	c10.viewUp = (0, 1, 0)
	c10.viewAngle = 10
	c10.parallelScale = 2660.43
	c10.nearPlane = -3000 #-5320.86
	c10.farPlane = 3000 #5320.86
	c10.imageZoom = 4
	c10.perspective = 1
 
	# Create a tuple of camera values and x values. The x values are weights
	# that help to determine where in [0,1] the control points occur.
	cpts = (c0,c1)
	x=[]
	for i in range(2):
		x = x + [float(i) / float(1.)]
	#first animation
	nsteps=100
	for i in range(nsteps):
		t = float(i) / float(nsteps - 1)
		c = EvalCubicSpline(t, x, cpts)
		c.nearPlane = -3000
		c.farPlane = 3000
		SetView3D(c)
		DrawPlots()
		SaveWindow()

	# Create a tuple of camera values and x values. The x values are weights
	# that help to determine where in [0,1] the control points occur.
	cpts = (c1,c2,c3,c4,c5,c6,c7,c8,c9)
	x=[]
	for i in range(9):
		x = x + [float(i) / float(8.)]
	#second animation
	nsteps=500
	for i in range(nsteps):
		t = float(i) / float(nsteps - 1)
		c = EvalCubicSpline(t, x, cpts)
		c.nearPlane = -3000
		c.farPlane = 3000
		SetView3D(c)
		if i>175:																																																				  
				MeshAtts = MeshAttributes()
				MeshAtts.opacity = (0.8-0.05)/25*(i-175)+0.05
				SetPlotOptions(MeshAtts)
		if i>200:
				MeshAtts = MeshAttributes()
				MeshAtts.opacity = 0.8
				SetPlotOptions(MeshAtts)

		if i>400:
				MeshAtts = MeshAttributes()
				MeshAtts.opacity = -(0.8-0.05)/25*(i-400)+0.8
				SetPlotOptions(MeshAtts)
		if i>425:
				MeshAtts = MeshAttributes()
				MeshAtts.opacity = 0.05
				SetPlotOptions(MeshAtts)

		DrawPlots()
		SaveWindow()
	
def rendering():

	window_options()

	# Flying around the initial data and showing rho and the mesh
	hdf5file = "outMatterU1LFp_000000.3d.hdf5"
	OpenDatabase("../old_data/" + hdf5file)
	DeleteAllPlots()
	annot_ini()
	rho()
	mesh()
	fly_ini()
	DeleteAllPlots()
	
	# Evolve reading next hdf5 files
	for i in range(begin,end,1):
		hdf5file = "outMatterU1LFp_000%03i.3d.hdf5"%i
		print("Analysing file " + hdf5file)
		OpenDatabase("../old_data/" + hdf5file)

		annot_evo()
		rho()
		bh()
		if i>140:
			time = 150 + (i-140)*8
			if i>165:
				gws(time)

			if i<191:
				meshgw(time,i)
				meshgwflat(time,i)
			if i>190:
				meshgw(time,191)
				meshgwflat(time,191)

		fly_evo(i, hdf5file)

		DrawPlots()
		SaveWindow()
		DeleteAllPlots()
	
if __visit_script_file__ == __visit_source_file__:
	rendering()
	os.remove("./visitlog.py")
