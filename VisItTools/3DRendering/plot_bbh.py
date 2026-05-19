# Code used for "Scalar fields around black hole binaries in LVK" movie https://www.youtube.com/watch?v=9yTvjDl9L9A
import os
import glob
import numpy as np

# file details
begin_file = 350
end_file = 400
file_step = 1
plt_prefix = "BinaryBHPlot_"
eta = 1
factor=1
G = "1"
path_to_hdf5_files = "../hdf5/"
# path_to_hdf5_files = "../final_runs/1e-3/2/hdf5/"

# plt_prefix2 = "BoostedBHp_"
# path_to_hdf5_files2 = "/cosma7/data/dp260/dc-aurr2/strings_GRDzhadzha/FlatEvolution/LocalStrings/hdf5/"
output_directory = "frames/"
plot_variable = "chi"
outputfile_prefix = plot_variable

G = float(G)
val = 1

percent = 75

# Create output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
    print(f"Created output directory: {output_directory}")

# Set output PNG format and check the last written PNG file
outputfile_pattern = os.path.join(output_directory, f"{outputfile_prefix}_*.png")
existing_files = glob.glob(outputfile_pattern)

# # If there are existing PNGs, find the last one
# if existing_files:
#     # Sort files to find the latest by comparing the number in the filename
#     existing_files.sort()
#     last_png = existing_files[-1]
    
#     # Extract the number from the filename (assumes format like "phi_000001.png")
#     last_file_number = int(last_png.split('_')[-1].split('.')[0])
    
#     # Restart from the next file
#     begin_file = last_file_number + file_step
# else:
#     # No PNGs found, start from the beginning
#     begin_file = 0


DefineScalarExpression("phi", "sqrt(0.5)*sqrt(((phi_Re)^2+(phi_Im)^2))")
DefineVectorExpression("Coords", "coord(Mesh)")
DefineScalarExpression("X", "Coords[0]-128")
DefineScalarExpression("Y", "Coords[1]-128")
DefineScalarExpression("Z", "Coords[2]-128")
DefineScalarExpression("R", "sqrt(X*X+Y*Y+Z*Z)")
DefineScalarExpression("rPsi4", "R*Weyl4_Re")

DefineScalarExpression("gws", "Weyl4_Re * Weyl4_Re + Weyl4_Im * Weyl4_Im")

    

def isosurface():
        # add pseudocolour plot
    AddPlot("Pseudocolor", (plot_variable), 1, 0)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    # limitsMode =  OriginalData, CurrentPlot
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData
    #PseudocolorAtts.minFlag = set_min_max
    #PseudocolorAtts.min = min_value
    #PseudocolorAtts.maxFlag = set_min_max
    #PseudocolorAtts.max = max_value
    PseudocolorAtts.colorTableName = "white"
    PseudocolorAtts.invertColorTable = 0
    # opacityType = ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque 
    PseudocolorAtts.smoothingLevel = 0
    # centering = Natural, Nodal, Zonal
    PseudocolorAtts.legendFlag = 0
    PseudocolorAtts.centering = PseudocolorAtts.Nodal
    SetPlotOptions(PseudocolorAtts)

    AddOperator("Reflect")
    ReflectAtts = ReflectAttributes()
    ReflectAtts.octant = ReflectAtts.NXNYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
    ReflectAtts.useXBoundary = 1
    ReflectAtts.specifiedX = 0
    ReflectAtts.useYBoundary = 1
    ReflectAtts.specifiedY = 0
    ReflectAtts.useZBoundary = 1
    ReflectAtts.specifiedZ = 0
    ReflectAtts.reflections = (1, 0, 0, 0, 1, 0, 0, 0)
    ReflectAtts.planePoint = (0, 0, 0)
    ReflectAtts.planeNormal = (0, 0, 1)
    ReflectAtts.reflectType = ReflectAtts.Axis  # Plane, Axis
    SetOperatorOptions(ReflectAtts, 0, 0)
    
    AddOperator("Isosurface", 0)
    IsosurfaceAtts = IsosurfaceAttributes()
    IsosurfaceAtts.contourNLevels = 10
    IsosurfaceAtts.contourValue = (0.4)
    IsosurfaceAtts.contourPercent = (percent)
    IsosurfaceAtts.contourMethod = IsosurfaceAtts.Value  # Level, Value, Percent
    IsosurfaceAtts.minFlag = 1
    IsosurfaceAtts.min = 0
    IsosurfaceAtts.maxFlag = 1
    IsosurfaceAtts.max = 0.4
    IsosurfaceAtts.scaling = IsosurfaceAtts.Linear  # Linear, Log
    IsosurfaceAtts.variable = "default"
    SetOperatorOptions(IsosurfaceAtts, 0, 1)

    AddOperator("Transform")
    TransformAtts = TransformAttributes()
    TransformAtts.doTranslate = 1
    TransformAtts.translateZ = 128
    SetOperatorOptions(TransformAtts)


ymax = 2e-2

def slice3(x_1, y_1, x_2, y_2):
    AddPlot("Pseudocolor", "rPsi4", 1, 0)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    # limitsMode =  OriginalData, CurrentPlot
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = -ymax
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = ymax
    PseudocolorAtts.colorTableName = "gws"
    PseudocolorAtts.invertColorTable = 0
    # opacityType = ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque 
    PseudocolorAtts.smoothingLevel = 0
    # centering = Natural, Nodal, Zonal
    PseudocolorAtts.legendFlag = 0
    PseudocolorAtts.lightingFlag = 0
    PseudocolorAtts.centering = PseudocolorAtts.Nodal
    SetPlotOptions(PseudocolorAtts)
    
    AddOperator("Clip")
    ClipAtts = ClipAttributes()
    ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
    ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
    ClipAtts.center = (256, 256, 0)
    ClipAtts.radius = 128
    ClipAtts.sphereInverse = 1
    ClipAtts.crinkleClip = 0
    SetOperatorOptions(ClipAtts, 0, 0)

    # AddOperator("Clip")
    # SliceAtts = SliceAttributes()
    # ClipAtts = ClipAttributes()
    # ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
    # ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
    # ClipAtts.center = (x_2, y_2, 0)
    # ClipAtts.radius = 5.5
    # ClipAtts.sphereInverse = 0
    # ClipAtts.crinkleClip = 0
    # SetOperatorOptions(ClipAtts)
    
    AddOperator("Clip")
    ClipAtts = ClipAttributes()
    ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
    ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
    ClipAtts.center = (256, 256, 0)
    ClipAtts.radius = max(8, 1.8*np.sqrt((x_1-256)**2 + (y_1-256)**2))
    ClipAtts.sphereInverse = 0
    ClipAtts.crinkleClip = 0
    SetOperatorOptions(ClipAtts, 1, 0)
    
    AddOperator("Slice")
    SliceAtts = SliceAttributes()
    SliceAtts.originType = SliceAtts.Percent
    SliceAtts.originPercent = 0
    SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.project2d = 1
    SetOperatorOptions(SliceAtts)
    
    AddOperator("Elevate")
    ElevateAtts = ElevateAttributes()
    ElevateAtts.useXYLimits = ElevateAtts.Always  # Never, Auto, Always
    ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
    ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
    ElevateAtts.skewFactor = 1
    ElevateAtts.minFlag = 1
    ElevateAtts.min = -20
    ElevateAtts.maxFlag = 1
    ElevateAtts.max = 20
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "default"
    SetOperatorOptions(ElevateAtts)


    AddOperator("Transform")
    TransformAtts = TransformAttributes()
    TransformAtts.doTranslate = 1
    TransformAtts.translateZ = -15
    SetOperatorOptions(TransformAtts)

 
def volume():
    AddPlot("Volume", "rho", 1, 0)
    VolumeAtts = VolumeAttributes()
    VolumeAtts.OSPRayEnabledFlag = 0
    VolumeAtts.OSPRayRenderType = VolumeAtts.SciVis  # SciVis, PathTracer
    VolumeAtts.OSPRayShadowsEnabledFlag = 0
    VolumeAtts.OSPRayUseGridAcceleratorFlag = 0
    VolumeAtts.OSPRayPreIntegrationFlag = 0
    VolumeAtts.OSPRaySingleShadeFlag = 0
    VolumeAtts.OSPRayOneSidedLightingFlag = 0
    VolumeAtts.OSPRayAOTransparencyEnabledFlag = 0
    VolumeAtts.OSPRaySPP = 50
    VolumeAtts.OSPRayAOSamples = 0
    VolumeAtts.OSPRayAODistance = 100000
    VolumeAtts.OSPRayMinContribution = 0.001
    VolumeAtts.OSPRayMaxContribution = 2
    VolumeAtts.legendFlag = 1
    VolumeAtts.lightingFlag = 0
    VolumeAtts.colorControlPoints.SetNumControlPoints(9)
    VolumeAtts.colorControlPoints.GetControlPoints(0).colors = (247, 252, 253, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(0).position = 0
    VolumeAtts.colorControlPoints.GetControlPoints(1).colors = (224, 236, 244, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(1).position = 0.125
    VolumeAtts.colorControlPoints.GetControlPoints(2).colors = (191, 211, 230, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(2).position = 0.25
    VolumeAtts.colorControlPoints.GetControlPoints(3).colors = (158, 188, 218, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(3).position = 0.375
    VolumeAtts.colorControlPoints.GetControlPoints(4).colors = (140, 150, 198, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(4).position = 0.5
    VolumeAtts.colorControlPoints.GetControlPoints(5).colors = (140, 107, 177, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(5).position = 0.625
    VolumeAtts.colorControlPoints.GetControlPoints(6).colors = (136, 65, 157, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(6).position = 0.75
    VolumeAtts.colorControlPoints.GetControlPoints(7).colors = (129, 15, 124, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(7).position = 0.875
    VolumeAtts.colorControlPoints.GetControlPoints(8).colors = (77, 0, 75, 255)
    VolumeAtts.colorControlPoints.GetControlPoints(8).position = 1
    VolumeAtts.colorControlPoints.smoothing = VolumeAtts.colorControlPoints.Linear  # NONE, Linear, CubicSpline
    VolumeAtts.colorControlPoints.equalSpacingFlag = 0
    VolumeAtts.colorControlPoints.discreteFlag = 0
    VolumeAtts.colorControlPoints.tagNames = ()
    VolumeAtts.opacityAttenuation = 1
    VolumeAtts.opacityMode = VolumeAtts.GaussianMode  # FreeformMode, GaussianMode, ColorTableMode
    VolumeAtts.opacityControlPoints.SetNumControlPoints(1)
    VolumeAtts.opacityControlPoints.GetControlPoints(0).x = 0.810448
    VolumeAtts.opacityControlPoints.GetControlPoints(0).height = 0.108108
    VolumeAtts.opacityControlPoints.GetControlPoints(0).width = 0.0940298
    VolumeAtts.opacityControlPoints.GetControlPoints(0).xBias = 0.0552239
    VolumeAtts.opacityControlPoints.GetControlPoints(0).yBias = 1.74286
    VolumeAtts.resampleType = VolumeAtts.SingleDomain  # NoResampling, OnlyIfRequired, SingleDomain, ParallelRedistribute, ParallelPerRank
    VolumeAtts.resampleTarget = 50000000
    VolumeAtts.resampleCentering = VolumeAtts.NativeCentering  # NativeCentering, NodalCentering, ZonalCentering
    VolumeAtts.opacityVariable = "default"
    VolumeAtts.freeformOpacity = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 26, 44, 57, 65, 70, 80, 88, 93, 104, 109, 114, 119, 124, 127, 127, 130, 132, 132, 132, 135, 135, 135, 135, 135, 135, 124, 119, 109, 104, 91, 78, 72, 62, 52, 41, 33, 26, 13, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    VolumeAtts.useColorVarMin = 0
    VolumeAtts.colorVarMin = -3
    VolumeAtts.useColorVarMax = 1
    VolumeAtts.colorVarMax = 2
    VolumeAtts.useOpacityVarMin = 0
    VolumeAtts.opacityVarMin = 0.0001
    VolumeAtts.useOpacityVarMax = 0
    VolumeAtts.opacityVarMax = 600
    VolumeAtts.smoothData = 0
    VolumeAtts.samplesPerRay = 500
    VolumeAtts.rendererType = VolumeAtts.Serial  # Serial, Parallel, Composite, Integration, SLIVR
    VolumeAtts.gradientType = VolumeAtts.SobelOperator  # CenteredDifferences, SobelOperator
    VolumeAtts.scaling = VolumeAtts.Log  # Linear, Log, Skew
    VolumeAtts.skewFactor = 1
    VolumeAtts.limitsMode = VolumeAtts.OriginalData  # OriginalData, CurrentPlot
    VolumeAtts.sampling = VolumeAtts.Rasterization  # KernelBased, Rasterization, Trilinear
    VolumeAtts.rendererSamples = 3
    VolumeAtts.lowGradientLightingReduction = VolumeAtts.Lower  # Off, Lowest, Lower, Low, Medium, High, Higher, Highest
    VolumeAtts.lowGradientLightingClampFlag = 0
    VolumeAtts.lowGradientLightingClampValue = 1
    VolumeAtts.materialProperties = (0.4, 0.75, 0, 15)
    SetPlotOptions(VolumeAtts)

    AddOperator("Reflect")
    ReflectAtts = ReflectAttributes()
    ReflectAtts.octant = ReflectAtts.NXNYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
    ReflectAtts.useXBoundary = 1
    ReflectAtts.specifiedX = 0
    ReflectAtts.useYBoundary = 1
    ReflectAtts.specifiedY = 0
    ReflectAtts.useZBoundary = 1
    ReflectAtts.specifiedZ = 0
    ReflectAtts.reflections = (1, 0, 0, 0, 1, 0, 0, 0)
    ReflectAtts.planePoint = (0, 0, 0)
    ReflectAtts.planeNormal = (0, 0, 1)
    ReflectAtts.reflectType = ReflectAtts.Axis  # Plane, Axis
    SetOperatorOptions(ReflectAtts, 0, 0)
    
    AddOperator("Transform")
    TransformAtts = TransformAttributes()
    TransformAtts.doTranslate = 1
    TransformAtts.translateZ = 128
    SetOperatorOptions(TransformAtts)
    
def mesh():
    AddPlot("Mesh", "phi", 1, 1)
    # SetActivePlots(0)
    # SetActivePlots(0)
    MeshAtts = MeshAttributes()
    MeshAtts.showInternal = 0
    MeshAtts.meshColorSource = MeshAtts.Foreground  # Foreground, MeshCustom
    MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
    MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
    # MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
    MeshAtts.opacity = 0.2
    SetPlotOptions(MeshAtts)

    AddOperator("Slice")
    SliceAtts = SliceAttributes()
    SliceAtts.originType = SliceAtts.Percent
    SliceAtts.originPercent = 50
    SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.project2d = 0
    SetOperatorOptions(SliceAtts)


# max and min values for colourbar
set_min_max = 1 # 1 for true, 0 for false
min_value = -0.01
max_value = 0.01
# slice origin and normal direction
origin_point_x = 256
origin_point_y = 256
origin_point_z = 20
normal_in_x = 0
normal_in_y = 0
normal_in_z = 1
# max and min coords in the sliced plane (e.g. x and y if normal to z)
# NB relative to origin as defined above
min_u = -10.0
max_u = 10.0
min_v = -10.0
max_v = 10.0

# --------------------------------------------------------------------------
# From here you only need to amend if you want to do something non standard
# --------------------------------------------------------------------------


def view2():
    ResetView()
    # Begin spontaneous state
    View3DAtts = View3DAttributes()
    View3DAtts.viewNormal = (0.5, -0.8, 0.214584) #(0.5, -0.5, 0.214584)
    View3DAtts.focus = (256, 256, -512)
    View3DAtts.viewUp = (0, 0, 1)
    View3DAtts.viewAngle = 30
    View3DAtts.parallelScale = 384
    View3DAtts.nearPlane = -768
    View3DAtts.farPlane = 768
    View3DAtts.imagePan = (0, 0)
    View3DAtts.imageZoom = 8
    View3DAtts.perspective = 1
    View3DAtts.eyeAngle = 2
    View3DAtts.centerOfRotationSet = 0
    View3DAtts.centerOfRotation = (256, 256, 0)
    View3DAtts.axis3DScaleFlag = 0
    View3DAtts.axis3DScales = (1, 1, 1)
    View3DAtts.shear = (0, 0, 1)
    View3DAtts.windowValid = 1
    SetView3D(View3DAtts)
    # End spontaneous state


x_1 = np.loadtxt("../data/punctures.dat", usecols = 1)
y_1 = np.loadtxt("../data/punctures.dat", usecols = 2)
x_2 = np.loadtxt("../data/punctures.dat", usecols = 4)
y_2 = np.loadtxt("../data/punctures.dat", usecols = 5)

# OpenComputeEngine("localhost", ("-l", "mpirun", "-np", "28", "-nn", "1","-p", "spr"))
# Loop through the files, and plot each slice
for i in range ( begin_file , end_file + file_step , file_step) :
    hdf5filename = plt_prefix + "%06i.3d.hdf5"%i
    full_file_path = path_to_hdf5_files + hdf5filename
    # hdf5filename2 = plt_prefix2 + "%06i.3d.hdf5"%i
    # full_file_path2 = path_to_hdf5_files2 + hdf5filename2

    # Check if the file exists
    if not os.path.exists(full_file_path):
        print(f"File {full_file_path} not found. Stopping.")
        break
    print("Plotting file " + hdf5filename)

    OpenDatabase(full_file_path)

    isosurface()
    volume()
    slice3(x_1[i], y_1[i], x_2[i], y_2[i])
    view2()


    # # save settings
    # saveAtts = SaveWindowAttributes()
    # saveAtts.format = saveAtts.PNG
    # saveAtts.resConstraint = saveAtts.NoConstraint
    # saveAtts.width = 1280
    # saveAtts.height = 620
    # SetSaveWindowAttributes(saveAtts)




    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.timeInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.axesArray.visible = 0
    # AnnotationAtts.backgroundColor = (0, 0, 0, 255)
    # AnnotationAtts.foregroundColor = (255, 255, 255, 255)
    AnnotationAtts.axes3D.bboxFlag = 0
    AnnotationAtts.axes3D.triadFlag = 0
    SetAnnotationAttributes(AnnotationAtts)

    # plot all levels
    # silr = SILRestriction()
    # silr.TurnOnAll()
    # SetPlotSILRestriction(silr ,1) #Do i need this? idk


    InvertBackgroundColor()


    # engage!
    DrawPlots()

    # Zoom in
    #View2DAtts = View2DAttributes()
    #View2DAtts.windowCoords = (min_u, max_u, min_v, max_v)
    #View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    #View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
    #View2DAtts.fullFrameAutoThreshold = 100
    #View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
    #View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
    #View2DAtts.windowValid = 1
    #SetView2D(View2DAtts)

    # Create a Matplotlib figure
    # plt.figure()
    # plt.plot(x, y)
    # plt.title("Matplotlib Plot Example")
    # plt.xlabel("X-axis")
    # plt.ylabel("Y-axis")

    
    # Save Me
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputDirectory = output_directory
    SaveWindowAtts.fileName = outputfile_prefix
    SaveWindowAtts.family = 0
    SaveWindowAtts.fileName = f"{output_directory}{outputfile_prefix}_{i:06d}"  # Ensures %06i format
        # format = BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, 
        # RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.format = SaveWindowAtts.PNG
    SaveWindowAtts.width = 8192
    SaveWindowAtts.height = 5120
    SaveWindowAtts.quality = 100
    # resConstraint = NoConstraint, EqualWidthHeight, ScreenProportions
    # SaveWindowAtts.resConstraint = SaveWindowAtts.EqualWidthHeight 
    SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint 
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()

    # clear plots ready for next timestep in loop
    DeleteAllPlots()
    CloseDatabase(full_file_path)


# loop ends
print ("I've finished!")
