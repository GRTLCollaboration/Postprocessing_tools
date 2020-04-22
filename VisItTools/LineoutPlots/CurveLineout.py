# This script gets the data of a lineout for different components

begin_file = 0
end_file = 20
file_step = 20
plt_prefix = "name_of_your_hdf5_files"
path_to_hdf5_files = "path_to_hdf5s"

# Choose components
components = ["phi", "lapse", "chi"]

# specify a begin and end point for the lineout, e.g. two (x,y) points if you sliced normal to z
[u_min, v_min] = [0, 0]
[u_max, v_max] = [3000, 0]

def rendering():
    for name in components:

        def lineout():

            AddPlot("Curve", "operators/Lineout/" + name, 1, 1)
            LineoutAtts = LineoutAttributes()
            LineoutAtts.point1 = (u_min, v_min, 0)
            LineoutAtts.point2 = (u_max, v_max, 0)
            SetOperatorOptions(LineoutAtts, 1)    

        def window_options():

            SaveWindowAtts = SaveWindowAttributes()
            SaveWindowAtts.outputToCurrentDirectory = 1
            SaveWindowAtts.fileName = name + "_"
            SaveWindowAtts.family = 1
            SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
            SetSaveWindowAttributes(SaveWindowAtts)


        print("Component: " + name)

        # Evolve reading next hdf5 files
        for i in range(begin_file ,end_file ,file_step):
            hdf5filename = plt_prefix + "%06i.3d.hdf5"%i
            print("Analysing file " + hdf5filename)
            OpenDatabase(path_to_hdf5_files + hdf5filename)

            window_options()
            lineout()
            DrawPlots()
            SaveWindow()
            DeleteAllPlots()

            CloseDatabase(path_to_hdf5_files + hdf5filename) 

    
if __visit_script_file__ == __visit_source_file__:
    rendering()
    os.remove("./visitlog.py")
