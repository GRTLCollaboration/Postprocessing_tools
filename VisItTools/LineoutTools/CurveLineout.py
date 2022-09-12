# This script gets the data of a lineout for different components

begin_file = 0
end_file = 20
file_step = 20
plt_prefix = "name_of_your_hdf5_files"
path_to_hdf5_files = "path_to_hdf5s"

# Choose components
components = ["phi", "lapse", "chi"]

# specify a begin and end point for the lineout, e.g. two (x,y) points if you sliced normal to z
[x_min, y_min, z_min] = [0, 0, 0]
[x_max, y_max, z_max] = [3000, 0, 0]

def rendering():
    
        def lineout(name):

            AddPlot("Curve", "operators/Lineout/" + name, 1, 1)
            LineoutAtts = LineoutAttributes()
            LineoutAtts.point1 = (x_min, y_min, z_min)
            LineoutAtts.point2 = (x_max, y_max, z_max)
            SetOperatorOptions(LineoutAtts, 1)    

        def window_options(name):

            SaveWindowAtts = SaveWindowAttributes()
            SaveWindowAtts.outputToCurrentDirectory = 1
            SaveWindowAtts.fileName = name + "_"
            SaveWindowAtts.family = 1
            SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
            SetSaveWindowAttributes(SaveWindowAtts)


        # Evolve reading next hdf5 files
        for i in range(begin_file ,end_file ,file_step):
            hdf5filename = plt_prefix + "%06i.3d.hdf5"%i
            print("Analysing file " + hdf5filename)
            OpenDatabase(path_to_hdf5_files + hdf5filename)

            for name in components:

                window_options(name)
                lineout(name)
                DrawPlots()
                SaveWindow()
                DeleteAllPlots()
                print("Component: " + name)

            CloseDatabase(path_to_hdf5_files + hdf5filename) 

    
if __visit_script_file__ == __visit_source_file__:
    rendering()
    os.remove("./visitlog.py")
