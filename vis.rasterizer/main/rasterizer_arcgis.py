'''
@author: Karl J. Smith
'''
#TODO: Find a style guide for file headers. Then use style guide for file headers.

# Import libraries.
import arcpy, os, time
from arcpy import env, Raster, GetCellValue_management #@UnusedImport
from arcpy.sa import * #@UnusedWildImport

# Import functions generally.
from functions import clock

# Import functions individually.
from functions.projective_array import projective_array
from functions.surface_intersect import intersect_3d_line_with_surface
from functions.surface_intersect import check_pickup_points

#===============================================================================
# SET PARAMETERS
# 
# WARNING: THERE ARE OTHER PARAMETERS IN THE MAIN BODY OF THE SCRIPT.
#===============================================================================

#TODO: Add option to let user enter parameters via console. Would need a function? Otherwise probably tedious.
# Fudge it for now.
input_to_close = False

# Set overwrite status.
overwrite = True

# Set path for input observe points.
in_obs_pts = "C:\\GIS\\GAO_2023\\Temp\\Adam_demo\\TD_01_clip.shp"

#TODO: Iterate over targets as well. Will require and extra loop and reconfiguration of output D:
in_target_grid = "C:\\GIS\\GAO_2023\\Temp\\Adam_demo\\u189grd.shp"
in_target_mask = "C:\\GIS\\GAO_2023\\Temp\\Adam_demo\\u189msk.shp"
    
# set path elevation surface.
elev_path = "C:\\GIS\\GAO_2023\\Temp\\Adam_demo\\TD_01_lidar.tif"

# Set folder for processing datasets.
timestep_processing_folder = "C:\\GIS\\GAO_2023\\Temp\\Adam_demo\\processing"

# Set array parameters.
#TODO: BRING THE REST OF THE PARAMETERS IN HERE. CURRENTLY BUFFERS AT 5% BY DEFAULT.

# Set the resolution of the array (in minutes).
ray_array_res = 10

# Set array maximum range (in map units).
ray_range_max = 20000

# Set array minimum range (in map units).
ray_range_min = 1

#===============================================================================
# MAIN LOOP
#===============================================================================

# Create dictionary to hold timing information.
clock_dict = {}

# Get the script time start (in seconds).
script_start_time = time.time()

#------------------------------------------------------------------------------ 

# Print starting message.
print("Started script 'rasterizer_arcgis' at {}.".format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(script_start_time))))
print("")

#------------------------------------------------------------------------------ 

#TODO: Check CRS for each input and compare. Via function? Might add validation file.

#------------------------------------------------------------------------------ 

# Convert input mask shapefile to geometry object.
target_mask_geometry = arcpy.CopyFeatures_management(in_target_mask, arcpy.Geometry())[0]

#------------------------------------------------------------------------------ 

# Create raster object for intersect.
SRTM_raster = Raster(elev_path)
# Generate raster mosaic for intersect.
mosaic_cellsize_x = arcpy.GetRasterProperties_management(elev_path, "CELLSIZEX").getOutput(0)
mosaic_cellsize_y = arcpy.GetRasterProperties_management(elev_path, "CELLSIZEY").getOutput(0)
mosaic_max_cellsize = float(max((mosaic_cellsize_x, mosaic_cellsize_y))) * 1.1 # Multiply cell size by 1.1 to avoid edge effects.

#------------------------------------------------------------------------------ 

# Add field for to hold visibility values to input points (if it doesn't exist already).
in_pts_desc = arcpy.Describe(in_obs_pts)
if "VIS3D" not in [f.name for f in in_pts_desc.fields]:
    arcpy.management.AddField(in_obs_pts, "VIS3D","FLOAT")

#------------------------------------------------------------------------------ 

pr_count = 0
time_elapsed_list = []

cursor_fields = ['SHAPE@XYZ','VIS3D']

if overwrite:
    where_clause = None
else:
    where_clause = '"VIS3D" = -1'

with arcpy.da.UpdateCursor(in_obs_pts, cursor_fields, where_clause) as count_cursor: #@UndefinedVariable
    
    obs_pt_count = 0
    
    for row in count_cursor:
        
        if overwrite == True or int(row[1]) == -1:
        
            obs_pt_count += 1
        
del count_cursor

print("Processing {} observation points...".format(obs_pt_count))

#------------------------------------------------------------------------------ 
    
with arcpy.da.UpdateCursor(in_obs_pts, cursor_fields, where_clause) as obs_cursor: #@UndefinedVariable
    
    for row in obs_cursor:
        
        row_start_time = time.time()
        
        obs_x, obs_y, obs_z = row[0]
        
        #------------------------------------------------------------------------------ 

        # Generate projective array using adapted 'Directed_Projective_Array' tool.
        
        # Get list of coordinates for target from shapefile.
        in_target_grid_list = [i[0] for i in arcpy.da.SearchCursor(in_target_grid, ["SHAPE@XYZ"])] #@UndefinedVariableFromImport
        
        # Set filepath for array.
        ray_output = "{}\\array_{}_{}.shp".format(timestep_processing_folder, str(int(obs_x)), str(int(obs_y)))
        
        # If 'ray_range_min' is negative, set to 0.
        if ray_range_min < 0:
            
            ray_range_min = 0
        
        # Create the array.
        #TODO: Configure function to accept optional inputs. Also to handle errors if it fails.
        pa_result, pa_time = projective_array((obs_x, obs_y, obs_z),
                                              in_target_grid_list, # Grid int filepath.
                                              ray_output,
                                              "",
                                              ray_array_res,
                                              ray_range_max,
                                              ray_range_min, 
                                              False,
                                              "INCREMENT",
                                              5,
                                              5,
                                              5,
                                              5,
                                              mosaic_max_cellsize,
                                              "")
        
        # Clock function.
        clock.update_clock_dict(clock_dict, "Directed_Projective_Array", pa_time)    
        
        #------------------------------------------------------------------------------ 
        
        # Set variables for 'Intersect_3D_Line_With_Surface' tool.
        
        # Set folder path for rays ('pickup lines').
        pickup_lines_file = "{}\\{}.shp".format(timestep_processing_folder, "junk_lines") # Line features are generated by the ArcGIS tool but not used by the script.
        
        # Set folder path and filename for points ('pickup points').
        pickup_points_filename = "pp_{}_{}.shp".format(str(int(obs_x)), str(int(obs_y)))
        pickup_points_file = "{}\\{}".format(timestep_processing_folder, pickup_points_filename)
        
        # Intersect array with surface.
        int3d_result, int3d_time = intersect_3d_line_with_surface(ray_output, SRTM_raster, pickup_lines_file, pickup_points_file)
        
        # If intersection function was successful, 'int3d_tool_result' will be 'True'.
        if int3d_result:
            
            # Create variables for selected pickup point file.
            selected_pickup_points_filename = "upp_{}_{}.shp".format(str(int(obs_x)), str(int(obs_y)))
            selected_pickup_points_file = "{}\\{}".format(timestep_processing_folder, selected_pickup_points_filename)

            # Select pickup points and return count of visible points.
            #TODO: Add error handling for 'pickup_points_result'.
            pickup_points_result, pickup_points_count, pickup_points_time = check_pickup_points(ray_output, pickup_points_file, target_mask_geometry, selected_pickup_points_file)
            
        else:
            
            # Update pr_count.
            pr_count += 1
            
            pickup_points_count = 0
   
        #------------------------------------------------------------------------------ 
        
        # Update table.
        row[1] = pickup_points_count
        
        obs_cursor.updateRow(row)
        
        # Update pr_count.
        pr_count += 1
        
        # Compute timings and print line to console.
        time_elapsed = (time.time() - script_start_time)
        row_time_elapsed = (time.time() - row_start_time)
        time_elapsed_str = time.strftime("%H:%M:%S", time.gmtime(time_elapsed))
        row_time_elapsed_str = time.strftime("%M:%S", time.gmtime(row_time_elapsed))
        
        time_elapsed_list.append(row_time_elapsed)
        time_elapsed_average = float(sum(time_elapsed_list[-50:]))/float(len(time_elapsed_list[-50:])) # Sample the last fifty elements of the list to get a moving average.
        
        time_left = round((time_elapsed_average * (obs_pt_count - pr_count))/3600,2) # Get an average in hours.
        
        print("{:<8} : Point {}/{} : {:>6} Pickup Pts. : Took {}, approx. {} hrs. to go...".format(time_elapsed_str, pr_count, obs_pt_count, pickup_points_count, row_time_elapsed_str, time_left))
        
        # Clean up.
        for filepath in [pickup_lines_file, pickup_points_file, selected_pickup_points_file]:
            if os.path.exists(filepath):
                arcpy.Delete_management(filepath)

#TODO: Clean up and structure output.
print("Script finished!")
print("")

# Print the clock dict, to get some idea of how long everything took to run.
clock.print_clock_dict(clock_dict)

# Make the user press a key to close (if enabled).
if input_to_close:
    print("")
    input("Press any key to close...")
