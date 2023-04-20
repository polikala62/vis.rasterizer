'''
Created on Mar 16, 2023

@author: kjsmi
'''

import arcpy, math, os, time
from arcpy import env, Raster, GetCellValue_management #@UnusedImport
from arcpy.sa import * #@UnusedWildImport

clock_dict = {}

def update_clock_dict(in_dict, func_name, func_start_time):
    
    if func_name not in in_dict.keys():
        in_dict[func_name] = [(time.time() - func_start_time)]
    else:
        in_dict[func_name].append((time.time() - func_start_time))
        
def print_clock_dict(in_dict):
    
    measured_time = sum([sum(in_dict[i]) for i in in_dict.keys()])
    
    def format_secs(secs):
        if type(secs) == int or type(secs) == float:
            if round(secs,2) > 0:
                return str(round(secs,2))
            else:
                return "< 0.01"
        else:
            return str(secs)
    
    fields_list = ["NAME","TIMES_CALLED","MAX_TIME","MIN_TIME","TOTAL_TIME","AVG_TIME","%_MEASURED"]
    print("{:<36} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12}".format(*fields_list))
    
    for func in sorted(in_dict.keys()):
        
        times_called = len(in_dict[func])
        max_time = max(in_dict[func])
        min_time = min(in_dict[func])
        ttl_time = sum(in_dict[func])
        avg_time = ttl_time/times_called
        prcnt_total = (ttl_time/measured_time)*100
        
        row_list = [format_secs(i) for i in [func, times_called, max_time, min_time, ttl_time, avg_time, prcnt_total]]
        
        print("{:<36} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12}".format(*row_list))
        
def print_prcnt_complete(pr_count, pr_total, increment):
    
    inc_list = [i/pr_total for i in range(0,100,increment)]
    
    if pr_count/pr_total in inc_list:
        
        print("Completed {}% of total.".format(str(int(pr_count/pr_total))))
        
    elif pr_count == pr_total:
        
        print("Completed 100% of total.")
        

def search_cell(in_search_pt, in_cell_pt, cell_width): # Defines a square based on a centrepoint and a width, and check for 2D intersection.
    
    cell_x, cell_y, cell_z = in_cell_pt #@UnusedVariable
    search_x, search_y, search_z = in_search_pt #@UnusedVariable
    
    cell_xmin = cell_x - (cell_width * 0.5)
    cell_xmax = cell_x + (cell_width * 0.5)
    cell_ymin = cell_y - (cell_width * 0.5)
    cell_ymax = cell_y + (cell_width * 0.5)
    
    if search_x >= cell_xmin and search_x <= cell_xmax and search_y >= cell_ymin and search_y <= cell_ymax:
        
        return True
    
    else:
        
        return False

#===============================================================================
# COORDINATE DISTANCE TO GEOMETRY.

# This function does not report errors.
#===============================================================================

def Coordinate_Distance_to_Geometry(coordinates, geometry):
    
    # Start benchmarking, if enabled.
    #btime = rec.Benchmark_Start()
                        
    coord_x, coord_y = coordinates
    
    coord_geometry = arcpy.PointGeometry(arcpy.Point(coord_x, coord_y), arcpy.SpatialReference(3035))
    
    distance_list = []
    
    '''
    # Loop through features in input geometry.
    for feature in geometry:
        
        distance_list.append(coord_geometry.distanceTo(feature))
    '''
    distance_list.append(coord_geometry.distanceTo(geometry))
    
    # Finish benchmarking, if enabled.
    #rec.Benchmark_End(btime, bdict,"nav.Coordinate_Distance_to_Geometry")
    
    del coord_geometry
    
    return min(distance_list)

#===========================================================================
# Directed_Projective_Array
#===========================================================================

'''
Note:   This code has been copied from the 'Directed Projected Array' script developed as part of my MST.
        The following code has been modified so that:
            - User interface has been removed.
            - Some dataset validation has been removed / deprecated (because said validation occurs in the main script).
            - Loops have been removed, so that the script does not iterate through observer-target combinations. UIDs are received as variables instead.
            - Inputs have changed:  'in_obs_feature' must be singlepoint shapefile with only one point.
                                    'in_tar_feature' has been replaced by 'tar_coordinate_list'.
'''

def Directed_Projective_Array(in_obs_feature, # Filepath for the shapefile / feature class that contains observer points. Must be Z-enabled single point features.
                              tar_coordinate_list, # List of X/Y/Z coordinates representing features.
                              out_fc, # Filepath for output array.
                              out_table, # Filepath for optional ouput table.
                              ):

    #===============================================================================
    # I. SCRIPT SETUP
    #===============================================================================
    
    # Enable overwrite.
    env.overwriteOutput = 1
    
    dpa_start_time = time.time()
    
    #------------------------------------------------------------------------------ 
    
    # Set 'obs_id' and 'tar_id' to the filenames of the input observer and target files.
    obsID = "OBS"
    tarID = "TAR"
    #obsID = os.path.basename(in_obs_feature).split(".")[0]
    #tarID = os.path.basename(in_tar_feature).split(".")[0]
    
    # Set 'arrayID' to a combination of 'obsID' and 'tarID'.
    arrayID = "0"
    
    #------------------------------------------------------------------------------ 

    #===============================================================================
    # III. CREATE OUTPUT DATASET(S)
    #===============================================================================
    
    # 1. Create output array features.
    
    # Set variables for output folder and filename.
    out_fc_folder = os.path.dirname(out_fc)
    out_fc_filename = os.path.basename(out_fc)
    
    # Create polyline feature class.
    arcpy.CreateFeatureclass_management(out_fc_folder, out_fc_filename, "POLYLINE", "", "DISABLED", "ENABLED", arcpy.SpatialReference(3035))
    
    # Add ID fields identifying unique observer-target combinations with an arbitrary value.
    arcpy.AddField_management(out_fc, "Array_ID", "LONG")
    
    # Add ID field identifying observer point.
    arcpy.AddField_management(out_fc, "From_Obs", "TEXT")
    
    # Add ID field identifying target feature.
    arcpy.AddField_management(out_fc, "To_Feature", "TEXT")
    
    # Add diagnostic fields to radial array.
    for x in ("Azimuth", "Altitude", "Length", "From_X", "From_Y", "From_Z", "To_X", "To_Y", "To_Z"):
        arcpy.AddField_management(out_fc, x, "FLOAT")
    
    # Print message to console.
    #NFX arcpy.AddMessage("{} Output array features created in '{}'".format(ETstamp(startTime), out_fc)) 
    
    #------------------------------------------------------------------------------ 
    
    # 1. Create output table.
    
    if out_table != "":
        # Set variables for output folder and filename.
        out_table_folder = os.path.dirname(out_table)
        out_table_filename = os.path.basename(out_table)
        
        # Create table.
        arcpy.CreateTable_management(out_table_folder, out_table_filename)
        
        # Add ID fields identifying unique observer-target combinations with an arbitrary value.
        arcpy.AddField_management(out_table, "Array_ID", "LONG")
    
        # Add ID field identifying observer point.
        arcpy.AddField_management(out_table, "From_Obs", "TEXT")
    
        # Add ID field identifying target feature.
        arcpy.AddField_management(out_table, "To_Feature", "TEXT")
        
        # Add diagnostic fields to radial array.
        for x in ("MIN_ALT", "MAX_ALT", "MIN_AZI", "MAX_AZI", "MAX_DIST", "COUNT_2TAR", "COUNT_TTL"):
            arcpy.AddField_management(out_table, x, "FLOAT")
            
        # Delete automatically generated 'Field1'
        arcpy.DeleteField_management(out_table, "Field1")
        
        # Print message to console.
        #NFX arcpy.AddMessage("{} Output table created in '{}'.".format(ETstamp(startTime), out_table)) 
    
    #===============================================================================
    # IV. DEFINE FUNCTIONS
    #===============================================================================
    
    # 1. Define a function for finding the 3D distance between an origin and a target point.
    
    # Define function variables 'origin' and 'target', where both are tuples containing x, y, and z coordinates.
    def EDist(origin, target):
        
        # Calculate differences between input x, y, and z coordinates.
        Xdiff, Ydiff, Zdiff = (target[0] - origin[0]), (target[1] - origin[1]), (target[2] - origin[2])
        
        # Calculate the XY distance between input points by finding the hypotenuse of a triangle with sides 'Xdiff' & 'Ydiff'.
        XY_dist = math.hypot(Xdiff, Ydiff)
        
        # Calculate and return the 3D distance by finding the hypotenuse of a triangle with sides 'XY_dist' & 'Zdiff'.
        return math.sqrt((XY_dist * XY_dist) + (Zdiff * Zdiff))
    
    #------------------------------------------------------------------------------ 
    
    # 2. Define a function for calculating the right-handed azimuth angle (in minutes) defined by an origin point, a target point, and an initial ray y = 0.
    
    # Define function variables 'origin' and 'target', where both are tuples containing x and y coordinates.
    def Define_AZI(origin, target):
        
        # Calculate differences between maximum and minimum X and Y values.
        Xdiff, Ydiff = (target[0] - origin[0]), (target[1] - origin[1])
        
        # Return the azimuth by calculating the arc tangent of X and Y distances to target point.
        return math.degrees(math.atan2(Ydiff,Xdiff)) * 60
    
    #------------------------------------------------------------------------------ 
    
    # 3. Define a function for calculating the altitude angle (in minutes) from an up vector defined by (x, y = 0).
    
    def Define_ALT(origin, target):
        
        # Calculate differences between maximum and minimum X, Y, and Z values.
        Xdiff, Ydiff, Zdiff = (target[0] - origin[0]), (target[1] - origin[1]), (target[2] - origin[2])
        
        # Return length of the line as a variable, by calculating the hypotenuse of a triangle with sides 'Xdiff' & 'Ydiff'.
        XY_dist = math.hypot(Xdiff, Ydiff)
        
        # Return altitude by calculating the arc tangent of the line length and elevation.
        return math.degrees(math.atan2(XY_dist,Zdiff)) * 60
    
    #------------------------------------------------------------------------------ 
    
    # 4. Define a function for calculating the angular boundaries for a feature represented by a list of azimuths.
    
    def BoundaryAngles(angle_list):
        
        # Sort input angle list from smallest value to largest value.
        sorted_angle_list = sorted(angle_list)
        
        # Add last entry in the smallest angle list to the angle list.
        adjusted_angle_list = [sorted_angle_list[-1] - 360] + sorted_angle_list
        
        # Create list of values where value is subtracted from succeeding value (diff_list).
        diff_list = [abs(y-x) for x, y in zip(adjusted_angle_list, adjusted_angle_list[1:])]
        
        # Get positions for maximum values in the diff_list.
        positions = [i for i, j in enumerate(diff_list) if j == max(diff_list)]
        
        # Return first value in the positions list as a variable.
        firstposition = positions[0]
    
        # If greatest difference is between first and last elements of the list, return end of list as second value.
        if firstposition == 0:
            return (sorted_angle_list[firstposition],  sorted_angle_list[-1])
        
        # Otherwise, return the value preceding the maximum difference value.
        else:
            return (sorted_angle_list[firstposition], sorted_angle_list[firstposition - 1] + 360)

        
    #===============================================================================
    # VI. GENERATE ARRAYS & COMPLETE SCRIPT
    #===============================================================================
    
    # 1. Iterate observer points.
    
    # Create 'array_count_total', and 'observer count' variables.
    array_count_total = 0
    #obs_count = 0
    
    # If output statistics is enabled, create tuple to hold statistics tuples.
    if out_table != "":
        tuple_stats = []
    '''
    # Create search cursor for observer feature class.
    with arcpy.da.SearchCursor(in_obs_feature, feature_fields) as obscursor: #@UndefinedVariableFromImport
    
        # Iterate through features in observer feature class.
        for OBSrow in obscursor:
            
            # Increment observer counter.
            obs_count += 1
            
            # Create variables for observer point.
            obsX, obsY, obsZ = OBSrow
    '''
    obsX, obsY, obsZ = in_obs_feature
    
    # Delete observer cursor.
    #del obscursor
    
    #------------------------------------------------------------------------------ 
    
    # Define a 'count' variable to record the number of lines created for the array.
    array_count = 0
    
    for TARrow in tar_coordinate_list:
            
        #------------------------------------------------------------------------------ 
        
        # 3. Calculate angular values for target point.
        
        # Return target XYZ values for target point as variables.
        tarX, tarY, tarZ = TARrow
        
        # Create insert cursor for observer feature class.
        insert_line_fields = ["SHAPE@", "Array_ID", "From_Obs", "To_Feature", "Azimuth", "Altitude", "Length", "From_X", "From_Y", "From_Z", "To_X", "To_Y", "To_Z"]
        insert_line_cursor = arcpy.da.InsertCursor(out_fc, insert_line_fields) #@UndefinedVariableFromImport
        
        ray_coords = [[obsX, obsY, obsZ], [tarX, tarY, tarZ]]
        
        maxDistance = math.sqrt((tarX - obsX)**2 + (tarY - obsY)**2 + (tarZ - obsZ)**2)
                
        # Insert target feature.
        insert_line_cursor.insertRow((ray_coords, arrayID, obsID, tarID, 0, 0, maxDistance, obsX, obsY, obsZ, tarX, tarY, tarZ))
        
        # Add to 'array_count'.
        array_count += 1
        
    # Increment 'array_count_total'.
    array_count_total += array_count
        
    # If statistics output is enabled, append tuple containing statistics to 'tuple_stats'.
    if out_table != "":
        tuple_stats.append([arrayID, obsID, tarID, 0, 0, 0, 0, maxDistance, array_count])

    #------------------------------------------------------------------------------ 
    
    # 6. Populate output statistics table (if enabled), and print 'finished' message to console.
    
    # Insert rows in output table for each entry in 'tuple_stats'.
    if out_table != "":
        tableFields = ["Array_ID", "From_Obs", "To_Feature", "MIN_ALT", "MAX_ALT", "MIN_AZI", "MAX_AZI", "MAX_DIST", "COUNT_2TAR", "COUNT_TTL"]
        tablecursor = arcpy.da.InsertCursor(out_table, tableFields) #@UndefinedVariableFromImport
        
        for x in tuple_stats:
            tablecursor.insertRow((x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], array_count_total))
        
    # Clock function.
    update_clock_dict(clock_dict, "Directed_Projective_Array", dpa_start_time)

#===============================================================================
# INTERSECT_3D_LINE WITH SURFACE
#
# This script produces reports.
#===============================================================================

# This script runs the 'Intersect3DLineWithSurface_3d' tool, and returns error messages if the tool fails.
def Intersect_3D_Line_With_Surface(ray_output, SRTM_raster, pl_file, pp_file):
    
    ils_start_time = time.time()
    
    #function = "Intersect_3d_Line_With_Surface"
    
    report = True
                        
    # Start benchmarking, if enabled.
    #btime = rec.Benchmark_Start()

    try:
        # Check out extension.
        arcpy.CheckOutExtension("3D")
        
        # Run tool.
        arcpy.Intersect3DLineWithSurface_3d(ray_output, SRTM_raster, pl_file, pp_file)
        
        # Check in extension.
        arcpy.CheckInExtension("3D")
        
        # Finish benchmarking, if enabled.
        #rec.Benchmark_End(btime, bdict, function) 
        
        # Clock function.
        update_clock_dict(clock_dict, "Intersect_3D_Line_With_Surface", ils_start_time)
        
        # Return report and 'True' value.
        return (report, True)
    
    except:
        
        # Get the number of messages returned by the tool.
        message_count = arcpy.GetMessageCount()
        
        # Create list to hold messages.
        message_list = []
        
        # Add messages to list.
        for i in range(0,message_count):
            message_list.append(arcpy.GetMessage(i))
            print(arcpy.GetMessage(i))
        
        # Finish benchmarking, if enabled.
        #rec.Benchmark_End(btime, bdict, function) 
        
        # File report.
        #report = exc.File_Report(function, "Unknown error.", False, messages=message_list, outputs=[pl_file, pp_file], inputs=[ray_output, SRTM_raster])
        
        # Clock function.
        update_clock_dict(clock_dict, "Intersect_3D_Line_With_Surface", ils_start_time)
        
        # Return report and empty value.
        return (report, None)
    
#===============================================================================
# Update_Pickup_Points
#
# Copied from 'Visibility Analysis' --> update pickup points.
# The script has been substantially rewritten, in order to improve performance.
#
#===============================================================================

def Check_Pickup_Points(in_array, in_points, out_fc, output_status, in_grid_resolution):
    
    cpp_start_time = time.time()
    
    # Allow overwrite.
    env.overwriteOutput = 1
    
    #------------------------------------------------------------------------------ 
    
    # Create output point features, if enabled.
    
    if output_status == True:
    
        # Set variables for output folder and filename.
        out_fc_folder = os.path.dirname(out_fc)
        out_fc_filename = os.path.basename(out_fc)
        
        # Create polyline feature class.
        arcpy.CreateFeatureclass_management(out_fc_folder, out_fc_filename, "POINT", "", "DISABLED", "ENABLED", in_array)
        
        # Add ID fields to output points.
        arcpy.AddField_management(out_fc, "OID_Line", "TEXT")
        arcpy.AddField_management(out_fc, "Array_ID", "LONG")
        arcpy.AddField_management(out_fc, "From_Obs", "TEXT")
        arcpy.AddField_management(out_fc, "To_Feature", "TEXT")
        
        # Add diagnostic fields to output points.
        for x in ("Azimuth", "Altitude", "Dist_Along", "PICKUP"):
            arcpy.AddField_management(out_fc, x, "FLOAT")
            
    #------------------------------------------------------------------------------ 
    # if output_status == True: #Not sure why this check is here.
    # Create 'arrayDict'.
    array_fields = ["OID@", "Array_ID", "From_Obs", "To_Feature", "Azimuth", "Altitude", "SHAPE@"]
    arrayDict = {}
    
    
    with arcpy.da.SearchCursor(in_array, array_fields) as arrayCursor: #@UndefinedVariableFromImport
        for row in arrayCursor:
            arrayDict[row[0]] = {}
            arrayDict[row[0]]["Array_ID"] = row[1]
            arrayDict[row[0]]["From_Obs"] = row[2]
            arrayDict[row[0]]["To_Feature"] = row[3]
            arrayDict[row[0]]["Azimuth"] = row[4]
            arrayDict[row[0]]["Altitude"] = row[5]
            
            # Add endpoint to dictionary.
            '''
            for part in row[6]:
                print(len(part))
                for pnt in part:
                    print(pnt, pnt.X, pnt.Y, pnt.Z)
            '''
            
            for part in row[6]:
                for pnt in part:
                    end_pt = pnt # The last coordinate of the vector should be the endpoint.
            
            arrayDict[row[0]]["Endpoint"] = [end_pt.X, end_pt.Y, end_pt.Z]
            
            #arrayDict[row[0]]["Endpoint"] = row[6].positionAlongLine(100, True)
            
            
            # Add ray length to dictionary.
            arrayDict[row[0]]["Length"] = row[6].getLength("PRESERVE_SHAPE")
    
    #------------------------------------------------------------------------------ 
    
    # Create 'pointDict'.
    pointDict = {}
    
    #------------------------------------------------------------------------------ 
    
    # Use search cursor to iterate through point features, and use insert cursor to populate output feature class.
    
    # Create insert cursor, if points are set to export.
    if output_status == True:
        # Create insert cursor for output features.
        insertfields = ["SHAPE@X", "SHAPE@Y", "SHAPE@Z", "OID_Line", "Array_ID", "From_Obs", "To_Feature", "Azimuth", "Altitude", "Dist_Along", "PICKUP"]
        insertCursor = arcpy.da.InsertCursor(out_fc, insertfields) #@UndefinedVariableFromImport
    
    # Set counter for pickup points.
    pickup_point_count = 0
    
    # Iterate.
    pointCursor_fields = ["OID_LINE", "DIST_ALONG", "SHAPE@X", "SHAPE@Y", "SHAPE@Z"]
    with arcpy.da.SearchCursor(in_points, pointCursor_fields) as pointCursor: #@UndefinedVariableFromImport
        for row in pointCursor:

            # Add 'distance along' to pointDict, if no value exists.
            #if pointDict.has_key(row[0]) == False:
            if row[0] not in pointDict.keys():
                pointDict[row[0]] = {} # Create a dictionary for the ray.
                # If distance along is less than the max distance, add shape from cursor.
                if row[1] < arrayDict[row[0]]["Length"]:
                    pointDict[row[0]]["DIST_ALONG"] = row[1]
                    pointDict[row[0]]["SHAPE@XYZ"] = [row[2],row[3],row[4]]
                
                # Otherwise add endpoint of ray.
                else:
                    pointDict[row[0]]["DIST_ALONG"] = arrayDict[row[0]]["Length"]
                    pointDict[row[0]]["SHAPE@XYZ"] = arrayDict[row[0]]["Endpoint"]
                
            # If a value exists, overwrite if new 'distance along' value is smaller than the previous value.
            else:
                if pointDict[row[0]]["DIST_ALONG"] > row[1]: # If entry already in dictionary is greater.
                    pointDict[row[0]]["DIST_ALONG"] = row[1]
                    pointDict[row[0]]["SHAPE@XYZ"] = [row[2],row[3],row[4]]
                    
    #------------------------------------------------------------------------------ 
    # Iterate through entries in 'pointDict'.
    for key in pointDict.keys():
        # For each point, attempt to retrieve values from the site mask.
        # LANDMARK MASK IS POSITIVE - IF THIS IS CHANGED, CHANGE BOOLEAN TO TRUE.
        
        ray_endpoint = arrayDict[key]["Endpoint"]
        
        check_point = pointDict[key]["SHAPE@XYZ"]
        
        if search_cell(check_point, ray_endpoint, in_grid_resolution): # Being generous, sample pts have 10m resolution.
            
            pickup_point_count += 1
            
            pointDict[key]["PICKUP"] = float(1)
            
        else:
            
            pointDict[key]["PICKUP"] = 0
        
        '''
        DEPRECATED:
        check_dist = math.sqrt((ray_x - check_x)**2 + (ray_y - check_y)**2 + (ray_z - check_z)**2)
        
        if check_dist <= 5: # Being generous, sample pts have 10m resolution.
            
            pickup_point_count += 1
            
            pointDict[key]["PICKUP"] = float(1)
            
        else:
            
            pointDict[key]["PICKUP"] = 0
        
        checkpoint = arcpy.PointGeometry(arcpy.Point(pointDict[key]["SHAPE@XYZ"][0], pointDict[key]["SHAPE@XYZ"][1]), arcpy.SpatialReference(3035)) # ASSUMES THAT POLYGONS ARE IN ETRS89.
        if landmark_mask_geometry.disjoint(checkpoint) == False:
            
            pickup_point_count += 1
        '''
        
        # Add point to output, if enabled.
        if output_status == True:
            
            # Return variables from arrayDict.
            shape_x, shape_y, shape_z = pointDict[key]["SHAPE@XYZ"]
            arrayID = arrayDict[key]["Array_ID"]
            from_obs = arrayDict[key]["From_Obs"]
            to_feature = arrayDict[key]["To_Feature"]
            azimuth = arrayDict[key]["Azimuth"]
            altitude = arrayDict[key]["Altitude"]
            dist_along = pointDict[key]["DIST_ALONG"]
            pickup_status = pointDict[key]["PICKUP"]
            
            # Insert feature.
            #print shape_x, shape_y, shape_z, str(key), arrayID, from_obs, to_feature, azimuth, altitude, dist_along
            #["SHAPE@X", "SHAPE@Y", "SHAPE@Z", "OID_Line", "Array_ID", "From_Obs", "To_Feature", "Azimuth", "Altitude", "Dist_Along"]
            insertCursor.insertRow((shape_x, shape_y, shape_z, str(key), arrayID, from_obs, to_feature, azimuth, altitude, dist_along, pickup_status))
                
        #except:
            #pass
                
    # Cleanup
    if output_status == True:
        del arrayDict
    del pointDict
    
    # Clock function.
    update_clock_dict(clock_dict, "Check_Pickup_Points", cpp_start_time)
    
    return pickup_point_count
    
    
#===============================================================================
# MAIN LOOP
#===============================================================================

# Set input.
overwrite = False
save_processing_datasets = False

#in_obs_pts = "C:\\GIS\\GAO_2023\\Data\\Temp\\test_fnet_label_z_select.shp"
in_obs_pts = "C:\\GIS\\GAO_2023\\Data\\Obs_Pts\\TD_01_ext.shp"
in_target_grid = "C:\\GIS\\GAO_2023\\Data\\Targets\\TD_10m_grid.shp"
in_grid_resolution = 10 # In metres, in this case.
in_target_mask = "C:\\GIS\\GAO_2023\\Script_Datasets\\Landmarks\\Unified\\masks\\u189.shp"

landmark_mask_geometry = arcpy.CopyFeatures_management(in_target_mask, arcpy.Geometry())[0]

in_pts_desc = arcpy.Describe(in_obs_pts)
if "VIS2D" not in [f.name for f in in_pts_desc.fields]:
    arcpy.management.AddField(in_obs_pts, "VIS2D","FLOAT")
    
# Read in elevation surface.
elev_path = "C:\\GIS\\GAO_2023\\Data\\LiDAR\\TD_02.tif"
#elev_path = "C:\\GIS\\GAO_2023\\Data\\Temp\\test_raster2.tif"

# Create raster object for intersect.
SRTM_raster = Raster(elev_path)
# Generate raster mosaic for intersect.
mosaic_cellsize_x = arcpy.GetRasterProperties_management(elev_path, "CELLSIZEX").getOutput(0)
mosaic_cellsize_y = arcpy.GetRasterProperties_management(elev_path, "CELLSIZEY").getOutput(0)
mosaic_max_cellsize = float(max((mosaic_cellsize_x, mosaic_cellsize_y))) * 1.1 # Multiply cell size by 1.1 to avoid edge effects.

# Set folder for processing datasets.
timestep_processing_folder = "C:\\GIS\\GAO_2023\\Processing"

#------------------------------------------------------------------------------ 

script_start_time = time.time()

#------------------------------------------------------------------------------ 

pr_count = 0
time_elapsed_list = []

cursor_fields = ['SHAPE@XYZ','VIS2D']

if overwrite:
    where_clause = None
else:
    where_clause = '"VIS2D" = -1'

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
    
        # Set parameters for 'Directed_Projective_Array' tool.
        
        # Get list of coordinates for target from shapefile.
        in_target_grid_list = [i[0] for i in arcpy.da.SearchCursor(in_target_grid, ["SHAPE@XYZ"])] #@UndefinedVariableFromImport
        
        # Set filepath for array.
        ray_output = "{}\\array_{}_{}.shp".format(timestep_processing_folder, str(int(obs_x)), str(int(obs_y)))
        
        Directed_Projective_Array((obs_x, obs_y, obs_z),
                                  in_target_grid_list, # Grid int filepath.
                                  ray_output,
                                  ""
                                  )           
        
        #------------------------------------------------------------------------------ 
        
        # Set variables for 'Intersect_3D_Line_With_Surface' tool.
        
        # Set folder path for rays ('pickup lines').
        #pickup_lines_folder = "{}\\{}".format(timestep_processing_folder, "pickup_lines")
        pickup_lines_file = "{}\\{}.shp".format(timestep_processing_folder, "junk_lines") # Line features are generated by the ArcGIS tool but not used by the script.
        
        # Set folder path and filename for points ('pickup points').
        #pickup_points_folder = "{}\\{}".format(timestep_processing_folder, "pickup_points")
        pickup_points_filename = "pp_{}_{}.shp".format(str(int(obs_x)), str(int(obs_y)))
        pickup_points_file = "{}\\{}".format(timestep_processing_folder, pickup_points_filename)
        
        # Intersect array with surface.
        int3d_tool_report, int3d_tool_result = Intersect_3D_Line_With_Surface(ray_output, SRTM_raster, pickup_lines_file, pickup_points_file)
        
        # If intersection function was successful, 'int3d_tool_result' will be 'True'.
        if int3d_tool_result == True:
            
            # Create variables for selected pickup point file.
            #selected_pickup_points_folder = "{}\\{}".format(timestep_processing_folder, "updated_pickup_points")
            selected_pickup_points_filename = "upp_{}_{}.shp".format(str(int(obs_x)), str(int(obs_y)))
            selected_pickup_points_file = "{}\\{}".format(timestep_processing_folder, selected_pickup_points_filename)

            # Select pickup points and return count of visible points.
            pickup_points_count = Check_Pickup_Points(ray_output, pickup_points_file, selected_pickup_points_file, save_processing_datasets, in_grid_resolution)
            
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
        time_elapsed_average = float(sum(time_elapsed_list[-100:]))/float(len(time_elapsed_list[-100:])) # Sample the last fifty elements of the list to get a moving average.
        
        #time_left = time.strftime("%D:%H:%M:%S", time.gmtime((time_elapsed_average * (obs_pt_count - pr_count))))
        time_left = round((time_elapsed_average * (obs_pt_count - pr_count))/3600,2) # Get an average in hours.
        
        print("{:<8} : Point {}/{} : {:>6} Pickup Pts. : Took {}, approx. {} hrs. to go...".format(time_elapsed_str, pr_count, obs_pt_count, pickup_points_count, row_time_elapsed_str, time_left))
        
        # Clean up.
        if save_processing_datasets == False:
            for filepath in [pickup_lines_file, pickup_points_file, selected_pickup_points_file]:
                if os.path.exists(filepath):
                    arcpy.Delete_management(filepath)   
        
print("Script finished!")
print("")
print_clock_dict(clock_dict)