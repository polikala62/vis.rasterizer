'''
Created on Mar 16, 2023

@author: kjsmi
'''

# Import external libraries.
import arcpy, math, os, time
from arcpy import env, Raster, GetCellValue_management #@UnusedImport
from arcpy.sa import * #@UnusedWildImport

# Create dictionary to hold timing information.
clock_dict = {}

#===============================================================================
# DEFINE SCRIPT FUNCTIONS
#===============================================================================

# Define function to update clock_dict with the function name and the time it takes to run the function.
def update_clock_dict(in_dict, func_name, func_start_time):
    
    if func_name not in in_dict.keys():
        in_dict[func_name] = [(time.time() - func_start_time)]
    else:
        in_dict[func_name].append((time.time() - func_start_time))

# Define function to print the clock_dict to the console.   
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

# Define function to print a 'percent complete' message to the console.
def print_prcnt_complete(pr_count, pr_total, increment):
    
    inc_list = [i/pr_total for i in range(0,100,increment)]
    
    if pr_count/pr_total in inc_list:
        
        print("Completed {}% of total.".format(str(int(pr_count/pr_total))))
        
    elif pr_count == pr_total:
        
        print("Completed 100% of total.")

#===============================================================================
# COORDINATE DISTANCE TO GEOMETRY.

# This function does not report errors.
#===============================================================================

def Coordinate_Distance_to_Geometry(coordinates, geometry):
    
    # Get coordinate variables from input tuple.               
    coord_x, coord_y = coordinates
    
    # Create Point Geometry Object from coordinates.
    coord_geometry = arcpy.PointGeometry(arcpy.Point(coord_x, coord_y), arcpy.SpatialReference(3035))
    
    # Create list to hold distances.
    distance_list = []
    
    # Add distances to list.
    distance_list.append(coord_geometry.distanceTo(geometry))
    
    # Delete geometry.
    del coord_geometry
    
    # Get the minimum distance.
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
                              array_res, # The resolution (in minutes) of the projective array.
                              array_range, # The maximum radial coordinate (in map units) for lines in the array.
                              array_range_min, #NFX The minimum radial coordinate. Used to truncate arrays.
                              array_range_override, # If 'true', uses the 'array_range' value to determine the length of lines in the array.
                              buffer_type, # Set the buffer type. Accepted inputs are 'NONE', 'INCREMENT', & 'PERCENT'.
                              
                              inc_upper, # If buffer type = 'INCREMENT', set increment (in minutes) that output array will extend beyond its upper and lower boundaries.
                              inc_lower,
                              
                              inc_left, # If buffer type = 'INCREMENT', set increment (in minutes) that output array will extend beyond its left and right boundaries.
                              inc_right, 
                              
                              inc_length, # If buffer type = 'INCREMENT', set increment (in map units) that lines in output array will extend beyond their calculated length.
                              
                              prcnt_buffer # If buffer type = 'PERCENT', angular values for target features are incremented by a specified percentage of that value.
                              ):

    #===============================================================================
    # I. SCRIPT SETUP
    #===============================================================================
    
    # Enable overwrite.
    env.overwriteOutput = 1
    
    # Get the start time for the function.
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
    
    # If buffer type = 'NONE', set increments to 0.
    if buffer_type == 'NONE':
        inc_upper, inc_lower, inc_left, inc_right, inc_length = 0, 0, 0, 0, 0
    
    #------------------------------------------------------------------------------ 
    
    # 5. Raise exception if buffer type is not supported, or if input increments are not compatible with buffer type.
    if buffer_type not in ["NONE", "INCREMENT", "PERCENT"]:
        arcpy.AddError('Buffer type must be "NONE", "INCREMENT", or "PERCENT".')

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
        
        # Get the azimuth by calculating the arc tangent of X and Y distances to target point.
        deg_angle = math.degrees(math.atan2(Ydiff,Xdiff))
        
        # Modify angle to fit within 0-360 range.
        if deg_angle < 0:
            deg_angle += 360
            
        elif deg_angle > 360:
            deg_angle -= 360
            
        # Convert to minutes and return.
        return deg_angle * 60
    
    #------------------------------------------------------------------------------ 
    
    # 3. Define a function for calculating the altitude angle (in minutes) from an up vector defined by (x, y = 0).
    
    def Define_ALT(origin, target):
        
        # Calculate differences between maximum and minimum X, Y, and Z values.
        Xdiff, Ydiff, Zdiff = (target[0] - origin[0]), (target[1] - origin[1]), (target[2] - origin[2])
        
        # Return length of the line as a variable, by calculating the hypotenuse of a triangle with sides 'Xdiff' & 'Ydiff'.
        XY_dist = math.hypot(Xdiff, Ydiff)
        
        # Get altitude by calculating the arc tangent of the line length and elevation.
        deg_angle = math.degrees(math.atan2(XY_dist,Zdiff))
        
        # Modify angle to fit within 0-360 range.
        if deg_angle < 0:
            deg_angle += 360
            
        elif deg_angle > 360:
            deg_angle -= 360
        
        # Convert to minutes and return.
        return deg_angle * 60
    
    #------------------------------------------------------------------------------ 
    
    # 4. Define a function for calculating the angular boundaries for a feature represented by a list of azimuths.
    
    def BoundaryAngles(angle_list, angle_limit): # angle_limit should be 360 for degrees, 21600 for minutes.
        
        # Sort input angle list from smallest value to largest value.
        sorted_angle_list = sorted(angle_list)
        
        # Add last entry in the smallest angle list to the angle list.
        adjusted_angle_list = [sorted_angle_list[-1] - angle_limit] + sorted_angle_list # CHANGED FROM 360 to 21600
        
        # Create list of values where value is subtracted from succeeding value (diff_list).
        diff_list = [abs(y-x) for x, y in zip(adjusted_angle_list, adjusted_angle_list[1:])]
        
        # Get positions for maximum values in the diff_list.
        positions = [i for i, j in enumerate(diff_list) if j == max(diff_list)]
        
        # Return first value in the positions list as a variable.
        firstposition = positions[0]
    
        # If greatest difference is between first and last elements of the list, return end of list as second value.
        if firstposition == 0:
            return ([sorted_angle_list[firstposition],  sorted_angle_list[-1]])
        
        # Otherwise, return the value preceding the maximum difference value.
        else:
            return ([sorted_angle_list[firstposition], sorted_angle_list[firstposition - 1] + angle_limit])

    
    
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
    
    obsX, obsY, obsZ = in_obs_feature
    
    #------------------------------------------------------------------------------ 
        
    # Create / reset lists for azimuths, altitudes, and 3D distances defined by observer point and target points.
    AZI_list, ALT_list, length_list = [], [], []
    
    # Define a 'count' variable to record the number of lines created for the array.
    array_count = 0
    
    for TARrow in tar_coordinate_list:
            
        #------------------------------------------------------------------------------ 
        
        # 3. Calculate angular values for target point.
        
        # Return target XYZ values for target point as variables.
        tarX, tarY, tarZ = TARrow
        
        # Calculate the 3D distance between observer point and target point (see IV-1).
        tar_EDist = EDist((obsX, obsY, obsZ), (tarX, tarY, tarZ))
        
        # Calculate angular values for target point if falls within the input radius. If input radius is 0, calculate for all points.
        if array_range == 0 or tar_EDist <= array_range:
            
            # Append 3D distance between observer point and target point to list of distances.
            length_list.append(tar_EDist)
            
            # Find the azimuth angle between observer point and target point (see IV-2), and append to list of azimuths.
            compXY = Define_AZI((obsX, obsY), (tarX, tarY))
            AZI_list.append(compXY)
            
            # Find the altitude angle between observer point and target point (seeIV-3), and append to list of altitudes.
            compZ = Define_ALT((obsX, obsY, obsZ), (tarX, tarY, tarZ))
            ALT_list.append(compZ)
    
    #------------------------------------------------------------------------------ 

    # 4. Find angular boundaries.
    
    # Define a function to determine whether angular boundaries can be calculated.
    def RangeCheck(in_list):
        
        # Return a 'true' value if the buffer type is 'NONE', and 'in_list' has at least one entry.
        if buffer_type == 'NONE' and len(in_list) > 0:
            return True
        
        # Return a 'true' value if the buffer type is 'INCREMENT', and 'in_list' has at least one entry.
        elif buffer_type == 'INCREMENT' and len(in_list) > 0:
            return True
        
        # Return a 'true' value if the buffer type is 'PERCENT', and 'in_list' has at least one entry.
        elif buffer_type == 'PERCENT' and len(in_list) > 1:
            return True
        
        # If no 'true' value is returned and statistics table is enabled, add dummy entry to output statistics table.
        elif out_table != "":
            tuple_stats.append([arrayID, obsID, tarID, -9999, -9999, -9999, -9999, -9999, 0])

    # Calculate angular boundaries for target features, if points representing those features are in range.
    if RangeCheck(AZI_list):
    
        # Calculate minimum and maximum azimuths using the 'BoudaryAngles' function, which returns a tuple.
        XY_tuple = BoundaryAngles(AZI_list, 21600)
        
        # Calculate minimum and maximum altitudes, and return a tuple.
        Z_tuple = (min(ALT_list), max(ALT_list))
        
        # If buffer type is 'PERCENT', and percent values are positive, calculate increments.
        if buffer_type == "PERCENT" and prcnt_buffer > 0:
            inc_upper = (Z_tuple[1] - Z_tuple[0]) * prcnt_buffer
            inc_lower = (Z_tuple[1] - Z_tuple[0]) * prcnt_buffer
            inc_left = (XY_tuple[1] - XY_tuple[0]) * prcnt_buffer
            inc_right = (XY_tuple[1] - XY_tuple[0]) * prcnt_buffer
            inc_length = (max(length_list) - min(length_list)) * prcnt_buffer

        # Round minimum altitude and azimuth values down to the nearest multiple of the input resolution, and subtract buffer increments.
        minAltitude = int(Z_tuple[0] - (Z_tuple[0] % array_res) - inc_upper)
        minAzimuth = int(XY_tuple[0] - (XY_tuple[0] % array_res) - inc_right)
        
        # Round maximum altitude and azimuth values up to the nearest multiple of the input resolution, and add buffer increments.                
        maxAltitude = int(Z_tuple[1] - (Z_tuple[1] % array_res) + array_res + inc_lower)
        maxAzimuth = int(XY_tuple[1] - (XY_tuple[1] % array_res) + array_res + inc_left)
        
        # If 'array_override' is enabled, set the max distance to the input array range.
        if array_range_override:
            maxDistance = array_range
        
        # If 'array_override' is disabled, use measured distances to set the array range.
        else:
            maxDistance = max(length_list) + inc_length
            
            minDistance = array_range_min - inc_length #NFX used to truncate arrays.

        #------------------------------------------------------------------------------ 
        
        # 5. Generate array.

        # Create insert cursor for observer feature class.
        insert_line_fields = ["SHAPE@", "Array_ID", "From_Obs", "To_Feature", "Azimuth", "Altitude", "Length", "From_X", "From_Y", "From_Z", "To_X", "To_Y", "To_Z"]
        insert_line_cursor = arcpy.da.InsertCursor(out_fc, insert_line_fields) #@UndefinedVariableFromImport
        
        # Loop azimuth and altitude.
        for AZI in range(minAzimuth, maxAzimuth, array_res):
            for ALT in range(minAltitude, maxAltitude, array_res):
                
                # Return azimuth and altitude values in radians as variables.
                rAZI = math.radians(float(AZI) / 60)
                rALT = math.radians(float(ALT) / 60)
                
                # Calculate XYZ coordinates for line endpoints based on input radian values.
                rayX = obsX + (maxDistance * math.cos(rAZI) * math.sin(rALT))
                rayY = obsY + (maxDistance * math.sin(rAZI) * math.sin(rALT))
                rayZ = obsZ + (maxDistance * math.cos(rALT))
                
                # Return line endpoint coordinates as a list item.
                ray_pnt = [rayX, rayY, rayZ]
                
                # Get ray coordinates.
                ray_start_X = obsX + (minDistance * math.cos(rAZI) * math.sin(rALT))
                ray_start_Y = obsY + (minDistance * math.sin(rAZI) * math.sin(rALT))
                ray_start_Z = obsZ + (minDistance * math.cos(rALT))
                
                # Get line vertices as a list item.
                ray_coords = [[ray_start_X, ray_start_Y, ray_start_Z], ray_pnt]
                
                # Insert target feature.
                insert_line_cursor.insertRow((ray_coords, arrayID, obsID, tarID, AZI, ALT, maxDistance, obsX, obsY, obsZ, rayX, rayY, rayZ))
                
                # Add to 'array_count'.
                array_count += 1
        
        # Increment 'array_count_total'.
        array_count_total += array_count
        
        # If statistics output is enabled, append tuple containing statistics to 'tuple_stats'.
        if out_table != "":
            tuple_stats.append([arrayID, obsID, tarID, minAltitude, maxAltitude, minAzimuth, maxAzimuth, maxDistance, array_count])

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
    
    # Get the start time for the function.
    ils_start_time = time.time()

    try:
        # Check out extension.
        arcpy.CheckOutExtension("3D")
        
        # Run tool.
        arcpy.Intersect3DLineWithSurface_3d(ray_output, SRTM_raster, pl_file, pp_file)
        
        # Check in extension.
        arcpy.CheckInExtension("3D")
        
        # Clock function.
        update_clock_dict(clock_dict, "Intersect_3D_Line_With_Surface", ils_start_time)
        
        # Return 'True' value.
        return True
    
    except:
        
        # Get the number of messages returned by the tool.
        message_count = arcpy.GetMessageCount()
        
        # Create list to hold messages.
        message_list = []
        
        # Add messages to list.
        for i in range(0,message_count):
            message_list.append(arcpy.GetMessage(i))
            print(arcpy.GetMessage(i))
        
        # Clock function.
        update_clock_dict(clock_dict, "Intersect_3D_Line_With_Surface", ils_start_time)
        
        # Return False.
        return False
    
#===============================================================================
# Update_Pickup_Points
#
# Copied from 'Visibility Analysis' --> update pickup points.
# The script has been substantially rewritten, in order to improve performance.
#
#===============================================================================

def Check_Pickup_Points(in_array, in_points, out_fc, landmark_mask_geometry, output_status):
    
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
        for x in ("Azimuth", "Altitude", "Dist_Along"):
            arcpy.AddField_management(out_fc, x, "FLOAT")
            
    #------------------------------------------------------------------------------ 
    if output_status == True:
        # Create 'arrayDict'.
        array_fields = ["OID@", "Array_ID", "From_Obs", "To_Feature", "Azimuth", "Altitude"]
        arrayDict = {}
        
        
        with arcpy.da.SearchCursor(in_array, array_fields) as arrayCursor: #@UndefinedVariableFromImport
            for row in arrayCursor:
                arrayDict[row[0]] = {}
                arrayDict[row[0]]["Array_ID"] = row[1]
                arrayDict[row[0]]["From_Obs"] = row[2]
                arrayDict[row[0]]["To_Feature"] = row[3]
                arrayDict[row[0]]["Azimuth"] = row[4]
                arrayDict[row[0]]["Altitude"] = row[5]
    
    #------------------------------------------------------------------------------ 
    
    # Create 'pointDict'.
    pointDict = {}
    
    #------------------------------------------------------------------------------ 
    
    # Use search cursor to iterate through point features, and use insert cursor to populate output feature class.
    
    # Create insert cursor, if points are set to export.
    if output_status == True:
        # Create insert cursor for output features.
        insertfields = ["SHAPE@X", "SHAPE@Y", "SHAPE@Z", "OID_Line", "Array_ID", "From_Obs", "To_Feature", "Azimuth", "Altitude", "Dist_Along"]
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
                pointDict[row[0]] = {}
                pointDict[row[0]]["DIST_ALONG"] = row[1]
                pointDict[row[0]]["SHAPE@XYZ"] = [row[2],row[3],row[4]]
            # If a value exists, overwrite if new 'distance along' value is smaller than the previous value.
            else:
                if pointDict[row[0]]["DIST_ALONG"] > row[1]:
                    pointDict[row[0]]["DIST_ALONG"] = row[1]
                    pointDict[row[0]]["SHAPE@XYZ"] = [row[2],row[3],row[4]]
                    
    #------------------------------------------------------------------------------ 
    # Iterate through entries in 'pointDict'.
    for key in pointDict.keys():    
        # For each point, attempt to retrieve values from the site mask.
        # LANDMARK MASK IS POSITIVE - IF THIS IS CHANGED, CHANGE BOOLEAN TO TRUE.
        
        checkpoint = arcpy.PointGeometry(arcpy.Point(pointDict[key]["SHAPE@XYZ"][0], pointDict[key]["SHAPE@XYZ"][1]), arcpy.SpatialReference(3035)) # ASSUMES THAT POLYGONS ARE IN ETRS89.
        if landmark_mask_geometry.disjoint(checkpoint) == False:
            
            pickup_point_count += 1

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
            
            # Insert feature.
            insertCursor.insertRow((shape_x, shape_y, shape_z, str(key), arrayID, from_obs, to_feature, azimuth, altitude, dist_along))
            
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

in_obs_pts = "C:\\GIS\\GAO_2023\\Data\\Obs_Pts\\TD_01_ext.shp"
in_target_grid = "C:\\GIS\\GAO_2023\\Script_Datasets\\Landmarks\\Unified\\grids\\u189.shp"
in_target_mask = "C:\\GIS\\GAO_2023\\Script_Datasets\\Landmarks\\Unified\\masks\\u189.shp"

landmark_mask_geometry = arcpy.CopyFeatures_management(in_target_mask, arcpy.Geometry())[0]

in_pts_desc = arcpy.Describe(in_obs_pts)
if "VIS3D" not in [f.name for f in in_pts_desc.fields]:
    arcpy.management.AddField(in_obs_pts, "VIS3D","FLOAT")
    
# Read in elevation surface.
elev_path = "C:\\GIS\\GAO_2023\\Data\\LiDAR\\TD_01_fudge.tif"

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
        
        # Set array resolution to 10 minutes.
        ray_array_res = 10
        
        # Set array maximum range equal to the maximum visible distance plus 1x the resolution of the raster.
        ray_range_max = 20000 # 20 km should be well outside the boundary.
        
        # Set array minimum range equal to the observer-shore distance minus 1x the resolution of the raster.
        ray_range_min = 1
        
        # If 'ray_range_min' is negative, set to 0.
        if ray_range_min < 0:
            
            ray_range_min = 0
        
        Directed_Projective_Array((obs_x, obs_y, obs_z),
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
        
        #------------------------------------------------------------------------------ 
        
        # Set variables for 'Intersect_3D_Line_With_Surface' tool.
        
        # Set folder path for rays ('pickup lines').
        pickup_lines_file = "{}\\{}.shp".format(timestep_processing_folder, "junk_lines") # Line features are generated by the ArcGIS tool but not used by the script.
        
        # Set folder path and filename for points ('pickup points').
        pickup_points_filename = "pp_{}_{}.shp".format(str(int(obs_x)), str(int(obs_y)))
        pickup_points_file = "{}\\{}".format(timestep_processing_folder, pickup_points_filename)
        
        # Intersect array with surface.
        int3d_tool_result = Intersect_3D_Line_With_Surface(ray_output, SRTM_raster, pickup_lines_file, pickup_points_file)
        
        # If intersection function was successful, 'int3d_tool_result' will be 'True'.
        if int3d_tool_result:
            
            # Create variables for selected pickup point file.
            selected_pickup_points_filename = "upp_{}_{}.shp".format(str(int(obs_x)), str(int(obs_y)))
            selected_pickup_points_file = "{}\\{}".format(timestep_processing_folder, selected_pickup_points_filename)

            # Select pickup points and return count of visible points.
            pickup_points_count = Check_Pickup_Points(ray_output, pickup_points_file, selected_pickup_points_file, landmark_mask_geometry, False)
            
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
            
print("Script finished!")
print("")
print_clock_dict(clock_dict)