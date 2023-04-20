'''
@author: Karl J. Smith
'''
#TODO: Find a style guide for file headers. Then use style guide for file headers.

# Import libraries.
import arcpy, math, os, time
from arcpy import env

#===========================================================================
# Directed_Projective_Array
#===========================================================================

'''
Note:   This code has been copied from the 'Directed Projected Array' script developed as part of my MST.
        It was further modified to be part of the 'Navigation Simulator' used for my DPhil.
        Some options are superfluous.
'''

def projective_array(in_obs_feature, # Filepath for the shapefile / feature class that contains observer points. Must be Z-enabled single point features.
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
    start_time = time.time()
    
    #------------------------------------------------------------------------------ 
    
    # Set 'obs_id' and 'tar_id' to the filenames of the input observer and target files.
    # I think this is inherited from an older version - might need to streamline and get rid of these variables.
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
    
    #------------------------------------------------------------------------------ 
    
    # Return True (the function completed), and the start time.
    return True, start_time
    #TODO: Add error handling, false conditions for empty arrays, bad parameters, etc.
    