'''
@author: Karl J. Smith
'''
#TODO: Find a style guide for file headers. Then use style guide for file headers.

# Import libraries.
import arcpy, os, time
from arcpy import env

#===============================================================================
# INTERSECT_3D_LINE WITH SURFACE
#===============================================================================

# This script runs the 'Intersect3DLineWithSurface_3d' tool, and returns error messages if the tool fails.
def intersect_3d_line_with_surface(in_array,    # Path to the input ray shapefile.
                                   in_raster,   # Path to the input raster elevation surface.
                                   pl_file,     # Path to the output polyline file (intersected rays, not used but required by arcgis).
                                   pp_file):    # Path to the output point file (intersection points).
    
    # Get the start time for the function.
    start_time = time.time()
    
    # ArcPro seems to be a bit more stable - in the ArcMap/Python 2.6 days this used to fail all the time.
    # I've kept the exception, it's good to have some error handling around the place.
    try:
        # Check out extension.
        arcpy.CheckOutExtension("3D")
        
        # Run tool.
        arcpy.Intersect3DLineWithSurface_3d(in_array, in_raster, pl_file, pp_file)
        
        # Check in extension.
        arcpy.CheckInExtension("3D")
        
        # Return True (the function completed), and the start time.
        return True, start_time
    
    except:
        
        # Print error to console.
        print("The tool 'Intersect 3D Line with Surface' encountered an error.")
        
        # Get the number of messages returned by the tool.
        message_count = arcpy.GetMessageCount()
        
        # Add messages to list.
        for i in [str(arcpy.GetMessage(j)) for j in range(0,message_count)]:
            print("  {}".format(i))
        
        # Return True (the function did not complete), and the amount of time it took to run the function.
        return False, start_time
    
#===============================================================================
# CHECK_PICKUP_POINTS
# 
# For each ray in the input, this function selects first point where the ray intersects a surface. It then tests
# whether that point intersects the 2D input features to determine whether that ray intersects the target features.
#
# If output is enabled, it will produce points representing the first intersection per ray.
#===============================================================================

def check_pickup_points(in_array,               # Path to the original, non-intersected rays file.
                        in_points,              # Path to the points of intersection produced by 3Dint.
                        target_mask_geometry,   # Path to output points (default is empty string, no points will be produced).
                        out_fc=""):             # Arcgis Geometry Object representing the target.
    
    # Get the start time for the function.
    start_time = time.time()
    
    # Allow overwrite.
    env.overwriteOutput = 1
    
    #TODO: Check that out_fc is a valid filepath, and delete if necessary. Could add overwrite as an input parameter?
    
    #------------------------------------------------------------------------------ 
    
    # Create output point features and add rows, if enabled.
    if out_fc != "":
    
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
    
    # Not sure why we're reading information from the points we just created??
    #TODO: Figure out what this section was supposed to do...
    if out_fc != "":
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
    if out_fc != "":
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
        if target_mask_geometry.disjoint(checkpoint) == False:
            
            pickup_point_count += 1

        # Add point to output, if enabled.
        if out_fc != "":
            
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
    if out_fc != "":
        del arrayDict
    del pointDict
    
    # Return True (the function completed), the number of pickup points found, and the start time.
    return True, pickup_point_count, start_time
    #TODO: Add error handling, false conditions for empty arrays, bad parameters, etc.
    