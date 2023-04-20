# vis.rasterizer
This project contains scripts that quantify the projective extent of landscape features using ArcGIS geoprocessing tools.

## Description
These scripts measure the projective extent (how 'large' an object appears) of an input landscape feature from a set of observer points in a GIS shapefile. It creates a '3DVIS' field in the input points and updates that field with the object's apparent size in square minutes of arc.

The approach used in this script is similar to that of rasterizers used in computer graphics - a regular projective array is constructed from the observer point to the limits of the target features, that array is intersected with a surface, and the intersection points are compared with the input features to determine which unobstructed rays intersect the target features.

These scripts have been written for specific analyses/papers, so they will probably fail spectacularly when used with other datasets! Be warned.

## Installation
These scripts use default Python 3 libraries, and [AcrPy](https://pro.arcgis.com/en/pro-app/latest/arcpy/get-started/what-is-arcpy-.htm). I use the [Eclipse IDE](https://eclipseide.org/) and [PyDev](https://www.pydev.org/) to create and manage the project.

## Usage
I run these scripts in Eclipse, there is no GUI or console option at present. For each python file there should be a clearly-marked 'parameters' section, usually near the start of the main loop.

You will need the following data to run the main script:

 - A polygon shapefile representing the extent of your target features.
 - A 3D point shapefile (i.e. with z-values) with spot heights within the target features.
 - A 3D point shapefile representing the 'eye' of one or more observers.
 - A raster digital elevation model, ideally a .tiff.

## Credits
ArcPy is a proprietary geoprocessing library made by [ESRI](https://www.esri.com/en-us/home).

This project originated in scripts I wrote for my masters and PhD research. More details and bibliographies can be found in:

 - Smith, K. 2016. "Maritime Perspectives on the Coastal Promontory Forts of Pembrokeshire." MSt thesis, University of Oxford.
 - Smith, KJ. 2020. “Modelling Seafaring in Iron Age Atlantic Europe.” PhD thesis, University of Oxford.

The latter is available via the [Oxford Research Archive](https://ora.ox.ac.uk/objects/uuid:6c266b3d-1cb4-4b43-9592-2a0db3cbe924).
