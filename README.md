Contains scripts and data files required for pre- and post-processing MODPATH simulations with the LACPGM-MF6

Overall philosophy: 
this workflow is built with the goal of removing ArcGIS from our lives 
requires a GDAL installation - try Anaconda: https://opensourceoptions.com/blog/how-to-install-gdal-with-anaconda/
the scripts begin with function declarations; use folding in your preferred IDE to collapse the ones you don't need
the scripts end with calls to those functions within for loops (I think traditionally these would be enclosed in an if __name__ == '__main__': statement)

Suggestions for Streamlining: many of these scripts contain the same geospatial functions (e.g. row_col_to_utm(); ReprojectCoords(); lat_lon_labels(); zone_boundary())
check out the method in zone_boundary()
these could probably be worked into a separate package e.g. "gridmaptools"; you could then do something like "import gridmaptools"
