Contains scripts and data files required for pre- and post-processing MODPATH simulations with the LACPGM-MF6

This workflow is built with the goal of removing ArcGIS from our lives 

Requires a GDAL installation - try Anaconda:

conda create --name pygdal

conda activate pygdal

conda install -c conda-forge gdal=2.4.4

The latest GDAL distribution seems to cause problems with conversions from UTMs to Lat/Lon; specify gdal=2.4.4
otherwise refer to:
https://opensourceoptions.com/blog/how-to-install-gdal-with-anaconda/

Scripts that render animations require an ffmpeg install:
pip install ffmpeg-python

Scripts generally begin with function declarations; use folding in your preferred IDE to collapse the ones you don't need

Scripts generally end with calls to those functions within for loops (I think traditionally these would be enclosed in an if __name__ == '__main__': statement)

Suggestions for streamlining: many of these scripts contain the same geospatial functions (e.g. row_col_to_utm(); ReprojectCoords(); lat_lon_labels(); zone_boundary() - check out this method) that could probably be worked into a separate package e.g. "gridmaptools"; you could then do something like "import gridmaptools"
