Sample scripts that demonstrate how to map gridded model outputs onto a basemap using GDAL instead of Arc products

Requires a GDAL installation - try Anaconda:

conda create --name pygdal

conda activate pygdal

conda install -c conda-forge gdal=2.4.4

The latest GDAL distribution seems to cause problems with conversions from UTMs to Lat/Lon; specify gdal=2.4.4
otherwise refer to:
https://opensourceoptions.com/blog/how-to-install-gdal-with-anaconda/

Scripts that render animations require an ffmpeg install:
conda install -c conda-forge ffmpeg
