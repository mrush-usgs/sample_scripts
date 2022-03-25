Overall philosophy: 
this workflow is built with the goal of removing ArcGIS from our lives 
requires a GDAL installation - try Anaconda: https://opensourceoptions.com/blog/how-to-install-gdal-with-anaconda/
the scripts begin with function declarations; use folding to collapse the ones you don't need
the scripts end with calls to those functions within for loops (I think traditionally these would be enclosed in an if __name__ == '__main__': statement)

Suggestions for Streamlining: many of these scripts contain the same geospatial functions (e.g. row_col_to_utm(); ReprojectCoords(); lat_lon_labels(); zone_boundary())
I am especially stoked about the method in zone_boundary()
These could probably be worked into a separate package e.g. "gridmaptools"; you could then do something like import gridmaptools

3d_animation: processes pathline file to create 3D animations of flow (not currently in article or slides but they are quite cool)
capture_time: processes time series file to create plots of time required for well capture (Figure 16)
combined_sensitivity: processes pathline file to create figure demonstrating sensitity of results (vertical and lateral flow) to model extension approach and porosity (Figure 6 and 7)
compute_likely_paths: processes pathline file to come up with a "representative" particle path for each water source; this was used to pass coordinates to DJ for cross-sections (Figure 12)
create_hydrograph_key: generic LA model script for putting together a legend that describes the layers corresponding to the "official" report colors
extend_model.py: reads in 2019_MF6_Model and writes EXT_1984_2004 model (extended through 2040 using stress period data 1984-2004)
extend_model_backwards.py: reads in 2019_MF6_Model and writes EXT_2019_1999 model (extended through 2040 using stress period data 2019-1999)
layer_flow.py: processes endpoint and pathline files to create plot of projected paths with layers called out as colors (Figure 13)
map_figure1.py: attempt at making a basemap using Python; not the neatest but it has some good stuff (Figure 1)
particle_combined.py: makes two-part figure with map of particle cells on left and particle configuration on right (Figure 5)
particle_configuration.py: makes right part of Figure 5 (not currently used in article)
particle_map.py: makes left part of Figure 5 (not currently used in article)
pathline_3d.py: makes 3D pathline plots; these are difficult to render into an image that makes sense; rather these are best used in "interactive" matplotlib mode that allow you to play with the object
plot_budget_and_percent_diff.py: compares LACPGM-MF6 and LACPGM-USG model budgets; taken from another directory so paths may need to be modified (Figure 3)
plot_ghb_layer_thickness.py: makes figure of layer thicknesses at the Whittier  Narrows GHB (in presentation, slide 5, not in article)
plot_head_map.py: makes map of 1971-1975 average heads in all layers in the Montebello Forebay (Figure 10)
plot_pumping_recharge.py: analyzes pumping and recharge and then plots the most "representative" years selected (Figure 4)
plot_spreading.py: makes bar chart of annual spreading grounds water (Figure 2)
process_endpoint.py: analyzes endpoint file and plots stacked bar of layer endpoints (Figure 15) and zone endpoints (not used)
process_pathline.py: the magnum opus, analyzes pathline file and plots pathlines on map-like visuals; lots of options to try out (Figure 9, Figure 14, slide 1, slide 29) 
process_timseries.py: a very slow and inefficient script, analyzes time series file and computes time to leave Montebello Forebay; to use this for the first time set first_time = True and run it overnight, then use first_time = False (Figure 17)
timeseries_animation.py: analyzes time series file and makes animations of particle movement (slides 16-17, not in article)
tracer.py: does linear regression between tracer concentrations and particle #s (Figure 8, top)
tracer_map.py: makes a map of particle locations for three time periods with layers called out as colors (Figure 8, bottom)
write_coords.py: makes a list of UTM coordinates from the shapefile corresponding to cross-section A-A' in LA model report (auxiliary script)
write_oc.py: changes output control (OC) package to output heads for all time steps needed for MODPATH (auxiliary script)
write_sloc.py: writes starting locations for particles (auxiliary script)
write_sloc_tracer.py: writes starting locations for particles in continuous-release simulation (auxiliary script)
xsect_pathline1.py: maps particle paths onto cross-section A-A' (Figure 11)
xsect_pathline2.py: maps particle paths onto "representative" cross-sections (Figure 12)
