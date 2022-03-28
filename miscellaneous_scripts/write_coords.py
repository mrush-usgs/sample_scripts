import arcpy
import numpy as np
import os


# as you can see, this script requires arcpy
# coord_list.txt already exists in ../data/
# accordingly, there shouldn't really be a need to run this
# if necessary, have to run with an ArcPy Python interpreter
input_feature_class = r'..\data\2017ReportSection\a_a.shp'
coord_list = []
with arcpy.da.SearchCursor(input_feature_class, ['SHAPE@']) as s_cur:
    for row in s_cur:
        polyline = row[0]
        for feature in polyline:
            for point in feature:
                new_point = (float(str(point).split()[0]), float(str(point).split()[1]))
                coord_list.append(new_point)
np.savetxt(os.path.join('..', 'data', 'coord_list.txt'), coord_list)

