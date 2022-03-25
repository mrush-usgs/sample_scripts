import arcpy
import numpy as np
import os

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

