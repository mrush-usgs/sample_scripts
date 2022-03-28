import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import pandas as pd
import matplotlib.ticker as plticker
from matplotlib import markers


def build_endpoints(path_):
    ep_file = open(os.path.join(path_, 'lab_model.ep'))
    ep_list = []
    end_header = False
    for line in ep_file:
        line_split = line.split()
        if end_header:
            new_line = [int(line_split[1]), int(line_split[2]), int(line_split[3]), int(line_split[7]), int(line_split[17]), int(line_split[16])]
            ep_list.append(new_line)
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    column_list = ['Particle Group', 'Particle ID', 'Status', 'Initial Layer', 'Final Layer', 'Final Cell Number']
    ep_df = pd.DataFrame(data=ep_list, columns=column_list)
    return ep_df


def ReprojectCoords(coords,src_srs,tgt_srs):
    """ Reproject a list of x,y coordinates. """
    trans_coords=[]
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords


def row_col_to_utm(row_, col_, centered=False):
    if centered:
        conv_row = dYCoordinate - row_ * dCellSize - dCellSize / 2
        conv_col = dXCoordinate + col_ * dCellSize + dCellSize / 2
    else:
        conv_row = dYCoordinate - row_ * dCellSize
        conv_col = dXCoordinate + col_ * dCellSize
    return conv_row, conv_col


def build_data_structure(path_):
    pl_file = open(os.path.join(path_, 'lab_model.pl'))
    particle_groups = {1: 'Rio Hondo\nSpreading\nGrounds',
                       2: 'San Gabriel River\nSpreading\nGrounds',
                       3: 'Whittier\nNarrows\nReservoir',
                       4: 'Whittier Narrows\nUnderflow\n(Shallow)',
                       5: 'Whittier Narrows\nUnderflow\n(Deep)'}
    pl_dict = {} # particle groups will be keys
    for gp in particle_groups.keys():
        pl_dict[particle_groups[gp]] = {}
    end_header = False
    for ind, line in enumerate(pl_file):
        line_split = line.split()
        if ind == 1:
            x_offset = float(line_split[2])
            y_offset = float(line_split[3]) - nrow * dCellSize + dCellSize/2
        if end_header:
            if len(line_split) == 4:
                particle_grp = int(line_split[1])
                particle_num = int(line_split[2])
                pl_dict[particle_groups[particle_grp]][particle_num] = np.zeros([int(line_split[3]), 3])
                sub_ind = 0
            else:
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 0] = float(x_offset) + float(line_split[1]) * 0.3048
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 1] = float(y_offset) + float(line_split[2]) * 0.3048
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 2] = float(line_split[3]) * 0.3048
                sub_ind += 1
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    return pl_dict, particle_groups


def compute_paths():
    mean_paths = np.nan * np.ones([1660, 10])  # five sites, x and y
    min_max_dict = {}
    for group_name in pl_data.keys():
        min_max_dict[group_name] = [1e10, 1e10, 0, 0]
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        num_ptcls = len(group_data.keys())
        path_x = np.nan * np.ones([1660, num_ptcls])
        path_y = np.nan * np.ones([1660, num_ptcls])
        for prtcl_ind, prtcl in enumerate(group_data.keys()):
            pdata = group_data[prtcl]
            path_x[0:len(pdata[:, 0]), prtcl_ind] = pdata[:, 0]
            path_y[0:len(pdata[:, 1]), prtcl_ind] = pdata[:, 1]
        mean_paths[0:len(np.nanmean(path_x, axis=0)), group_ind*2] = np.nanmean(path_x, axis=0)
        mean_paths[0:len(np.nanmean(path_y, axis=0)), group_ind*2 + 1] = np.nanmean(path_y, axis=0)
        if np.nanmin(path_x) < min_max_dict[group_name][0]:
            min_max_dict[group_name][0] = np.nanmin(path_x)
        if np.nanmin(path_y) < min_max_dict[group_name][1]:
            min_max_dict[group_name][1] = np.nanmin(path_y)
        if np.nanmax(path_x) > min_max_dict[group_name][2]:
            min_max_dict[group_name][2] = np.nanmax(path_x)
        if np.nanmax(path_y) > min_max_dict[group_name][3]:
            min_max_dict[group_name][3] = np.nanmax(path_y)
    representative_paths = np.nan * np.ones([1660, 10])  # five sites, x and y
    for group_ind, group_name in enumerate(pl_data.keys()):
        max_frac = 0
        ind_to_ptcl = {}
        group_data = pl_data[group_name]
        num_ptcls = len(group_data.keys())
        path_distance = np.nan * np.ones([1660, num_ptcls])
        max_distance = np.sqrt((min_max_dict[group_name][2] - min_max_dict[group_name][0]) ** 2 + (min_max_dict[group_name][3] - min_max_dict[group_name][1]) ** 2)
        for prtcl_ind, prtcl in enumerate(group_data.keys()):
            pdata = group_data[prtcl]
            ind_to_ptcl[prtcl_ind] = prtcl
            prtcl_max_distance = np.sqrt((np.nanmax(pdata[:, 0]) - np.nanmin(pdata[:, 0])) ** 2 + (np.nanmax(pdata[:, 1]) - np.nanmin(pdata[:, 1])) ** 2)
            #prtcl_max_distance = np.sqrt((pdata[-1, 0] -pdata[0, 0]) ** 2 + (pdata[-1, 1] -pdata[0, 1]) ** 2)
            fraction = np.abs(prtcl_max_distance / max_distance)
            if fraction > max_frac:
                max_frac = fraction
            if fraction > 0.60:
                path_distance[0:len(pdata[:, 0]), prtcl_ind] = np.sqrt((pdata[:, 0] - mean_paths[0:len(pdata[:, 0]), group_ind * 2]) ** 2 + (pdata[:, 1] - mean_paths[0:len(pdata[:, 1]), group_ind * 2 + 1]) ** 2)
            else:
                path_distance[:, prtcl_ind] = np.nan * np.ones(len(path_distance))
        # print(max_frac)
        min_ind = ind_to_ptcl[np.nanargmin(np.nanmean(path_distance, axis=1))]  # should be column number of most representative path
        representative_paths[0:len(group_data[min_ind][:, 0]), group_ind * 2] = group_data[min_ind][:, 0]
        representative_paths[0:len(group_data[min_ind][:, 1]), group_ind * 2 + 1] = group_data[min_ind][:, 1]
    for group_ind, group_name in enumerate(pl_data.keys()):
        rep_x = representative_paths[:, group_ind * 2]
        rep_y = representative_paths[:, group_ind * 2 + 1]
        out_arr = np.zeros([np.sum(rep_x != np.nan), 2])
        out_arr[:, 0] = rep_x
        out_arr[:, 1] = rep_y
        np.savetxt('xsec_' + str(group_ind + 1) + '.csv', out_arr, delimiter=',', fmt='%.10f')
    return


youngest_on_top = False
dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
x, y = np.meshgrid(np.linspace(dXCoordinate, dXCoordinate + ncol * dCellSize, ncol + 1), np.linspace(dYCoordinate, dYCoordinate - nrow * dCellSize, nrow + 1))
fname = os.path.join('..', 'basemaps', 'basemap.tif')
img = plt.imread(os.path.join(fname))
# gdal
ds = gdal.Open(fname)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper\nWilmington A', '8: Upper\nWilmington B', '9: Lower\nWilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C\nand BC']
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
run = 'EXT_1984_2004'
porosity = '01'
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
dis = md.dis
idomain = dis.idomain.array
figure_path = os.path.join('Figures', 'POROSITY' + porosity, run)
sim_path = os.path.join('.', 'MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
pl_data, p_groups = build_data_structure(sim_path)
compute_paths()
