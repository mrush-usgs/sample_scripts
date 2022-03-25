import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import pandas as pd
from shapely.geometry import LineString
import matplotlib.ticker as plticker
from matplotlib import markers
import matplotlib.lines as mlines

def convert_linecoords_to_figurecoords(s_, z_, ptcl_conv_list, fig_conv_list):  # this works with arrays
    # line cooordinates (s, z) to figure coordinates (x, y)
    p_x0, p_y0, p_xr, p_yr = ptcl_conv_list
    f_x0, f_y0, f_xr, f_yr = fig_conv_list
    x_slope = (f_xr - f_x0) / (p_xr - p_x0)
    y_slope = (f_yr - f_y0) / (p_yr - p_y0)
    x_intercept = p_xr - f_xr * x_slope
    y_intercept = p_yr - f_yr * y_slope
    x_ = s_ * x_slope + x_intercept
    y_ = z_ * y_slope + y_intercept
    return x_, y_


def convert_figurecoords_to_linecoords(x_, y_, ptcl_conv_list, fig_conv_list):  # this works with arrays
    # figure coordinates (x, y) to line cooordinates (s, z)
    p_x0, p_y0, p_xr, p_yr, p_ymin = ptcl_conv_list
    f_x0, f_y0, f_xr, f_yr = fig_conv_list
    x_slope = (p_xr - p_x0) / (f_xr - f_x0)
    y_slope = (p_yr - p_ymin) / (f_yr - f_y0)  # note that "slope" is derived from (p_yr - p_ymin)
    x_intercept = p_xr - f_xr * x_slope
    y_intercept = p_yr - f_yr * y_slope
    s_ = (x_ - x_intercept) / x_slope
    z_ = (y_ - y_intercept) / y_slope
    return s_, z_


def convert_utms_to_linecoords(x_, y_, s_array, distance_thresh=1000):  # this probably doesn't work with arrays
    # UTM particle cooordinates (x, y) to line coordinates (s)
    distances = np.sqrt((x_ - s_array[:, 0])**2 + (y_ - s_array[:, 1])**2)
    if np.amin(distances) < distance_thresh:
        min_ind = np.argmin(distances)
        s_ = s_array[min_ind, 2]
    else:
        s_ = np.nan
    # min_ind = np.argmin(distances)
    # s_ = s_array[min_ind, 2]
    return s_


def make_s_array(coord_list, n=int(1e4)):
    # make array that relates X, Y along line to distance S along line
    line = LineString(coord_list)
    distances = np.linspace(0, line.length, n)
    even_points = [line.interpolate(distance) for distance in distances]
    s_array = np.zeros([len(even_points) - 1, 3]) # x, y, s
    for step in range(len(s_array)):
        point = even_points[step]
        next_point = even_points[step + 1]
        s_array[step, 0] = next_point.x
        s_array[step, 1] = next_point.y
        s_array[step, 2] = np.sqrt((point.x - next_point.x)**2 + (point.y - next_point.y)**2) / 0.3048
    s_array[:, 2] = np.nancumsum(s_array[:, 2])
    return s_array


def test_fig():
    ptcl_conv_dict_ = {
        0: [160, 1220, 3760, 20, 1220],  # x_origin, y_origin, x_r, y_r, y_min (sometimes different from origin)
        1: [160, 1050, 3760, 60, 1020],
        2: [133, 438, 435, 18, 438],
        3: [160, 770, 3462, 20, 740],
        4: [133, 438, 735, 18, 438]
    }
    fig_conv_dict_ = {
        0: [0, -3000, 12000, 1000],  # x_origin, y_origin, x_r, y_r
        1: [0, -2400, 12000, 800],
        2: [0, -400, 1000, 1000],
        3: [0, -1400, 11000, 1000],
        4: [0, -400, 2000, 1000]
    }
    for group_ind, group_name in enumerate(pl_data.keys()):
        fname = os.path.join('..', 'data', 'xsec_outputs_chopped_reppeto', 'xsec_' + str(group_ind + 1) + '.jpg')
        img = plt.imread(os.path.join(fname))
        fig, axes = plt.subplots()
        axes.imshow(img)
        axes.axis('off')
        axes.xaxis.set_visible(False)
        axes.yaxis.set_visible(False)
        coords = ptcl_conv_dict_[group_ind]
        axes.scatter(coords[0], coords[1], color='b', s=0.2, linewidths=None)  # origin (always same as xmin)
        axes.scatter(coords[0], coords[4], color='c', s=0.2, linewidths=None)  # y-min (sometimes different from origin)
        axes.scatter(coords[2], coords[1], color='g', s=0.2, linewidths=None)  # xmax
        axes.scatter(coords[0], coords[3], color='r', s=0.2, linewidths=None)  # ymax
        fig.savefig(os.path.join(figure_path, group_name.replace('\n', ' ') + '_test.png'), dpi=500)
    return ptcl_conv_dict_, fig_conv_dict_


def get_xsect_line(xsect_coords_):
    line = LineString(xsect_coords_)
    n = int(1e3)
    distances = np.linspace(0, line.length, n)
    even_points = [line.interpolate(distance) for distance in distances]
    processed_list = []
    for point in even_points:
        processed_list.append((point.x, point.y))
    processed_array = np.array(processed_list)
    return processed_array


def identify_xsect_paths(xsect_coords, group_name_, distance_thresh=1000, num_thresh=1000):
    xsect_line = get_xsect_line(xsect_coords)
    prtcl_list = []
    path_count = 0
    group_data = pl_data[group_name_]
    for prtcl in group_data.keys():
        pdata = group_data[prtcl]  # x, y in columns
        distance_arr = np.array([np.amin(np.sqrt((pdata[point, 0] - xsect_line[:, 0])**2 + (pdata[point, 1] - xsect_line[:, 1])**2)) for point in range(len(pdata))])
        if np.sum(distance_arr < distance_thresh) > num_thresh:
            prtcl_list.append(prtcl)
            path_count += 1
    print(str(path_count) + ' paths identified...')
    return prtcl_list


def plot_xsect_paths():
    lims = [
        (0, 7.9, -4000, 100),
        (0, 6.45, -2000, 100),
        (0.05, 0.65, -1000, 300),
        (0, 7.1, -2000, 100),
        (0, 1.4, -1000, 300)
    ]
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        fname = os.path.join('..', 'data', 'xsec_outputs_chopped_reppeto', 'xsec_' + str(group_ind + 1) + '.jpg')
        img = plt.imread(os.path.join(fname))
        fig, axes = plt.subplots()
        ptcl_conv_dict, fig_conv_dict = test_fig()
        x_min, y_max = convert_figurecoords_to_linecoords(0, 0, ptcl_conv_dict[group_ind], fig_conv_dict[group_ind])  # upper left
        x_max, y_min = convert_figurecoords_to_linecoords(img.shape[1], img.shape[0], ptcl_conv_dict[group_ind], fig_conv_dict[group_ind])  # lower right
        new_extent = (x_min * 3.281 / 5280, x_max * 3.281 / 5280, y_min * 3.281, y_max * 3.281)
        axes.imshow(img, extent=new_extent)
        # axes.axis('off')
        # axes.xaxis.set_visible(False)
        # axes.yaxis.set_visible(False)
        if group_ind in path_dict:
            print('plotting ' + group_name.replace('\n', ' ') + '...')
            xsect_xyz = path_dict[group_ind]
            x_ptcl = xsect_xyz[:, 0]
            y_ptcl = xsect_xyz[:, 1]
            z_ptcl = xsect_xyz[:, 2]
            s_ptcl = np.zeros(len(xsect_xyz))
            first_nan = np.where(np.isnan(xsect_xyz[:, 0:2]))[0][0]
            clip_nans = xsect_xyz[0:first_nan, 0:2]
            s_array = make_s_array(clip_nans)
            for p in range(len(xsect_xyz)):
                s_ptcl[p] = convert_utms_to_linecoords(x_ptcl[p], y_ptcl[p], s_array=s_array)
            plot_s = s_ptcl / 5280
            plot_z = z_ptcl
            axes.plot(plot_s, plot_z, color='k', linewidth=1)
            # Plot Nearby Paths:
            # additional_paths = identify_xsect_paths(clip_nans, group_name)
            # for prtcl in group_data.keys():
            #     if prtcl in additional_paths:
            #         pdata = group_data[prtcl]
            #         x_ptcl = pdata[:, 0]
            #         y_ptcl = pdata[:, 1]
            #         z_ptcl = pdata[:, 2]
            #         s_ptcl = np.zeros(len(pdata))
            #         for p in range(len(pdata)):
            #             s_ptcl[p] = convert_utms_to_linecoords(x_ptcl[p], y_ptcl[p], s_array=s_array)
            #         plot_s = s_ptcl / 5280
            #         plot_z = z_ptcl
            #         axes.plot(plot_s, plot_z, color='k', linewidth=0.5)
        axes.set_xlim([lims[group_ind][0], lims[group_ind][1]])
        axes.set_ylim([lims[group_ind][2], lims[group_ind][3]])
        axes.set_ylabel('Elevation (ft)', fontsize=8)
        axes.set_xlabel('Distance along Cross-Section (miles)', fontsize=8)
        axes.set_aspect('auto')
        # axes.tick_params(labelsize=8)
        # fig.legend(fontsize=5, ncol=4, loc=8)
        # fig.tight_layout()
        # fig.subplots_adjust(bottom=0.15)
        handles = []
        # labels = []
        for lay in range(12):
            line = mlines.Line2D([], [], color=scotts_colors[lay], label=layer_list[lay], linewidth=5)
            handles.append(line)
        fig.legend(handles=handles, loc=7, fontsize=6)
        fig.subplots_adjust(right=0.8)
        fig.suptitle('Representative Flow Path: ' + group_name.replace('\n', ' '), fontsize=10)
        fig.savefig(os.path.join(figure_path, 'xsect_path_' + group_name.replace('\n', ' ') + '.png'), dpi=500)


def plot_sg_xsect_paths():
    lims = [
        (0, 7.9, -4000, 100),
        (0, 6.45, -2000, 100),
        (0.05, 0.65, -1000, 300),
        (0, 7.1, -2000, 100),
        (0, 1.4, -1000, 300)
    ]
    fig, axes = plt.subplots(2, 1)
    for group_ind, group_name in enumerate(pl_data.keys()):
        if group_ind < 2:
            group_data = pl_data[group_name]
            fname = os.path.join('..', 'data', 'xsec_outputs_chopped_reppeto', 'xsec_' + str(group_ind + 1) + '.jpg')
            img = plt.imread(os.path.join(fname))
            ptcl_conv_dict, fig_conv_dict = test_fig()
            x_min, y_max = convert_figurecoords_to_linecoords(0, 0, ptcl_conv_dict[group_ind], fig_conv_dict[group_ind]) # upper left
            x_max, y_min = convert_figurecoords_to_linecoords(img.shape[1], img.shape[0], ptcl_conv_dict[group_ind], fig_conv_dict[group_ind])  # lower right
            new_extent = (x_min * 3.281 / 5280, x_max * 3.281 / 5280, y_min * 3.281, y_max * 3.281)
            axes[group_ind].imshow(img, extent=new_extent)
            # axes.axis('off')
            # axes.xaxis.set_visible(False)
            # axes.yaxis.set_visible(False)
            if group_ind in path_dict:
                print('plotting ' + group_name.replace('\n', ' ') + '...')
                xsect_xyz = path_dict[group_ind]
                x_ptcl = xsect_xyz[:, 0]
                y_ptcl = xsect_xyz[:, 1]
                z_ptcl = xsect_xyz[:, 2]
                s_ptcl = np.zeros(len(xsect_xyz))
                first_nan = np.where(np.isnan(xsect_xyz[:, 0:2]))[0][0]
                clip_nans = xsect_xyz[0:first_nan, 0:2]
                s_array = make_s_array(clip_nans)
                for p in range(len(xsect_xyz)):
                    s_ptcl[p] = convert_utms_to_linecoords(x_ptcl[p], y_ptcl[p], s_array=s_array)
                plot_s = s_ptcl / 5280
                plot_z = z_ptcl
                axes[group_ind].plot(plot_s, plot_z, color='k', linewidth=1)
                # Plot Nearby Paths:
                # additional_paths = identify_xsect_paths(clip_nans, group_name)
                # for prtcl in group_data.keys():
                #     if prtcl in additional_paths:
                #         pdata = group_data[prtcl]
                #         x_ptcl = pdata[:, 0]
                #         y_ptcl = pdata[:, 1]
                #         z_ptcl = pdata[:, 2]
                #         s_ptcl = np.zeros(len(pdata))
                #         for p in range(len(pdata)):
                #             s_ptcl[p] = convert_utms_to_linecoords(x_ptcl[p], y_ptcl[p], s_array=s_array)
                #         plot_s = s_ptcl / 5280
                #         plot_z = z_ptcl
                #         axes.plot(plot_s, plot_z, color='k', linewidth=0.5)
                axes[group_ind].set_xlim([lims[group_ind][0], lims[group_ind][1]])
                axes[group_ind].set_ylim([lims[group_ind][2], lims[group_ind][3]])
                axes[group_ind].set_ylabel('Elevation (ft)', fontsize=6)
                axes[group_ind].set_xlabel('Distance along Cross-Section (miles)', fontsize=6)
                axes[group_ind].set_aspect('auto')
                axes[group_ind].set_title(group_name.replace('\n', ' '), fontsize=8)
                axes[group_ind].tick_params(labelsize=6)
    handles = [mlines.Line2D([], [], color='k', label='Particle Path', linewidth=1)]
    for lay in range(12):
        line = mlines.Line2D([], [], color=scotts_colors[lay], label=layer_list[lay], linewidth=5)
        handles.append(line)
    fig.legend(handles=handles, loc=7, fontsize=6)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(top=0.9)
    fig.suptitle('Representative Particle Paths for Montebello Forebay Spreading Grounds', fontsize=10)
    fig.savefig(os.path.join(figure_path, 'xsect_path_sg.png'), dpi=500)


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
    representative_paths = np.nan * np.ones([1660, 15])  # five sites, x y and z
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
        representative_paths[0:len(group_data[min_ind][:, 0]), group_ind * 3] = group_data[min_ind][:, 0]
        representative_paths[0:len(group_data[min_ind][:, 1]), group_ind * 3 + 1] = group_data[min_ind][:, 1]
        representative_paths[0:len(group_data[min_ind][:, 1]), group_ind * 3 + 2] = group_data[min_ind][:, 2]
    path_dict_ = {}
    for group_ind, group_name in enumerate(pl_data.keys()):
        rep_x = representative_paths[:, group_ind * 3]
        rep_y = representative_paths[:, group_ind * 3 + 1]
        rep_z = representative_paths[:, group_ind * 3 + 2]
        out_arr = np.zeros([np.sum(rep_x != np.nan), 3])
        out_arr[:, 0] = rep_x
        out_arr[:, 1] = rep_y
        out_arr[:, 2] = rep_z
        path_dict_[group_ind] = out_arr
        # np.savetxt('xsec_' + str(group_ind + 1) + '.csv', out_arr, delimiter=',', fmt='%.10f')
    return path_dict_


youngest_on_top = False
dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
x, y = np.meshgrid(np.linspace(dXCoordinate, dXCoordinate + ncol * dCellSize, ncol + 1), np.linspace(dYCoordinate, dYCoordinate - nrow * dCellSize, nrow + 1))
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
figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
pl_data, p_groups = build_data_structure(sim_path)
path_dict = compute_paths()
#plot_xsect_paths()
plot_sg_xsect_paths()