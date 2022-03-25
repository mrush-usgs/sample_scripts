import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import pandas as pd
from mpl_toolkits import mplot3d


def build_data_structure(path_):
    pl_file = open(os.path.join(path_, 'lab_model.pl'))
    particle_groups = {1: 'Rio Hondo\nSpreading\nGrounds',
                       2: 'San Gabriel River\nSpreading\nGrounds',
                       3: 'Whittier\nNarrows\nReservoir',
                       4: 'Whittier Narrows\nUnderflow',
                       5: 'Whittier Narrows\nUnderflow'}
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
                pl_dict[particle_groups[particle_grp]][particle_num] = np.zeros([int(line_split[3]), 4])
                sub_ind = 0
            else:
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 0] = float(x_offset) + float(line_split[1]) * 0.3048
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 1] = float(y_offset) + float(line_split[2]) * 0.3048
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 2] = float(line_split[3]) * 0.3048
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 3] = int(line_split[8])  # layer yo
                sub_ind += 1
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    return pl_dict, particle_groups


def plot_paths_by_group():
    lw = 0.5
    fig = plt.figure(figsize=(8, 8))
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        axes = fig.add_subplot(2, 2, group_ind+1, projection='3d')
        for ind, prtcl in enumerate(group_data.keys()):
            pdata = group_data[prtcl]
            x_adj = -(pdata[:, 1] - pdata[0, 1]) / 5280
            y_adj = -(pdata[:, 0] - pdata[0, 0]) / 5280
            z_adj = pdata[:, 2]
            if ind == 0:
                axes.plot3D(x_adj, y_adj, z_adj, color=color_list[group_ind], linewidth=lw, label=group_name)
            else:
                axes.plot3D(x_adj, y_adj, z_adj, color=color_list[group_ind], linewidth=lw)
        axes.set_title(group_name.replace('\n', ' '), fontsize=10)
        axes.set_xlabel('Southerly Displacement (miles)', fontsize=5)
        axes.set_ylabel('Westerly Displacement (miles)', fontsize=5)
        axes.set_zlabel('Elevation (ft)', fontsize=5)
        axes.set_ylim([3, -0.5])
        axes.set_xlim([-0.5, 3])
        axes.set_zlim([-600, 50])
        axes.tick_params(labelsize=5)
    fig.suptitle('Montebello Forebay Flowpaths')
    fig.savefig(os.path.join(figure_path, '3d_paths_by_group.png'), dpi=500)


def plot_spreading():
    lw = 0.5
    fig = plt.figure()
    axes = fig.add_subplot(projection='3d')
    for group_ind, group_name in enumerate(pl_data.keys()):
        if group_ind < 3:
            group_data = pl_data[group_name]
            for ind, prtcl in enumerate(group_data.keys()):
                pdata = group_data[prtcl]
                x_adj = -(pdata[:, 1] - pdata[0, 1]) / 5280
                y_adj = -(pdata[:, 0] - pdata[0, 0]) / 5280
                z_adj = pdata[:, 2]
                if ind == 0:
                    axes.plot3D(x_adj, y_adj, z_adj, color=color_list[group_ind], linewidth=lw, label=group_name)
                else:
                    axes.plot3D(x_adj, y_adj, z_adj, color=color_list[group_ind], linewidth=lw)
    axes.set_xlabel('Southerly Displacement (miles)')
    axes.set_ylabel('Westerly Displacement (miles)')
    axes.set_zlabel('Elevation (ft)')
    axes.set_ylim([1.7, -0.5])
    axes.set_xlim([-0.5, 2.1])
    axes.set_zlim([-600, 0])
    fig.legend()
    fig.suptitle('Montebello Forebay\nFlowpaths', ma='left', ha='right')
    fig.savefig(os.path.join(figure_path, '3d_spreading_paths.png'), dpi=500)


def plot_projected_paths_by_group():
    lw = 0.5
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        min_z = 0
        for prtcl in group_data.keys():
            pdata = group_data[prtcl]
            xy = np.sqrt((pdata[:, 0] - pdata[0, 0])**2 + (pdata[:, 1]-pdata[0, 1])**2) / 5280
            z = pdata[:, 2]
            axes[group_ind].plot(xy, z, color=color_list[group_ind], linewidth=lw)
            axes[group_ind].set_title(group_name.replace('\n', ' '), fontsize=8)
            if group_ind > 1:
                axes[group_ind].set_xlabel('Horizontal Displacement (miles)', fontsize=8)
            if (group_ind == 0) or (group_ind == 2):
                axes[group_ind].set_ylabel('Elevation\n(ft)', ha='right', va='center', ma='left', rotation=0, fontsize=8)
            if np.nanmin(z) < min_z:
                min_z = np.nanmin(z)
        print(min_z)
        axes[group_ind].set_ylim([-650, 100])
        axes[group_ind].set_xlim([0, 3.1])
        axes[group_ind].tick_params(labelsize=8)
    fig.suptitle('Vertical and Lateral Flow Paths')
    fig.tight_layout()
    fig.savefig(os.path.join(figure_path, 'projected_paths_by_group.png'), dpi=500)


def plot_projected_paths():
    lw = 0.5
    fig, axes = plt.subplots()
    zorders = [2, 3, 4, 1]
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        for ind, prtcl in enumerate(group_data.keys()):
            pdata = group_data[prtcl]
            xy = np.sqrt((pdata[:, 0] - pdata[0, 0])**2 + (pdata[:, 1]-pdata[0, 1])**2) / 5280
            z = pdata[:, 2]
            if ind == 0:
                axes.plot(xy, z, color=color_list[group_ind], linewidth=lw, zorder=zorders[group_ind], label=group_name)
            else:
                axes.plot(xy, z, color=color_list[group_ind], linewidth=lw, zorder=zorders[group_ind])
        axes.set_xlabel('Horizontal Displacement (miles)', fontsize=10)
        axes.set_ylabel('Elevation (ft)', fontsize=10)
        axes.set_ylim([-650, 100])
        axes.set_xlim([0, 3.1])
        if group_ind == 3:
            axes.legend()
    fig.suptitle('Flow Paths: Montebello Forebay')

    fig.tight_layout()
    fig.savefig(os.path.join(figure_path, 'projected_paths_all.png'), dpi=500)


def plot_layer_paths_by_group():
    lw = 0.5
    fig = plt.figure(figsize=(8, 8))
    for group_ind, group_name in enumerate(pl_data.keys()):
        if group_ind == 0:
            group_data = pl_data[group_name]
            if group_ind == 4:
                plot_ind = 3
            else:
                plot_ind = group_ind
            axes = fig.add_subplot(2, 2, plot_ind+1, projection='3d')
            for ind, prtcl in enumerate(group_data.keys()):
                pdata = group_data[prtcl]
                for lay in range(12):
                    lay_subset = pdata[:, 3] == (lay + 1)
                    bool_to_int = [int(bool) for bool in lay_subset]
                    lay_diff = np.diff(bool_to_int)
                    if np.sum(lay_diff != 0) > 2:  # switches from False to True or vice versa more than twice
                        firsts = np.where(lay_diff == 1)[0].tolist()  # goes from False to True
                        lasts = np.where(lay_diff == -1)[0].tolist()  # goes from True to False
                        for chunk_ind, first in enumerate(firsts):
                            new_subset = np.copy(lay_subset)
                            new_subset[:first] = False
                            if chunk_ind != len(lasts):
                                last = lasts[chunk_ind]
                                new_subset[last + 1:] = False
                            x_adj = -(pdata[new_subset, 1] - pdata[0, 1]) / 5280
                            y_adj = -(pdata[new_subset, 0] - pdata[0, 0]) / 5280
                            z_adj = pdata[new_subset, 2]
                            if (ind == 0) and (group_ind == 0):
                                axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, label=layer_list[lay], zorder=(12 - lay),
                                                    color=scotts_colors[lay])
                            else:
                                axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, zorder=(12 - lay), color=scotts_colors[lay])
                    else:
                        x_adj = -(pdata[lay_subset, 1] - pdata[0, 1]) / 5280
                        y_adj = -(pdata[lay_subset, 0] - pdata[0, 0]) / 5280
                        z_adj = pdata[lay_subset, 2]
                        if (ind == 0) and (group_ind == 0):
                            axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, label=layer_list[lay], zorder=(12 - lay),
                                                color=scotts_colors[lay])
                        else:
                            axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, zorder=(12 - lay), color=scotts_colors[lay])
            axes.set_title(group_name.replace('\n', ' '), fontsize=10)
            axes.set_xlabel('Southerly Displacement (miles)', fontsize=5)
            axes.set_ylabel('Westerly Displacement (miles)', fontsize=5)
            axes.set_zlabel('Elevation (ft)', fontsize=5)
            axes.set_ylim([3, 0.0])
            axes.set_xlim([0.0, 3])
            axes.set_zlim([-700, 50])
            axes.tick_params(labelsize=5)
    fig.legend()
    fig.suptitle('Montebello Forebay Flowpaths')
    fig.savefig(os.path.join(figure_path, '3d_layer_paths_by_group.png'), dpi=500)
    return


def plot_layer_paths_individually():
    lw = 0.5
    for group_ind, group_name in enumerate(pl_data.keys()):
        sub_axes = plt.axes(projection='3d')
        group_data = pl_data[group_name]
        for ind, prtcl in enumerate(group_data.keys()):
            pdata = group_data[prtcl]
            for lay in range(12):
                lay_subset = pdata[:, 3] == (lay + 1)
                bool_to_int = [int(bool) for bool in lay_subset]
                lay_diff = np.diff(bool_to_int)
                if np.sum(lay_diff != 0) > 2:  # switches from False to True or vice versa more than twice
                    firsts = np.where(lay_diff == 1)[0].tolist()  # goes from False to True
                    lasts = np.where(lay_diff == -1)[0].tolist()  # goes from True to False
                    for chunk_ind, first in enumerate(firsts):
                        new_subset = np.copy(lay_subset)
                        new_subset[:first] = False
                        if chunk_ind != len(lasts):
                            last = lasts[chunk_ind]
                            new_subset[last + 1:] = False
                        x_adj = -(pdata[new_subset, 1] - pdata[0, 1]) / 5280
                        y_adj = -(pdata[new_subset, 0] - pdata[0, 0]) / 5280
                        z_adj = pdata[new_subset, 2]
                        if ind == 0:
                            sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, label=layer_list[lay], zorder=(12 - lay),
                                            color=scotts_colors[lay])
                        else:
                            sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, zorder=(12 - lay), color=scotts_colors[lay])
                else:
                    x_adj = -(pdata[lay_subset, 1] - pdata[0, 1]) / 5280
                    y_adj = -(pdata[lay_subset, 0] - pdata[0, 0]) / 5280
                    z_adj = pdata[lay_subset, 2]
                    if ind == 0:
                        sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, label=layer_list[lay], zorder=(12 - lay),
                                        color=scotts_colors[lay])
                    else:
                        sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, zorder=(12 - lay), color=scotts_colors[lay])
        sub_axes.set_title(group_name.replace('\n', ' ') + ' Flow Paths', fontsize=10)
        sub_axes.set_xlabel('Southerly Displacement (miles)', fontsize=5)
        sub_axes.set_ylabel('Westerly Displacement (miles)', fontsize=5)
        sub_axes.set_zlabel('Elevation (ft)', fontsize=5)
        if group_ind < 3:
            sub_axes.set_ylim([3, 0])
        else:
            sub_axes.set_ylim([1, 0])
        plt.show()
        # plt.legend(fontsize=5, loc=2)
        # plt.savefig(os.path.join(figure_path, group_name.replace('\n','_') + '_3d.png'), dpi=500)
    return


dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper\nWilmington A', '8: Upper\nWilmington B', '9: Lower\nWilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C\nand BC']
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
runs = ['EXT_1984_2004']
porosities = ['01']
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        pl_data, p_groups = build_data_structure(sim_path)
        # plot_paths_by_group()
        # plot_spreading()
        # plot_layer_paths_by_group()
        plot_layer_paths_individually()
        # plot_projected_paths_by_group()
        # plot_projected_paths()
