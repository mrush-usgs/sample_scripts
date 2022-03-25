import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import matplotlib.animation as manimation


def build_data_structure(path_):
    ts_file = open(os.path.join(path_, 'lab_model.ts'))
    ts_list = []
    end_header = False
    for ind, line in enumerate(ts_file):
        line_split = line.split()
        if ind == 1:
            x_offset = float(line_split[2])
            y_offset = float(line_split[3]) - nrow * dCellSize
        if end_header:
            new_line = [float(line_split[2]),
                        int(line_split[4]),
                        int(line_split[5]),
                        0.3048*float(line_split[10]) + x_offset,
                        0.3048*float(line_split[11]) + y_offset,
                        0.3048*float(line_split[12]),
                        int(line_split[13])]
            ts_list.append(new_line)
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    column_list = ['Time', 'Particle Group', 'Particle ID', 'X', 'Y', 'Z', 'Layer']
    ts_df = pd.DataFrame(ts_list, columns=column_list)
    return ts_df


def plot_individually():
    lw = 0.5
    for group in range(5):
        lw = 0.5
        fig = plt.figure()
        sub_axes = fig.add_subplot(projection='3d')
        tsteps = np.unique(ts_data['Time'])
        last_ts = tsteps[0]
        with writer.saving(fig, os.path.join(figure_path, particle_groups[group + 1].replace('\n', '_') + '_3d.mp4'), dpi=500):
            for ts in tsteps:
                print('processing timestep ' + str(ts))
                time_group_subset = ((ts_data['Time'] == ts) | (ts_data['Time'] == last_ts)) & (ts_data['Particle Group'] == int(group + 1))
                particles = np.unique(ts_data['Particle ID'].loc[time_group_subset])
                for ind, particle in enumerate(particles):
                    init_subset = (ts_data['Time'] == tsteps[0]) & (ts_data['Particle ID'] == particle) & (ts_data['Particle Group'] == int(group + 1))
                    for lay in range(12):
                        lay_subset = time_group_subset & (ts_data['Particle ID'] == particle) & (ts_data['Layer'] == (lay + 1))
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
                                x_adj = -(ts_data['X'].loc[new_subset] - float(ts_data['X'].loc[init_subset])) / 5280
                                y_adj = -(ts_data['Y'].loc[new_subset] - float(ts_data['Y'].loc[init_subset])) / 5280
                                z_adj = ts_data['Z'].loc[new_subset]
                                if ind == 0:
                                    sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, label=layer_list[lay],
                                                    zorder=(12 - lay), color=scotts_colors[lay])
                                else:
                                    sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, zorder=(12 - lay),
                                                    color=scotts_colors[lay])
                        elif np.nansum(lay_subset) > 0:
                            x_adj = -(ts_data['X'].loc[lay_subset] - float(ts_data['X'].loc[init_subset])) / 5280
                            y_adj = -(ts_data['Y'].loc[lay_subset] - float(ts_data['Y'].loc[init_subset])) / 5280
                            z_adj = ts_data['Z'].loc[lay_subset]
                            if ind == 0:
                                sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, label=layer_list[lay],
                                                zorder=(12 - lay),
                                                color=scotts_colors[lay])
                            else:
                                sub_axes.plot3D(x_adj, y_adj, z_adj, linewidth=lw, zorder=(12 - lay),
                                                color=scotts_colors[lay])
                sub_axes.set_title(particle_groups[group + 1].replace('\n', ' ') + ' Flow Paths', fontsize=10)
                sub_axes.set_xlabel('Southerly Displacement (miles)', fontsize=5)
                sub_axes.set_ylabel('Westerly Displacement (miles)', fontsize=5)
                sub_axes.set_zlabel('Elevation (ft)', fontsize=5)
                sub_axes.set_ylim([1.7, -0.5])
                sub_axes.set_xlim([-0.5, 2.1])
                sub_axes.set_zlim([-600, 0])
                writer.grab_frame()
                last_ts = ts
    return


def plot_spreading_grounds():
    lw = 0.5
    fig = plt.figure()
    axes = fig.add_subplot(projection='3d')
    tsteps = np.unique(ts_data['Time'])
    last_ts = tsteps[0]
    with writer.saving(fig, os.path.join(figure_path, 'spreading_3d.mp4'), dpi=500):
        for ts in tsteps:
            print('processing timestep ' + str(ts))
            for group in range(3):
                time_group_subset = ((ts_data['Time'] == ts) | (ts_data['Time'] == last_ts)) & (ts_data['Particle Group'] == int(group + 1))
                particles = np.unique(ts_data['Particle ID'].loc[time_group_subset])
                for ind, particle in enumerate(particles):
                    particle_subset = time_group_subset & (ts_data['Particle ID'] == particle)
                    init_subset = (ts_data['Time'] == tsteps[0]) & (ts_data['Particle ID'] == particle) & (ts_data['Particle Group'] == int(group + 1))
                    x_series = -(ts_data['X'].loc[particle_subset] - float(ts_data['X'].loc[init_subset])) / 5280
                    y_series = -(ts_data['Y'].loc[particle_subset] - float(ts_data['Y'].loc[init_subset])) / 5280
                    z_series = ts_data['Z'].loc[particle_subset]
                    axes.plot3D(y_series, x_series, z_series, color=color_list[group], linewidth=lw)
                axes.set_xlabel('Southerly Displacement (miles)', fontsize=5)
                axes.set_ylabel('Westerly Displacement (miles)', fontsize=5)
                axes.set_zlabel('Elevation (ft)', fontsize=5)
                axes.set_ylim([1.7, -0.5])
                axes.set_xlim([-0.5, 2.1])
                axes.set_zlim([-600, 0])
            axes.set_title('Montebello Forebay Flowpaths: ' + str(int(1971 + np.round(ts / 365, 0))))
            writer.grab_frame()
            last_ts = ts


def plot_underflow():
    lw = 0.5
    fig = plt.figure()
    axes = fig.add_subplot(projection='3d')
    tsteps = np.unique(ts_data['Time'])
    last_ts = tsteps[0]
    with writer.saving(fig, os.path.join(figure_path, 'underflow_3d.mp4'), dpi=500):
        for ts in tsteps:
            print('processing timestep ' + str(ts))
            for group in range(3, 5):
                time_group_subset = ((ts_data['Time'] == ts) | (ts_data['Time'] == last_ts)) & (ts_data['Particle Group'] == int(group + 1))
                particles = np.unique(ts_data['Particle ID'].loc[time_group_subset])
                for ind, particle in enumerate(particles):
                    particle_subset = time_group_subset & (ts_data['Particle ID'] == particle)
                    init_subset = (ts_data['Time'] == tsteps[0]) & (ts_data['Particle ID'] == particle) & (ts_data['Particle Group'] == int(group + 1))
                    x_series = -(ts_data['X'].loc[particle_subset] - float(ts_data['X'].loc[init_subset])) / 5280
                    y_series = -(ts_data['Y'].loc[particle_subset] - float(ts_data['Y'].loc[init_subset])) / 5280
                    z_series = ts_data['Z'].loc[particle_subset]
                    axes.plot3D(y_series, x_series, z_series, color=color_list[group], linewidth=lw)
                axes.set_xlabel('Southerly Displacement (miles)', fontsize=5)
                axes.set_ylabel('Westerly Displacement (miles)', fontsize=5)
                axes.set_zlabel('Elevation (ft)', fontsize=5)
                axes.set_ylim([3, -0.5])
                axes.set_xlim([-0.5, 3])
                axes.set_zlim([-600, 0])
            axes.set_title('Whittier Narrows Underflow: ' + str(int(1971 + np.round(ts / 365, 0))))
            writer.grab_frame()
            last_ts = ts


dCellSize = 201.168
dXCoordinate = 355046.48
dYCoordinate = 3777028.88
nrow = 256
ncol = 312
color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
particle_groups = {1: 'Rio Hondo\nSpreading Grounds',
                   2: 'San Gabriel River\nSpreading Grounds',
                   3: 'Whittier Narrows\nReservoir',
                   4: 'Whittier Narrows\nUnderflow (Shallow)',
                   5: 'Whittier Narrows\nUnderflow (Deep)'}
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper\nWilmington A', '8: Upper\nWilmington B', '9: Lower\nWilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C\nand BC']
num_particle = {0: 1000,
                1: 375,
                2: 200,
                3: 525}
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
runs = ['EXT_1984_2004']
porosities = ['01']
FFMpegWriter = manimation.writers['ffmpeg']
writer = FFMpegWriter(fps=5)
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        ts_data = build_data_structure(sim_path)
        #plot_individually()
        plot_spreading_grounds()
        plot_underflow()
