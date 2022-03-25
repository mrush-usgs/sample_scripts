import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.lines as mlines


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


def build_data_structure(path_):
    pl_file = open(os.path.join(path_, 'lab_model.pl'))
    particle_groups = {1: 'Rio Hondo\nSpreading Grounds',
                       2: 'San Gabriel River\nSpreading Grounds',
                       3: 'Whittier Narrows\nReservoir',
                       4: 'Whittier Narrows\nUnderflow (Shallow)',
                       5: 'Whittier Narrows\nUnderflow (Deep)'}
    pl_dict = {} # particle groups will be keys
    for gp in particle_groups.keys():
        pl_dict[particle_groups[gp]] = {}
    end_header = False
    for ind, line in enumerate(pl_file):
        line_split = line.split()
        if ind == 1:
            x_offset = float(line_split[2])
            y_offset = float(line_split[3]) - nrow * dCellSize
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
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 3] = float(line_split[8])  # layer yo
                sub_ind += 1
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    return pl_dict, particle_groups


def plot_projected_paths_by_group():
    lw = 0.5
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    cm = plt.get_cmap('gist_earth_r')
    rhsg_lay = []
    labeled_scatter = False
    for group_ind, group_name in enumerate(pl_data.keys()):
        if group_ind == 4:
            plot_ind = 3
        else:
            plot_ind = group_ind
        group_data = pl_data[group_name]
        min_z = 0
        for p_ind, prtcl in enumerate(group_data.keys()):
            pdata = group_data[prtcl]
            xy = np.sqrt((pdata[:, 0] - pdata[0, 0])**2 + (pdata[:, 1]-pdata[0, 1])**2) / 0.3048 / 5280
            z = pdata[:, 2] / 0.3048
            if group_ind == 2:
                unique_lays = np.unique(pdata[:, 3])
                for layers in unique_lays:
                    if layers not in rhsg_lay:
                        rhsg_lay.append(layers)
            for lay in range(12):
                if youngest_on_top:
                    zorders = 12 - lay
                else:
                    zorders = lay + 1
                lay_subset = pdata[:, 3] == (lay + 1)
                bool_to_int = [int(bool) for bool in lay_subset]
                lay_diff = np.diff(bool_to_int)
                if np.sum(lay_diff != 0) > 2: # switches from False to True or vice versa more than twice
                    firsts = np.where(lay_diff == 1)[0].tolist()  # goes from False to True
                    lasts = np.where(lay_diff == -1)[0].tolist()  # goes from True to False
                    for chunk_ind, first in enumerate(firsts):
                        new_subset = np.copy(lay_subset)
                        new_subset[:first] = False
                        if chunk_ind != len(lasts):
                            last = lasts[chunk_ind]
                            new_subset[last+1:] = False
                        if (p_ind == 0) and (group_ind == 0):
                            axes[plot_ind].plot(xy[new_subset], z[new_subset], linewidth=lw, label=layer_list[lay],
                                                 color=scotts_colors[lay], zorder=zorders)
                        else:
                            axes[plot_ind].plot(xy[new_subset], z[new_subset], linewidth=lw, color=scotts_colors[lay],
                                                zorder=zorders)
                else:
                    if (p_ind == 0) and (plot_ind == 0):
                        axes[plot_ind].plot(xy[lay_subset], z[lay_subset], linewidth=lw, label=layer_list[lay],
                                            color=scotts_colors[lay], zorder=zorders)
                    else:
                        axes[plot_ind].plot(xy[lay_subset], z[lay_subset], linewidth=lw, color=scotts_colors[lay],
                                            zorder=zorders)
                # ep_subset = (ep_data['Particle Group'] == (group_ind + 1)) & (ep_data['Particle ID'] == prtcl)
                # status = ep_data['Status'].loc[ep_subset].item()
                # if status == 5:
                #     if not labeled_scatter:
                #         axes[plot_ind].scatter(xy[-1:], z[-1:], color='k', zorder=15, s=2, linewidths=0,
                #                                label='Pathline\nTerminates\nat Well', marker='o')
                #         labeled_scatter = True
                #     else:
                #         axes[plot_ind].scatter(xy[-1:], z[-1:], color='k', zorder=15, s=2, marker='o', linewidths=0)
            if plot_ind == 3:
                axes[plot_ind].set_title('Whittier Narrows Underflow', fontsize=8)
            else:
                axes[plot_ind].set_title(group_name.replace('\n', ' '), fontsize=8)
            if plot_ind > 1:
                axes[plot_ind].set_xlabel('Horizontal Displacement (miles)', fontsize=8)
            if (plot_ind == 0) or (plot_ind == 2):
                axes[plot_ind].set_ylabel('Elevation\n(ft)', ha='right', va='center', ma='left', rotation=0, fontsize=8)
            if np.nanmin(z) < min_z:
                min_z = np.nanmin(z)
        print(min_z)
        axes[plot_ind].set_ylim([-2400, 250])
        axes[plot_ind].set_xlim([0, 13])
        axes[plot_ind].tick_params(labelsize=8)
    print(rhsg_lay)
    leg = fig.legend(loc=7, fontsize=6)
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    fig.suptitle('Vertical and Lateral Flow Paths')
    fig.tight_layout()
    fig.subplots_adjust(right=0.83)
    fig.savefig(os.path.join(figure_path, 'projected_paths_by_group.png'), dpi=500)


youngest_on_top=False
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
dCellSize = 201.168
dXCoordinate = 355046.48
dYCoordinate = 3777028.88
nrow = 256
ncol = 312
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red']
#runs = ['BASE', 'EXT_1981_2001', 'EXT_2019_1999']
runs = ['EXT_1984_2004']
# porosities = ['01', '02', '03']
porosities = ['01']
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper\nWilmington A', '8: Upper\nWilmington B', '9: Lower\nWilmington', '10: Long\nBeach A',
              '11: Long\nBeach B', '12: Long\nBeach C\nand BC']
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        pl_data, p_groups = build_data_structure(sim_path)
        ep_data = build_endpoints(sim_path)
        plot_projected_paths_by_group()
