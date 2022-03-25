import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import flopy


def build_basin_dict():
    basins = pd.read_csv(os.path.join('..', 'data', 'gwsubbasins.csv'))
    basin_dict_ = {}
    subbasin_dict_ = {}
    basin_names = ['West Coast', 'Central', 'Santa Monica', 'Hollywood', 'Orange']
    subbasin_names = ['Montebello Forebay', 'Central Basin Pressure Area', 'Los Angeles Forebay', 'Whittier Area']
    name_dict = {}
    num_dict_ = {}
    for num, name in enumerate(basin_names):
        basin_dict_[name] = np.zeros([256, 312], dtype=int)
        num_dict_[name] = num+1
    sb_num_dict_ = {}
    for num, name in enumerate(subbasin_names):
        subbasin_dict_[name] = np.zeros([256, 312], dtype=int)
        sb_num_dict_[name] = num+1
    basin_array_ = np.zeros([256, 312], dtype=int)
    subbasin_array_ = np.zeros([256, 312], dtype=int)
    name_dict['CENTRAL BASIN, PRESSURE AREA, WRD'] = 'Central Basin Pressure Area'
    name_dict['CENTRAL BASIN PRESSURE AREA, NON-WRD'] = 'Central Basin Pressure Area'
    name_dict['LOS ANGELES FOREBAY, NON-PRESSURE AREA, NON-WRD'] = 'Los Angeles Forebay'
    name_dict['LOS ANGELES FOREBAY, NON-PRESSURE AREA, WRD'] = 'Los Angeles Forebay'
    name_dict['WHITTIER AREA'] = 'Whittier Area'
    name_dict['MONTEBELLO FOREBAY, NON-PRESSURE AREA'] = 'Montebello Forebay'
    for dummy, line in basins.iterrows():
        basin = line['GW_Basin']
        if basin in basin_dict_.keys():
            row = line['Row'] - 1
            col = line['Column_'] - 1
            basin_dict_[basin][row, col] = int(1)
            basin_num = num_dict_[basin]
            basin_array_[row, col] = basin_num
        if line['subbasin'] in name_dict.keys():
            subbasin = name_dict[line['subbasin']]
            row = line['Row'] - 1
            col = line['Column_'] - 1
            subbasin_dict_[subbasin][row, col] = int(1)
            subbasin_num = sb_num_dict_[subbasin]
            subbasin_array_[row, col] = subbasin_num
    return basin_dict_, subbasin_dict_, basin_array_, subbasin_array_


def make_node_dict(nlay, nrow, ncol):
    n_nodes = nlay * nrow * ncol
    node_list = np.arange(1, n_nodes + 1)
    node_list.shape = (nlay, nrow, ncol)
    node_dict_ = {}
    for lay in range(nlay):
        for row in range(nrow):
            for col in range(ncol):
                node = node_list[lay, row, col]
                node_dict_[node] = (lay, row, col)
    return node_dict_


def build_well_array():
    well_array_ = np.zeros([12, 256, 312])
    for line in connection_data:
        well_array_[line['cellid'][0], line['cellid'][1], line['cellid'][2]] = 1
    return well_array_


def build_data_structure(path_):
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
    return ep_df, particle_groups


def plot_zones_by_group():
    subbasin_names = ['Montebello\nForebay', 'Central Basin\nPressure Area', 'Los Angeles\nForebay', 'Whittier\nArea']
    basin_dict, subbasin_dict, basin_array, subbasin_array = build_basin_dict()
    data_arr = np.zeros([4, 4])
    fig, axes = plt.subplots(figsize=(15, 10))
    # cm = plt.get_cmap('viridis')
    total_out = 0
    for group_num in range(4):
        zone_totals = np.zeros(4)
        group_subset = (ep_data['Particle Group'] == int(group_num + 1))
        group_data = ep_data['Final Cell Number'].loc[group_subset]
        for ptcl in group_data:
            row = node_dict[ptcl][1]
            col = node_dict[ptcl][2]
            zone_no = subbasin_array[row, col] - 1
            zone_totals[zone_no] = zone_totals[zone_no] + 1
            if basin_array[row, col] != 2:
                total_out += 1
        zone_totals = 100 * zone_totals / len(group_data)
        data_arr[:, group_num] = zone_totals
    print(total_out)
    pos_cum_sum = np.nancumsum(data_arr.clip(min=0), axis=0)
    cum_pos = np.zeros(np.shape(data_arr))
    cum_pos[1:] = pos_cum_sum[:-1]
    new_color_list = ['tab:purple', 'tab:brown', 'tab:pink', 'tab:olive']
    for zone in range(4):
        axes.bar(np.arange(4), data_arr[zone], label=subbasin_names[zone], bottom=cum_pos[zone], color=new_color_list[zone]) #color=cm((zone+1)/4))
        # handles, labels = axes.get_legend_handles_labels()
    axes.set_ylabel('Percent of\nParticles', ha='right', va='center', ma='left', rotation=0, fontsize=15)
    axes.set_xlabel('Particle Group', fontsize=15)
    axes.get_xticks()
    axes.set_xticks(np.arange(4))
    axes.set_xticklabels([p_groups[gp] for gp in p_groups.keys()], fontsize=15)
    axes.set_ylim([0, 105])
    axes.tick_params(labelsize=15)
    fig.legend(loc=7, fontsize=15)
    fig.suptitle('Particle Endpoints', fontsize=20)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    fig.savefig(os.path.join(figure_path, 'zone_endpoints.png'), dpi=500)


def plot_layers_by_group():
    data_arr = np.zeros([12, 4])
    fig, axes = plt.subplots(figsize=(15, 10))
    bool_lay = np.full(12, False)
    cm = plt.get_cmap('gist_earth_r')
    for group_num in range(4):
        if group_num == 3:
            group_subset = (ep_data['Particle Group'] == int(group_num + 1)) | (ep_data['Particle Group'] == int(group_num + 2))
            group_data = ep_data['Final Layer'].loc[group_subset]
            for lay in range(12):
                subset = (group_data == int(lay + 1))
                data_arr[lay, group_num] = 100 * np.sum(subset) / len(group_data)
                if data_arr[lay, group_num] != 0:
                    bool_lay[lay] = True
        else:
            group_subset = (ep_data['Particle Group'] == int(group_num + 1))
            group_data = ep_data['Final Layer'].loc[group_subset]
            for lay in range(12):
                subset = (group_data == int(lay + 1))
                data_arr[lay, group_num] = 100 * np.sum(subset) / len(group_data)
                if data_arr[lay, group_num] != 0:
                    bool_lay[lay] = True
    pos_cum_sum = np.nancumsum(data_arr.clip(min=0), axis=0)
    cum_pos = np.zeros(np.shape(data_arr))
    cum_pos[1:] = pos_cum_sum[:-1]
    for lay in range(12):
        axes.bar(np.arange(4), data_arr[lay], label=layer_list[lay], bottom=cum_pos[lay], color=scotts_colors[lay]) #, color=cm((lay + 1)/12))
        handles, labels = axes.get_legend_handles_labels()
    axes.set_ylabel('Percent of\nParticles', ha='right', va='center', ma='left', rotation=0, fontsize=15)
    axes.set_xlabel('Particle Group', fontsize=15)
    axes.get_xticks()
    axes.set_xticks(np.arange(4))
    axes.set_xticklabels([p_groups[gp] for gp in p_groups.keys()], fontsize=15)
    axes.set_ylim([0, 100])
    axes.tick_params(labelsize=15)
    axes.invert_yaxis()
    filt_handles = [i for (i, v) in zip(handles, bool_lay) if v]
    filt_labels = [i for (i, v) in zip(labels, bool_lay) if v]
    fig.legend(filt_handles, filt_labels, loc=7, fontsize=15)
    fig.suptitle('Particle Endpoints', fontsize=20)
    fig.tight_layout()
    # fig.subplots_adjust(top=0.92)
    fig.subplots_adjust(right=0.78)
    # fig.subplots_adjust(bottom=0.18)
    fig.savefig(os.path.join(figure_path, 'layer_endpoints.png'), dpi=500)


layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper Wilmington A', '8: Upper Wilmington B', '9: Lower Wilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach\nC and BC']
particle_groups = {1: 'Rio Hondo\nSpreading Grounds',
                   2: 'San Gabriel River\nSpreading Grounds',
                   3: 'Whittier Narrows\nReservoir',
                   4: 'Whittier Narrows\nUnderflow'}
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
runs = ['EXT_1984_2004']
porosities = ['01']
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red']
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
maw = md.maw
connection_data = maw.connectiondata.get_data()
node_dict = make_node_dict(12, 256, 312)
well_array = build_well_array()
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        ep_data, p_groups = build_data_structure(sim_path)
        plot_layers_by_group()
        plot_zones_by_group()