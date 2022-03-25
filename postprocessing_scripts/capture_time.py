import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import flopy
import matplotlib.lines as mlines
import pickle


def save_pickles(dict_, filename_):
    with open(os.path.join(filename_+'.pickle'), 'wb') as handle:
        pickle.dump(dict_, handle)
    return


def open_pickles(filename_):
    with open(os.path.join(filename_+'.pickle'), 'rb') as handle:
        df_ = pickle.load(handle)
    return df_


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


def build_ep_data_structure(path_):
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


def build_ts_data_structure(path_):
    ts_file = open(os.path.join(path_, 'lab_model.ts'))
    ts_list = []
    end_header = False
    node_dict = make_node_dict(12, 256, 312)
    for ind, line in enumerate(ts_file):
        line_split = line.split()
        if end_header:
            row = node_dict[int(line_split[6])][1]
            col = node_dict[int(line_split[6])][2]
            new_line = [float(line_split[2]), int(line_split[4]), int(line_split[5]), int(line_split[13]), row, col]
            ts_list.append(new_line)
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    column_list = ['Time', 'Particle Group', 'Particle ID', 'Layer', 'Row', 'Column']
    ts_df = pd.DataFrame(ts_list, columns=column_list)
    return ts_df


def build_captured_dict(porosity_dict_):
    captured_dict_ = {}
    for phi in porosity_dict_.keys():
        total_captured = 0
        total_active = 0
        group_dict = {}
        ep_data_ = porosity_dict_[phi]
        for group_num in range(5):
            group_captured = 0
            group_subset = (ep_data_['Particle Group'] == int(group_num + 1))
            group_data = ep_data_.loc[group_subset]
            particle_id_list = []
            num_active = np.sum((ep_data_['Particle Group'] == int(group_num + 1)) & (ep_data_['Status'] != 2) & (ep_data_['Status'] != 7))
            for dummy, line in group_data.iterrows():
                cell_no = line['Final Cell Number']
                status = line['Status']
                if status == 5:
                    if well_array[node_dict[cell_no][0], node_dict[cell_no][1], node_dict[cell_no][2]] == 1:
                        particle_id_list.append(line['Particle ID'])
                        group_captured += 1
            group_dict[group_num] = (particle_id_list, num_active)
            # print(100 * group_captured/num_active)
            total_captured += group_captured
            total_active += num_active
        # print(100 * total_captured / total_active)
        captured_dict_[phi] = group_dict
    return captured_dict_


def plot_well_particles(captured_dict_, porosity_ts_dict_, porosity_ep_dict_):
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    processed_dict = {}
    for phi in captured_dict_.keys():
        captured_data = captured_dict_[phi]
        time_data = porosity_ts_dict_[phi]
        ep_data = porosity_ep_dict_[phi]
        tsteps = np.unique(time_data['Time'])
        num_ts = len(tsteps)
        p_dict = {}
        for group in range(5):
            particle_id_list = captured_data[group][0]
            max_ts_array = np.zeros(len(particle_id_list))
            ts_array = np.zeros([num_ts, 2])
            for particle_ind, particle in enumerate(particle_id_list):
                subset = (time_data['Particle Group'] == int(group + 1)) & \
                         (time_data['Particle ID'] == int(particle))
                max_ts_array[particle_ind] = np.amax(time_data['Time'].loc[subset])
                # ep check:
                ep_subset = (ep_data['Particle ID'] == int(particle)) & (ep_data['Particle Group'] == int(group + 1))
                if int(ep_data['Status'].loc[ep_subset].item()) != 5:
                    print('error')
            for ts_ind, ts in enumerate(tsteps):
                ts_array[ts_ind, 0] = ts
                ts_array[ts_ind, 1] = 100 * np.sum(max_ts_array <= ts) / captured_data[group][1]
            p_dict[group] = (ts_array, max_ts_array)
        processed_dict[phi] = p_dict
    for group in range(5):
        if group == 4:
            plt_gp = 3
        else:
            plt_gp = group
        max_arr = np.zeros([len(processed_dict['01'][group][0][:, 0] / 365), len(phi_list)])
        for ind, phi in enumerate(phi_list):
            time = processed_dict[phi][group][0][:, 0] / 365
            pct = processed_dict[phi][group][0][:, 1]
            max_arr[:, ind] = pct
            if ind == 1:
                axes[plt_gp].plot(time, pct, color=color_list[group])
        left_line = np.amax(max_arr, axis=1)
        right_line = np.amin(max_arr, axis=1)
        if plt_gp > 2:
            if group == 3:
                axes[plt_gp].set_title('Whittier Narrows Underflow', fontsize=8)
                axes[plt_gp].fill_between(time, left_line, right_line, color=color_list[group], alpha=0.5, label='Shallow')
            else:
                axes[plt_gp].fill_between(time, left_line, right_line, color=color_list[group], alpha=0.5, label='Deep')
                axes[plt_gp].legend(fontsize=8)
        else:
            axes[plt_gp].set_title(particle_groups[plt_gp + 1].replace('\n', ' '), fontsize=8)
            axes[plt_gp].fill_between(time, left_line, right_line, color=color_list[group], alpha=0.5)
        if run == 'EXT_1984_2004':
            print('\n')
            print(left_line[-1])
            print(max_arr[-1, 0])
            print(right_line[-1])
        if (plt_gp == 0) or (plt_gp == 2):
            axes[group].set_ylabel('Percent of\nParticles', ha='right', va='center', ma='left', rotation=0, fontsize=8)
        if plt_gp > 1:
            axes[plt_gp].set_xlabel('Number of Years', fontsize=8)
        axes[plt_gp].set_ylim([0, 60])
        axes[plt_gp].tick_params(labelsize=8)
        # if run == 'EXT_1984_2004':
        #     for phis in phi_list:
        #         if phis == '01':
        #             print('\nAt ' + particle_groups[plt_gp + 1] + ' with porosity=' + phis)
        #             prtcl_times = processed_dict[phis][group][1] / 365
        #             print('Mean time to be pumped is ' + str(np.nanmean(prtcl_times)) + ' years')
        #             print('Median time to be pumped is ' + str(np.nanmedian(prtcl_times)) + ' years')
        #             print('with standard deviation ' + str(np.nanstd(prtcl_times)) + ' years')
    fig.suptitle('Travel Time to Wells')
    fig.tight_layout()
    fig.savefig(os.path.join(figure_path, '..', '..', 'well_particles' + run+ '.png'), dpi=500)


layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper Wilmington A', '8: Upper Wilmington B', '9: Lower Wilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach\nC and BC']
particle_groups = {1: 'Rio Hondo\nSpreading Grounds',
                   2: 'San Gabriel River\nSpreading Grounds',
                   3: 'Whittier Narrows\nReservoir',
                   4: 'Whittier Narrows\nUnderflow (Shallow)',
                   5: 'Whittier Narrows\nUnderflow (Deep)'}
runs = ['EXT_1984_2004']
porosities = ['005', '01', '015', '02', '025', '03']
phi_list = porosities
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
any_sim_path = os.path.join('..', '2019_MF6_MODEL', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
maw = md.maw
connection_data = maw.connectiondata.get_data()
node_dict = make_node_dict(12, 256, 312)
well_array = build_well_array()
for run in runs:
    porosity_ep_dict = {}
    porosity_ts_dict = {}
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        ep_data = build_ep_data_structure(sim_path)
        ts_data = build_ts_data_structure(sim_path)
        porosity_ep_dict[porosity] = ep_data
        porosity_ts_dict[porosity] = ts_data
    captured_dict = build_captured_dict(porosity_ep_dict)
    plot_well_particles(captured_dict, porosity_ts_dict, porosity_ep_dict)
