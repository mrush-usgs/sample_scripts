import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.lines as mlines
import flopy


def save_pickles(dict_, filename_):
    with open(os.path.join(filename_+'.pickle'), 'wb') as handle:
        pickle.dump(dict_, handle)
    return


def open_pickles(filename_):
    with open(os.path.join(filename_+'.pickle'), 'rb') as handle:
        df_ = pickle.load(handle)
    return df_


def build_well_array():
    well_array_ = np.zeros([12, 256, 312])
    for line in connection_data:
        well_array_[line['cellid'][0], line['cellid'][1], line['cellid'][2]] = 1
    return well_array_


def build_basin_dict():
    basins = pd.read_csv(os.path.join('..', 'data', 'gwsubbasins.csv'))
    basin_dict_ = {}
    subbasin_dict_ = {}
    basin_names = ['West Coast', 'Central', 'Santa Monica', 'Hollywood', 'Orange']
    subbasin_names = ['Central Basin Pressure Area', 'Los Angeles Forebay', 'Whittier Area', 'Montebello Forebay']
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


def build_data_structure(sim_path):
    ts_file = open(os.path.join(sim_path, 'lab_model.ts'))
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


def build_ep_structure(path_):
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


def plot_layers_by_group():
    fig, axes = plt.subplots(3, 3) #, figsize=(13.33, 7.5))
    axes = axes.flat
    tsteps = np.unique(ts_data['Time'])
    num_ts = len(tsteps)
    ts_array = np.zeros([num_ts, 4, 9])
    for ts_ind, ts in enumerate(tsteps):
        for group in range(4):
            for lay in range(9):
                subset = (ts == ts_data['Time']) & (ts_data['Particle Group'] == int(group + 1)) & (ts_data['Layer'] == int(lay + 1))
                ts_array[ts_ind, group, lay] = 100 * np.sum(subset) / num_particle[group]
    for lay in range(9):
        for group in range(4):
            axes[lay].plot(1971 + tsteps/365, ts_array[:, group, lay], color=color_list[group], label=particle_groups[int(group+1)], linewidth=1)
        if lay == 8:
            handles, labels = axes[lay].get_legend_handles_labels()
        if np.remainder(lay, 3) == 0:
            axes[lay].set_ylabel('Percent of\nParticles\nin Layer', ha='right', va='center', ma='left', rotation=0, fontsize=5)
        if lay > 5:
            axes[lay].set_xlabel('Time (years)', fontsize=5)
        # axes[lay].set_ylim([0, 60])
        axes[lay].set_title(layer_list[lay], fontsize=5)
        axes[lay].tick_params(labelsize=5)
    fig.legend(handles, labels, loc=8, ncol=4, fontsize=8)
    fig.suptitle('Particle Distribution with Depth', fontsize=10)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.15)
    fig.savefig(os.path.join(figure_path, 'breakthrough_curve.png'), dpi=500)


def compute_travel_times():
    basin_dict, subbasin_dict, basin_array, subbasin_array = build_basin_dict()
    tsteps = np.unique(ts_data['Time'])
    p_dict = {}
    for group in range(5):
        group_subset = (ts_data['Particle Group'] == int(group + 1))
        num_ptcl = len(np.unique(ts_data['Particle ID'].loc[group_subset]))
        p_dict[group] = np.nan * np.ones(num_ptcl)
    for ts_ind, ts in enumerate(tsteps):
        for group in range(5):
            time_group_subset = (ts_data['Time'] == ts) & (ts_data['Particle Group'] == int(group + 1))
            particle_ids = ts_data['Particle ID'].loc[time_group_subset]
            for particle_ind, particle in enumerate(particle_ids):
                particle_subset = time_group_subset & (ts_data['Particle ID'] == particle)
                p_row = ts_data['Row'].loc[particle_subset].values.item()
                p_col = ts_data['Column'].loc[particle_subset].values.item()
                if (subbasin_array[p_row, p_col] != 4) and np.isnan(p_dict[group][particle_ind]):  # leaving Montebello Forebay
                    p_dict[group][particle_ind] = ts_data['Time'].loc[particle_subset].values.item()
    for group in range(5):
        mean = np.nanmean(p_dict[group]) / 365
        stdev = np.nanstd(p_dict[group]) / 365
        print('\nAt ' + particle_groups[group+1])
        print('particles take ' + str(mean) + ' +- ' + str(stdev) + ' years to leave the Montebello Forebay')
    return


def plot_mtb_particles(porosity_dict_, ep_dict_):
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    basin_dict, subbasin_dict, basin_array, subbasin_array = build_basin_dict()
    first_time_2 = True
    if first_time_2:
        processed_dict = {}
        for phi in porosity_dict_.keys():
            phi_data = porosity_dict[phi]
            ep_data_ = ep_dict_[phi]
            tsteps = np.unique(phi_data['Time'])
            num_ts = len(tsteps)
            ts_array = np.zeros([num_ts, 6])
            p_dict = {}
            for group in range(5):
                group_subset = (phi_data['Particle Group'] == int(group + 1))
                num_active = np.sum(ep_data_['Particle Group'] == int(group+1)) # & (ep_data_['Status'] != 2) & (ep_data_['Status'] != 7))
                num_ptcl = len(np.unique(phi_data['Particle ID'].loc[group_subset]))
                p_dict[group] = np.nan * np.ones(num_ptcl)
                for ts_ind, ts in enumerate(tsteps):
                    ts_array[ts_ind, 0] = ts
                    time_group_subset = (phi_data['Time'] == ts) & (phi_data['Particle Group'] == int(group + 1))
                    particle_ids = phi_data['Particle ID'].loc[time_group_subset]
                    # num_active = np.sum(time_group_subset)
                    for particle_ind, particle in enumerate(particle_ids):
                        particle_subset = time_group_subset & (phi_data['Particle ID'] == particle)
                        p_row = phi_data['Row'].loc[particle_subset].values.item()
                        p_col = phi_data['Column'].loc[particle_subset].values.item()
                        if (subbasin_array[p_row, p_col] != 4) and np.isnan(p_dict[group][particle_ind]):  # leaving Montebello Forebay
                            p_dict[group][particle_ind] = phi_data['Time'].loc[particle_subset].values.item()
                    ts_array[ts_ind, group + 1] = 100 * np.sum(~np.isnan(p_dict[group])) / num_active
            processed_dict[phi] = (ts_array, p_dict)
        save_pickles(processed_dict, run + '_processed_dict')
    else:
        processed_dict = open_pickles(run + '_processed_dict')
    for group in range(5):
        if group >= 3:
            plt_gp = 3
        else:
            plt_gp = group
        time = processed_dict['01'][0][:, 0] / 365
        max_arr = np.zeros([len(processed_dict['01'][0][:, 0] / 365), len(phi_list)])
        for ind, phi in enumerate(phi_list):
            time = processed_dict[phi][0][:, 0] / 365
            pct = processed_dict[phi][0][:, group + 1]
            max_arr[:, ind] = pct
            if ind == 1:
                axes[plt_gp].plot(time, pct, color=color_list[group])
        if group == 0:
            frac_1 = num_act_particle[0] / 3075
            frac_2 = num_act_particle[1] / 3075
            frac_3 = num_act_particle[2] / 3075
            frac_4 = num_act_particle[3] / 3075
            frac_5 = num_act_particle[4] / 3075
            phi_1 = frac_1 * processed_dict['005'][0][-1, 1] + frac_2 * processed_dict['005'][0][-1, 2] + frac_3 * processed_dict['005'][0][-1, 3] + frac_4 * processed_dict['005'][0][-1, 4] + frac_5 * processed_dict['005'][0][-1, 5]
            phi_2 = frac_1 * processed_dict['01'][0][-1, 1] + frac_2 * processed_dict['01'][0][-1, 2] + frac_3 * processed_dict['01'][0][-1, 3] + frac_4 * processed_dict['01'][0][-1, 4] + frac_5 * processed_dict['01'][0][-1, 5]
            phi_3 = frac_1 * processed_dict['015'][0][-1, 1] + frac_2 * processed_dict['015'][0][-1, 2] + frac_3 * processed_dict['015'][0][-1, 3] + frac_4 * processed_dict['015'][0][-1, 4] + frac_5 * processed_dict['015'][0][-1, 5]
            print('\n')
            print(phi_1)
            print(phi_2)
            print(phi_3)
            print('\n')
        left_line = np.amax(max_arr, axis=1)
        right_line = np.amin(max_arr, axis=1)
        if run == 'EXT_1984_2004':
            print('\n')
            print(left_line[-1])
            print(right_line[-1])
        if plt_gp > 2:
            if group == 3:
                axes[plt_gp].set_title('Whittier Narrows Underflow', fontsize=8)
                axes[plt_gp].fill_between(time, left_line, right_line, color=color_list[group], alpha=0.5, label='Shallow')
            else:
                axes[plt_gp].fill_between(time, left_line, right_line, color=color_list[group], alpha=0.5, label='Deep')
                axes[plt_gp].legend(fontsize=8, loc=7)
        else:
            axes[plt_gp].set_title(particle_groups[plt_gp + 1].replace('\n', ' '), fontsize=8)
            axes[plt_gp].fill_between(time, left_line, right_line, color=color_list[plt_gp], alpha=0.5)
        if (plt_gp == 0) or (plt_gp == 2):
            axes[plt_gp].set_ylabel('Percent of\nParticles', ha='right', va='center', ma='left', rotation=0, fontsize=8)
        if plt_gp > 1:
            axes[plt_gp].set_xlabel('Number of Years', fontsize=8)
        axes[plt_gp].set_ylim([0, 100])
        axes[plt_gp].tick_params(labelsize=8)
        if run == 'EXT_1984_2004':
            for phis in phi_list:
                if phis == '01':
                    print('\nAt ' + particle_groups[group + 1] + ' with porosity=' + phis)
                    prtcl_times = processed_dict[phis][1][group] / 365
                    print('Mean time to leave Montebello Forebay is ' + str(np.nanmean(prtcl_times)) + ' years')
                    print('Median time to leave Montebello Forebay is ' + str(np.nanmedian(prtcl_times)) + ' years')
                    print('with standard deviation ' + str(np.nanstd(prtcl_times)) + ' years')
    fig.suptitle('Time to Flow out of the Montebello Forebay')
    fig.tight_layout()
    if first_time:
        fig.savefig(os.path.join(figure_path, '..', '..', 'mtb_particles' + run+ '.png'), dpi=500)
    else:
        fig.savefig(os.path.join(figure_path, 'mtb_particles' + run + '.png'), dpi=500)


layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper Wilmington A', '8: Upper Wilmington B', '9: Lower Wilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C and BC']
particle_groups = {1: 'Rio Hondo\nSpreading Grounds',
                   2: 'San Gabriel River\nSpreading Grounds',
                   3: 'Whittier Narrows\nReservoir',
                   4: 'Whittier Narrows\nUnderflow (Shallow)',
                   5: 'Whittier Narrows\nUnderflow (Deep)'}
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
# runs = ['BASE', 'EXT_1981_2001', 'EXT_2019_1999']
num_act_particle = {
    0: 1000,
    1: 375,
    2: 200,
    3: 675,
    4: 825
}
runs = ['EXT_1984_2004']
porosities = ['005', '01', '015', '02', '025', '03']
phi_list = porosities
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
maw = md.maw
connection_data = maw.connectiondata.get_data()
well_array = False
first_time = True
for run in runs:
    if first_time:
        porosity_dict = {}
        ep_dict = {}
        for porosity in porosities:
            figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
            sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
            ts_data = build_data_structure(sim_path)
            ep_data, p_groups = build_ep_structure(sim_path)
            porosity_dict[porosity] = ts_data
            ep_dict[porosity] = ep_data
            # plot_layers_by_group()
            compute_travel_times()
        save_pickles(porosity_dict, 'porosity_dict')
        save_pickles(ep_dict, 'ep_dict')
        plot_mtb_particles(porosity_dict, ep_dict)
    else:
        porosity_dict = open_pickles('porosity_dict')
        ep_dict = open_pickles('ep_dict')
        figure_path = os.path.join('Figures')
        plot_mtb_particles(porosity_dict, ep_dict)
