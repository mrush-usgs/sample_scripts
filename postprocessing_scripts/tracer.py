import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.lines as mlines
import flopy
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import linregress


def save_pickles(dict_, filename_):
    with open(os.path.join(filename_+'.pickle'), 'wb') as handle:
        pickle.dump(dict_, handle)
    return


def open_pickles(filename_):
    with open(os.path.join(filename_+'.pickle'), 'rb') as handle:
        df_ = pickle.load(handle)
    return df_


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
            row = node_dict[int(line_split[6])][1] + 1
            col = node_dict[int(line_split[6])][2] + 1
            new_line = [float(line_split[2]), int(line_split[4]), int(line_split[5]), int(line_split[13]), row, col]
            ts_list.append(new_line)
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    column_list = ['Time', 'Particle Group', 'Particle ID', 'Layer', 'Row', 'Column']
    ts_df = pd.DataFrame(ts_list, columns=column_list)
    return ts_df


def plot_concentrations(ts_data_):
    sp_list = [177, 183, 189, 177, 183, 189]  # 1-indexed: end of Q1 2015, March 2016, September 2016
    rad_range = 10
    ptcl_gp_list = [1, 2]
    with PdfPages(os.path.join(figure_path, 'sep_concentration_porosity_' + porosity + '.pdf')) as pdf:
        for radius in range(1, rad_range):
            fig, axarr = plt.subplots()
            label1 = False
            label2 = False
            s_point_list = []
            s_particle_list = []
            a_point_list = []
            a_particle_list = []
            for ind, line in well_df.iterrows():
                subset = line.iloc[5:]
                for point_ind, point in enumerate(subset):
                    if ~np.isnan(point):
                        if sp_list[point_ind] <= 180:
                            elapsed = sp_list[point_ind] * 91.25
                        else:
                            elapsed = 180 * 91.25 + (sp_list[point_ind] - 180) * 91.25/3
                        row_list = []
                        col_list = []
                        for adj in range(-radius, radius):
                            row_list.append(line['Row'] + radius)
                            col_list.append(line['Column'] + radius)
                        test_subset = ts_data_['Time'] == elapsed
                        # if ind == 0:
                        #     print(np.sum(test_subset))
                        # layer_list = [line['Layer'] - 1, line['Layer'], line['Layer'] + 1]
                        layer_list = [line['Layer']]
                        time_subset = (ts_data_['Time'] == elapsed) & (ts_data_['Layer'].isin(layer_list)) & \
                                      (ts_data_['Row'].isin(row_list)) & (ts_data_['Column'].isin(col_list)) & \
                                      (ts_data['Particle Group'].isin(ptcl_gp_list))
                        # time_subset = (ts_data_['Time'] == elapsed) & (ts_data_['Row'].isin(row_list)) & \
                        #               (ts_data_['Column'].isin(col_list)) & \
                        #               (ts_data['Particle Group'].isin(ptcl_gp_list))
                        particles_in_cell = np.sum(time_subset)

                        if particles_in_cell != 0:
                            if point_ind > 2:
                                if not label2:
                                    axarr.scatter(point, particles_in_cell, color='g', label='Sucralose')
                                    label2=True
                                else:
                                    axarr.scatter(point, particles_in_cell, color='g')
                                s_point_list.append(point)
                                s_particle_list.append(particles_in_cell)
                            else:
                                if not label1:
                                    axarr.scatter(point, particles_in_cell, color='b', label='Acesulfame')
                                    label1 = True
                                else:
                                    axarr.scatter(point, particles_in_cell, color='b')
                                a_point_list.append(point)
                                a_particle_list.append(particles_in_cell)
            # acesulfame
            a_point_array = np.array(a_point_list)
            a_particle_array = np.array(a_particle_list)
            res = linregress(a_point_array, a_particle_array)
            x_range = np.linspace(np.min(a_point_array), np.max(a_point_array))
            y = res.intercept + res.slope * x_range
            stat_lab = 'Acesulfame:\n$R^2$ = ' + str(np.round(res.rvalue**2, 2)) + '\n' + 'p = ' + '{:.1e}'.format(res.pvalue)
            axarr.plot(x_range, y, color='k', label=stat_lab)
            # sucralose
            s_point_array = np.array(s_point_list)
            s_particle_array = np.array(s_particle_list)
            res = linregress(s_point_array, s_particle_array)
            x_range = np.linspace(np.min(s_point_array), np.max(s_point_array))
            y = res.intercept + res.slope * x_range
            stat_lab = 'Sucralose:\n$R^2$ = ' + str(np.round(res.rvalue**2, 2)) + '\n' + 'p = ' + '{:.1e}'.format(res.pvalue)
            axarr.plot(x_range, y, color='k', label=stat_lab)
            axarr.legend()
            axarr.set_ylabel('Number of Particles')
            axarr.set_xlabel('Tracer Concentration')
            fig.suptitle('Porosity = ' + porosity[-1] + '0%' + ', Cell Radius = ' + str(radius))
            fig.savefig(pdf, format='pdf', dpi=200)
    return


def plot_wells(ts_data_):
    sp_list = [177, 183, 189, 177, 183, 189]  # 1-indexed: end of Q1 2015, March 2016, September 2016
    ptcl_gp_list = [1, 2]
    gp_list = ['Bell Gardens 1', 'Cerritos 2', 'Downey 1', 'Lakewood 2', 'Montebello 1', 'Norwalk 2', 'Pico 1', 'Pico 2', 'Rio Hondo 1', 'South Gate 1', 'South Gate 2']
    with PdfPages(os.path.join(figure_path, 'porosity_' + porosity + '.pdf')) as pdf:
        for gp_ind, group in enumerate(gp_list):
            fig, axarr = plt.subplots()
            s_point_list = []
            s_particle_list = []
            s_layer_list = []
            a_point_list = []
            a_particle_list = []
            a_layer_list = []
            for ind, line in well_df.iterrows():
                if group == line['Well']:
                    subset = line.iloc[5:]
                    for point_ind, point in enumerate(subset):
                        if ~np.isnan(point):
                            if sp_list[point_ind] <= 180:
                                elapsed = sp_list[point_ind] * 91.25
                            else:
                                elapsed = 180 * 91.25 + (sp_list[point_ind] - 180) * 91.25/3
                            row_list = []
                            col_list = []
                            row_list.append(line['Row'])
                            col_list.append(line['Column'])
                            layer_list = [line['Layer']]
                            time_subset = (ts_data_['Time'] == elapsed) & (ts_data_['Layer'].isin(layer_list)) & \
                                          (ts_data_['Row'].isin(row_list)) & (ts_data_['Column'].isin(col_list)) & \
                                          (ts_data_['Particle Group'].isin(ptcl_gp_list))
                            particles_in_cell = np.sum(time_subset)
                            if particles_in_cell != 0:
                                if point_ind > 2:
                                    s_point_list.append(point)
                                    s_particle_list.append(particles_in_cell)
                                    s_layer_list.append(str(line['Layer']))
                                else:
                                    a_point_list.append(point)
                                    a_particle_list.append(particles_in_cell)
                                    a_layer_list.append(str(line['Layer']))
            # acesulfame
            axarr.scatter(a_point_list, a_particle_list, color='b', label='Acesulfame')
            for label in range(len(a_layer_list)):
                axarr.annotate(a_layer_list[label], (a_point_list[label], a_particle_list[label]))
            if len(a_point_list) > 0:
                a_point_array = np.array(a_point_list)
                a_particle_array = np.array(a_particle_list)
                res = linregress(a_point_array, a_particle_array)
                x_range = np.linspace(np.min(a_point_array), np.max(a_point_array))
                y = res.intercept + res.slope * x_range
                stat_lab = 'Acesulfame:\n$R^2$ = ' + str(np.round(res.rvalue**2, 2)) + '\n' + 'p = ' + '{:.1e}'.format(res.pvalue)
                axarr.plot(x_range, y, color='k', label=stat_lab)
            # sucralose
            axarr.scatter(s_point_list, s_particle_list, color='g', label='Sucralose')
            for label in range(len(s_layer_list)):
                axarr.annotate(s_layer_list[label], (s_point_list[label], s_particle_list[label]))
            if len(s_point_list) > 0:
                s_point_array = np.array(s_point_list)
                s_particle_array = np.array(s_particle_list)
                res = linregress(s_point_array, s_particle_array)
                x_range = np.linspace(np.min(s_point_array), np.max(s_point_array))
                y = res.intercept + res.slope * x_range
                stat_lab = 'Sucralose:\n$R^2$ = ' + str(np.round(res.rvalue**2, 2)) + '\n' + 'p = ' + '{:.1e}'.format(res.pvalue)
                axarr.plot(x_range, y, color='k', label=stat_lab)
            axarr.legend()
            axarr.set_ylabel('Number of Particles')
            axarr.set_xlabel('Tracer Concentration (ng/L)')
            fig.suptitle(group)
            if group == 'Bell Gardens 1':
                fig.savefig(os.path.join(figure_path, 'bellgardens_test_' + porosity + '.png'), dpi=200)
            fig.savefig(pdf, format='pdf', dpi=200)
    return


def plot_bell_gardens(porosity_dict_):
    sp_list = [177, 183, 189, 177, 183, 189]  # 1-indexed: end of Q1 2015, March 2016, September 2016
    ptcl_gp_list = [1, 2]
    group = 'Bell Gardens 1'
    fig, axarr = plt.subplots(2, 3, figsize=(14, 5))
    axarr = axarr.flat
    for phi_ind, phi in enumerate(porosities):
        ts_data_ = porosity_dict_[phi]
        a_point_list = []
        a_particle_list = []
        a_layer_list = []
        for ind, line in well_df.iterrows():
            if group == line['Well']:
                subset = line.iloc[5:]
                for point_ind, point in enumerate(subset):
                    if ~np.isnan(point):
                        if sp_list[point_ind] <= 180:
                            elapsed = sp_list[point_ind] * 91.25
                        else:
                            elapsed = 180 * 91.25 + (sp_list[point_ind] - 180) * 91.25/3
                        row_list = []
                        col_list = []
                        row_list.append(line['Row'])
                        col_list.append(line['Column'])
                        layer_list = [line['Layer']]
                        time_subset = (ts_data_['Time'] == elapsed) & (ts_data_['Layer'].isin(layer_list)) & \
                                      (ts_data_['Row'].isin(row_list)) & (ts_data_['Column'].isin(col_list)) & \
                                      (ts_data_['Particle Group'].isin(ptcl_gp_list))
                        particles_in_cell = np.sum(time_subset)
                        if particles_in_cell != 0:
                            if point_ind <= 2:
                                a_point_list.append(point)
                                a_particle_list.append(particles_in_cell)
                                a_layer_list.append(str(line['Layer']))
        # acesulfame
        axarr[phi_ind].scatter(a_point_list, a_particle_list, color='k')
        for label in range(len(a_layer_list)):
            axarr[phi_ind].annotate(a_layer_list[label], (a_point_list[label], a_particle_list[label]))
        a_point_array = np.array(a_point_list)
        a_particle_array = np.array(a_particle_list)
        res = linregress(a_point_array, a_particle_array)
        x_range = np.linspace(np.min(a_point_array), np.max(a_point_array))
        y = res.intercept + res.slope * x_range
        stat_lab = 'Acesulfame:\n$R^2$ = ' + str(np.round(res.rvalue**2, 2)) + '\n' + 'p = ' + '{:.1e}'.format(res.pvalue)
        axarr[phi_ind].plot(x_range, y, color='k', label=stat_lab)
        axarr[phi_ind].legend()
        axarr[phi_ind].set_ylabel('Number of Particles')
        axarr[phi_ind].set_xlabel('Tracer Concentration (ng/L)')
        axarr[phi_ind].set_title('Porosity = ' + title_list[phi_ind])
        fig.suptitle(group)
        fig.tight_layout()
        fig.savefig(os.path.join(figure_path, 'bellgardens_test.png'), dpi=200)
    return


def plot_all_wells(porosity_dict_):
    sp_list = [177, 183, 189, 177, 183, 189]  # 1-indexed: end of Q1 2015, March 2016, September 2016
    ptcl_gp_list = [1, 2]
    gp_list = ['Bell Gardens 1', 'Cerritos 2', 'Downey 1', 'Lakewood 2', 'Montebello 1', 'Norwalk 2', 'Pico 1', 'Pico 2', 'Rio Hondo 1', 'South Gate 1', 'South Gate 2']
    for gp_ind, group in enumerate(gp_list):
        fig, axarr = plt.subplots(2, 3, figsize=(10, 5))
        axarr = axarr.flat
        for phi_ind, phi in enumerate(porosities):
            ts_data_ = porosity_dict_[phi]
            a_point_list = []
            a_particle_list = []
            a_layer_list = []
            for ind, line in well_df.iterrows():
                if group == line['Well']:
                    subset = line.iloc[5:]
                    for point_ind, point in enumerate(subset):
                        if ~np.isnan(point):
                            if sp_list[point_ind] <= 180:
                                elapsed = sp_list[point_ind] * 91.25
                            else:
                                elapsed = 180 * 91.25 + (sp_list[point_ind] - 180) * 91.25/3
                            row_list = []
                            col_list = []
                            row_list.append(line['Row'])
                            col_list.append(line['Column'])
                            layer_list = [line['Layer']]
                            time_subset = (ts_data_['Time'] == elapsed) & (ts_data_['Layer'].isin(layer_list)) & \
                                          (ts_data_['Row'].isin(row_list)) & (ts_data_['Column'].isin(col_list)) & \
                                          (ts_data_['Particle Group'].isin(ptcl_gp_list))
                            particles_in_cell = np.sum(time_subset)
                            if particles_in_cell != 0:
                                if point_ind <= 2:
                                    a_point_list.append(point)
                                    a_particle_list.append(particles_in_cell)
                                    a_layer_list.append(str(line['Layer']))
            # acesulfame
            if len(a_particle_list) > 0:
                axarr[phi_ind].scatter(a_point_list, a_particle_list, color='k')
                # for label in range(len(a_layer_list)):
                #     axarr[phi_ind].annotate(a_layer_list[label], (a_point_list[label], a_particle_list[label]))
                a_point_array = np.array(a_point_list)
                a_particle_array = np.array(a_particle_list)
                res = linregress(a_point_array, a_particle_array)
                x_range = np.linspace(np.min(a_point_array), np.max(a_point_array))
                y = res.intercept + res.slope * x_range
                stat_lab = '$R^2$ = ' + str(np.round(res.rvalue**2, 2)) + '\n' + 'p = ' + '{:.1e}'.format(res.pvalue)
                axarr[phi_ind].plot(x_range, y, color='k', label=stat_lab)
                axarr[phi_ind].legend()
                axarr[phi_ind].set_ylabel('Number of Particles')
                axarr[phi_ind].set_xlabel('Tracer Concentration (ng/L)')
                axarr[phi_ind].set_title('Porosity = ' + title_list[phi_ind])
        fig.suptitle(group)
        fig.tight_layout()
        fig.savefig(os.path.join(figure_path, group + '_test.png'), dpi=200)
    return


layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper Wilmington A', '8: Upper Wilmington B', '9: Lower Wilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C and BC']
particle_groups = {1: 'Rio Hondo\nSpreading Grounds',
                   2: 'San Gabriel River\nSpreading Grounds',
                   3: 'Whittier Narrows\nReservoir',
                   4: 'Whittier Narrows\nUnderflow'}
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:red']
num_act_particle = {
    0: 3000,
    1: 1125,
    2: 600,
    3: 1500
}
run = 'BASE'
porosities = ['005', '01', '015', '02', '025', '03']
title_list = ['0.05', '0.10', '0.15', '0.20', '0.25', '0.30']
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
maw = md.maw
connection_data = maw.connectiondata.get_data()
well_df = pd.read_csv(os.path.join('..', 'data', 'Concentrations.csv'))
first_time = True
if first_time:
    porosity_dict = {}
    ep_dict = {}
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures')
        sim_path = os.path.join('MODPATH', 'TRACER_SHORT3', 'POROSITY_' + porosity, run, 'lab_model')
        ts_data = build_data_structure(sim_path)
        porosity_dict[porosity] = ts_data
        # plot_concentrations(ts_data)
        # plot_wells(ts_data)
    # plot_bell_gardens(porosity_dict)
    plot_all_wells(porosity_dict)
    save_pickles(porosity_dict, 'porosity_dict')
    save_pickles(ep_dict, 'ep_dict')
else:
    porosity_dict = open_pickles('porosity_dict')
    figure_path = os.path.join('Figures')
    for porosity in porosities:
        ts_data = porosity_dict[porosity]
        plot_concentrations(ts_data)
        # plot_wells(ts_data)
    # plot_all_wells(porosity_dict)
    # plot_bell_gardens(porosity_dict)
