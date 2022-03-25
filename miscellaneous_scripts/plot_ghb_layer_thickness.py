import os
import flopy
import numpy as np
import matplotlib.pyplot as plt


def plot_layers_by_group():
    dim_ghb = len(list_of_ghb_cells)
    thick_arr = np.zeros([12, dim_ghb])
    bot_arr = np.zeros([12, dim_ghb])
    fig, axes = plt.subplots(figsize=(15, 10))
    bool_lay = np.full(12, False)
    cm = plt.get_cmap('gist_earth_r')
    for group_num in range(dim_ghb):
        well_coords = list_of_ghb_cells[group_num]
        for lay in range(12):
            # if idomain[lay, well_coords[0], well_coords[1]] != 1:
            #     thick_arr[lay, group_num] = 0
            #     bot_arr[lay, group_num] = 0
            # print('Cell ' + str(lay) + ', ' + str(well_coords[0]) + ', ' +  str(well_coords[1]) + ' is inactive')
            if lay == 0:
                thick_arr[lay, group_num] = top[well_coords[0] - 1, well_coords[1] - 1] - botm[lay, well_coords[0] - 1, well_coords[1] - 1]
            else:
                thick_arr[lay, group_num] = botm[lay-1, well_coords[0] - 1, well_coords[1] - 1] - botm[lay, well_coords[0] - 1, well_coords[1] - 1]
            bot_arr[lay, group_num] = botm[lay, well_coords[0] - 1, well_coords[1] - 1]
            if thick_arr[lay, group_num] != 0:
                bool_lay[lay] = True
    for lay in range(12):
        axes.bar(np.arange(dim_ghb), thick_arr[lay], label=layer_list[lay], bottom=bot_arr[lay], color=cm((lay + 1)/12))
        handles, labels = axes.get_legend_handles_labels()
    axes.set_ylabel('Elevation\n(ft)', ha='right', va='center', ma='left', rotation=0, fontsize=15)
    axes.set_xlabel('Cell', fontsize=15)
    axes.get_xticks()
    axes.set_xticks(np.arange(dim_ghb))
    axes.tick_params(labelsize=15)
    axes.set_xticklabels([str(list_of_ghb_cells[gp]) for gp in range(dim_ghb)], fontsize=5)
    axes.set_ylim([-700, 300])
    filt_handles = [i for (i, v) in zip(handles, bool_lay) if v]
    filt_labels = [i for (i, v) in zip(labels, bool_lay) if v]
    fig.legend(filt_handles, filt_labels, loc=7, fontsize=15)
    fig.suptitle('Layer Elevations at Whittier Narrows General Head Boundary (GHB)', fontsize=20)
    fig.tight_layout()
    fig.subplots_adjust(right=0.78)
    fig.savefig(os.path.join('..', 'figures', 'ghb_layer_thickness.png'), dpi=500)


list_of_ghb_cells = [
    (56, 229),
    (57, 230),
    (58, 230),
    (59, 231),
    (60, 232),
    (61, 232),
    (62, 233),
    (63, 233),
    (64, 234),
    (65, 234),
    (66, 235),

]


layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper Wilmington A', '8: Upper Wilmington B', '9: Lower Wilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach\nC and BC']
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
dis = md.dis
idomain = dis.idomain.array
top = dis.top.array
botm = dis.botm.array

plot_layers_by_group()