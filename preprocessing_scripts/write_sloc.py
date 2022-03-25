import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from flopy.mf6 import MFSimulation
import flopy.utils.binaryfile as bf
import datetime
import pandas as pd
from osgeo import gdal, osr


def particle_cells():
    locations = pd.read_csv(os.path.join('..', 'data', 'new_locations.csv'), header=0)
    particles_ = {}
    for num, line in locations.iterrows():
        key = line['SITE']
        if key in particles_:
            existing = particles_[key]
            new_line = (line['LAYER'], line['ROW'], line['COLUMN'])
            existing.append(new_line)
            particles_[key] = existing
        else:
            new_line = (line['LAYER'], line['ROW'], line['COLUMN'])
            particles_[key] = [new_line]
    return particles_


def write_sloc():
    filename_list = ['sgrsg.sloc',
                     'rhsg.sloc',
                     'wndsg.sloc',
                     'wnghb_dom.sloc',
                     'wnghb_lbc.sloc']
    name_dict = {0: 'qsgrsg',
                 1: 'qrhsg',
                 2: 'qwnd',
                 3: 'ghbwht_dom',
                 4: 'ghbwht_lbc'}
    path = os.path.join('.')
    time_offset = 0
    drape = 1
    for name_ind, name in enumerate(filename_list):
        loc_file = open(os.path.join(path, filename_list[name_ind]), 'w')
        loc_file.write('1\n1\n')
        data = particles[name_dict[name_ind]]
        if name_ind >= 3:
            particle_count = len(data) * 25 * 3
        else:
            particle_count = len(data) * 25
        loc_file.write(str(particle_count) + ', 0\n')
        for data_line in data:
            if name_ind >= 3:
                for zrange in range(3):
                    for xrange in range(5):
                        for yrange in range(5):
                            local_x = (66 + xrange * 132) / 660
                            local_y = (66 + yrange * 132) / 660
                            local_z = 0.25 + zrange*0.25
                            line = str(data_line[0]) + ', ' + str(data_line[1]) + ', ' + str(data_line[2]) + ', ' + str(
                                local_x) \
                                   + ', ' + str(local_y) + ', ' + str(local_z) + ', ' + str(time_offset) + ', ' + str(
                                drape) + '\n'
                            loc_file.write(line)
            else:
                initial_head = heads[data_line[0] - 1, data_line[1] - 1, data_line[2] - 1]
                bottom_n = bottom[data_line[0] - 1, data_line[1] - 1, data_line[2] - 1]
                top_n = top[data_line[0] - 1, data_line[1] - 1, data_line[2] - 1]
                if bottom_n <= initial_head <= top_n:
                    local_z = (initial_head - bottom_n) / (top_n - bottom_n)
                elif initial_head > top_n:
                    local_z = 0.9
                    print('warning, particle is being placed near the top of the cell at ' + name)
                for xrange in range(5):
                    for yrange in range(5):
                        local_x = (66 + xrange*132) / 660
                        local_y = (66 + yrange*132) / 660
                        line = str(data_line[0]) + ', ' + str(data_line[1]) + ', ' + str(data_line[2]) + ', ' + str(local_x) \
                               + ', ' + str(local_y) + ', ' + str(local_z) + ', ' + str(time_offset) + ', ' + str(drape) + '\n'
                        loc_file.write(line)

            # for zrange in range(3):
            #     for xrange in range(5):
            #         for yrange in range(5):
            #             local_x = (66 + xrange * 132) / 660
            #             local_y = (66 + yrange * 132) / 660
            #             local_z = 0.25 + zrange*0.25
            #             line = str(data_line[0]) + ', ' + str(data_line[1]) + ', ' + str(data_line[2]) + ', ' + str(
            #                 local_x) \
            #                    + ', ' + str(local_y) + ', ' + str(local_z) + ', ' + str(time_offset) + ', ' + str(
            #                 drape) + '\n'
            #             loc_file.write(line)
            # if name == 'wnghb.sloc':
            #     bottom_n = bottom[data_line[0] - 1, data_line[1] - 1, data_line[2] - 1]
            #     top_n = top[data_line[0] - 1, data_line[1] - 1, data_line[2] - 1]
            #     if data_line[0] == 1:
            #         print('\n')
            #         print(top_n)
            #         print(bottom_n)
    return


def get_initial_heads():
    model_pth = os.path.join('..', '2019_MF6_Model', 'lab_model')
    sim_ = MFSimulation.load(sim_name='lab_model', version='mf6', sim_ws=model_pth, verbosity_level=0)
    gwf_ = sim_.get_model('lab_model')
    top_1_ = gwf_.dis.top.array  # array of top elevation of top layer
    botm_ = gwf_.dis.botm.array  # array of bottom elevation of all layers
    top_all = np.zeros([12, 256, 312])
    for lay in range(12):
        if lay == 0:
            top_all[lay] = top_1_
        else:
            top_all[lay] = botm_[lay-1]
    ic_ = gwf_.ic
    heads_ = ic_.strt.array  # starting heads
    idomain_ = gwf_.dis.idomain.array
    return heads_, botm_, top_all, idomain_


particles = particle_cells()
heads, bottom, top, idomain = get_initial_heads()
write_sloc()