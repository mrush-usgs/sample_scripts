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
        if '_' in line['SITE']:
            key = line['SITE'].split('_')[0]
        else:
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
    filename_list = ['sgrsg_tracer.sloc',
                     'rhsg_tracer.sloc',
                     'wndsg_tracer.sloc']
    name_dict = {0: 'qsgrsg',
                 1: 'qrhsg',
                 2: 'qwnd'}
    path = os.path.join('.')
    drape = 1
    start_sp = 0 # 80
    for name_ind, name in enumerate(filename_list):
        loc_file = open(os.path.join(path, filename_list[name_ind]), 'w')
        loc_file.write('1\n1\n')
        data = particles[name_dict[name_ind]]
        # particle_count = len(data) * 25 * (180 * 3 + 12)
        particle_count = len(data) * 25 * ((180 - start_sp) * 3 + 12)
        loc_file.write(str(particle_count) + ', 0\n')
        for sp in range(start_sp, 192):
            if sp < 180:
                for mon in range(3):
                    time_offset = (sp * 91.25) + (mon * 30.42)
                    for data_line in data:
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
            else:
                time_offset = 180 * 91.25 + (sp - 180) * 30.42
                for data_line in data:
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
                            local_x = (66 + xrange * 132) / 660
                            local_y = (66 + yrange * 132) / 660
                            line = str(data_line[0]) + ', ' + str(data_line[1]) + ', ' + str(data_line[2]) + ', ' + str(
                                local_x) \
                                   + ', ' + str(local_y) + ', ' + str(local_z) + ', ' + str(time_offset) + ', ' + str(
                                drape) + '\n'
                            loc_file.write(line)
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