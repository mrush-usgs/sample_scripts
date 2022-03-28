from shapely.geometry import LineString
import numpy as np
import os
import matplotlib.pyplot as plt
import flopy


def build_data_structure(path_):
    pl_file = open(os.path.join(path_, 'lab_model.pl'))
    particle_groups = {1: 'Rio Hondo\nSpreading\nGrounds',
                       2: 'San Gabriel River\nSpreading\nGrounds',
                       3: 'Whittier\nNarrows\nReservoir',
                       4: 'Whittier Narrows\nUnderflow\n(Shallow)',
                       5: 'Whittier Narrows\nUnderflow\n(Deep)'}
    pl_dict = {} # particle groups will be keys
    for gp in particle_groups.keys():
        pl_dict[particle_groups[gp]] = {}
    end_header = False
    lay_list = []
    for ind, line in enumerate(pl_file):
        line_split = line.split()
        if ind == 1:
            x_offset = float(line_split[2])
            y_offset = float(line_split[3]) - nrow * dCellSize + dCellSize/2
        if end_header:
            if len(line_split) == 4:
                particle_grp = int(line_split[1])
                particle_num = int(line_split[2])
                pl_dict[particle_groups[particle_grp]][particle_num] = np.zeros([int(line_split[3]), 3])
                sub_ind = 0
            else:
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 0] = float(x_offset) + float(line_split[1]) * 0.3048
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 1] = float(y_offset) + float(line_split[2]) * 0.3048
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 2] = float(line_split[3]) * 0.3048
                if particle_grp == 5:
                    if float(line_split[8]) not in lay_list:
                        lay_list.append(float(line_split[8]))
                sub_ind += 1
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    # print(lay_list)
    return pl_dict, particle_groups


def get_xsect_coords():
    coord_list = np.loadtxt(os.path.join('..', 'data', 'coord_list.txt'))
    line = LineString(coord_list)
    n = int(1e3)
    distances = np.linspace(0, line.length, n)
    even_points = [line.interpolate(distance) for distance in distances]
    processed_list = []
    for point in even_points:
        processed_list.append((point.x, point.y))
    processed_array = np.array(processed_list)
    return processed_array


def identify_xsect_paths(distance_thresh=1000, num_thresh=1000):
    xsect_path_dict = {}
    path_count = 0
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        for prtcl in group_data.keys():
            pdata = group_data[prtcl]  # x, y in columns
            distance_arr = np.array([np.amin(np.sqrt((pdata[point, 0] - xsect_coords[:, 0])**2 + (pdata[point, 1] - xsect_coords[:, 1])**2)) for point in range(len(pdata))])
            if np.sum(distance_arr < distance_thresh) > num_thresh:
                if group_name not in xsect_path_dict:
                    xsect_path_dict[group_name] = [prtcl]
                else:
                    existing = xsect_path_dict[group_name]
                    existing.append(prtcl)
                    xsect_path_dict[group_name] = existing
                path_count += 1
    print(str(path_count) + ' paths identified...')
    return xsect_path_dict


def convert_linecoords_to_figurecoords(s_, z_):  # this works with arrays
    # line cooordinates (s, z) to figure coordinates (x, y)
    x_slope = 9820 / 220000  # s range 220000 corresponds to  (10150 - 330 = 9820)
    y_slope = - 4120 / 9000  # z range (1000 - -8000 = 9000) corresponds to y range (4310 - 190 = 4120)
    x_intercept = 10150 - 220000 * x_slope  # s = 220000 corresponds to x = 10150
    y_intercept = 190 - 1000 * y_slope  # z = 1000 corresponds to y = 190
    x_ = s_ * x_slope + x_intercept
    y_ = z_ * y_slope + y_intercept
    return x_, y_


def convert_figurecoords_to_linecoords(x_, y_):  # this works with arrays
    # figure coordinates (x, y) to line cooordinates (s, z)
    x_slope = 9820 / (220000 / 4)  # s range 220000/4 corresponds to  (10150 - 330 = 9820)
    y_slope = - 4120 / 9000  # z range (1000 - -8000 = 9000) corresponds to y range (4310 - 190 = 4120)
    x_intercept = 10150 - (220000 / 4) * x_slope  # s = 220000 corresponds to x = 10150
    y_intercept = 190 - 1000 * y_slope  # z = 1000 corresponds to y = 190
    s_ = (x_ - x_intercept) / x_slope
    z_ = (y_ - y_intercept) / y_slope
    return s_, z_


def convert_utms_to_linecoords(x_, y_, s_array, distance_thresh=1000):  # this probably doesn't work with arrays
    # UTM particle cooordinates (x, y) to line coordinates (s)
    distances = np.sqrt((x_ - s_array[:, 0])**2 + (y_ - s_array[:, 1])**2)
    # if np.amin(distances) < distance_thresh:
    #     min_ind = np.argmin(distances)
    #     s_ = s_array[min_ind, 2]
    # else:
    #     s_ = np.nan
    min_ind = np.argmin(distances)
    s_ = s_array[min_ind, 2]
    # min_ind = np.argmin(distances)
    # s_ = s_array[min_ind, 2]
    return s_


def make_s_array(n=int(1e4)):
    # make array that relates X, Y along line to distance S along line
    coord_list = np.flip(np.loadtxt('..', 'data', 'coord_list.txt'), axis=0)
    line = LineString(coord_list)
    distances = np.linspace(0, line.length, n)
    even_points = [line.interpolate(distance) for distance in distances]
    s_array = np.zeros([len(even_points) - 1, 3]) # x, y, s
    for step in range(len(s_array)):
        point = even_points[step]
        next_point = even_points[step + 1]
        s_array[step, 0] = next_point.x
        s_array[step, 1] = next_point.y
        s_array[step, 2] = np.sqrt((point.x - next_point.x)**2 + (point.y - next_point.y)**2) / 0.3048
    s_array[:, 2] = np.nancumsum(s_array[:, 2])
    return s_array


def plot_xsect_paths():
    lab_arr = [False, False, False, False, False]
    lw = 0.5
    fig, axes = plt.subplots()
    x_min, y_max = convert_figurecoords_to_linecoords(0, 0)  # upper left
    x_max, y_min = convert_figurecoords_to_linecoords(img.shape[1], img.shape[0])  # lower right
    new_extent = (x_min/5280, x_max/5280, y_min, y_max)
    axes.imshow(img, extent=new_extent)
    # axes.axis('off')
    # axes.xaxis.set_visible(False)
    # axes.yaxis.set_visible(False)
    s_array = make_s_array()
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        if group_name in xsects:
            print('plotting ' + group_name.replace('\n', ' ') + '...')
            xsect_list = xsects[group_name]
            for prtcl in group_data.keys():
                if prtcl in xsect_list:
                    pdata = group_data[prtcl]
                    x_ptcl = pdata[:, 0]
                    y_ptcl = pdata[:, 1]
                    z_ptcl = pdata[:, 2]
                    s_ptcl = np.zeros(len(pdata))
                    for p in range(len(pdata)):
                        s_ptcl[p] = convert_utms_to_linecoords(x_ptcl[p], y_ptcl[p], s_array=s_array)
                    # plot_s, plot_z = convert_linecoords_to_figurecoords(s_ptcl, z_ptcl / 0.3048)
                    plot_s = s_ptcl / 5280
                    plot_z = z_ptcl / 0.3048
                    if not lab_arr[group_ind]:
                        axes.plot(plot_s, plot_z, color=color_list[group_ind], linewidth=lw, label=group_name)
                        lab_arr[group_ind] = True
                    else:
                        axes.plot(plot_s, plot_z, color=color_list[group_ind], linewidth=lw)
    # axes.scatter(330, 190, color='k', s=.1, linewidths=None)
    # axes.scatter(330, 4310, color='k', s=.1, linewidths=None)
    # axes.scatter(10150, 4400, color='k', s=.1, linewidths=None)
    #axes.set_ylim([-2500, 200])
    axes.set_ylim([-3000, 200])
    axes.set_xlim([0.1, 13])
    axes.set_ylabel('Elevation (ft)', fontsize=8)
    axes.set_xlabel('Distance along Cross-Section (miles)', fontsize=8)
    axes.set_aspect('auto')
    axes.tick_params(labelsize=8)
    fig.legend(fontsize=5, ncol=4, loc=8)
    fig.suptitle("Particle Paths along Cross-Section A-A'", fontsize=10)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.15)
    fig.savefig(os.path.join(figure_path, 'xsect_paths.png'), dpi=500)


xsect_coords = get_xsect_coords()
fname = os.path.join('..', 'basemaps', 'base_mickey.jpg')
img = plt.imread(os.path.join(fname))
dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
runs = ['EXT_1984_2004']
porosities = ['01']
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
dis = md.dis
idomain = dis.idomain.array
active = np.amax(idomain, axis=0)
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        pl_data, p_groups = build_data_structure(sim_path)
        xsects = identify_xsect_paths()
        plot_xsect_paths()
