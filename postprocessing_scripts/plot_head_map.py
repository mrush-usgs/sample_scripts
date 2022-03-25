import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import pandas as pd
import matplotlib.ticker as plticker
from matplotlib import markers


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


def zones():
    zones = pd.read_csv(os.path.join('..', 'data', 'LA_BUDGET_ZONES_UPDATE.csv'), header=0)
    zone_array_ = np.empty([256, 312])
    for num, line in zones.iterrows():
        zone_array_[line['RowVal'] - 1, line['ColVal_'] - 1] = line['BudZone']
    return zone_array_


def get_all_wells():
    connection_data = md.maw.connectiondata.get_data()
    well_lrc_list = []
    for line in connection_data:
        well_lrc_list.append(line['cellid'])
    return well_lrc_list


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


def build_basin_dict():
    basins = pd.read_csv(os.path.join('..', 'data', 'subbasins.csv'))
    basin_names = ['West Coast', 'Central', 'Santa Monica', 'Hollywood', 'Orange', 'Ocean']
    subbasin_names = ['Pressure Area', 'Los Angeles Forebay', 'Whittier Area', 'Montebello Forebay']
    name_dict = {}
    sb_name_dict = {}
    num_dict = {}
    sb_num_dict = {}
    for num, name in enumerate(basin_names):
        num_dict[name] = num+1
    for num, name in enumerate(subbasin_names):
        sb_num_dict[name] = num+1
    basin_array_ = np.zeros([256, 312], dtype=int)
    subbasin_array_ = np.zeros([256, 312], dtype=int)
    name_dict['CENTRAL BASIN, PRESSURE AREA, WRD'] = 'Central'
    name_dict['CENTRAL BASIN PRESSURE AREA, NON-WRD'] = 'Central'
    name_dict['LOS ANGELES FOREBAY, NON-PRESSURE AREA, NON-WRD'] = 'Central'
    name_dict['LOS ANGELES FOREBAY, NON-PRESSURE AREA, WRD'] = 'Central'
    name_dict['MONTEBELLO FOREBAY, NON-PRESSURE AREA'] = 'Central'
    name_dict['WHITTIER AREA'] = 'Central'
    name_dict['ORANGE COUNTY CENTRAL'] = 'Orange'
    name_dict['ORANGE COUNTY NORTH'] = 'Orange'
    name_dict['ORANGE COUNTY SOUTH'] = 'Orange'
    name_dict['WEST COAST'] = 'West Coast'
    name_dict['SANTA MONICA'] = 'Santa Monica'
    name_dict['HOLLYWOOD'] = 'Hollywood'
    name_dict['OCEAN'] = 'Ocean'
    sb_name_dict['CENTRAL BASIN, PRESSURE AREA, WRD'] = 'Pressure Area'
    sb_name_dict['CENTRAL BASIN PRESSURE AREA, NON-WRD'] = 'Pressure Area'
    sb_name_dict['LOS ANGELES FOREBAY, NON-PRESSURE AREA, NON-WRD'] = 'Los Angeles Forebay'
    sb_name_dict['LOS ANGELES FOREBAY, NON-PRESSURE AREA, WRD'] = 'Los Angeles Forebay'
    sb_name_dict['WHITTIER AREA'] = 'Whittier Area'
    sb_name_dict['MONTEBELLO FOREBAY, NON-PRESSURE AREA'] = 'Montebello Forebay'
    for dummy, line in basins.iterrows():
        if line['subbasin'] in name_dict.keys():
            row = line['Row'] - 1
            col = line['Column_'] - 1
            basin_array_[row, col] = num_dict[name_dict[line['subbasin']]]
        if line['subbasin'] in sb_name_dict.keys():
            row = line['Row'] - 1
            col = line['Column_'] - 1
            subbasin_array_[row, col] = sb_num_dict[sb_name_dict[line['subbasin']]]
    return basin_array_, subbasin_array_


def zone_boundary(zone_):
    basin_array, subbasin_array = build_basin_dict()
    if zone_ == 'basin':
        zone_array = basin_array
    elif zone_ == 'subbasin':
        zone_array = subbasin_array
    ver_seg = np.where(zone_array[:, 1:] != zone_array[:, :-1])
    hor_seg = np.where(zone_array[1:, :] != zone_array[:-1, :])
    l = []
    for p in zip(*hor_seg):
        l.append((p[1], p[0] + 1))
        l.append((p[1] + 1, p[0] + 1))
        l.append((np.nan, np.nan))
    for p in zip(*ver_seg):
        l.append((p[1] + 1, p[0]))
        l.append((p[1] + 1, p[0] + 1))
        l.append((np.nan, np.nan))
    segments = np.array(l)
    bound_xcoords = segments[:, 0]
    bound_ycoords = segments[:, 1]
    bound_y_utms, bound_x_utms = row_col_to_utm(bound_ycoords, bound_xcoords)
    return bound_x_utms, bound_y_utms


def plot_basin_boundaries(axes_, fs=5, legend=True, lw=0.5, label=False):
    axes_.plot(x_basin, y_basin, color='w', linewidth=lw, zorder=14)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    if label:
        axes_.text(x=0.25*xrange + xlim[0], y=0.75*yrange + ylim[0], s='Santa\nMonica', ha='center', fontsize=fs, color='w')
        if legend:
            axes_.text(x=0.47*xrange + xlim[0], y=0.32*yrange + ylim[0], s='West\nCoast', ha='center', fontsize=fs, color='w')
        else:
            axes_.text(x=0.4*xrange + xlim[0], y=0.4*yrange + ylim[0], s='West\nCoast', ha='center', fontsize=fs, color='w')
        axes_.text(x=0.6*xrange + xlim[0], y=0.38*yrange + ylim[0], s='Central', ha='center', fontsize=fs, color='w')
        axes_.text(x=0.74*xrange + xlim[0], y=0.35*yrange + ylim[0], s='Orange\nCounty', ha='center', fontsize=fs, color='w')
        if fs<8:
            axes_.text(x=0.35*xrange + xlim[0], y=0.85*yrange + ylim[0], s='Hollywood', ha='center', fontsize=fs, rotation=10, color='w')
        else:
            axes_.text(x=0.35*xrange + xlim[0], y=0.86*yrange + ylim[0], s='Hollywood', ha='center', fontsize=fs, rotation=10, color='w')
    return axes_


def plot_subbasin_boundaries(axes_, lw=0.5):
    axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw, zorder=13)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    # axes_.text(x=0.5*xrange + xlim[0], y=0.73*yrange + ylim[0], s='Los Angeles\nForebay', ha='center', fontsize=fs, color='w', zorder=14)
    # axes_.text(x=0.6*xrange + xlim[0], y=0.42*yrange + ylim[0], s='Pressure\nArea', ha='center', fontsize=fs, color='w', zorder=14)
    # axes_.text(x=0.73*xrange + xlim[0], y=0.58*yrange + ylim[0], s='Whittier\nArea', ha='center', fontsize=fs, color='w', zorder=14)
    # axes_.text(x=0.62*xrange + xlim[0], y=0.75*yrange + ylim[0], s='Montebello\nForebay', ha='center', fontsize=fs, color='w', zorder=14)
    return axes_


def ReprojectCoords(coords,src_srs,tgt_srs):
    """ Reproject a list of x,y coordinates. """
    trans_coords=[]
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords


def crop_edges(axes_, xmin, xmax, ymin, ymax, x_spacing=5000, y_spacing=5000):
    orig_ylim = axes_.get_ylim()
    new_ymin = orig_ylim[0] + ymin * (orig_ylim[1] - orig_ylim[0])
    new_ymax = orig_ylim[1] - ymax * (orig_ylim[1] - orig_ylim[0])
    axes_.set_ylim([new_ymin, new_ymax])
    orig_xlim = axes_.get_xlim()
    new_xmin = orig_xlim[0] + xmin * (orig_xlim[1] - orig_xlim[0])
    new_xmax = orig_xlim[1] - xmax * (orig_xlim[1] - orig_xlim[0])
    axes_.set_xlim([new_xmin, new_xmax])
    axes_.set_aspect('equal')
    x_loc = plticker.MultipleLocator(base=x_spacing)
    y_loc = plticker.MultipleLocator(base=y_spacing)
    axes_.xaxis.set_major_locator(x_loc)
    axes_.yaxis.set_major_locator(y_loc)
    return axes_


def lat_lon_labels(axes_, fs=10, xlabel=True, ylabel=True):
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds.GetProjection())
    tgt_srs = src_srs.CloneGeogCS()
    xlabels = axes_.get_xticks()
    ylabels = axes_.get_yticks()
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xcoords = [(xlab, ylim[0]) for xlab in xlabels]
    ycoords = [(xlim[0], ylab) for ylab in ylabels]
    new_xcoords = ReprojectCoords(xcoords, src_srs, tgt_srs)
    new_ycoords = ReprojectCoords(ycoords, src_srs, tgt_srs)
    xlons = [-np.round(lon[0], 2) for lon in new_xcoords]
    ylats = [np.round(lat[1], 2) for lat in new_ycoords]
    axes_.set_xticks(axes_.get_xticks()[1:-1])
    axes_.set_yticks(axes_.get_yticks()[1:-1])
    axes_.set_xticklabels(xlons[1:-1], rotation=-45, ha='left')
    axes_.set_yticklabels(ylats[1:-1])
    if ylabel:
        axes_.set_ylabel('Latitude ($^\circ$N)', fontsize=fs)
    if xlabel:
        axes_.set_xlabel('Longitude ($^\circ$W)', fontsize=fs)
    axes_.tick_params(labelsize=fs)
    return axes_


def row_col_to_utm(row_, col_, centered=False):
    if centered:
        conv_row = dYCoordinate - row_ * dCellSize - dCellSize / 2
        conv_col = dXCoordinate + col_ * dCellSize + dCellSize / 2
    else:
        conv_row = dYCoordinate - row_ * dCellSize
        conv_col = dXCoordinate + col_ * dCellSize
    return conv_row, conv_col


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
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 2] = float(line_split[8])  # layer yo
                sub_ind += 1
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    return pl_dict, particle_groups


def plot_heads():
    fig, axes = plt.subplots(3, 4, figsize=(5.4, 5.8))
    axes = axes.flat
    model_head_path = os.path.join('EXT_1984_2004', 'lab_model')
    head_sim = flopy.mf6.MFSimulation.load(sim_ws=model_head_path, verbosity_level=0)
    head_file = head_sim.simulation_data.mfdata['lab_model', 'HDS', 'HEAD']
    head_sp_array = np.zeros([20, 12, 256, 312])
    basin_array, subbasin_array = build_basin_dict()
    for sp in range(20):
        head_sp = head_file[sp]
        head_sp[head_sp == 1e30] = np.nan
        head_sp_array[sp] = head_sp
    head_array = np.nanmean(head_sp_array, axis=0)
    label_list = []
    for lay in range(12):
        axes[lay].imshow(img, extent=extent, origin='upper')
        plot_array = np.copy(head_array[lay])
        plot_array[subbasin_array != 4] = np.nan
        im = axes[lay].pcolormesh(x, y, plot_array, cmap='viridis', vmin=0, vmax=150)
        u = head_array[lay, :, 1:] - head_array[lay, :, :-1]
        v = head_array[lay, 1:, :] - head_array[lay, :-1, :]
        # axes[lay].streamplot(x[1:-1, 1:-1], y[1:-1, 1:-1], -u[:-1, :], v[:, :-1], density=5, linewidth=0.5, arrowsize=0.3, color='k')
        # axes[lay].quiver(x[1:-1, 1:-1], y[1:-1, 1:-1], -u[:-1, :], v[:, :-1], color='k')
        # levs = np.arange(round(np.nanmin(head_array[lay]) / 10) * 10,
        #                  round(np.nanmax(head_array[lay]) / 10) * 10, 10)
        # cs = axes[lay].contour(x_center, y_center, head_array[lay], levels=levs, colors='k', linewidths=0.5)
        # axes[lay].clabel(cs, cs.levels, inline=True, fontsize=3, fmt='%1.0f')
        print('\n')
        print(lay+1)
        print(np.nanmin(head_array[lay]))
        print(np.nanmax(head_array[lay]))
        for name_ind in name_dict.keys():
            data = particles[name_dict[name_ind]]
            for particle in data:
                row, col = row_col_to_utm(particle[1] - 1, particle[2] - 1, centered=True)
                color = color_list[name_ind]
                plot_name = name_list[name_ind]
                if particle[0] == (lay + 1):
                    if plot_name not in label_list:
                        axes[lay].scatter(col, row, color=color, s=0.8, marker='s', zorder=20, linewidths=0.1, label=plot_name, edgecolors='k')
                        label_list.append(plot_name)
                    else:
                        axes[lay].scatter(col, row, color=color, s=0.8, marker='s', zorder=20, linewidths=0.1, edgecolors='k')
                else:
                    axes[lay].scatter(col, row, color=color, s=0.8, marker='s', zorder=20, linewidths=0.1, edgecolors='k', facecolors='none')
        axes[lay] = plot_subbasin_boundaries(axes[lay])
        axes[lay] = plot_basin_boundaries(axes[lay], fs=5)
        axes[lay] = crop_edges(axes[lay], 0.58, 0.30, 0.55, 0.18, x_spacing=1500, y_spacing=1500)
        if np.remainder(lay, 4) == 0:
            axes[lay] = lat_lon_labels(axes[lay], fs=4, xlabel=False)
        if lay > 7:
            axes[lay] = lat_lon_labels(axes[lay], fs=4, ylabel=False)
        else:
            axes[lay] = lat_lon_labels(axes[lay], fs=4, xlabel=False, ylabel=False)
        axes[lay].set_title(layer_list[lay].replace('\n', ' '), fontsize=6)
        if np.remainder(lay, 4) != 0:
            axes[lay].set_yticks([])
            axes[lay].yaxis.set_ticklabels([])
        if lay < 8:
            axes[lay].set_xticks([])
            axes[lay].xaxis.set_ticklabels([])
    fig.suptitle('Montebello Forebay Hydraulic Heads: 1971-1975', fontsize=8)
    fig.tight_layout(pad=0.01, h_pad=0.01, w_pad=0.01)
    cbar_ax = fig.add_axes([0.9, 0.4, 0.03, 0.5])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=5)
    cbar_ax.set_title('Hydraulic\nHead\n(ft)\n', fontsize=6)
    fig.subplots_adjust(right=0.85)
    fig.subplots_adjust(top=0.92)
    fig.legend(loc=(0.85, 0.12), fontsize=4, title='Simulated\nGroundwater\nSource\nLocations:', title_fontsize=6, markerscale=2)
    fig.savefig(os.path.join(figure_path, 'heads.png'), dpi=700)

# model_grid = flopy.utils.MfGrdFile('lab_model.dis.grb')
# print(model_grid.get_modelgrid())

dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
x, y = np.meshgrid(np.linspace(dXCoordinate, dXCoordinate + ncol * dCellSize, ncol + 1), np.linspace(dYCoordinate, dYCoordinate - nrow * dCellSize, nrow + 1))
# x_center, y_center = np.meshgrid(np.linspace(dXCoordinate + dCellSize/2, dXCoordinate + ncol * dCellSize, ncol), np.linspace(dYCoordinate - dCellSize/2, dYCoordinate - nrow * dCellSize, nrow))
fname = os.path.join('..', 'basemaps', 'basemap.tif')
img = plt.imread(os.path.join(fname))
zone_array = zones()
# gdal
ds = gdal.Open(fname)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])
color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
name_list = ['San Gabriel River\nSpreading Grounds',
             'Rio Hondo\nSpreading Grounds',
             'Whittier Narrows\nReservoir',
             'Whittier Narrows\nUnderflow\n(Shallow)',
             'Whittier Narrows\nUnderflow\n(Deep)']
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper\nWilmington A', '8: Upper\nWilmington B', '9: Lower\nWilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C\nand BC']
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
name_dict = {0: 'qsgrsg',
             1: 'qrhsg',
             2: 'qwnd',
             3: 'ghbwht_dom',
             4: 'ghbwht_lbc'}
runs = ['EXT_1984_2004']
porosities = ['01']
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
dis = md.dis
idomain = dis.idomain.array
active = np.amax(idomain, axis=0)
x_basin, y_basin = zone_boundary('basin')
x_subbasin, y_subbasin = zone_boundary('subbasin')
particles = particle_cells()
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        pl_data, p_groups = build_data_structure(sim_path)
        ep_data = build_endpoints(sim_path)
        plot_heads()