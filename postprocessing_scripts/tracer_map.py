import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import matplotlib.animation as manimation
import matplotlib.ticker as plticker


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


def zone_boundary(zone_):
    basin_dict, subbasin_dict, basin_array, subbasin_array = build_basin_dict()
    if zone_ == 'basin':
        zone_array = basin_array
    elif zone_ == 'subbasin':
        zone_array = subbasin_array
    elif zone_ == 'model':
        active_mod = np.copy(active)
        active_mod[33, 99] = 1
        zone_array = active_mod
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
    axes_.plot(x_basin, y_basin, color='w', linewidth=lw, zorder=2)
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


def plot_subbasin_boundaries(axes_, fs=5, lw=0.5):
    axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw, zorder=1)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    axes_.text(x=0.45*xrange + xlim[0], y=0.76*yrange + ylim[0], s='Los\nAngeles\nForebay', ha='center', fontsize=fs, color='k',  fontfamily='serif')
    axes_.text(x=0.47*xrange + xlim[0], y=0.48*yrange + ylim[0], s='Pressure\nArea', ha='center', fontsize=fs, color='k',  fontfamily='serif')
    axes_.text(x=0.74*xrange + xlim[0], y=0.58*yrange + ylim[0], s='Whittier\nArea', ha='center', fontsize=fs, color='k',  fontfamily='serif')
    axes_.text(x=0.63*xrange + xlim[0], y=0.76*yrange + ylim[0], s='Montebello\nForebay', ha='center', fontsize=fs, color='k',  fontfamily='serif')
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
    axes_.set_xticklabels(xlons[1:-1])
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
    ts_file = open(os.path.join(path_, 'lab_model.ts'))
    ts_list = []
    end_header = False
    for ind, line in enumerate(ts_file):
        line_split = line.split()
        if ind == 1:
            x_offset = float(line_split[2])
            y_offset = float(line_split[3]) - nrow * dCellSize + dCellSize/2
        if end_header:
            new_line = [float(line_split[2]),
                        int(line_split[4]),
                        int(line_split[5]),
                        0.3048*float(line_split[10]) + x_offset,
                        0.3048*float(line_split[11]) + y_offset,
                        int(line_split[13])]
            ts_list.append(new_line)
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    column_list = ['Time', 'Particle Group', 'Particle ID', 'X', 'Y', 'Layer']
    ts_df = pd.DataFrame(ts_list, columns=column_list)
    return ts_df


def plot_spreading_grounds(porosity_):
    sp_list = [177, 183, 189]  # 1-indexed: end of Q1 2015, March 2016, September 2016
    sp_time_list = ['Spring 2015', 'Spring 2016', 'Fall 2016']
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ind, sp in enumerate(sp_list):
        axes[ind].imshow(img, extent=extent, origin='upper')
        axes[ind] = plot_subbasin_boundaries(axes[ind], fs=10, lw=1)
        axes[ind] = plot_basin_boundaries(axes[ind], fs=10, lw=1)
        axes[ind] = crop_edges(axes[ind], .5, .22, .5, .22)
        axes[ind] = lat_lon_labels(axes[ind], fs=10)
        if sp <= 180:
            elapsed = sp * 91.25
        else:
            elapsed = 180 * 91.25 + (sp - 180) * 91.25 / 3
        plotted_leg = False
        for group in range(2):
            time_group_subset = (elapsed == ts_data['Time']) & (ts_data['Particle Group'] == int(group + 1)) & \
                                (extent[0] <= ts_data['X']) & (ts_data['X'] <= extent[1]) & (extent[2] <= ts_data['Y']) & (ts_data['Y'] <= extent[3])
            for lay in range(12):
                lay_subset = time_group_subset & (ts_data['Layer'] == (lay+1))
                x_series = ts_data['X'].loc[lay_subset]
                y_series = ts_data['Y'].loc[lay_subset]
                if (ind == 0) and (group == 0):
                    axes[ind].scatter(x_series, y_series, color=scotts_colors[lay], label=layer_list[lay], s=1, linewidths=0)
                else:
                    axes[ind].scatter(x_series, y_series, color=scotts_colors[lay], s=1, linewidths=0)
        bg_rowcol = (94, 193)
        bg_y, bg_x = row_col_to_utm(bg_rowcol[0], bg_rowcol[1])
        axes[ind].scatter(bg_x, bg_y, color='w', s=50)
        axes[ind].annotate(text='Bell\nGardens\n1', xy=(bg_x+500, bg_y), color='w', va='center', fontfamily='serif')
        if not plotted_leg:
            fig.legend(loc=8, fontsize=10, markerscale=5, ncol=6)
            plotted_leg = True
        axes[ind].set_title(sp_time_list[ind], fontsize=10)
    fig.suptitle('Locations of Particles Representing Acesulfame Tracer')
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(os.path.join(figure_path, str(porosity_) + '_tracer_map.png'), dpi=500)


def plot_spreading_grounds_rotate(porosity_):
    sp_list = [177, 183, 189]  # 1-indexed: end of Q1 2015, March 2016, September 2016
    sp_time_list = ['Spring 2015', 'Spring 2016', 'Fall 2016']
    fig, axes = plt.subplots(3, 1, figsize=(8, 15))
    for ind, sp in enumerate(sp_list):
        axes[ind].imshow(img, extent=extent, origin='upper')
        axes[ind] = plot_subbasin_boundaries(axes[ind], fs=10, lw=1)
        axes[ind] = plot_basin_boundaries(axes[ind], fs=10, lw=1)
        axes[ind] = crop_edges(axes[ind], .5, .22, .5, .22)
        axes[ind] = lat_lon_labels(axes[ind], fs=10)
        if sp <= 180:
            elapsed = sp * 91.25
        else:
            elapsed = 180 * 91.25 + (sp - 180) * 91.25 / 3
        plotted_leg = False
        for group in range(2):
            time_group_subset = (elapsed == ts_data['Time']) & (ts_data['Particle Group'] == int(group + 1)) & \
                                (extent[0] <= ts_data['X']) & (ts_data['X'] <= extent[1]) & (extent[2] <= ts_data['Y']) & (ts_data['Y'] <= extent[3])
            for lay in range(12):
                lay_subset = time_group_subset & (ts_data['Layer'] == (lay+1))
                x_series = ts_data['X'].loc[lay_subset]
                y_series = ts_data['Y'].loc[lay_subset]
                if (ind == 0) and (group == 0):
                    axes[ind].scatter(x_series, y_series, color=scotts_colors[lay], label=layer_list[lay], s=1, linewidths=0)
                else:
                    axes[ind].scatter(x_series, y_series, color=scotts_colors[lay], s=1, linewidths=0)
        bg_rowcol = (94, 193)
        bg_y, bg_x = row_col_to_utm(bg_rowcol[0], bg_rowcol[1])
        axes[ind].scatter(bg_x, bg_y, color='w', s=50)
        axes[ind].annotate(text='Bell\nGardens\n1', xy=(bg_x+500, bg_y), color='w', va='center', fontfamily='serif')
        if not plotted_leg:
            fig.legend(loc=5, fontsize=9, markerscale=5)
            plotted_leg = True
        axes[ind].set_title(sp_time_list[ind], fontsize=10)
    fig.suptitle('Locations of Particles Representing Acesulfame Tracer')
    fig.tight_layout()
    fig.subplots_adjust(right=0.9)
    fig.savefig(os.path.join(figure_path, str(porosity_) + '_tracer_map_sideways.png'), dpi=500)


dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
x, y = np.meshgrid(np.linspace(dXCoordinate, dXCoordinate + ncol * dCellSize, ncol + 1), np.linspace(dYCoordinate, dYCoordinate - nrow * dCellSize, nrow + 1))
fname = os.path.join('..', 'basemaps', 'tracer.tif')
img = plt.imread(os.path.join(fname))
# gdal
ds = gdal.Open(fname)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
particle_groups = {1: 'Rio Hondo\nSpreading Grounds',
                   2: 'San Gabriel River\nSpreading Grounds',
                   3: 'Whittier Narrows\nReservoir',
                   4: 'Whittier Narrows\nUnderflow (Shallow)',
                   5: 'Whittier Narrows\nUnderflow (Deep)'}
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper\nWilmington A', '8: Upper\nWilmington B', '9: Lower\nWilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C\nand BC']
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
run = 'BASE'
porosities = ['01']
#porosities = ['03']
any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
dis = md.dis
idomain = dis.idomain.array
active = np.amax(idomain, axis=0)
x_model, y_model = zone_boundary('model')
x_basin, y_basin = zone_boundary('basin')
x_subbasin, y_subbasin = zone_boundary('subbasin')
for porosity in porosities:
    figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
    sim_path = os.path.join('MODPATH', 'TRACER_SHORT3', 'POROSITY_' + porosity, run, 'lab_model')
    ts_data = build_data_structure(sim_path)
    plot_spreading_grounds(porosity)
    plot_spreading_grounds_rotate(porosity)
