import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import matplotlib.animation as manimation
import matplotlib.ticker as plticker


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
    axes_.plot(x_basin, y_basin, color='w', linewidth=lw)
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
    axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    axes_.text(x=0.5*xrange + xlim[0], y=0.73*yrange + ylim[0], s='Los Angeles\nForebay', ha='center', fontsize=fs, color='k', fontfamily='serif')
    axes_.text(x=0.6*xrange + xlim[0], y=0.42*yrange + ylim[0], s='Pressure\nArea', ha='center', fontsize=fs, color='k', fontfamily='serif')
    axes_.text(x=0.73*xrange + xlim[0], y=0.58*yrange + ylim[0], s='Whittier\nArea', ha='center', fontsize=fs, color='k', fontfamily='serif')
    axes_.text(x=0.63*xrange + xlim[0], y=0.76*yrange + ylim[0], s='Montebello\nForebay', ha='center', fontsize=fs, color='k', fontfamily='serif')
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
                        0.3048*float(line_split[11]) + y_offset]
            ts_list.append(new_line)
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    column_list = ['Time', 'Particle Group', 'Particle ID', 'X', 'Y']
    ts_df = pd.DataFrame(ts_list, columns=column_list)
    return ts_df


def plot_spreading_grounds(run_, porosity_):
    lw = 0.5
    fig, axes = plt.subplots()
    axes.imshow(img, extent=extent, origin='upper')
    axes = plot_subbasin_boundaries(axes, fs=8, lw=1)
    axes = plot_basin_boundaries(axes, fs=8, lw=1)
    axes = crop_edges(axes, .43, .23, .35, .18)
    axes = lat_lon_labels(axes, fs=8)
    tsteps = np.unique(ts_data['Time'])
    plotted_leg = False
    last_ts = tsteps[0]
    with writer.saving(fig, os.path.join(figure_path, 'spreading_' + run_ + '_' + porosity_ +'.mp4'), dpi=500):
        for ts in tsteps:
            print('processing timestep ' + str(ts))
            for group in range(3):
                time_group_subset = ((ts_data['Time'] == ts) | (ts_data['Time'] == last_ts)) & (ts_data['Particle Group'] == int(group + 1))
                particles = np.unique(ts_data['Particle ID'].loc[time_group_subset])
                for ind, particle in enumerate(particles):
                    particle_subset = time_group_subset & (ts_data['Particle ID'] == particle)
                    x_series = ts_data['X'].loc[particle_subset]
                    y_series = ts_data['Y'].loc[particle_subset]
                    if ind == 0:
                        axes.plot(x_series, y_series, color=color_list[group], linewidth=lw, label=particle_groups[int(group+1)])
                    else:
                        axes.plot(x_series, y_series, color=color_list[group], linewidth=lw)
            if not plotted_leg:
                fig.legend(loc=7, fontsize=8)
                plotted_leg = True
                fig.subplots_adjust(right=0.75)
            fig.suptitle('Montebello Forebay Spreading Grounds Flowpaths: ' + str(int(1971 + np.round(ts / 365, 0))))
            writer.grab_frame()
            last_ts = ts


def plot_underflow(run_, porosity_):
    lw = 0.5
    fig, axes = plt.subplots()
    axes.imshow(img, extent=extent, origin='upper')
    axes = plot_subbasin_boundaries(axes, fs=8, lw=1)
    axes = plot_basin_boundaries(axes, fs=8, lw=1)
    axes = crop_edges(axes, .43, .23, .35, .18)
    axes = lat_lon_labels(axes, fs=8)
    tsteps = np.unique(ts_data['Time'])
    plotted_leg = False
    last_ts = tsteps[0]
    with writer.saving(fig, os.path.join(figure_path, 'underflow_' + run_ + '_' + porosity_ +'.mp4'), dpi=500):
        for ts in tsteps:
            print('processing timestep ' + str(ts))
            for group in range(3, 5):
                time_group_subset = ((ts_data['Time'] == ts) | (ts_data['Time'] == last_ts)) & (ts_data['Particle Group'] == int(group + 1))
                particles = np.unique(ts_data['Particle ID'].loc[time_group_subset])
                for ind, particle in enumerate(particles):
                    particle_subset = time_group_subset & (ts_data['Particle ID'] == particle)
                    x_series = ts_data['X'].loc[particle_subset]
                    y_series = ts_data['Y'].loc[particle_subset]
                    if ind == 0:
                        axes.plot(x_series, y_series, color=color_list[group], linewidth=lw, label=particle_groups[int(group+1)])
                    else:
                        axes.plot(x_series, y_series, color=color_list[group], linewidth=lw)
            if not plotted_leg:
                fig.legend(loc=7, fontsize=8)
                plotted_leg = True
                fig.subplots_adjust(right=0.75)
            fig.suptitle('Whittier Narrows Underflow Flowpaths: ' + str(int(1971 + np.round(ts / 365, 0))))
            writer.grab_frame()
            last_ts = ts


dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
x, y = np.meshgrid(np.linspace(dXCoordinate, dXCoordinate + ncol * dCellSize, ncol + 1), np.linspace(dYCoordinate, dYCoordinate - nrow * dCellSize, nrow + 1))
fname = os.path.join('..', 'basemaps', 'basemap.tif')
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
FFMpegWriter = manimation.writers['ffmpeg']
writer = FFMpegWriter(fps=5)
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        ts_data = build_data_structure(sim_path)
        plot_spreading_grounds(run, porosity)
        plot_underflow(run, porosity)
