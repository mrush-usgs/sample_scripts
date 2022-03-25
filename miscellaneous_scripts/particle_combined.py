import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from flopy.mf6 import MFSimulation
import flopy.utils.binaryfile as bf
import datetime
import pandas as pd
from osgeo import gdal, osr


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


def plot_basin_boundaries(axes_, lw=0.5):
    axes_.plot(x_basin, y_basin, color='w', linewidth=lw)
    return axes_


def particle_cells():
    locations = pd.read_csv(os.path.join('..', 'data', 'new_locations.csv'), header=0)
    particles_ = []
    sites_ = np.zeros([nrow, ncol], dtype=int)
    for num, line in locations.iterrows():
        particles_.append((line['ROW']-1, line['COLUMN']-1))
        if line['SITE'] == 'qsgrsg':
            sites_[line['ROW'] - 1, line['COLUMN'] - 1] = 1
        elif line['SITE'] == 'qrhsg':
            sites_[line['ROW'] - 1, line['COLUMN'] - 1] = 2
        elif line['SITE'] == 'qwnd':
            sites_[line['ROW'] - 1, line['COLUMN'] - 1] = 3
        elif 'ghbwht' in line['SITE']:
            sites_[line['ROW'] - 1, line['COLUMN'] - 1] = 4
    return particles_, sites_


def ReprojectCoords(coords,src_srs,tgt_srs):
    """ Reproject a list of x,y coordinates. """
    trans_coords=[]
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords


def row_col_to_utm(row_, col_, centered=False):
    if centered:
        conv_row = dYCoordinate - row_ * dCellSize - dCellSize / 2
        conv_col = dXCoordinate + col_ * dCellSize + dCellSize / 2
    else:
        conv_row = dYCoordinate - row_ * dCellSize
        conv_col = dXCoordinate + col_ * dCellSize
    return conv_row, conv_col


def lat_lon_labels(axes_, crop_edges=True, decimals=1):
    # a couple unrelated axes adjustments
    if crop_edges:
        orig_ylim = axes_.get_ylim()
        new_ymin = orig_ylim[0] + 0.56 * (orig_ylim[1] - orig_ylim[0])
        new_ymax = orig_ylim[1] - 0.12 * (orig_ylim[1] - orig_ylim[0])
        axes_.set_ylim([new_ymin, new_ymax])
        orig_xlim = axes_.get_xlim()
        new_xmin = orig_xlim[0] + 0.1 * (orig_xlim[1] - orig_xlim[0])
        new_xmax = orig_xlim[1] - 0.2 * (orig_xlim[1] - orig_xlim[0])
        axes_.set_xlim([new_xmin, new_xmax])
    axes_.set_aspect('equal')
    # now onto the lat lon business
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
    xlons = [-np.round(lon[0], decimals) for lon in new_xcoords]
    ylats = [np.round(lat[1], decimals) for lat in new_ycoords]
    axes_.set_xticklabels(xlons)
    axes_.set_yticklabels(ylats)
    axes_.set_ylabel('Latitude ($^\circ$N)', fontsize=15)
    axes_.set_xlabel('Longitude ($^\circ$W)', fontsize=15)
    axes_.tick_params(labelsize=15)
    return axes_


def particle_map(axes_):
    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    name_list = ['San Gabriel River\nSpreading Grounds',
                 'Rio Hondo\nSpreading Grounds',
                 'Whittier Narrows\nReservoir',
                 'Whittier Narrows\nGeneral Head\nBoundary']

    im_base = axes_[0].imshow(img, extent=extent, origin='upper')
    label_ind = 0
    for particle in particles:
        row, col = row_col_to_utm(particle[0], particle[1], centered=True)
        color = color_list[sites[particle[0], particle[1]] - 1]
        axes_[0].scatter(col, row, color=color, s=50, marker='s')
        plot_name = name_list[sites[particle[0], particle[1]] - 1]
        if label_ind == (sites[particle[0], particle[1]] - 1):
            if (label_ind == 0) or (label_ind == 3):
                ha = 'left'
                adj = 1e3
            else:
                ha= 'right'
                adj = -1e3
            axes_[0].text(x=col+adj, y=row, s=plot_name, color=color, fontsize=15, backgroundcolor='k', ha=ha, ma='center', fontfamily='serif')
            label_ind += 1
    axes_[0] = plot_basin_boundaries(axes_[0], lw=2)
    axes_[0].set_xlim([extent[0], extent[1]])
    axes_[0].set_ylim([extent[2], extent[3] - 0.05 * (extent[3] - extent[2])])
    axes_[0] = lat_lon_labels(axes_[0], crop_edges=False, decimals=2)
    axes_[0].set_title('Model Cells Representing Spreading Grounds\nand Whittier Narrows Underflow', fontsize=20)
    axes_[0].set_aspect('equal')
    return axes_


def particle_configuration(axes_):
    xrange = np.arange(0.1, 1, 0.2) * 660
    yrange = np.arange(0.1, 1, 0.2) * 660
    for row in xrange:
        for col in yrange:
            axes_[1].scatter(row, col, color='k', s=80)
    axes_[1].set_xlim([0, 660])
    axes_[1].set_ylim([0, 660])
    axes_[1].tick_params(labelsize=15)
    axes_[1].set_ylabel('Cell Length (Feet)', fontsize=15)
    axes_[1].set_xlabel('Cell Width (Feet)', fontsize=15)
    axes_[1].set_title('Initial Particle Configuration\n(Plan View)', fontsize=20)
    axes_[1].set_aspect('equal')
    return axes_


# basemap:
# numbers to map to UTM coordinates
dCellSize = 201.15818348
# dXCoordinate = 355047.067638 + (dCellSize / 2)
# dYCoordinate = 3776931.47329
# THIS IS THE UPPER LEFT CORNER
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nlay = 12
nrow = 256
ncol = 312
fname = os.path.join(os.path.join('basemaps', 'basemap_local.tif'))
img = plt.imread(os.path.join(fname))
# gdal
ds = gdal.Open(fname)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])

x_basin, y_basin = zone_boundary('basin')
x_subbasin, y_subbasin = zone_boundary('subbasin')

particles, sites = particle_cells()

fig, axes = plt.subplots(1, 2, figsize=(18, 9))
axes = particle_map(axes)
axes = particle_configuration(axes)
fig.tight_layout()
# fig.subplots_adjust(right=0.9)
# fig.subplots_adjust(left=0.1)
fig.savefig(os.path.join('..', 'figures', 'particles_combined.png'), dpi=300)