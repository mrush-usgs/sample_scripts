import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import pandas as pd
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
    zone_array = np.zeros([258, 314])  # extra cell all aroond
    if zone_ == 'basin':
        zone_array[1:-1, 1:-1] = basin_array
    elif zone_ == 'subbasin':
        zone_array[1:-1, 1:-1] = subbasin_array
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
    bound_xcoords = segments[:, 0] - 1
    bound_ycoords = segments[:, 1] - 1
    bound_y_utms, bound_x_utms = row_col_to_utm(bound_ycoords, bound_xcoords)
    return bound_x_utms, bound_y_utms


def plot_basin_boundaries(axes_, fs=5, legend=True, lw=0.5, label=False, legend_label=True):
    if not legend_label:
        axes_.plot(x_basin, y_basin, color='w', linewidth=lw, zorder=14, label='Groundwater\nBasin\nBoundaries')
    else:
        axes_.plot(x_basin, y_basin, color='w', linewidth=lw, zorder=14)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    if label:
        axes_.text(x=0.19*xrange + xlim[0], y=0.60*yrange + ylim[0], s='Santa\nMonica\nBay', ha='center', fontsize=fs, color='w', zorder=14, fontfamily='serif')
        axes_.text(x=0.25*xrange + xlim[0], y=0.75*yrange + ylim[0], s='Santa\nMonica\nBasin', ha='center', fontsize=fs, color='w', zorder=14, fontfamily='serif')
        if legend:
            axes_.text(x=0.47*xrange + xlim[0], y=0.32*yrange + ylim[0], s='West\nCoast\nBasin', ha='center', fontsize=fs, color='w', zorder=14, fontfamily='serif')
        else:
            axes_.text(x=0.35*xrange + xlim[0], y=0.4*yrange + ylim[0], s='West\nCoast\nBasin', ha='center', fontsize=fs, color='w', zorder=14, fontfamily='serif')
        axes_.text(x=0.58*xrange + xlim[0], y=0.13*yrange + ylim[0], s='San\nPedro\nBay', ha='center', fontsize=fs, color='w', zorder=14, fontfamily='serif')
        axes_.text(x=0.5*xrange + xlim[0], y=0.5*yrange + ylim[0], s='Central\nBasin', ha='center', fontsize=fs, color='w', zorder=14, fontfamily='serif')
        axes_.text(x=0.73*xrange + xlim[0], y=0.35*yrange + ylim[0], s='Orange\nCounty', ha='center', fontsize=fs, color='w', zorder=14, fontfamily='serif')
        if fs<8:
            axes_.text(x=0.37*xrange + xlim[0], y=0.855*yrange + ylim[0], s='Hollywood\nBasin', ha='center', fontsize=fs, rotation=10, color='w', zorder=14, fontfamily='serif')
        else:
            axes_.text(x=0.37*xrange + xlim[0], y=0.84*yrange + ylim[0], s='Hollywood\nBasin', ha='center', fontsize=fs, rotation=10, color='w', zorder=14, fontfamily='serif')
    return axes_


def plot_subbasin_boundaries(axes_, fs=5, lw=0.5, legend_label=True, zoomed=False):
    if not legend_label:
        axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw, zorder=13, label='Central Basin\nSub-Area\nBoundaries')
    else:
        axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw, zorder=13)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    if not zoomed:
        axes_.text(x=0.45*xrange + xlim[0], y=0.73*yrange + ylim[0], s='Los Angeles\nForebay\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
        axes_.text(x=0.6*xrange + xlim[0], y=0.48*yrange + ylim[0], s='Pressure\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
        axes_.text(x=0.75*xrange + xlim[0], y=0.58*yrange + ylim[0], s='Whittier\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
        axes_.text(x=0.62*xrange + xlim[0], y=0.725*yrange + ylim[0], s='Montebello\nForebay\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
    else:
        axes_.text(x=0.5*xrange + xlim[0], y=0.73*yrange + ylim[0], s='Los Angeles\nForebay\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
        axes_.text(x=0.6*xrange + xlim[0], y=0.48*yrange + ylim[0], s='Pressure\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
        axes_.text(x=0.74*xrange + xlim[0], y=0.58*yrange + ylim[0], s='Whittier\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
        axes_.text(x=0.62*xrange + xlim[0], y=0.725*yrange + ylim[0], s='Montebello\nForebay\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
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


def lat_lon_labels(axes_, fs=10, ts=10, xlabel=True, ylabel=True):
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
    axes_.tick_params(labelsize=ts)
    return axes_


def row_col_to_utm(row_, col_, centered=False):
    if centered:
        conv_row = dYCoordinate - row_ * dCellSize - dCellSize / 2
        conv_col = dXCoordinate + col_ * dCellSize + dCellSize / 2
    else:
        conv_row = dYCoordinate - row_ * dCellSize
        conv_col = dXCoordinate + col_ * dCellSize
    return conv_row, conv_col


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


def build_data_structure(porosity_, run_):
    pl_file = open(os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity_, run_, 'lab_model', 'lab_model.pl'))
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
                pl_dict[particle_groups[particle_grp]][particle_num][sub_ind, 2] = float(line_split[3]) * 0.3048
                sub_ind += 1
        if line_split[0] == 'END' and line_split[1] == 'HEADER':
            end_header = True
    return pl_dict, particle_groups


def plot_lateral_sensitivity():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    zorders = [2, 4, 5, 1, 3]
    for run_ind, run_name in enumerate(run_dict.keys()):
        if run_ind == 4:
            plot_ind = 3
        else:
            plot_ind = run_ind
        run_data, ep_data_ = run_dict[run_name]
        axes[plot_ind].imshow(img, extent=extent, origin='upper')
        for group_ind, group_name in enumerate(run_data.keys()):
            group_data = run_data[group_name]
            for prtcl in group_data.keys():
                pdata = group_data[prtcl]
                if (prtcl == 1) and (plot_ind == 0):
                    axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, label=group_name, zorder=zorders[group_ind])
                else:
                    axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, zorder=zorders[group_ind])
                ep_subset = (ep_data_['Particle Group'] == (group_ind + 1)) & (ep_data_['Particle ID'] == prtcl)
                status = ep_data_['Status'].loc[ep_subset].item()
        axes[plot_ind] = plot_basin_boundaries(axes[run_ind], fs=6, label=False, legend_label=basin_label)
        axes[plot_ind] = plot_subbasin_boundaries(axes[run_ind], fs=4, legend_label=subbasin_label, zoomed=True)
        basin_label = True
        subbasin_label = True
        axes[plot_ind] = crop_edges(axes[run_ind], .43, .23, .33, .17)
        fs = 6
        axes[plot_ind] = lat_lon_labels(axes[run_ind], fs=fs, ts=5)
        axes[plot_ind].set_title(run_titles[run_ind], fontsize=7)
    fig.legend(loc=7, fontsize=7, facecolor='tab:gray', framealpha=0.2)
    fig.suptitle('Sensitivity of Lateral Flow to Porosity and Model Extension Approach', fontsize=10)
    fig.tight_layout(pad=0.5, h_pad=0.1, w_pad=0.1)
    fig.subplots_adjust(right=0.82)
    fig.subplots_adjust(left=0.1)
    fig.savefig(os.path.join(figure_path, 'combined_lateral_sensitivity.png'), dpi=500)


def plot_vertical_sensitivity():
    lw = 0.5
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    zorders = [2, 4, 5, 1, 3]
    for run_ind, run_name in enumerate(run_dict.keys()):
        if run_ind == 4:
            plot_ind = 3
        else:
            plot_ind = run_ind
        run_data, ep_data_ = run_dict[run_name]
        for group_ind, group_name in enumerate(run_data.keys()):
            group_data = run_data[group_name]
            max_xy = 0
            max_z = 0
            min_z = 0
            for ind, prtcl in enumerate(group_data.keys()):
                pdata = group_data[prtcl]
                xy = np.sqrt((pdata[:, 0] - pdata[0, 0])**2 + (pdata[:, 1]-pdata[0, 1])**2) / 0.3048 / 5280
                z = pdata[:, 2] / 0.3048
                if (ind == 0) and (run_ind == 0):
                    axes[plot_ind].plot(xy, z, color=color_list[group_ind], linewidth=lw, zorder=zorders[group_ind], label=group_name)
                else:
                    axes[plot_ind].plot(xy, z, color=color_list[group_ind], linewidth=lw, zorder=zorders[group_ind])
                if np.amax(xy) > max_xy:
                    max_xy = np.amax(xy)
                if np.amin(z) < min_z:
                    min_z = np.amin(z)
                if np.amax(z) > max_z:
                    max_z = np.amax(z)
                ep_subset = (ep_data_['Particle Group'] == (group_ind + 1)) & (ep_data_['Particle ID'] == prtcl)
                status = ep_data_['Status'].loc[ep_subset].item()
                # if status == 5:
                #     axes[plot_ind].plot(xy[-2:], z[-2:], color='k', linewidth=lw, zorder=15)
            print('\n')
            print(run_name)
            print(group_name)
            print(max_xy)
            print(min_z)
            # print(max_z)
        if plot_ind > 1:
            axes[plot_ind].set_xlabel('Horizontal Displacement (miles)', fontsize=8)
        if (run_ind == 0) or (run_ind == 2):
            axes[plot_ind].set_ylabel('Elevation\n(ft)', ha='right', va='center', ma='left', rotation=0, fontsize=8)
        axes[plot_ind].set_title(run_titles[run_ind], fontsize=8)
        axes[plot_ind].set_ylim([-3500, 250])
        axes[plot_ind].set_xlim([0, 14.5])
        axes[plot_ind].tick_params(labelsize=8)
    fig.legend(ncol=5, loc=8, fontsize=6)
    fig.suptitle('Sensitivity of Vertical and Lateral Flow to Porosity and Model Extension Approach', fontsize=10)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(os.path.join(figure_path, 'combined_vertical_sensitivity.png'), dpi=500)


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
porosities = ['005', '03'] # choose 01, 02, 03
runs = ['EXT_1984_2004', 'EXT_2019_1999']
run_titles = ['Porosity = 0.05\nModel Extended with 1984-2004 Data',
              'Porosity = 0.05\nModel Extended with 2019-1999 Data',
              'Porosity = 0.30\nModel Extended with 1984-2004 Data',
              'Porosity = 0.30\nModel Extended with 2019-1999 Data',
              ]
figure_path = os.path.join('..', 'figures')
any_sim_path = os.path.join('..', '2019_MF6_MODEL', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
dis = md.dis
idomain = dis.idomain.array
active = np.amax(idomain, axis=0)
x_basin, y_basin = zone_boundary('basin')
x_subbasin, y_subbasin = zone_boundary('subbasin')

run_dict = {}
for porosity in porosities:
    for run in runs:
        pl_data, p_groups = build_data_structure(porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        ep_data = build_endpoints(sim_path)
        key = porosity + '_' + run
        run_dict[key] = (pl_data, ep_data)

plot_lateral_sensitivity()
plot_vertical_sensitivity()