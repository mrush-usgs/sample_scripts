import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import pandas as pd
import matplotlib.ticker as plticker
from matplotlib import markers
from matplotlib.lines import Line2D


def get_all_wells():
    stress_data = md.maw.perioddata.get_data()
    pumping_list = []
    for sp in stress_data.keys():
        data = stress_data[sp]
        for line in data:
            if (line['mawsetting'] == 'rate') and (line['mawsetting_data'] < 0):
                if line['wellno'] not in pumping_list:
                    pumping_list.append(line['wellno'])
    connection_data = md.maw.connectiondata.get_data()
    well_lrc_list = []
    for line in connection_data:
        if line['wellno'] in pumping_list:
            if line['cellid'] not in well_lrc_list:
                well_lrc_list.append(line['cellid'])
    return well_lrc_list


def get_faults():
    hfb = md.get_package('hfb')
    hfb_data_ = hfb.stress_period_data.get_data()
    hfb_spd = hfb_data_[0]
    lay_seg_dict = {}
    for lay in range(12):
        lay_seg_dict[lay] = {}
    for line in hfb_spd:
        cell1 = line[0]
        cell2 = line[1]
        if cell1[0] == cell2[0]:
            if cell1[0] in lay_seg_dict:
                if cell1[1] == cell2[1]:  # same row, vertical segment
                    if 'ver_seg' in lay_seg_dict[cell1[0]]:
                        existing = lay_seg_dict[cell1[0]]['ver_seg']
                        lay_seg_dict[cell1[0]]['ver_seg'] = (np.append(existing[0], np.minimum(cell1[1], cell2[1])), np.append(existing[1], np.minimum(cell1[2], cell2[2])))
                    else:
                        lay_seg_dict[cell1[0]]['ver_seg'] = (np.array([np.minimum(cell1[1], cell2[1])]), np.array([np.minimum(cell1[2], cell2[2])]))
                if cell1[2] == cell2[2]:  # same column, horizontal segment
                    if 'hor_seg' in lay_seg_dict[cell1[0]]:
                        existing = lay_seg_dict[cell1[0]]['hor_seg']
                        lay_seg_dict[cell1[0]]['hor_seg'] = (np.append(existing[0], np.minimum(cell1[1], cell2[1])), np.append(existing[1], np.minimum(cell1[2], cell2[2])))
                    else:
                        lay_seg_dict[cell1[0]]['hor_seg'] = (np.array([np.minimum(cell1[1], cell2[1])]), np.array([np.minimum(cell1[2], cell2[2])]))
    fault_dict = {}
    for lay in range(12):
        if ('ver_seg' in lay_seg_dict[lay]) and ('hor_seg' in lay_seg_dict[lay]):
            ver_seg = lay_seg_dict[lay]['ver_seg']
            hor_seg = lay_seg_dict[lay]['hor_seg']
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
            fault_dict[lay] = (bound_x_utms, bound_y_utms)
        else:
            fault_dict[lay] = None
    return fault_dict


def get_plotting_streams():
    stream_coords = pd.read_csv(os.path.join('..', 'data', 'plotting_streams.csv'))
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds.GetProjection())
    src_srs = src_srs.CloneGeogCS()
    tgt_srs=osr.SpatialReference()
    tgt_srs.ImportFromEPSG(4269)
    coord_dict = {}
    for dummy, line in stream_coords.iterrows():
        if line['GNIS_Name'] not in coord_dict:
            coord_dict[line['GNIS_Name']] = [(line['POINT_X'], line['POINT_Y'])]
        elif line['GNIS_Name'] in coord_dict:
            existing = coord_dict[line['GNIS_Name']]
            existing.append((line['POINT_X'], line['POINT_Y']))
            coord_dict[line['GNIS_Name']] = existing
    # converted_coord_dict = {}
    # for new_key in coord_dict.keys():
    #     converted_coord_dict[new_key] = ReprojectCoords(coord_dict[new_key], src_srs, tgt_srs)
    #     print(converted_coord_dict[new_key])
    return coord_dict


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


def plot_basin_boundaries(axes_, fs=5, legend=True, lw=0.5):
    axes_.plot(x_basin, y_basin, color='w', linewidth=lw, zorder=14, label='Groundwater\nBasin\nBoundaries')
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
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


def plot_subbasin_boundaries(axes_, fs=5, lw=0.5):
    axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw, zorder=13, label='Central Basin\nSub-Area\nBoundaries')
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    axes_.text(x=0.45*xrange + xlim[0], y=0.73*yrange + ylim[0], s='Los Angeles\nForebay\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
    axes_.text(x=0.6*xrange + xlim[0], y=0.48*yrange + ylim[0], s='Pressure\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
    axes_.text(x=0.75*xrange + xlim[0], y=0.58*yrange + ylim[0], s='Whittier\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
    axes_.text(x=0.62*xrange + xlim[0], y=0.725*yrange + ylim[0], s='Montebello\nForebay\nArea', ha='center', fontsize=fs, color='k', zorder=13, fontfamily='serif')
    return axes_


def plot_xsect_coords(axes_):
    coord_list = np.loadtxt(os.path.join('..', 'data', 'coord_list.txt'))
    axes_.plot(coord_list[:, 0], coord_list[:, 1], label="Cross-Section\nA-A'", color='r', linewidth=1)
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


def miscellaneous_annotations(axes_, fs=5):
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    fill = dict(facecolor='tab:gray', edgecolor='none', alpha=0.0)
    axes_.text(x=0.35*xrange + xlim[0], y=0.25*yrange + ylim[0], s='Palos\nVerdes\nHills', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    axes_.text(x=0.15*xrange + xlim[0], y=0.85*yrange + ylim[0], s='Santa\nMonica\nMountains', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    axes_.text(x=0.38*xrange + xlim[0], y=0.93*yrange + ylim[0], s='Hollywood\nHills', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    axes_.text(x=0.47*xrange + xlim[0], y=0.85*yrange + ylim[0], s='Elysian\nHills', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill,  fontfamily='serif')
    axes_.text(x=0.55*xrange + xlim[0], y=0.86*yrange + ylim[0], s='Repetto\nHills', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    axes_.text(x=0.59*xrange + xlim[0], y=0.8*yrange + ylim[0], s='Merced\nHills', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    axes_.text(x=0.72*xrange + xlim[0], y=0.76*yrange + ylim[0], s='Whittier\nNarrows', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    axes_.text(x=0.75*xrange + xlim[0], y=0.68*yrange + ylim[0], s='Puente\nHills', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    axes_.text(x=0.49*xrange + xlim[0], y=0.89*yrange + ylim[0], s='Los\nAngeles\nNarrows', ha='center', fontsize=fs, color='tab:olive', zorder=14, bbox=fill, fontfamily='serif')
    return axes_


def plot_faults(axes_, lw=1):
    faults = get_faults()
    leg = False
    for index_ in range(12):
        if faults[index_]:
            if not leg:
                axes_.plot(faults[index_][0], faults[index_][1], color='tab:cyan', linewidth=lw, label='Horizontal\nFlow\nBarriers\n(Faults)')
                leg = True
            else:
                axes_.plot(faults[index_][0], faults[index_][1], color='tab:cyan', linewidth=lw)
    return axes_


def plot_streams(axes_, lw=1):
    streams = get_plotting_streams()
    leg = False
    for stream in streams.keys():
        stream_array = np.array(streams[stream])
        stream_array = stream_array[stream_array[:, 1].argsort()]
        stream_array[stream_array[:, 1] > extent[3], :] = np.nan
        stream_array[stream_array[:, 1] < extent[2], :] = np.nan
        stream_array[stream_array[:, 0] > extent[1], :] = np.nan
        if stream == 'Rio Hondo':
            diff = stream_array[1:, 1] - stream_array[:-1, 1]
            max_diff_ind = np.argmax(np.abs(diff))
            stream_array = np.insert(arr=stream_array, obj=max_diff_ind+1, values=np.array([np.nan, np.nan]), axis=0)
        if not leg:
            axes_.plot(stream_array[:, 0], stream_array[:, 1], color='tab:blue', linewidth=lw, label='Stream\nChannels')
            leg = True
        else:
            axes_.plot(stream_array[:, 0], stream_array[:, 1], color='tab:blue', linewidth=lw)
    return axes_


def plot_wells(axes_):
    labeled_wells = False
    well_list = get_all_wells()
    for lrc in well_list:
        y_utm, x_utm = row_col_to_utm(lrc[1], lrc[2])
        if not labeled_wells:
            axes_.scatter(x_utm, y_utm, zorder=13, color='tab:pink', s=1, label='Pumping\nCells',
                                   linewidths=0)
            labeled_wells = True
        else:
            axes_.scatter(x_utm, y_utm, zorder=13, color='tab:pink', s=1, linewidths=0)
    return axes_


def plot_map():
    fig, axes = plt.subplots()
    axes.imshow(img, extent=extent, origin='upper')
    axes = plot_wells(axes)
    axes = plot_faults(axes, lw=0.5)
    axes = plot_streams(axes)
    axes = miscellaneous_annotations(axes)
    axes = plot_basin_boundaries(axes, fs=6, legend=False)
    axes = plot_subbasin_boundaries(axes, fs=6)
    axes = crop_edges(axes, .1, .1, 0.05, 0.02, x_spacing=10000, y_spacing=10000)
    axes = lat_lon_labels(axes, fs=6)
    axes = plot_xsect_coords(axes)
    handles, labels = axes.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color='tab:olive', lw=1))
    labels.append('Upland and\nCanyon\nAreas')
    handles = [handles[2], handles[3], handles[0], handles[1], handles[4], handles[5], handles[6]]
    labels = [labels[2], labels[3], labels[0], labels[1], labels[4], labels[5], labels[6]]
    fig.suptitle('Los Angeles Coastal Plain Groundwater Model (LACPGM) Area')
    fig.legend(handles=handles, labels=labels, loc=5, fontsize=5, facecolor='tab:gray', framealpha=0.2, markerscale=2)
    fig.tight_layout(pad=0.5, h_pad=0.1, w_pad=0.1)
    fig.subplots_adjust(right=0.92)
    fig.savefig(os.path.join('..', 'figures', 'map_figure1.png'), dpi=600)


dCellSize = 201.168
dXCoordinate = 355047.067638
dYCoordinate = 3777031.579092
nrow = 256
ncol = 312
x, y = np.meshgrid(np.linspace(dXCoordinate, dXCoordinate + ncol * dCellSize, ncol + 1), np.linspace(dYCoordinate, dYCoordinate - nrow * dCellSize, nrow + 1))
fname = os.path.join(os.path.join('..', 'basemaps', 'basemap.tif'))
img = plt.imread(os.path.join(fname))
# gdal
ds = gdal.Open(fname)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])
color_list = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper\nWilmington A', '8: Upper\nWilmington B', '9: Lower\nWilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach C\nand BC']
scotts_colors = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                 '#A020F0', '#FFA500', '#BDB76B']
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
plot_map()
