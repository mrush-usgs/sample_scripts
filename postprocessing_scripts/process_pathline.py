import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from osgeo import gdal, osr
import pandas as pd
import matplotlib.ticker as plticker
from matplotlib import markers


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


def plot_basin_boundaries(axes_, fs=5, legend=True, lw=0.5, label=False, legend_label=True, h_font=5):
    if not legend_label:
        axes_.plot(x_basin, y_basin, color='w', linewidth=lw, zorder=15, label='Groundwater\nBasin\nBoundaries')
    else:
        axes_.plot(x_basin, y_basin, color='w', linewidth=lw, zorder=15)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    if label:
        axes_.text(x=0.19*xrange + xlim[0], y=0.60*yrange + ylim[0], s='Santa\nMonica\nBay', ha='center', fontsize=fs, color='w', fontfamily='serif')
        axes_.text(x=0.25*xrange + xlim[0], y=0.75*yrange + ylim[0], s='Santa\nMonica\nBasin', ha='center', fontsize=fs, color='w', fontfamily='serif')
        if legend:
            axes_.text(x=0.47*xrange + xlim[0], y=0.32*yrange + ylim[0], s='West\nCoast\nBasin', ha='center', fontsize=fs, color='w',  fontfamily='serif', zorder=15)
        else:
            axes_.text(x=0.35*xrange + xlim[0], y=0.4*yrange + ylim[0], s='West\nCoast\nBasin', ha='center', fontsize=fs, color='w',  fontfamily='serif', zorder=15)
        axes_.text(x=0.58*xrange + xlim[0], y=0.13*yrange + ylim[0], s='San\nPedro\nBay', ha='center', fontsize=fs, color='w',fontfamily='serif')
        axes_.text(x=0.5*xrange + xlim[0], y=0.5*yrange + ylim[0], s='Central\nBasin', ha='center', fontsize=fs, color='w', fontfamily='serif', zorder=15)
        axes_.text(x=0.73*xrange + xlim[0], y=0.35*yrange + ylim[0], s='Orange\nCounty', ha='center', fontsize=fs, color='w',  fontfamily='serif', zorder=15)
        if fs<8:
            axes_.text(x=0.37*xrange + xlim[0], y=0.86*yrange + ylim[0], s='Hollywood\nBasin', ha='center', fontsize=h_font, rotation=10, color='w',fontfamily='serif')
        else:
            axes_.text(x=0.37*xrange + xlim[0], y=0.84*yrange + ylim[0], s='Hollywood\nBasin', ha='center', fontsize=h_font, rotation=10, color='w',  fontfamily='serif')
    return axes_


def plot_subbasin_boundaries(axes_, fs=5, lw=0.5, label=True, legend_label=True, zoomed=False):
    if not legend_label:
        axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw, zorder=13, label='Central Basin\nSub-Area\nBoundaries')
    else:
        axes_.plot(x_subbasin, y_subbasin, color='k', linewidth=lw, zorder=13)
    xlim = axes_.get_xlim()
    ylim = axes_.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[1] - ylim[0]
    if label:
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


def plot_zoomed_paths_by_group():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    for group_ind, group_name in enumerate(pl_data.keys()):
        if group_ind == 4:
            plot_ind = 3
        else:
            plot_ind = group_ind
        group_data = pl_data[group_name]
        axes[plot_ind].imshow(img, extent=extent, origin='upper')
        for prtcl in group_data.keys():
            pdata = group_data[prtcl]
            if prtcl == 1:
                axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, label=group_name)
            else:
                axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw)
        if group_ind < 4:
            axes[plot_ind] = plot_basin_boundaries(axes[plot_ind], fs=5, legend_label=basin_label)
            axes[plot_ind] = plot_subbasin_boundaries(axes[plot_ind], fs=5, legend_label=subbasin_label, zoomed=True)
            basin_label = True
            subbasin_label = True
            axes[plot_ind] = crop_edges(axes[plot_ind],  .43, .23, .4, .17)
            axes[plot_ind] = lat_lon_labels(axes[plot_ind], fs=5)
    handles = []
    labels = []
    for ax in fig.axes:
        axLine, axLabel = ax.get_legend_handles_labels()
        handles.extend(axLine)
        labels.extend(axLabel)
    handles = [handles[0], handles[3], handles[4], handles[5], handles[6], handles[1], handles[2]]
    labels = [labels[0], labels[3], labels[4], labels[5], labels[6], labels[1], labels[2]]
    fig.legend(handles=handles, labels=labels, loc=7, fontsize=6, facecolor='tab:gray', framealpha=0.2)
    fig.suptitle('Montebello Forebay Flow Paths')
    fig.subplots_adjust(right=0.8)
    fig.savefig(os.path.join(figure_path, 'paths_by_group.png'), dpi=500)


def plot_zoomed_layer_paths_by_group():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    cm = plt.get_cmap('gist_earth_r')
    labeled_scatter = False
    for group_ind, group_name in enumerate(pl_data.keys()):
        if group_ind == 4:
            plot_ind = 3
        else:
            plot_ind = group_ind
        group_data = pl_data[group_name]
        axes[plot_ind].imshow(img, extent=extent, origin='upper')
        for p_ind, prtcl in enumerate(group_data.keys()):
            pdata = group_data[prtcl]
            for lay in range(12):
                if youngest_on_top:
                    zorders = (12 - lay)
                else:
                    zorders = lay + 1
                lay_subset = pdata[:, 2] == (lay + 1)
                bool_to_int = [int(bool) for bool in lay_subset]
                lay_diff = np.diff(bool_to_int)
                if np.sum(lay_diff != 0) > 2: # switches from False to True or vice versa more than twice
                    firsts = np.where(lay_diff == 1)[0].tolist()  # goes from False to True
                    lasts = np.where(lay_diff == -1)[0].tolist()  # goes from True to False
                    for chunk_ind, first in enumerate(firsts):
                        new_subset = np.copy(lay_subset)
                        new_subset[:first] = False
                        if chunk_ind != len(lasts):
                            last = lasts[chunk_ind]
                            new_subset[last+1:] = False
                        if (p_ind == 0) and (group_ind == 0):
                            axes[plot_ind].plot(pdata[new_subset, 0], pdata[new_subset, 1], linewidth=lw, label=layer_list[lay], zorder=zorders, color=scotts_colors[lay]) # color=cm((lay + 1) / 12)
                        else:
                            axes[plot_ind].plot(pdata[new_subset, 0], pdata[new_subset, 1], linewidth=lw, zorder=zorders,
                                                color=scotts_colors[lay]) # color=cm((lay + 1) / 12)
                else:
                    if (p_ind == 0) and (group_ind == 0):
                        axes[plot_ind].plot(pdata[lay_subset, 0], pdata[lay_subset, 1], linewidth=lw,
                                            label=layer_list[lay], zorder=(12 - lay), color=scotts_colors[lay])
                    else:
                        axes[plot_ind].plot(pdata[lay_subset, 0], pdata[lay_subset, 1], linewidth=lw, zorder=zorders,
                                            color=scotts_colors[lay])
                # ep_subset = (ep_data['Particle Group'] == (group_ind + 1)) & (ep_data['Particle ID'] == prtcl)
                # status = ep_data['Status'].loc[ep_subset].item()
                # if status == 5:
                #     if not labeled_scatter:
                #         axes[plot_ind].scatter(pdata[-1:, 0], pdata[-1:, 1], color='k', zorder=15, s=2, linewidths=0,
                #                             label='Pathline\nTerminates\nat Well')
                #         labeled_scatter = True
                #     else:
                #         axes[plot_ind].scatter(pdata[-1:, 0], pdata[-1:, 1], color='k', zorder=15, s=2, linewidths=0)
        if group_ind < 4:
            axes[plot_ind] = plot_basin_boundaries(axes[plot_ind], fs=6, legend_label=basin_label)
            axes[plot_ind] = plot_subbasin_boundaries(axes[plot_ind], fs=5, legend_label=subbasin_label, zoomed=True)
            basin_label = True
            subbasin_label = True
            axes[plot_ind] = crop_edges(axes[plot_ind],  .43, .22, .35, .17)
            axes[plot_ind] = lat_lon_labels(axes[plot_ind], fs=5)
            if group_ind == 3:
                axes[plot_ind].set_title('Whittier Narrows Underflow', fontsize=8)
            else:
                axes[plot_ind].set_title(group_name.replace('\n', ' '), fontsize=8)
    leg = fig.legend(loc=7, fontsize=6, facecolor='tab:gray', framealpha=0.2)
    for ind, line in enumerate(leg.get_lines()):
        if ind < 12:
            line.set_linewidth(4.0)
    fig.suptitle('Montebello Forebay Flow Paths')
    fig.tight_layout(pad=0.5, h_pad=0.1, w_pad=0.1)
    fig.subplots_adjust(right=0.86)
    fig.savefig(os.path.join(figure_path, 'layer_paths_by_group.png'), dpi=500)


def plot_all_zoomed_paths():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots()
    axes.imshow(img, extent=extent, origin='upper')
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        for prtcl in group_data.keys():
            pdata = group_data[prtcl]
            if prtcl == 1:
                axes.plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw,
                              label=group_name)
            else:
                axes.plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw)
    axes = plot_basin_boundaries(axes, fs=6, legend_label=basin_label)
    axes = plot_subbasin_boundaries(axes, fs=8, legend_label=subbasin_label, zoomed=True)
    basin_label = True
    subbasin_label = True
    axes = crop_edges(axes, .43, .23, .4, .17)
    axes = lat_lon_labels(axes, fs=8)
    fig.legend(ncol=8, loc=8, fontsize=5, facecolor='tab:gray', framealpha=0.2)
    fig.subplots_adjust(bottom=0.2)
    fig.savefig(os.path.join(figure_path, 'all_zoomed_paths.png'), dpi=500)


def plot_zoomed_paths():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.5))
    axes = axes.flat
    for pl in range(2):
        axes[pl].imshow(img, extent=extent, origin='upper')
        if pl == 0:
            for group_ind, group_name in enumerate(pl_data.keys()):
                if group_ind < 3:
                    group_data = pl_data[group_name]
                    for prtcl in group_data.keys():
                        pdata = group_data[prtcl]
                        if prtcl == 1:
                            axes[pl].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw,
                                          label=group_name)
                        else:
                            axes[pl].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw)
        else:
            for group_ind, group_name in enumerate(pl_data.keys()):
                if group_ind >= 3:
                    group_data = pl_data[group_name]
                    for prtcl in group_data.keys():
                        pdata = group_data[prtcl]
                        if prtcl == 1:
                            axes[pl].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw,
                                          label=group_name)
                        else:
                            axes[pl].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw)
        axes[pl] = plot_basin_boundaries(axes[pl], fs=6, legend_label=basin_label)
        axes[pl] = plot_subbasin_boundaries(axes[pl], fs=6, legend_label=subbasin_label, zoomed=True)
        basin_label = True
        subbasin_label = True
        axes[pl] = crop_edges(axes[pl],  .43, .23, .4, .17)
        axes[pl] = lat_lon_labels(axes[pl], fs=5)
    handles = []
    labels = []
    for ax in fig.axes:
        axLine, axLabel = ax.get_legend_handles_labels()
        handles.extend(axLine)
        labels.extend(axLabel)
    handles = [handles[0], handles[1], handles[2], handles[5], handles[6], handles[3], handles[4]]
    labels = [labels[0], labels[1], labels[2], labels[5], labels[6], labels[3], labels[4]]
    fig.legend(handles=handles, labels=labels, ncol=8, loc=8, fontsize=5, facecolor='tab:gray', framealpha=0.2)
    fig.suptitle('Montebello Forebay Flow Paths')
    fig.subplots_adjust(bottom=0.12)
    fig.savefig(os.path.join(figure_path, 'zoomed_paths.png'), dpi=500)


def plot_zoomed_paths_combo():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots(2, 2)
    axes = axes.flat
    labeled_wells = False
    for pl in range(5):
        if pl == 4:
            plot_ind = 3
        else:
            plot_ind = pl
        axes[plot_ind].imshow(img, extent=extent, origin='upper')
        if (pl == 0) or (pl == 2):
            for group_ind, group_name in enumerate(pl_data.keys()):
                if group_ind < 3:
                    group_data = pl_data[group_name]
                    for prtcl in group_data.keys():
                        pdata = group_data[prtcl]
                        if (prtcl == 1) and (pl == 0):
                            axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, label=group_name)
                        else:
                            axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw)
                        # ep_subset = (ep_data['Particle Group'] == (group_ind + 1)) & (ep_data['Particle ID'] == prtcl)
                        # status = ep_data['Status'].loc[ep_subset].item()
                        # if status == 5:
                        #      axes[plot_ind].plot(pdata[-3:, 0], pdata[-3:, 1], color='k', zorder=15, linewidth=lw)
        else:
            for group_ind, group_name in enumerate(pl_data.keys()):
                if group_ind >= 3:
                    group_data = pl_data[group_name]
                    for prtcl in group_data.keys():
                        pdata = group_data[prtcl]
                        if (prtcl == 1) and (pl == 1):
                            axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, label=group_name)
                        else:
                            axes[plot_ind].plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw)
                        # ep_subset = (ep_data['Particle Group'] == (group_ind + 1)) & (ep_data['Particle ID'] == prtcl)
                        # status = ep_data['Status'].loc[ep_subset].item()
                        # if status == 5:
                        #      axes[plot_ind].plot(pdata[-3:, 0], pdata[-3:, 1], color='k', zorder=15, linewidth=lw)
        if pl <= 1:
            axes[plot_ind] = plot_basin_boundaries(axes[plot_ind], fs=4, legend=False, label=True, legend_label=basin_label, lw=0.3, h_font=2.5)
            axes[plot_ind] = plot_subbasin_boundaries(axes[plot_ind], fs=4, label=False, legend_label=subbasin_label, lw=0.3)
            basin_label = True
            subbasin_label = True
            axes[plot_ind] = crop_edges(axes[plot_ind], .13, .1, .1, .07, x_spacing=10000, y_spacing=10000)
        elif pl < 4:
            axes[plot_ind] = plot_basin_boundaries(axes[plot_ind], fs=4, legend_label=basin_label)
            axes[plot_ind] = plot_subbasin_boundaries(axes[plot_ind], fs=4, legend_label=subbasin_label, zoomed=True)
            basin_label = True
            subbasin_label = True
            axes[plot_ind] = crop_edges(axes[plot_ind], .43, .23, .35, .17)
        if pl < 4:
            axes[plot_ind] = lat_lon_labels(axes[plot_ind], fs=5)
            # well_list = get_all_wells()
            # for lrc in well_list:
            #     y_utm, x_utm = row_col_to_utm(lrc[1], lrc[2])
            #     if not labeled_wells:
            #         axes[plot_ind].scatter(x_utm, y_utm, zorder=13, color='k', s=1, label='Pumping\nCells',
            #                                linewidths=0)
            #         labeled_wells = True
            #     else:
            #         axes[plot_ind].scatter(x_utm, y_utm, zorder=13, color='k', s=1, linewidths=0)
    handles = []
    labels = []
    for ax in fig.axes:
        axLine, axLabel = ax.get_legend_handles_labels()
        handles.extend(axLine)
        labels.extend(axLabel)
    handles = [handles[0], handles[1], handles[2], handles[5], handles[6], handles[3], handles[4]]
    labels = [labels[0], labels[1], labels[2], labels[5], labels[6], labels[3], labels[4]]
    fig.legend(handles=handles, labels=labels, ncol=8, loc=8, fontsize=5, facecolor='tab:gray', framealpha=0.2)
    fig.suptitle('Montebello Forebay Flow Paths')
    fig.subplots_adjust(bottom=0.15)
    fig.savefig(os.path.join(figure_path, 'zoomed_paths_combo.png'), dpi=500)


def plot_all_paths():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots()
    axes.imshow(img, extent=extent, origin='upper')
    for group_ind, group_name in enumerate(pl_data.keys()):
        group_data = pl_data[group_name]
        if group_ind < 4:
            z_order = 4-group_ind
        else:
            z_order = 5
        for prtcl in group_data.keys():
            pdata = group_data[prtcl]
            if prtcl == 1:
                axes.plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, label=group_name, zorder=z_order)
            else:
                axes.plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, zorder=z_order)
        # if group_ind == 4:
    axes = plot_basin_boundaries(axes, fs=6, legend_label=basin_label, label=True, h_font=6)
    axes = plot_subbasin_boundaries(axes, fs=6, legend_label=subbasin_label, zoomed=True)
    l = axes.legend(loc=6, facecolor='tab:gray', framealpha=0.8, fontsize=8)
    basin_label = True
    subbasin_label = True
    l.set_zorder(20)
    axes = crop_edges(axes, .13, .1, .1, .07, x_spacing=10000, y_spacing=10000)
    axes = lat_lon_labels(axes)
    fig.suptitle('Montebello Forebay Flow Paths')
    fig.tight_layout()
    fig.savefig(os.path.join(figure_path, 'all_paths.png'), dpi=500)


def plot_spreading_grounds():
    basin_label = False
    subbasin_label = False
    lw = 0.5
    fig, axes = plt.subplots()
    axes.imshow(img, extent=extent, origin='upper')
    for group_ind, group_name in enumerate(pl_data.keys()):
        if group_ind < 3:
            group_data = pl_data[group_name]
            for prtcl in group_data.keys():
                pdata = group_data[prtcl]
                if prtcl == 1:
                    axes.plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, label=group_name, zorder=12)
                else:
                    axes.plot(pdata[:, 0], pdata[:, 1], color=color_list[group_ind], linewidth=lw, zorder=12)
            # if group_ind == 2:
    axes = plot_subbasin_boundaries(axes, fs=6, legend_label=subbasin_label, zoomed=True)
    axes = plot_basin_boundaries(axes, fs=6, legend_label=basin_label, label=True, h_font=6)
    l = axes.legend(loc=6, facecolor='tab:gray', framealpha=0.8, fontsize=8)
    l.set_zorder(20)
    basin_label = True
    subbasin_label = True
    axes = crop_edges(axes, .13, .1, .1, .07, x_spacing=10000, y_spacing=10000)
    axes = lat_lon_labels(axes)
    fig.suptitle('Montebello Forebay Spreading Grounds Flow Paths')
    fig.tight_layout()
    fig.savefig(os.path.join(figure_path, 'spreading_grounds_paths.png'), dpi=500)


youngest_on_top = False
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
for run in runs:
    for porosity in porosities:
        figure_path = os.path.join('..', 'figures', 'POROSITY' + porosity, run)
        sim_path = os.path.join('MODPATH', 'POROSITY', 'POROSITY_' + porosity, run, 'lab_model')
        pl_data, p_groups = build_data_structure(sim_path)
        ep_data = build_endpoints(sim_path)
        # plot_zoomed_paths_by_group()
        # plot_all_zoomed_paths()
        # plot_zoomed_paths()
        # plot_zoomed_paths_combo()
        # plot_all_paths()
        # plot_spreading_grounds()
        plot_zoomed_layer_paths_by_group()