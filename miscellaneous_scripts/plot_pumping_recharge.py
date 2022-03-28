import os
import flopy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


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


def year_from_sp(sp):
    if sp < 180:
        year_ind = np.floor(sp / 4)
        remain = np.remainder(sp + 1, 4)
        days = 91.25
    else:
        year_ind = 45 + np.floor((sp - 180) / 12)
        remain = np.remainder(sp - 180 + 1, 12)
        days = 91.25 / 3
    year = 1971 + year_ind
    return year, int(year_ind), remain, days


def get_yearly_pumping():
    basin_dict, subbasin_dict, basin_array, subbasin_array = build_basin_dict()
    stress_data = maw.perioddata.get_data()
    pumping_ts = np.zeros(len(years))
    cb_pumping = np.zeros(len(years))
    stress_list = []
    cb_stress_list = []
    cb_well_list = []
    conn_data = maw.connectiondata.get_data()
    for line in conn_data:
        if basin_array[line['cellid'][1], line['cellid'][2]] == 2:
            cb_well_list.append(line['wellno'])
    for sp in stress_data.keys():
        for line in stress_data[sp]:
            if (line[1] == 'rate') and (line[2] < 0):
                stress_list.append(line[2])
            if (line[1] == 'rate') and (line[2] < 0) and (line['wellno'] in cb_well_list):
                cb_stress_list.append(line[2])
        year, year_ind, remain, days = year_from_sp(sp)
        if remain == 0:
            pumping_ts[year_ind] = np.sum(stress_list) * days * 2.29569e-5 / 1000 # cubic feet per day
            stress_list = []
            cb_pumping[year_ind] = np.sum(cb_stress_list) * days * 2.29569e-5 / 1000 # cubic feet per day
            cb_stress_list = []
    return pumping_ts, cb_pumping


def get_yearly_recharge():
    stress_data = rch.recharge.get_data()
    recharge_ts = np.zeros(len(years))
    stress_list = []
    for sp in stress_data.keys():
        stress_list.append(np.sum(stress_data[sp]))
        year, year_ind, remain, days = year_from_sp(sp)
        if remain == 0:
            recharge_ts[year_ind] = np.sum(stress_list)  * 660 * 660 * days * 2.29569e-5 / 1000 # feet per day
            stress_list = []
    return recharge_ts


def plot_pumping_recharge(highlight=False, plot_ind=0):
    # pumping_label = 'PUMPING:\nMean ' + str(np.round(np.mean(pumping), decimals=0)) + '\nMedian ' + \
    #                 str(int(np.round(np.median(pumping), 0)))
    # recharge_label = 'RECHARGE:\nMean ' + str(np.round(np.mean(recharge), decimals=0)) + '\nMedian ' + \
    #                  str(int(np.round(np.median(recharge), 0)))
    if plot_ind == 0:
        ncol = 2
        axes[plot_ind].bar(years, pumping, color='tab:pink', label='Pumping')
        axes[plot_ind].bar(years, recharge, color='navy', label='Recharge')
        # print(np.mean(pumping))
        # print(np.mean(recharge))
        axes[plot_ind].plot(years, np.ones(len(years))*np.mean(pumping), color='tab:pink', label='Mean Annual Pumping', linestyle=':')
        axes[plot_ind].plot(years, np.ones(len(years))*np.mean(recharge), color='navy', label='Mean Annual Recharge', linestyle=':')
        axes[plot_ind].set_ylim([-470, 220])
    else:
        ncol = 1
        axes[plot_ind].bar(years, pumping, color='tab:pink')
        axes[plot_ind].bar(years, recharge, color='navy')

    if highlight:
        new_years = np.arange(2020, 2041)
        highlight_pumping = np.copy(pumping)
        highlight_pumping[0:13]= np.nan
        highlight_pumping[34:] = np.nan
        highlight_recharge = np.copy(recharge)
        highlight_recharge[0:13] = np.nan
        highlight_recharge[34:] = np.nan
        axes[plot_ind].bar(years, highlight_recharge, color='navy', label='Recharge Selected to Represent 2020-2040', edgecolor='k', linewidth=1)
        axes[plot_ind].bar(years, highlight_pumping, color='tab:pink', label='Pumping Selected to Represent 2020-2040', edgecolor='k', linewidth=1)
        axes[plot_ind].bar(new_years, highlight_recharge[13:34], color='navy', edgecolor='k', linewidth=1)
        axes[plot_ind].bar(new_years, highlight_pumping[13:34], color='tab:pink', edgecolor='k', linewidth=1)
        axes[plot_ind].set_ylim([-470, 220])
        axes[plot_ind].set_title('LACPGM Pumping and Recharge: 1971-2040', fontsize=15)
    else:
        axes[plot_ind].set_title('LACPGM Pumping and Recharge: 1971-2019', fontsize=15)
    axes[plot_ind].set_ylabel('Thousand\nAcre Feet', ha='right', va='center', ma='left', rotation=0, fontsize=12)
    axes[plot_ind].set_xlabel('Year', fontsize=12)
    axes[plot_ind].legend(ncol=ncol, loc=8, fontsize=12)
    axes[plot_ind].tick_params(labelsize=12)


def analyze_years():
    num_years = 21
    starts = np.arange(1971, 2000)
    pumping_devs = np.zeros(len(starts))
    recharge_devs = np.zeros(len(starts))
    pumping_std = np.std(pumping)
    recharge_std = np.std(recharge)
    for ind, start in enumerate(starts):
        pumping_devs[ind] = (np.mean(pumping[ind:ind+num_years]) - np.mean(pumping)) / pumping_std
        recharge_devs[ind] = (np.mean(recharge[ind:ind+num_years]) - np.mean(recharge)) / recharge_std
        # print(starts[ind])
        # # print(pumping_devs[ind])
        # # print(recharge_devs[ind])
        # print(np.abs(pumping_devs[ind]) + np.abs(recharge_devs[ind]))
        if start == 1999:
            print(ind)
    axes[1].bar(starts, np.abs(pumping_devs), color='tab:pink', label='Pumping')
    axes[1].bar(starts, np.abs(recharge_devs), bottom=np.abs(pumping_devs), color='navy', label='Recharge')
    axes[1].set_ylabel('Number of\nStandard\nDeviations', ha='right', va='center', ma='left', rotation=0, fontsize=12)
    axes[1].set_title('LACPGM Deviations from Average', fontsize=15)
    axes[1].set_xlabel('First Year of 21-Year Range', fontsize=12)
    axes[1].legend(fontsize=12, loc=9)
    axes[1].tick_params(labelsize=12)

    print(pumping_devs[28])
    print(recharge_devs[28])


years = np.arange(1971, 2020)
sim = flopy.mf6.MFSimulation.load(sim_ws=os.path.join('..', '2019_MF6_Model', 'lab_model'), verbosity_level=0)
md = sim.get_model()
maw = md.maw
rch = md.get_package('rch')

pumping, cb_pumping = get_yearly_pumping()
recharge = get_yearly_recharge()

# print(np.nansum(cb_pumping))
# print(np.nanmean(cb_pumping))

fig, axes = plt.subplots(3, 1, figsize=(12, 15))
plot_pumping_recharge()
analyze_years()
plot_pumping_recharge(highlight=True, plot_ind=2)
fig.suptitle('Selection of Representative Years for Model Extension', fontsize=20)
fig.tight_layout()
fig.subplots_adjust(top=0.92)
fig.savefig(os.path.join('..', 'figures', 'pumping_recharge.png'), dpi=500)
