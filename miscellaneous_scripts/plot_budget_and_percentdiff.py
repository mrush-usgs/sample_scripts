import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import calendar
import matplotlib.lines as mlines
import os

usg = pd.read_csv(os.path.join('..', 'data', 'USG.csv'))
mf6 = pd.read_csv(os.path.join('..', 'data', 'LA_Water_Budget.csv'))
fig, axes = plt.subplots(2, 1, figsize=(13.33, 10))
bud_list = ['WATER SPREADING',
            'PUMPING AND INJECTION',
            'RECHARGE',
            'STORAGE',
            'DRAINS',
            'GHB']
exclude_list = []
color_list = [
    'tab:purple',
    'tab:pink',
    'navy',
    'tab:gray',
    'tab:olive',
    'tab:cyan']
mf6_years = np.arange(1971, 2020)
usg_years = np.arange(1971, 2016)
data_arr_pct = np.zeros([len(bud_list), len(mf6_years)])
data_arr_bud = np.zeros([len(bud_list), len(mf6_years)])
year_total = np.sum(np.abs(mf6)/2, axis=1)[:-4]
for ind, bud in enumerate(bud_list):
    if bud not in exclude_list:
        if bud == 'STORAGE':
            mf6_plot = mf6['SPECIFIC STORAGE'] + mf6['SPECIFIC YIELD']
            usg_plot = usg[bud]
        elif bud == 'GHB':
            mf6_plot = mf6['GHB OCEAN'] + mf6['GHB OC'] + mf6['GHB WHT'] + mf6['GHB LAN']
            usg_plot = usg['GHB OCEAN'] + usg['GHB OC'] + usg['GHB WHT'] + usg['GHB LAN']
        elif bud == 'PUMPING AND INJECTION':
            mf6_plot = mf6['PUMPING']
            usg_plot = usg['INJECTION'] + usg['PUMPING']
        else:
            mf6_plot = mf6[bud]
            usg_plot = usg[bud]
        data_arr_pct[ind, :-4] = 100 * (mf6_plot[:-4] - usg_plot) / year_total
        data_arr_pct[ind, -4:] = 0
        data_arr_bud[ind, :] = mf6_plot

for subplot in range(2):
    if subplot == 0:
        plot_arr = data_arr_bud
        plot_title = '(a) Budget: LACPGM-MF6'
        ylab = 'Acre Feet\nper Year'
    else:
        plot_title = '(b) Budget Comparison by Percent: LACPGM-MF6 - LACPGM-USG'
        plot_arr = data_arr_pct
        ylab = 'Percent\nDifference'
    neg_cum_sum = np.nancumsum(plot_arr.clip(max=0), axis=0)
    cum_neg = np.zeros(np.shape(plot_arr))
    cum_neg[1:] = neg_cum_sum[:-1]
    pos_cum_sum = np.nancumsum(plot_arr.clip(min=0), axis=0)
    cum_pos = np.zeros(np.shape(plot_arr))
    cum_pos[1:] = pos_cum_sum[:-1]
    row_mask = plot_arr < 0
    cum_pos[row_mask] = cum_neg[row_mask]
    if subplot == 0:
        second_y = axes[subplot].twinx()
        second_y.tick_params(labelsize=15)
    for ind, bud in enumerate(bud_list):
        if bud not in exclude_list:
            if bud == 'GHB':
                axes[subplot].bar(mf6_years, plot_arr[ind], color=color_list[ind], label='Head\nBoundaries', bottom=cum_pos[ind])
            elif bud == 'WATER SPREADING':
                axes[subplot].bar(mf6_years, plot_arr[ind], color=color_list[ind], label='Water\nSpreading', bottom=cum_pos[ind])
            else:
                axes[subplot].bar(mf6_years, plot_arr[ind], color=color_list[ind], label=bud.capitalize(), bottom=cum_pos[ind])
            if subplot == 0:
                second_y.plot(mf6_years, np.nancumsum(plot_arr[ind])/1e6, color=color_list[ind])
                second_y.set_ylabel('Million\nAcre Feet\n(Cumulative)', ha='left', va='center', ma='left', rotation=0, fontsize=15)

            axes[subplot].set_ylabel(ylab, ha='right', va='center', ma='left', rotation=0, fontsize=15)
            axes[subplot].tick_params(labelsize=15)
            if subplot == 0:
                handles, labels = axes[subplot].get_legend_handles_labels()
    axes[subplot].set_title(plot_title, fontsize=15)

fig.legend(handles, labels, fontsize=15, loc=8, ncol=7)
fig.tight_layout()
# fig.subplots_adjust(right=0.84)
fig.subplots_adjust(bottom=0.1)
fig.savefig(os.path.join('..', 'figures', 'budget_and_percentdiff.png'), dpi=500)

