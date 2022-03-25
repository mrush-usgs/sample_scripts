import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import flopy


def wtryr_from_sp(sp_):
    if sp_ < 180:
        yr_ = 1971 + np.floor(((sp_ + 1) / 4))
    elif sp_ == 180:
        yr_ = 2016
    else:
        yr_ = 2016 + np.floor(((sp_ - 180 + 3) / 12))
    return yr_


any_sim_path = os.path.join('..', '2019_MF6_Model', 'lab_model')
sim = flopy.mf6.MFSimulation.load(sim_ws=any_sim_path, verbosity_level=0)
md = sim.get_model()
well_spd = md.wel.stress_period_data.get_data()

sgs = ['RHSG', 'SGRSG', 'WND']
title_list = ['Rio Hondo Spreading Grounds', 'San Gabriel River Spreading Grounds', 'Whittier Narrows Reservoir']
color_list = ['tab:orange', 'tab:blue', 'tab:green']

# begin_year = 2015 - 1971
# begin_year = 2010 - 1971
begin_year = 0
end_year = 2020 - 1971

fig, axes = plt.subplots(3, 1, figsize=(10, 15))
for sg_ind, sg in enumerate(sgs):
    sp_array = np.zeros([len(well_spd.keys()), 2])
    for sp in well_spd.keys():
        sp_data = well_spd[sp]
        sp_total = []
        if sp < 180:
            factor = 2.29569e-5 * 91.25  # convert from cubic feet per day to acre feet per stress period
        else:
            factor = 2.29569e-5 * 91.25 / 3 # convert from cubic feet per day to acre feet per stress period
        for line in sp_data:
            lrc = line[0]
            rate = line[1] # volumetric (cubic feet per day)
            id = line[2]
            if id == 'q' + sg.lower():
                sp_total.append(rate * factor)
        yr = wtryr_from_sp(sp)
        sp_array[sp, 0] = yr
        sp_array[sp, 1] = np.sum(sp_total)

    yr_array = np.zeros([49, 2])
    yr_array[:, 0] = np.arange(1971, 2020)
    for ind in range(len(yr_array)):
        wtryr_subset = (sp_array[:, 0] == yr_array[ind, 0])
        yr_array[ind, 1] = np.nansum(sp_array[wtryr_subset, 1])

    yr_array[0, 1] = 0

    axes[sg_ind].bar(yr_array[begin_year:end_year, 0], yr_array[begin_year:end_year, 1], color=color_list[sg_ind])
    axes[sg_ind].set_ylabel('Acre Feet\nper Year', ha='right', va='center', ma='left', rotation=0, fontsize=15)
    axes[sg_ind].set_title(title_list[sg_ind], fontsize=20)
    axes[sg_ind].set_ylim([0, 105000])
    axes[sg_ind].tick_params(labelsize=15)
fig.tight_layout()
fig.savefig(os.path.join('..', 'figures', 'spreading.png'), dpi=500)

