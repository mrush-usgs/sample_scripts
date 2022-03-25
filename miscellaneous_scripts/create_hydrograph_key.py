import matplotlib.pyplot as plt
import numpy as np
import os

color_list = ['#FFFF00', '#FFD700', '#40E0D0', '#FF0000', '#00FF00', '#0000FF', '#859B5F', '#556B2F', '#FFC0CB',
                    '#A020F0', '#FFA500', '#BDB76B']
lw_list = [1,1,1,2,2,2,2,2,4,4,4,4]
layer_list = ['1: Dominguez', '2: Mesa', '3: Pacific A', '4: Pacific', '5: Harbor', '6: Bent Spring',
              '7: Upper Wilmington A', '8: Upper Wilmington B', '9: Lower Wilmington', '10: Long Beach A',
              '11: Long Beach B', '12: Long Beach\nC and BC']
dummy_fig, dummy_axes = plt.subplots()
legend_fig = plt.figure(figsize=(12, 15))

for lay in range(12):
    dummy_axes.plot(np.arange(100), np.arange(100), color=color_list[lay], linewidth=lw_list[lay]*5, label=layer_list[lay])

handles, labels = dummy_axes.get_legend_handles_labels()
legend_fig.legend(handles, labels, fontsize=50)
# legend_fig.tight_layout()
legend_fig.savefig(os.path.join('..', 'figures', 'hydrograph_legend.pdf'))