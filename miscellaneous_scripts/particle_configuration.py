import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from flopy.mf6 import MFSimulation
import flopy.utils.binaryfile as bf
import datetime
import pandas as pd
from osgeo import gdal, osr



def particle_map():
    fig, axes = plt.subplots(figsize=(10, 10))
    xrange = np.arange(0.1, 1, 0.2) * 660
    yrange = np.arange(0.1, 1, 0.2) * 660
    print(0.2*660)
    for row in xrange:
        for col in yrange:
            axes.scatter(row, col, color='k', s=80)

    axes.set_xlim([0, 660])
    axes.set_ylim([0, 660])
    # axes.set_xticklabels([])
    # axes.set_yticklabels([])
    # axes.set_xticks([])
    # axes.set_yticks([])
    axes.tick_params(labelsize=15)
    axes.set_ylabel('Cell Length (Feet)', fontsize=15)
    axes.set_xlabel('Cell Width (Feet)', fontsize=15)
    axes.set_title('Initial Particle Configuration (Plan View)', fontsize=20)
    fig.tight_layout()
    fig.savefig(os.path.join('..', 'figures', 'particle_configuration.png'), dpi=300)
    return


particle_map()