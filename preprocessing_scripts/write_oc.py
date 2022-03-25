import os
import flopy
import numpy as np
import pandas


# load the lacpgm
print('loading model...')
sim = flopy.mf6.MFSimulation.load(sim_ws=os.path.join('..', '2019_MF6_Model', 'lab_model'), verbosity_level=0)
model = sim.get_model()


# oc package
print('processing oc...')
oc = model.get_package('oc')
sr = oc.saverecord.get_data()
new_sr = {}
for sp in sr.keys():
    sr_sp = sr[sp]
    sr_sp['ocsetting'] = 'all'
    new_sr[sp] = sr_sp

oc.saverecord.set_data(new_sr)
oc.filename = 'lab_model.oc'
oc.write()

