import os
import flopy
import numpy as np
import pandas


def extend_stress_period_data(stress_period_data):  # the dictionary version
    new_spd = {}
    for index in range(180):
        stress_period = stress_period_data[index]
        if start_sp <= index <= end_sp:
            new_spd[(index - start_sp) + 180] = stress_period
        new_spd[index] = stress_period
    return new_spd


def extend_record(stress_period_data):  # the recarray version
    new_spd = {}
    if len(stress_period_data) == 1:
        new_spd = stress_period_data  # only one stress period record, does not need to be duplicated
    else:
        for index, stress_period in stress_period_data.items():
            if start_sp <= index <= end_sp:
                new_spd[(index - start_sp) + 180] = stress_period
            new_spd[index] = stress_period
    return new_spd


def extend_time(pd_, new_sp_):
    for sp in range(new_sp_):
        pd_.append((91.25, 5, 1.0))
    return pd_


def extend_rch(stress_period_data):
    new_spd = {}
    new_irch = {}
    for index in range(180):
        stress_period = stress_period_data[index]
        rch_file = 'rch_{:03d}.asc'.format(index + 1)
        if start_sp <= index <= end_sp:
            new_irch[(index - start_sp) + 180] = {'filename': os.path.join('..', 'external_files', 'arrays', 'irch.asc')}
            new_spd[(index - start_sp) + 180] = {'filename': os.path.join('..', 'external_files', 'arrays', 'rch', rch_file), 'data': stress_period}
        new_irch[index] = {'filename': os.path.join('..', 'external_files', 'arrays', 'irch.asc')}
        new_spd[index] = {'filename': os.path.join('..', 'external_files', 'arrays', 'rch', rch_file), 'data': stress_period}
    return new_spd, new_irch


def make_directories(new_path_):
    if not os.path.isdir(new_path_):
        os.mkdir(new_path_)
    if not os.path.isdir(os.path.join(new_path_, 'lab_model')):
        os.mkdir(os.path.join(new_path_, 'lab_model'))
    if not os.path.isdir(os.path.join(new_path_, 'external_files')):
        os.mkdir(os.path.join(new_path_, 'external_files'))
    if not os.path.isdir(os.path.join(new_path_, 'external_files', 'arrays')):
        os.mkdir(os.path.join(new_path_, 'external_files', 'arrays'))
    if not os.path.isdir(os.path.join(new_path_, 'external_files', 'arrays', 'rch')):
        os.mkdir(os.path.join(new_path_, 'external_files', 'arrays', 'rch'))
    return


# some time conversion numbers:
new_sp = 21 * 4  # 21 years with quarterly stress periods
# start_sp = 40 # 1981
start_sp = 52 # 1984
end_sp = start_sp + new_sp # 2001

print('\nmaking directories...')
new_path = os.path.join('EXT_1984_2004')
make_directories(new_path)

# load the lacpgm
print('loading model...')
sim = flopy.mf6.MFSimulation.load(sim_ws=os.path.join('..', '2019_MF6_Model', 'lab_model'), verbosity_level=0)
model = sim.get_model()

# tdis package
print('processing tdis...')
tdis = sim.get_package('tdis')
pd = tdis.perioddata.get_data().tolist()
pd_new = extend_time(pd, new_sp)
tdis.perioddata.set_data(pd_new)
tdis.nper = 228 + new_sp

# oc package
print('processing oc...')
oc = model.get_package('oc')
pr = extend_record(oc.printrecord.get_data())
sr = extend_record(oc.saverecord.get_data())
oc.printrecord.set_data(pr)
oc.saverecord.set_data(sr)

# ghb packages
print('processing ghbtr...')
ghb_tr = model.get_package('ghbtr')
spd = extend_stress_period_data(ghb_tr.stress_period_data.get_data())
ghb_tr.stress_period_data.set_data(spd)
ghb_tr.observations = None

print('processing la_ghb...')
ghb_la = model.get_package('la_ghb')
spd = extend_stress_period_data(ghb_la.stress_period_data.get_data())
ghb_la.stress_period_data.set_data(spd)
ghb_la.observations = None

# wel package
print('processing wel...')
wel = model.get_package('wel')
spd = extend_stress_period_data(wel.stress_period_data.get_data())
wel.observations = None
wel.stress_period_data.set_data(spd)

# maw package
print('processing maw...')
maw = model.get_package('maw')
spd = extend_stress_period_data(maw.perioddata.get_data())
maw.observations = None
maw.perioddata.set_data(spd)

# recharge package
print('getting recharge...')
rch = model.get_package('rch')
rch_data = rch.recharge.get_data()

# now set new simulation path
print('setting new simulation path...')
sim.simulation_data.mfpath.set_sim_path(os.path.join(new_path, 'lab_model'))

# recharge package
print('processing recharge...')
recharge, irecharge = extend_rch(rch_data)
rch.recharge.set_data(recharge)
rch.irch.set_data(irecharge)
rch.observations = None

# write out updated simulation
print('writing simulation...')
sim.write_simulation(ext_file_action=flopy.mf6.mfbase.ExtFileAction.copy_all)
