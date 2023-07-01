import os,pickle

for feedback in ['SB','BW']:
    SimInfo = pickle.load(open(f'Code/SimulationInfo.{feedback}.pickle','rb'))
    for sim in SimInfo:
        for halo in SimInfo[sim]['halos']:
            os.system(f'scp Images/{sim}.{feedback}/{halo}/*.tar.gz glaurung:/myhome2/users/vannest/ObservedShapes/Images/{sim}.{feedback}/{halo}/')