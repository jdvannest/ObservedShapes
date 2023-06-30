import os,pickle

loop = True
while loop:
    i_dict = input('Pull Image Dictionary? (y/n): ')
    if i_dict in ['y','n']:
        loop = False 
if i_dict=='y': os.system('scp glaurung:/myhome2/users/vannest/ObservedShapes/Data/*.Images.pickle Data/')
loop = True
while loop:
    i_tar = input('Pull Image tarballs? (y/n): ')
    if i_tar in ['y','n']:
        loop = False
if i_tar=='y':
    for feedback in ['SB','BW']:
        SimInfo = pickle.load(open(f'SimulationInfo.{feedback}.pickle','rb'))
        for sim in SimInfo:
            for halo in SimInfo[sim]['halos']:
                os.system(f'scp glaurung:/myhome2/users/vannest/ObservedShapes/Images/{sim}.{feedback}/{halo}/*.tar.gz Images/{sim}.{feedback}/{halo}')