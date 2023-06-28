import os,pickle

print('Pulling Image Dictionaries...')
os.system('scp glaurung:/myhome2/users/vannest/ObservedShapes/Data/*.Images.pickle Data/')
print('Pulling Images...')
for feedback in ['SB','BW']:
    SimInfo = pickle.load(open(f'SimulationInfo.{feedback}.pickle','rb'))
    for sim in SimInfo:
        for halo in SimInfo[sim]['halos']:
            os.system(f'scp -r glaurung:/myhome2/users/vannest/ObservedShapes/Images/{sim}{halo}/*.tar.gz Images/{sim}{halo}')