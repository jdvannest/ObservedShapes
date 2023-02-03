import os,pickle

gitdir = os.getcwd().rstrip('Code')
loop,erase = True,False
while loop:
    query = input('Delete .png files after compressing? (y/n): ')
    if query in ['y','Y']:
        erase = True
        loop = False
    elif query in ['n','N']:
        erase = False
        loop = False

for feedback in ['BW','SB']:
    os.chdir(f'{gitdir}Code')
    SimInfo = pickle.load(open(f'SimulationInfo.{feedback}.pickle','rb'))
    for sim in SimInfo:
        for halo in SimInfo[sim]['halos']:
            os.chdir(f'{gitdir}Images/{sim}.{feedback}/{halo}/')
            os.system(f'tar -czf {halo}.Images.tar.gz *0.png')
            if erase: os.system('rm *0.png')
            os.system(f'tar -czf {halo}.Isophotes.tar.gz *.Isophote.png')
            if erase: os.system('rm *.Isophote.png')
print('Done')