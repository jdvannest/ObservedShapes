import os,pickle

loop,rewrite = True,False
while loop:
    query = input('Remake Image Directory? (y/n): ')
    if query in ['y','Y']:
        rewrite=True
        loop=False
    elif query in ['n','N']:
        loop=False
    else: print('Invalid input.')


config = {
    #Path to directory where datafiles are written
    'output_path' : 'Data/',
    #Path to prefered python executable
    'python_path' : '/myhome2/users/vannest/anaconda3/bin/python'
}

out = open('Config.pickle','wb') 
pickle.dump(config,out)
out.close()


marvel_path = '/myhome2/users/munshi/dwarf_volumes/'
dcjl_path = '/myhome2/users/munshi/e12gals/'

Sims = {
    'cptmarvel' : {
        'path' : marvel_path+'cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,2,3,5,6,7,10,11,13,14,24]
    },
    'elektra' : {
        'path' : marvel_path+'elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,2,3,4,5,8,9,10,11,12,17,36,64]
    },
    'storm' : {
        'path' : marvel_path+'storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,2,3,4,5,6,7,8,10,11,12,14,15,22,23,31,37,44,48,55,118]
    },
    'rogue' : {
        'path' : marvel_path+'rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,3,7,8,10,11,12,15,16,17,28,31,37,58,116]
    },
    'h148' : {
        'path' : dcjl_path+'h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096',
        'halos' : [2,3,4,6,7,11,12,13,15,20,23,27,28,29,33,34,37,38,41,43,51,59,65,75,86,94,109,114,122]
    },
    'h229' : {
        'path' : dcjl_path+'h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096',
        'halos' : [2,3,6,14,15,18,20,22,25,33,47,48,49,52,57,62,89,92,127]
    },
    'h242' : {
        'path' : dcjl_path+'h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096',
        'halos' : [8,10,21,26,30,34,38,42,44,45,63,70,81,138]
    },
    'h329' : {
        'path' : dcjl_path+'h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096',
        'halos' : [7,29,30,37,53,92,115,117,127]
    }
}

out = open('SimulationInfo.BW.pickle','wb') 
pickle.dump(Sims,out)
out.close()

#Create Image Directory
if rewrite:
    os.system('rm -rf ../Images/*')
    for s in Sims:
        os.system(f'mkdir ../Images/{s}.BW')
        for h in Sims[s]['halos']: os.system(f'mkdir ../Images/{s}.BW/{h}')
        os.system(f'rm ../Data/{s}.BW.Images.pickle')
        os.system(f'rm ../Data/{s}.BW.Profiles.pickle')
        data = {}
        for h in Sims[s]['halos']: data[str(h)] = {}
        pickle.dump(data,open(f'../Data/{s}.BW.Images.pickle','wb'))
        pickle.dump(data,open(f'../Data/{s}.BW.Profiles.pickle','wb'))

Sims = {
    'storm' : {
        'path' : marvel_path+'storm.cosmo25cmb.4096g1HsbBH/storm.cosmo25cmb.4096g1HsbBH.004096',
        #'halos' : [1,2,3,4,5,6,7,8,10,11,12,14,15,18,21,23,35,42,48,49,61,88,125,262,272,273,472,541]
        'halos' : [1,2,3,4,5,6,7,8,10,11,12,14,15,18,21,23,35,48,49,61,88,125,133,175,235,262,272,300,541]#42 is weird, idk why
    }
}

out = open('SimulationInfo.SB.pickle','wb') 
pickle.dump(Sims,out)
out.close()

#Create Image Directory
if rewrite:
    for s in Sims:
        os.system(f'mkdir ../Images/{s}.SB')
        for h in Sims[s]['halos']: os.system(f'mkdir ../Images/{s}.SB/{h}')
        os.system(f'rm ../Data/{s}.SB.Images.pickle')
        os.system(f'rm ../Data/{s}.SB.Profiles.pickle')
        data = {}
        for h in Sims[s]['halos']: data[str(h)] = {}
        pickle.dump(data,open(f'../Data/{s}.SB.Images.pickle','wb'))
        pickle.dump(data,open(f'../Data/{s}.SB.Profiles.pickle','wb'))