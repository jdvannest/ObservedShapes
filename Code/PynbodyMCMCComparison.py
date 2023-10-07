import pickle
import numpy as np
import matplotlib.pylab as plt

pB,pC,mB,mC = [],[],[],[]

SimInfo = pickle.load(open(f'SimulationInfo.BW.pickle','rb'))

for sim in SimInfo:
    prof = pickle.load(open(f'../Data/{sim}.BW.Profiles.pickle','rb'))
    pynb = pickle.load(open(f'../Data/{sim}.BW.3DShapes.pickle','rb'))
    mcmc = pickle.load(open(f'../Data/{sim}.BW.MCMC.pickle','rb'))

    for halo in SimInfo[sim]['goodhalos']:
        Reff = prof[str(halo)]['x000y000']['Rhalf']
        try:    
            ind_eff = np.argmin(abs(pynb[str(halo)]['rbins']-Reff))
            pB.append(pynb[str(halo)]['ba'][ind_eff])
            pC.append(pynb[str(halo)]['ca'][ind_eff])
        except:
            pB.append(np.NaN)
            pC.append(np.NaN)
        try:
            mB.append(mcmc[str(halo)]['Isophote']['alpha_50'][0])
            mC.append(mcmc[str(halo)]['Isophote']['alpha_50'][1])
        except:
            mB.append(np.NaN)
            mC.append(np.NaN)

f,ax = plt.subplot_mosaic([['B','B','B','.'],
                            ['S','S','S','C'],
                            ['S','S','S','C'],
                            ['S','S','S','C']],figsize=(7,7))
plt.subplots_adjust(hspace=0,wspace=0)
ax['S'].set_xlabel('B/A',fontsize=15)
ax['S'].set_ylabel('C/A',fontsize=15)
ax['S'].tick_params(which='both',labelsize=10)
ax['S'].plot([0,1],[0,1],c='k',linestyle='--')
for p in ['B','C','S']:
    ax[p].set_xlim([0,1])
    ax[p].set_ylim([0,1])
    if not p=='S':
        ax[p].set_xticks([])
        ax[p].set_yticks([])
ax['B'].set_ylim([0,4])
ax['C'].set_xlim([0,6])

ax['S'].scatter(pB,pC,c='b',label='Pynbody')
ax['S'].scatter(mB,mC,c='r',label='MCMC')

bins = np.linspace(0,1,26)
ax['B'].hist(pB,bins,histtype='step',edgecolor='b',density=True)
ax['B'].hist(mB,bins,histtype='step',edgecolor='r',density=True)
ax['C'].hist(pC,bins,histtype='step',edgecolor='b',density=True,orientation="horizontal")
ax['C'].hist(mC,bins,histtype='step',edgecolor='r',density=True,orientation="horizontal")

ax['S'].legend(loc='upper left',prop={'size':12})
f.savefig(f'../Images/MCMC/PynbodyMCMCComparison.png',bbox_inches='tight',pad_inches=.1)