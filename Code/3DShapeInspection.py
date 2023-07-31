import os,pickle
import matplotlib.pylab as plt

SimInfo = pickle.load(open('SimulationInfo.BW.pickle','rb'))

# loop,remake=True,False
# while loop:
#     rem = input('Remake Image Directory: Images/3DShapes/- (y/n): ')
#     if rem in ['y','n']:
#         loop = False
#         if rem=='y': remake = True
# if remake: 
#     os.system('rmdir -f ../Images/3DShapes')
#     os.system('mkdir ../Images/3DShapes')
#     for sim in SimInfo:
#         os.system(f'mkdir ../Images/3DShapes/{sim}.BW/')
    

for sim in SimInfo:
    Shapes = pickle.load(open(f'../Data/{sim}.BW.3DShapes.pickle','rb'))
    Profiles = pickle.load(open(f'../Data/{sim}.BW.Profiles.pickle','rb'))

    for hid in SimInfo[sim]['halos']:
        rbins,ba,ca=Shapes[str(hid)]['rbins'],Shapes[str(hid)]['ba'], Shapes[str(hid)]['ca']
        reffs = []
        if len(rbins)>0:
            for angle in Profiles[str(hid)]:
                reffs.append(Profiles[str(hid)][angle]['Reff'])
            
            f,ax=plt.subplots(1,1,figsize=(6,2))
            ax.set_xlim([0,max(rbins)])
            ax.set_ylim([0,1])
            ax.set_xlabel('R [kpc]',fontsize=15)
            ax.set_ylabel('Axis Ratio',fontsize=15)
            ax.tick_params(which='both',labelsize=10)

            ax.plot(rbins,ba,c='k',label='B/A')
            ax.plot(rbins,ca,c='k',linestyle='--',label='C/A')
            #ax.axvspan(min(reffs),max(reffs),color='k',alpha=0.3,label=r'R$_{eff}$')
            for reff in reffs:
                ax.axvline(reff,c='k',alpha=.05)
            ax.axvline(reffs[0],c='k',alpha=.05,label=r'R$_{eff}$')

            ax.legend(loc='lower right',prop={'size':12})
            f.savefig(f'../Images/3DShapes/{sim}.BW/{hid}.png',bbox_inches='tight',pad_inches=.1)
            plt.close()