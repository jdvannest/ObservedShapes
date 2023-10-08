import pickle
import numpy as np
import matplotlib.pylab as plt

SimInfo = pickle.load(open('SimulationInfo.BW.pickle','rb'))
mass_data = pickle.load(open('../Data/BasicData/Marvel_DCJL.Masses.pickle','rb'))
with open('../Data/BasicData/HaloTypes.txt') as f:
    halotype = f.readlines()
    del halotype[0]
types = {}
for line in halotype:
    l = line.split('\t')
    if l[0] not in types: types[l[0]] = {}
    types[l[0]][l[1]] = l[-2]



def T(ba,ca):
    return( (1-ba**2)/(1-ca**2) )
T_s,T_d = [],[]
B_s,C_s,B_d,C_d = [],[],[],[]
masses,htype = [],[]

for sim in SimInfo:
    StShapes = pickle.load(open(f'../Data/{sim}.BW.3DShapes.pickle','rb'))
    DMShapes = pickle.load(open(f'../Data/{sim}.BW.DMShapes.pickle','rb'))
    Profiles = pickle.load(open(f'../Data/{sim}.BW.Profiles.pickle','rb'))

    for hid in SimInfo[sim]['goodhalos']:
        try:
            rbins,rd,ba_s,ca_s,ba_d,ca_d = [StShapes[str(hid)]['rbins'],DMShapes[str(hid)]['rbins'],
                                            StShapes[str(hid)]['ba_smooth'],StShapes[str(hid)]['ca_smooth'],
                                            DMShapes[str(hid)]['ba_smooth'],DMShapes[str(hid)]['ca_smooth']]
            Reff = Profiles[str(hid)]['x000y000']['Reff']
            #indeff = np.argmin(np.abs(rd-Reff))
            B_s.append(ba_s(Reff))
            C_s.append(ca_s(Reff))
            T_s.append(T(ba_s(Reff),ca_s(Reff)))
            B_d.append(ba_d(Reff))
            C_d.append(ca_d(Reff))
            T_d.append(T(ba_d(Reff),ca_d(Reff)))
            masses.append(np.log10(mass_data[sim][str(hid)]['Mstar']))
        except:
            continue
        if sim not in types:
            htype.append('o')
        elif str(hid) not in types[sim]:
            htype.append('o')
        else:
            if types[sim][str(hid)] in ['Central','Backsplash']: htype.append('o')
            else: htype.append('v')
        
        if False:
            f,ax=plt.subplots(1,1,figsize=(6,2))
            ax.set_xlim([0,max(rbins)])
            ax.set_ylim([0,1])
            ax.set_xlabel('R [kpc]',fontsize=15)
            ax.set_ylabel('Axis Ratio',fontsize=15)
            ax.tick_params(which='both',labelsize=10)
            ax.plot([-1,-1],[-1,-1],c='.5',label='B/A')
            ax.plot([-1,-1],[-1,-1],c='.5',linestyle='--',label='C/A')

            ax.plot(rbins,ba_d,c='k',label='Dark Matter')
            ax.plot(rbins,ca_d,c='k',linestyle='--')

            ax.plot(rbins,ba_s(rbins),c='r',label='Stars')
            ax.plot(rbins,ca_s(rbins),c='r',linestyle='--')

            ax.legend(loc='lower right',prop={'size':12},ncol=2)
            f.savefig(f'../Images/3DShapes/{sim}.BW/{hid}.Dark.png',bbox_inches='tight',pad_inches=.1)
            plt.close()

#Convert lists to arrays for htype indexing
T_d,T_s,masses,htype = np.array(T_d),np.array(T_s),np.array(masses),np.array(htype)
B_s,C_s,B_d,C_d = np.array(B_s),np.array(C_s),np.array(B_d),np.array(C_d)

#T* vs Tdm
f,ax = plt.subplots(1,1,figsize=(5,5))
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.fill_between([0,1],[-.1,.9],[.1,1.1],color='0.75',alpha=.3)
ax.plot([0,1],[0,1],c='0.5',linestyle='--')
ax.set_ylabel(r'T$_*$',fontsize=15)
ax.set_xlabel(r'T$_{DM}$',fontsize=15)

norm = plt.Normalize(int(min(masses)),int(max(masses))+.1)
p = ax.scatter(T_d[htype=='o'],T_s[htype=='o'],marker='o',c=masses[htype=='o'],cmap='viridis',norm=norm)
ax.scatter(T_d[htype=='v'],T_s[htype=='v'],marker='v',c=masses[htype=='v'],cmap='viridis',norm=norm,label='Satellites')
cbar = f.colorbar(p,cax=f.add_axes([.91,.11,.03,.77]))
cbar.set_label(r'Log(M$_*$/M$_\odot$)',fontsize=15)

ax.legend(loc='lower left',prop={'size':12})
f.savefig(f'../Images/3DShapes/T_Comparison.png',bbox_inches='tight',pad_inches=.1)


#T vs Mstar
f,ax=plt.subplots(1,1,figsize=(12,4))
ax.set_xlim([5.8,9.5])
ax.set_ylim([0,1])
ax.set_yticks([])
ax.set_xlabel(r'Log(M$_*$/M$_\odot$)',fontsize=25)
ax.set_ylabel('T',fontsize=25)
ax.plot([4,9.5],[1/3,1/3],c='.75',linestyle='--',zorder=0)
ax.plot([4,9.5],[2/3,2/3],c='.75',linestyle='--',zorder=0)
ax.tick_params(which='both',labelsize=15)
ax.text(5.83,1/6,'Oblate',fontsize=18,rotation='vertical',verticalalignment='center',c='.5')
ax.text(5.83,3/6,'Triaxial',fontsize=18,rotation='vertical',verticalalignment='center',c='.5')
ax.text(5.83,5/6,'Prolate',fontsize=18,rotation='vertical',verticalalignment='center',c='.5')

for i in np.arange(len(masses)):
    ax.axvline(masses[i],ymin=min([T_d[i],T_s[i]]),ymax=max([T_d[i],T_s[i]]),c='.5',zorder=0)
ax.scatter(masses[htype=='o'],T_d[htype=='o'],c='k',label='Dark Matter',marker='o')
ax.scatter(masses[htype=='o'],T_s[htype=='o'],c='r',label='Stellar',marker='o')
ax.scatter(masses[htype=='v'],T_d[htype=='v'],c='k',marker='v')
ax.scatter(masses[htype=='v'],T_s[htype=='v'],c='r',marker='v')

ax.scatter(0,0,c='.5',marker='v',label='Satellites')
ax.legend(prop={'size':15},ncol=3,loc='center left', bbox_to_anchor=(0.02,0.1))
f.savefig(f'../Images/3DShapes/TvsMstar.png',bbox_inches='tight',pad_inches=.1)



#B/A vs C/A links
f,ax = plt.subplots(1,1,figsize=(5,5))
ax.set_xlim([0,1])
ax.set_ylim([0,1])
#ax.fill_between([0,1],[-.1,.9],[.1,1.1],color='0.75',alpha=.3)
ax.plot([0,1],[0,1],c='0.5',linestyle='--')
ax.set_xlabel(r'$S$',fontsize=20)
ax.set_ylabel(r'$Q$',fontsize=20)
ax.tick_params(which='both',labelsize=15)
ax.scatter(-1,-1)

for i in np.arange(len(B_s)):
    ax.plot([B_s[i],B_d[i]],[C_s[i],C_d[i]],c='.5',zorder=0)

f.savefig(f'../Images/3DShapes/CvB.LinksOnly.png',bbox_inches='tight',pad_inches=.1)

ax.scatter(B_d[htype=='o'],C_d[htype=='o'],c='k',label='Dark Matter')
ax.scatter(B_d[htype=='v'],C_d[htype=='v'],c='k',marker='v')
ax.scatter(B_s[htype=='o'],C_s[htype=='o'],c='r',label='Stellar')
ax.scatter(B_s[htype=='v'],C_s[htype=='v'],c='r',marker='v')
ax.scatter(-1,-1,c='.5',marker='v',label='Satellites')

ax.legend(loc='upper left',prop={'size':15})
f.savefig(f'../Images/3DShapes/CvB.Links.png',bbox_inches='tight',pad_inches=.1)