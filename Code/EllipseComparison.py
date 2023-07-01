import pickle
import numpy as np
import matplotlib.pylab as plt

mass_data = pickle.load(open('../Data/BasicData/Marvel_DCJL.Masses.pickle','rb'))

for xrot in np.arange(0,180,30):
    for yrot in np.arange(0,360,30):

        projected,observed,col = [],[],[]
        f,ax = plt.subplots(1,1,figsize=(5,5))
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        ax.fill_between([0,1],[-.1,.9],[.1,1.1],color='0.75',alpha=.3)
        ax.plot([0,1],[0,1],c='0.5',linestyle='--')
        ax.set_xlabel('Projected b/a',fontsize=15)
        ax.set_ylabel('Isophote b/a',fontsize=15)

        for fb in ['BW']:#,'SB']:
            for sim in ['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329']:
                obs_data = pickle.load(open(f'../Data/{sim}.{fb}.ShapeData.pickle','rb'))
                #obs_mask = pickle.load(open(f'../Data/{sim}.{fb}.Masking.pickle','rb'))
                prj_data = pickle.load(open(f'../Data/{sim}.{fb}.ProjectedData.pickle','rb'))

                for halo in obs_data:
                    col.append(np.log10(mass_data[sim][halo]['Mstar']))
                    pr = prj_data[halo][f'x{xrot:03d}y{yrot:03d}']['b/a']
                    if pr<=.1: pr = np.NaN
                    projected.append(pr)
                    try:
                        observed.append(obs_data[halo][f'x{xrot:03d}y{yrot:03d}']['b/a'])
                    except:
                        observed.append(np.NaN)
                        print(f'{sim} - {halo}')

        norm = plt.Normalize(int(min(col)),int(max(col))+1)
        p = ax.scatter(projected,observed,c=col,cmap='viridis',norm=norm)
        cbar = f.colorbar(p,cax=f.add_axes([.91,.11,.03,.77]))
        cbar.set_label(r'Log(M$_*$/M$_\odot$)',fontsize=15)
        f.savefig(f'../Images/EllipseComparison/EllipseComparison.x{xrot:03d}.y{yrot:03d}.png',bbox_inches='tight',pad_inches=.1)
        plt.close()