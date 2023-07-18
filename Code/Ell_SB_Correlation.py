import pickle,sys,warnings
import numpy as np
from pandas import *
import scipy.stats as stats
import matplotlib.pylab as plt
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
warnings.filterwarnings("ignore")

def correlation(ell,sb):
    if not isinstance(ell,np.ndarray): ell = np.array(ell)
    if not isinstance(sb,np.ndarray): sb = np.array(sb)
    n = np.nansum( (ell-np.nanmean(ell)) * (sb-np.nanmean(sb)) )
    d1 = np.nansum( (ell-np.nanmean(ell))**2 )
    d2 = np.nansum( (sb-np.nanmean(sb))**2 )
    return n/np.sqrt(d1*d2)

def myround(x, base=5):
    return base * round(x/base)

def randomorientation(base=30):
    phi = np.degrees(np.random.uniform(0,np.pi*2))
    y = int(myround(phi,base=base))
    theta = np.degrees(np.arccos(np.random.uniform(-1,1)))
    x = int(myround(theta,base=base))
    if x>165:
        x = x%180
        y = (180-y)
        if y<0:
            y = y+360
    if y == 360:
        y = 0
    return(f'x{x:03d}y{y:03d}')


#Xu&Randall2020 (https://ui.adsabs.harvard.edu/abs/2020ApJ...900...69X/abstract)
LGall = [-0.287,0.058,'#7B0B01','All LG Dsphs']
LGbright = [-0.353,0.059,'#FFCC99','LG Dsphs\n'+r'$M/L<100M_\odot/L_\odot$']
LGdim = [0.509,0.309,'#339AFE','LG Dsphs\n'+r'$M/L>100M_\odot/L_\odot$']
FIRE = [-0.322,0.213,'#411B52','FIRE Simulations']
FIRE_Data = read_csv('../Data/BasicData/FIRE_Data.csv')
FIRE_x = FIRE_Data['X'].tolist()
FIRE_y = FIRE_Data['Y'].tolist()


ProfData,ProjData,MorphData,ShapeData,nhalo = {},{},{},{},0
for sim in ['h148','h229','h242','h329','cptmarvel','elektra','storm','rogue']:
    ProfData[sim]=pickle.load(open(f'../Data/{sim}.BW.Profiles.pickle','rb'))
    ProjData[sim]=pickle.load(open(f'../Data/{sim}.BW.ProjectedData.pickle','rb'))
    MorphData[sim]=pickle.load(open(f'../Data/{sim}.BW.3DShapes.pickle','rb'))
    ShapeData[sim]=pickle.load(open(f'../Data/{sim}.BW.ShapeData.pickle','rb'))
    nhalo+=len(ProfData[sim])
Masses = pickle.load(open('../Data/BasicData/Marvel_DCJL.Masses.pickle','rb'))
Magnitudes = pickle.load(open('../Data/BasicData/Marvel_DCJL.Magnitudes.pickle','rb'))



n_itr = int(1e4)
x = np.linspace(-1,1,500)
Xu = True
StellarMass = True
Morphology = True
Scatter,Nscat = True,10
PlotDict = ShapeData



#Xu Comparison
if Xu:
    r_dist,ind_r = np.zeros(n_itr),np.zeros(nhalo)
    print('Running: 0.00%')
    for n in np.arange(n_itr):
        ell,sb,i = np.zeros(nhalo),np.zeros(nhalo),0
        for sim in ProfData:
            for halo in ProfData[sim]:
                angle = randomorientation()
                ell[i] = 1-PlotDict[sim][halo][angle]['b/a']
                sb[i] = ProfData[sim][halo][angle]['lum_den'][0]
                if n==n_itr-1:
                    ell_ind,sb_ind = [],[]
                    for a in ProfData[sim][halo]:
                        ell_ind.append(1-PlotDict[sim][halo][a]['b/a'])
                        sb_ind.append(ProfData[sim][halo][a]['lum_den'][0])
                    ind_r[i] = correlation(ell_ind,sb_ind)
                i+=1

        r_dist[n] = correlation(ell,sb)
        myprint(f'Running: {round((n+1)/n_itr*100,2)}%',clear=True)


    #Faceon
    ell_f,sb_f,ell_s,sb_s,i = np.zeros(nhalo),np.zeros(nhalo),np.zeros(nhalo),np.zeros(nhalo),0
    for sim in ProfData:
        for halo in ProfData[sim]:
            ell_f[i] = 1-PlotDict[sim][halo]['x000y000']['b/a']
            sb_f[i] = ProfData[sim][halo]['x000y000']['lum_den'][0]
            ell_s[i] = 1-PlotDict[sim][halo]['x090y000']['b/a']
            sb_s[i] = ProfData[sim][halo]['x090y000']['lum_den'][0]
            i+=1
    r_face = correlation(ell_f,sb_f)
    r_side = correlation(ell_s,sb_s)


    f,ax = plt.subplots(1,1,figsize=(8,6))
    ax.set_xlabel(r'r$_{\epsilon\Sigma_*}$',fontsize=25)
    ax.tick_params(which='both',labelsize=15)
    ax.set_xlim([-1,1])
    ind_ax = ax.twiny()
    ind_ax.set_xlim([-1,1])
    ind_ax.set_xticks(np.ndarray.tolist(ind_r[np.isfinite(ind_r)]))
    ind_ax.set_xticklabels([])
    ind_ax.tick_params(axis='x',length=8,direction='in')

    ax.axvspan(FIRE[0]-FIRE[1],FIRE[0]+FIRE[1],color=FIRE[2],alpha=.3)
    ax.plot(FIRE_x,FIRE_y,color=FIRE[2],linewidth=2,label=FIRE[3])
    ax.fill_between(FIRE_x,[0]*len(FIRE_x),FIRE_y,color=FIRE[2],alpha=.5)
    for obs in [LGall,LGbright,LGdim]:
        ax.axvline(obs[0],color=obs[2],linewidth=2,label=obs[3])
        ax.axvspan(obs[0]-obs[1],obs[0]+obs[1],color=obs[2],alpha=.3)

    density = stats.gaussian_kde(r_dist)
    ax.plot(x,density(x),c='k',linewidth=3,label='Marvel+DCJL Simulations')
    ax.set_ylim(bottom=0)
    ax.legend(loc='lower right',prop={'size':12})
    #plt.show()
    f.savefig('../Images/CorrelationTesting/Correlation.png',bbox_inches='tight',pad_inches=.1)
    plt.close()



#PARAMETER TESTING

#Stellar Mass
if StellarMass:
    print('Testing Stellar Mass...')
    limits = ['6.0','6.5','7.0','7.5','8.0']
else:
    limits = []
for lim in limits:
    r_dist,ind_r = np.zeros(n_itr),np.zeros(nhalo)
    r_dist_low,r_dist_high = np.zeros(n_itr),np.zeros(n_itr)
    ind_r_low,ind_r_high = np.zeros(nhalo),np.zeros(nhalo)
    print(f'\t{lim}: 0.00%')
    for n in np.arange(n_itr):
        ell,sb,i = np.zeros(nhalo),np.zeros(nhalo),0
        ell_l,sb_l,ell_h,sb_h = np.zeros(nhalo),np.zeros(nhalo),np.zeros(nhalo),np.zeros(nhalo)
        for sim in ProfData:
            for halo in ProfData[sim]:
                angle = randomorientation()
                ell[i] = 1-PlotDict[sim][halo][angle]['b/a']
                sb[i] = ProfData[sim][halo][angle]['lum_den'][0]
                if np.log10(Masses[sim][halo]['Mstar'])<float(lim):
                    ell_l[i] = 1-PlotDict[sim][halo][angle]['b/a']
                    sb_l[i] = ProfData[sim][halo][angle]['lum_den'][0]
                    ell_h[i] = np.NaN
                    sb_h[i] = np.NaN
                else:
                    ell_l[i] = np.NaN
                    sb_l[i] = np.NaN
                    ell_h[i] = 1-PlotDict[sim][halo][angle]['b/a']
                    sb_h[i] = ProfData[sim][halo][angle]['lum_den'][0]
                if n==n_itr-1:
                    ell_ind,sb_ind = [],[]
                    ell_ind_l,sb_ind_l,ell_ind_h,sb_ind_h = [],[],[],[]
                    for a in ProfData[sim][halo]:
                        ell_ind.append(1-PlotDict[sim][halo][a]['b/a'])
                        sb_ind.append(ProfData[sim][halo][a]['lum_den'][0])
                        if np.log10(Masses[sim][halo]['Mstar'])<float(lim):
                            ell_ind_l.append(1-PlotDict[sim][halo][a]['b/a'])
                            sb_ind_l.append(ProfData[sim][halo][a]['lum_den'][0])
                            ell_ind_h.append(np.NaN)
                            sb_ind_h.append(np.NaN)
                        else:
                            ell_ind_l.append(np.NaN)
                            sb_ind_l.append(np.NaN)
                            ell_ind_h.append(1-PlotDict[sim][halo][a]['b/a'])
                            sb_ind_h.append(ProfData[sim][halo][a]['lum_den'][0])
                    ind_r[i]=correlation(ell_ind,sb_ind)
                    ind_r_low[i]=correlation(ell_ind_l,sb_ind_l)
                    ind_r_high[i]=correlation(ell_ind_h,sb_ind_h)
                i+=1

        r_dist[n] = correlation(ell,sb)
        r_dist_low[n] = correlation(ell_l,sb_l)
        r_dist_high[n] = correlation(ell_h,sb_h)
        myprint(f'\t{lim}: {round((n+1)/n_itr*100,2)}%',clear=True)
    ind_r = np.array(ind_r)


    f,ax = plt.subplots(1,1,figsize=(8,6))
    ax.set_xlabel(r'r$_{\epsilon\Sigma_*}$',fontsize=25)
    ax.tick_params(which='both',labelsize=15)
    ax.set_xlim([-1,1])
    ind_ax = ax.twiny()
    ind_ax.set_xlim([-1,1])
    ind_ax.set_xticks(np.ndarray.tolist(ind_r[np.isfinite(ind_r)]))
    ind_ax.set_xticklabels([])
    ind_ax.tick_params(axis='x',length=8,direction='in')
    density = stats.gaussian_kde(r_dist)
    ax.plot(x,density(x),c='k',linewidth=3)
    #ax.hist(r_dist,x,histtype='step',edgecolor='k',linewidth=3,density=True)
    ax.scatter(ind_r_low,[max(density(x))+.1]*nhalo,c='r',s=2**2)
    ax.scatter(ind_r_high,[max(density(x))+.1]*nhalo,c='b',s=2**2)
    density = stats.gaussian_kde(r_dist_low)
    ax.plot(x,density(x),c='r',label=r'log(M$_*$/M$_\odot$)$<$'+lim)
    #ax.hist(r_dist_low,x,histtype='step',edgecolor='r',label=r'log(M$_*$/M$_\odot$)$<$'+lim,density=True)
    density = stats.gaussian_kde(r_dist_high)
    ax.plot(x,density(x),c='b',label=r'log(M$_*$/M$_\odot$)$>$'+lim)
    #ax.hist(r_dist_high,x,histtype='step',edgecolor='b',label=r'log(M$_*$/M$_\odot$)$<$'+lim,density=True)

    ax.legend(loc='lower right',ncol=1,prop={'size':15})
    ax.set_ylim(bottom=0)
    f.savefig(f'../Images/CorrelationTesting/Correlation.Mass_{lim}_.png',bbox_inches='tight',pad_inches=.1)
    plt.close()


#3D Morphology
def Oblate(ba,ca):
    return True if ba>.5 and ca<.5 else False
if Morphology:
    print('Testing Morphology: 0.00%')
    r_dist,ind_r = np.zeros(n_itr),np.zeros(nhalo)
    r_dist_low,r_dist_high = np.zeros(n_itr),np.zeros(n_itr)
    ind_r_low,ind_r_high = np.zeros(nhalo),np.zeros(nhalo)
    for n in np.arange(n_itr):
        ell,sb,i = np.zeros(nhalo),np.zeros(nhalo),0
        ell_l,sb_l,ell_h,sb_h = np.zeros(nhalo),np.zeros(nhalo),np.zeros(nhalo),np.zeros(nhalo)
        for sim in ProfData:
            for halo in ProfData[sim]:
                angle = randomorientation()
                ell[i] = 1-PlotDict[sim][halo][angle]['b/a']
                sb[i] = ProfData[sim][halo][angle]['lum_den'][0]
                try:
                    Reff = ProfData[sim][halo][angle]['Reff']
                    ind_eff = np.argmin(abs(MorphData[sim][halo]['rbins']-Reff))
                    ba,ca = MorphData[sim][halo]['ba'][ind_eff],MorphData[sim][halo]['ca'][ind_eff]
                except:
                    ba,ca = np.NaN,np.NaN
                    #if n==n_itr-1: print(f'{sim} - {halo}')
                if not Oblate(ba,ca):
                    ell_l[i] = 1-PlotDict[sim][halo][angle]['b/a']
                    sb_l[i] = ProfData[sim][halo][angle]['lum_den'][0]
                    ell_h[i] = np.NaN
                    sb_h[i] = np.NaN
                else:
                    ell_l[i] = np.NaN
                    sb_l[i] = np.NaN
                    ell_h[i] = 1-PlotDict[sim][halo][angle]['b/a']
                    sb_h[i] = ProfData[sim][halo][angle]['lum_den'][0]
                if n==n_itr-1:
                    ell_ind,sb_ind = [],[]
                    ell_ind_l,sb_ind_l,ell_ind_h,sb_ind_h = [],[],[],[]
                    for a in ProfData[sim][halo]:
                        ell_ind.append(1-PlotDict[sim][halo][a]['b/a'])
                        sb_ind.append(ProfData[sim][halo][a]['lum_den'][0])
                        if not Oblate(ba,ca):
                            ell_ind_l.append(1-PlotDict[sim][halo][a]['b/a'])
                            sb_ind_l.append(ProfData[sim][halo][a]['lum_den'][0])
                            ell_ind_h.append(np.NaN)
                            sb_ind_h.append(np.NaN)
                        else:
                            ell_ind_l.append(np.NaN)
                            sb_ind_l.append(np.NaN)
                            ell_ind_h.append(1-PlotDict[sim][halo][a]['b/a'])
                            sb_ind_h.append(ProfData[sim][halo][a]['lum_den'][0])
                    ind_r[i]=correlation(ell_ind,sb_ind)
                    ind_r_low[i]=correlation(ell_ind_l,sb_ind_l)
                    ind_r_high[i]=correlation(ell_ind_h,sb_ind_h)
                i+=1

        r_dist[n] = correlation(ell,sb)
        r_dist_low[n] = correlation(ell_l,sb_l)
        r_dist_high[n] = correlation(ell_h,sb_h)
        myprint(f'Testing Morphology: {round((n+1)/n_itr*100,2)}%',clear=True)
    ind_r = np.array(ind_r)


    f,ax = plt.subplots(1,1,figsize=(8,6))
    ax.set_xlabel(r'r$_{\epsilon\Sigma_*}$',fontsize=25)
    ax.tick_params(which='both',labelsize=15)
    ax.set_xlim([-1,1])
    ind_ax = ax.twiny()
    ind_ax.set_xlim([-1,1])
    ind_ax.set_xticks(np.ndarray.tolist(ind_r[np.isfinite(ind_r)]))
    ind_ax.set_xticklabels([])
    ind_ax.tick_params(axis='x',length=8,direction='in')
    density = stats.gaussian_kde(r_dist)
    ax.plot(x,density(x),c='k',linewidth=3)
    #ax.hist(r_dist,x,histtype='step',edgecolor='k',linewidth=3,density=True)
    ax.scatter(ind_r_low,[max(density(x))+.1]*nhalo,c='r',s=2**2)
    ax.scatter(ind_r_high,[max(density(x))+.1]*nhalo,c='b',s=2**2)
    density = stats.gaussian_kde(r_dist_low)
    ax.plot(x,density(x),c='r',label='Non-Oblate')
    #ax.hist(r_dist_low,x,histtype='step',edgecolor='r',label=r'log(M$_*$/M$_\odot$)$<$'+lim,density=True)
    density = stats.gaussian_kde(r_dist_high)
    ax.plot(x,density(x),c='b',label='Oblate')
    #ax.hist(r_dist_high,x,histtype='step',edgecolor='b',label=r'log(M$_*$/M$_\odot$)$<$'+lim,density=True)

    ax.legend(loc='lower right',ncol=1,prop={'size':15})
    ax.set_ylim(bottom=0)
    f.savefig(f'../Images/CorrelationTesting/Correlation.Morphology.png',bbox_inches='tight',pad_inches=.1)
    plt.close()


if Scatter:
    IndR = {}
    for sim in ProfData:
        IndR[sim] = {}
        for halo in ProfData[sim]:
            ell,sb = [],[]
            for a in ProfData[sim][halo]:
                ell.append(1-PlotDict[sim][halo][a]['b/a'])
                sb.append(ProfData[sim][halo][a]['lum_den'][0])
            IndR[sim][halo] = correlation(ell,sb)
    for i in np.arange(Nscat):
        f,ax = plt.subplots(1,1,figsize=(8,6))
        ax.set_xlabel(r'Surface Brightness [log(L$_\odot$/pc$^2$)]',fontsize=25)
        ax.set_ylabel(r'Ellipticity [1-$b/a$]',fontsize=25)
        ax.tick_params(which='both',labelsize=15)
        ax.set_ylim([0,1])
        ax.set_xlim([-4,1.2])

        ell,sb,j = np.zeros(nhalo),np.zeros(nhalo),0
        for sim in ProfData:
            for halo in ProfData[sim]:
                angle = randomorientation()
                E = 1-PlotDict[sim][halo][angle]['b/a']
                S = ProfData[sim][halo][angle]['lum_den'][0]
                ell[j],sb[j] = E,S
                R = IndR[sim][halo] 
                t = np.arctan(R)
                xdev = (.15-.1*abs(R))*np.cos(t)
                ydev = (.15-.1*abs(R))*np.sin(t)
                ax.plot([np.log10(S)-xdev,np.log10(S)+xdev],[E-ydev,E+ydev],c='.5',zorder=0)
                ax.scatter(np.log10(S),E,c='k',zorder=2)
                j+=1

        #Wheel
        cm = plt.get_cmap('bwr_r')
        for k in np.linspace(-1,1,100):
            xdev = (.15-.1212*abs(k))*np.cos(np.arctan(k))
            ydev = (.15-.1212*abs(k))*np.sin(np.arctan(k))
            ax.plot([.9,.9+xdev],[.9,.9+ydev],c=cm((k+1)/2),zorder=0)
        ax.scatter(.9,.9,c='k',zorder=5)
        ax.text(.85,.86,'prolate',ha='right',c='r',fontsize=15)
        ax.text(.85,.92,'oblate',ha='right',c='b',fontsize=15)

        r_dist = correlation(ell,sb)
        ax.set_title(r'r$_{\epsilon\Sigma_*}$= '+str(round(r_dist,3)),fontsize=25)
        f.savefig(f'../Images/CorrelationTesting/Scatter.{i+1}.png',bbox_inches='tight',pad_inches=.1)
        plt.close()