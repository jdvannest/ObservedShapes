import os,pickle,sys
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.interpolate import RegularGridInterpolator as RGI
sys.path.append('../')
from scipy.interpolate import griddata
from ShapeFunctions.Functions import RegularGridInterpolator as RGI
from ShapeFunctions.Functions import Ellipsoid,Project,kpc2pix,pix2kpc,Project_OLD
from matplotlib.patches import Rectangle as Rec
from matplotlib.collections import PatchCollection

def PadData(data):
    '''
    RegularGridInterpolator has hard boundar cuts at data edges.
    Data covers (theta,phi) from (0,0) to (150,330) but we want interpolations
    over (0,0) to (180,360).
    This Function pads the data with the edge points through symmetry Phi(0)=Phi(360), etc.
    '''
    xdim,ydim = data.shape
    padded = np.zeros((xdim+1,ydim+1))
    for x in np.arange(xdim+1):
        for y in np.arange(ydim):
            if x==xdim:
                if y<=ydim/2:
                    padded[x][y] = data[0][int(ydim/2-y)]
                else:
                    padded[x][y] = data[0][int(ydim-(y-ydim/2))]
            else:
                padded[x][y] = data[x][y]
    padded[:,-1] = padded[:,0]
    return padded

def FillNaN(data):
    '''
    RegularGridInterpolator doesn't accept NaN's, so a more shallow interpolator 
    will fill in the Major Grid points for galaxies with a few holes in data
    '''
    good,bad,fit = [],[],[]
    xdim,ydim = data.shape
    for x in np.arange(xdim):
        for y in np.arange(ydim):
            if np.isnan(data[x][y]):
                bad.append([x*30,y*30])
            else:
                good.append([x*30,y*30])
                fit.append(data[x][y])
    fill = griddata(good,fit,bad,method='cubic')
    for i in np.arange(len(bad)):
        x,y = int(bad[i][0]/30),int(bad[i][1]/30)
        data[x][y] = fill[i]
    return data,bad,fill


norm = plt.Normalize(0,1)
norm2 = plt.Normalize(-1,1)
norm3 = plt.Normalize(0,.75)
#parameters for raw data plots
xang = np.arange(0,180,30)
yang = np.arange(0,360,30)
Yedge = np.arange(-15,185,30)
Xedge = np.arange(-15,375,30)
#paraameters for padded data plots
xang_pad = np.arange(0,210,30)
yang_pad = np.arange(0,390,30)
Yedge_pad = np.arange(-15,225,30)
Xedge_pad = np.arange(-15,405,30)
#paarameters for interpolated data plots
da = 10
xang_int = np.arange(0,180+da,da)
yang_int = np.arange(0,360+da,da)
Yedge_int = np.arange(-da/2,180+da*1.5,da)
Xedge_int = np.arange(-da/2,360+da*1.5,da)

Plot = False
smooth = {
    'cptmarvel':[1,2,3,5,6,7,10,11,13],
    'elektra':[2,3,4,5,9,10,12,17,36],
    'storm':[2,3,4,5,6,7,8,14,15,23,31,44],
    'rogue':[7,8,10,11,12,15,16,17,28,31],
    'h148':[3,4,6,7,11,12,13,15,20,23,27,28,33,34,38,65,114],
    'h229':[2,3,6,14,18,22,49,92],
    'h242':[8,10,21,26,30,34,38],
    'h329':[7,29,115,117]
}

SimInfo = pickle.load(open('SimulationInfo.BW.pickle','rb'))
for sim in SimInfo:
    Shapes = pickle.load(open(f'../Data/{sim}.BW.ShapeData.pickle','rb'))
    Profiles = pickle.load(open(f'../Data/{sim}.BW.Profiles.pickle','rb'))
    Shapes3d = pickle.load(open(f'../Data/{sim}.BW.3DShapes.pickle','rb'))
    InterpFuncs = {}
    for hid in Shapes:
        InterpFuncs[hid] = np.NaN
        plot = np.zeros((len(xang),len(yang)))
        for x in np.arange(len(xang)):
            for y in np.arange(len(yang)):
                plot[x][y] = Shapes[hid][f'x{x*30:03d}y{y*30:03d}']['b/a']

        #Raw Data Plot
        f,ax = plt.subplots(1,1,figsize=(12,8))
        ax.tick_params(labelsize=15)
        ax.set_xticks(np.arange(0,375,30))
        ax.set_yticks(np.arange(0,195,30))
        ax.set_xlim([-15,345])
        ax.set_ylim([-15,165])
        ax.set_xlabel(r'$\phi$-rotation [$^o$]',fontsize=25)
        ax.set_ylabel(r'$\theta$-rotation [$^o$]',fontsize=25)

        p = ax.pcolor(Xedge,Yedge,plot,cmap='viridis',edgecolors='k',norm=norm)
        c = f.colorbar(p,ax=ax,pad=0.01,aspect=15)
        c.set_label('b/a',fontsize=25)
        c.ax.tick_params(labelsize=15)
        c.ax.plot([0,1],[np.amin(plot),np.amin(plot)],c='w') 
        c.ax.plot([0,1],[np.amax(plot),np.amax(plot)],c='w') 

        f.savefig(f'../Images/Interpolation/{sim}.BW/{hid}.Data.png',bbox_inches='tight',pad_inches=.1)
        plt.close(f)

        #Filled and Padded Data Plot
        f,ax = plt.subplots(1,1,figsize=(12,8))
        ax.tick_params(labelsize=15)
        ax.set_xticks(np.arange(0,405,30))
        ax.set_yticks(np.arange(0,225,30))
        ax.set_xlim([-15,375])
        ax.set_ylim([-15,195])
        ax.set_xlabel(r'$\phi$-rotation [$^o$]',fontsize=25)
        ax.set_ylabel(r'$\theta$-rotation [$^o$]',fontsize=25)

        plot_fill1,bad1,fill1 = FillNaN(plot)
        plot_pad = PadData(plot_fill1)
        plot_fill2,bad2,fill2 = FillNaN(plot_pad)

        p = ax.pcolor(Xedge_pad,Yedge_pad,plot_fill2,cmap='viridis',edgecolors='k',norm=norm)
        c = f.colorbar(p,ax=ax,pad=0.01,aspect=15)
        c.set_label('b/a',fontsize=25)
        c.ax.tick_params(labelsize=15)
        c.ax.plot([0,1],[np.amin(plot),np.amin(plot)],c='w') 
        c.ax.plot([0,1],[np.amax(plot),np.amax(plot)],c='w') 

        if len(bad1)>0:
            for i in np.arange(len(bad1)):
                ax.add_patch(Rec((bad1[i][1]-5,bad1[i][0]-5),10,10,color='w'))
        if len(bad2)>0:
            for i in np.arange(len(bad2)):
                ax.add_patch(Rec((bad2[i][1]-5,bad2[i][0]-5),10,10,color='w'))
        
        f.savefig(f'../Images/Interpolation/{sim}.BW/{hid}.Filled.png',bbox_inches='tight',pad_inches=.1)
        plt.close(f)

        if not True in np.isnan(plot_fill2):
            try:
                interpolate = RGI(points=(xang_pad,yang_pad),values=plot_fill2,method='cubic')
                InterpFuncs[hid] = interpolate

                if not Plot: continue

                plot_dif,plot_sig = np.zeros((len(xang_int),len(yang_int))),np.zeros((len(xang_int),len(yang_int)))
                for xrot in np.arange(len(xang_int)):
                    for yrot in np.arange(len(yang_int)):
                        #Get Projected b/a
                        reffx = 30*(xrot//30) if xrot<180 else 0
                        reffy = 30*(yrot//30) if yrot<360 else 0
                        reffkey = f'x{reffx:03d}y{reffy:03d}'
                        Rhalf = Profiles[str(hid)][reffkey]['Reff']
                        if np.isnan(Rhalf): Rhalf = Profiles[str(hid)][reffkey]['Rhalf']
                        rbins = Shapes3d[str(hid)]['rbins']
                        ind_eff = np.argmin(abs(rbins-Rhalf))
                        a,ba,ca,Es = [Shapes3d[str(hid)]['a'][ind_eff],Shapes3d[str(hid)]['ba'][ind_eff],
                                        Shapes3d[str(hid)]['ca'][ind_eff],Shapes3d[str(hid)]['Es'][ind_eff]]
                        if int(hid) in smooth[sim]:
                            ba = float(Shapes3d[str(hid)]['ba_smooth'](Rhalf))
                            ca = float(Shapes3d[str(hid)]['ca_smooth'](Rhalf))
                        x,y,z = Ellipsoid(a,ba,ca,Es,xrot*da,yrot*da)
                        x = -x
                        ap,bp,cen,phi = Project_OLD(x,y,z)
                        atrue,btrue = max([ap,bp]),min([ap,bp])

                        #Plot array for interpolated - projected
                        plot_dif[xrot][yrot] = interpolate((xrot*da,yrot*da)) - btrue/atrue
                        #Plot array for distance from 1-1 line
                        plot_sig[xrot][yrot] = abs(interpolate((xrot*da,yrot*da)) - btrue/atrue)/np.sqrt(2)

                f,ax = plt.subplots(1,1,figsize=(12,8))
                ax.tick_params(labelsize=15)
                ax.set_xticks(np.arange(0,360+2*da,2*da))
                ax.set_yticks(np.arange(0,180+2*da,2*da))
                ax.set_xlim([-da/2,360+da/2])
                ax.set_ylim([-da/2,180+da/2])
                ax.set_xlabel(r'$\phi$-rotation [$^o$]',fontsize=25)
                ax.set_ylabel(r'$\theta$-rotation [$^o$]',fontsize=25)

                p = ax.pcolor(Xedge_int,Yedge_int,plot_dif,cmap='seismic',edgecolors='k',norm=norm2)
                c = f.colorbar(p,ax=ax,pad=0.01,aspect=15)
                c.set_label(r'$\Delta$(Interp. - Proj.)',fontsize=25)
                c.ax.tick_params(labelsize=15)
                c.ax.plot([-1,1],[np.amin(plot_dif),np.amin(plot_dif)],c='gold') 
                c.ax.plot([-1,1],[np.amax(plot_dif),np.amax(plot_dif)],c='gold') 

                f.savefig(f'../Images/Interpolation/{sim}.BW/{hid}.Interpolated.Dif.png',bbox_inches='tight',pad_inches=.1)
                plt.close(f)


                f,ax = plt.subplots(1,1,figsize=(12,8))
                ax.tick_params(labelsize=15)
                ax.set_xticks(np.arange(0,360+2*da,2*da))
                ax.set_yticks(np.arange(0,180+2*da,2*da))
                ax.set_xlim([-da/2,360+da/2])
                ax.set_ylim([-da/2,180+da/2])
                ax.set_xlabel(r'$\phi$-rotation [$^o$]',fontsize=25)
                ax.set_ylabel(r'$\theta$-rotation [$^o$]',fontsize=25)

                p = ax.pcolor(Xedge_int,Yedge_int,plot_sig,cmap='viridis',edgecolors='k',norm=norm3)
                c = f.colorbar(p,ax=ax,pad=0.01,aspect=15)
                c.set_label(r'$\sigma_{1-1}$',fontsize=25)
                c.ax.tick_params(labelsize=15)
                c.ax.plot([-1,1],[np.amin(plot_sig),np.amin(plot_sig)],c='w') 
                c.ax.plot([-1,1],[np.amax(plot_sig),np.amax(plot_sig)],c='w') 

                f.savefig(f'../Images/Interpolation/{sim}.BW/{hid}.Interpolated.Sig.png',bbox_inches='tight',pad_inches=.1)
                plt.close(f)
            except:
                continue
    pickle.dump(InterpFuncs,open(f'../Data/{sim}.BW.InterpolationFunctions.pickle','wb'))