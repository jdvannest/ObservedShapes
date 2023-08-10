import os,pickle,sys
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.interpolate import RegularGridInterpolator as RGI
sys.path.append('../')
from scipy.interpolate import griddata
from ShapeFunctions.Functions import RegularGridInterpolator as RGI
from matplotlib.patches import Rectangle as Rec
from matplotlib.collections import PatchCollection

SimInfo = pickle.load(open('../SimulationInfo.BW.pickle','rb'))

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
    
LowResFull = pickle.load(open(f'../../Data/storm.BW.ShapeData.pickle','rb'))
LowRes = LowResFull['3']
HighRes = pickle.load(open('Storm.3.ShapeData.pickle','rb'))
xang_h,yang_h = np.arange(90,160,10),np.arange(150,220,10)
xang_l,yang_l = np.arange(0,180,30),np.arange(0,360,30)

plot_h = np.zeros((len(xang_h),len(yang_h)))
plot_l = np.zeros((len(xang_l),len(yang_l)))
for x in np.arange(len(xang_l)):
    for y in np.arange(len(yang_l)):
        plot_l[x][y] = LowRes[f'x{x*30:03d}y{y*30:03d}']['b/a']
for x in np.arange(len(xang_h)):
    for y in np.arange(len(yang_h)):
        xa,ya = xang_h[x],yang_h[y]
        plot_h[x][y] = HighRes[f'x{xa:03d}y{ya:03d}']['b/a']

inter = RGI(points=(xang_l,yang_l),values=plot_l,method='cubic')

# xang,yang = np.arange(90,151),np.arange(150,211)
# X,Y = np.meshgrid(xang,yang)
# pts = zip(X.ravel(),Y.ravel())

diffplot = np.zeros((len(xang_h),len(yang_h)))
err = []
for x in np.arange(len(xang_h)):
    for y in np.arange(len(yang_h)):
        diffplot[x][y] = plot_h[x][y] - inter((xang_h[x],yang_h[y]))
        err.append(plot_h[x][y] - inter((xang_h[x],yang_h[y])))
Xedge,Yedge = np.arange(145,225,10),np.arange(85,165,10)

f,ax = plt.subplots(1,1,figsize=(12,8))
ax.tick_params(labelsize=15)
#ax.set_xticks(xang_h)
#ax.set_yticks(yang_h)
ax.set_xlim([Xedge[0],Xedge[-1]])
ax.set_ylim([Yedge[0],Yedge[-1]])
ax.set_xlabel(r'$\phi$-rotation [$^o$]',fontsize=25)
ax.set_ylabel(r'$\theta$-rotation [$^o$]',fontsize=25)

norm = plt.Normalize(-1,1)
p = ax.pcolor(Xedge,Yedge,diffplot,cmap='seismic',edgecolors='k',norm=norm)
c = f.colorbar(p,ax=ax,pad=0.01,aspect=15)
c.set_label(r'$\Delta$(Data-Interpolation)',fontsize=25)
c.ax.tick_params(labelsize=15)
c.ax.plot([-1,1],[np.amin(diffplot),np.amin(diffplot)],c='k') 
c.ax.plot([-1,1],[np.amax(diffplot),np.amax(diffplot)],c='k') 

err = np.array(err)
rms = np.sqrt(np.mean(err**2))
ax.set_title(f'RMS: {round(rms,3)}',fontsize=25)

f.savefig(f'DeltaGrid.png',bbox_inches='tight',pad_inches=.1)
plt.close(f)