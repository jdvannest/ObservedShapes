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

# loop,remake=True,False
# while loop:
#     rem = input('Remake Image Directory: Images/Interpolation/- (y/n): ')
#     if rem in ['y','n']:
#         loop = False
#         if rem=='y': remake = True
# if remake: 
#     os.system('rmdir -f ../Images/Interpolation')
#     os.system('mkdir ../Images/Interpolation')
#     for sim in SimInfo:
#         os.system(f'mkdir ../Images/Interpolation/{sim}.BW/')

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
xang = np.arange(0,180,30)
yang = np.arange(0,360,30)
Yedge = np.arange(-15,185,30)
Xedge = np.arange(-15,375,30)
Xgrid,Ygrid = np.meshgrid(yang,xang)
Xgrid_inter,Ygrid_inter = np.meshgrid(np.arange(0,390,30),np.arange(0,210,30))
ptszip = zip(Xgrid.ravel(),Ygrid.ravel())
xang_fill = np.arange(0,210,30)
yang_fill = np.arange(0,390,30)
Yedge_fill = np.arange(-15,215,30)
Xedge_fill = np.arange(-15,405,30)



for sim in SimInfo:
    Shapes = pickle.load(open(f'../../Data/{sim}.BW.ShapeData.pickle','rb'))
    for hid in SimInfo[sim]['halos']:
        plot = np.zeros((len(xang),len(yang)))
        Z3 = np.zeros(Xgrid.shape)
        for x in np.arange(len(xang)):
            for y in np.arange(len(yang)):
                plot[x][y] = Shapes[str(hid)][f'x{x*30:03d}y{y*30:03d}']['b/a']
                Z3[x][y] = Shapes[str(hid)][f'x{Ygrid[:,0][x]:03d}y{Xgrid[0][y]:03d}']['b/a']
        
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

        f.savefig(f'../../Images/Interpolation/{sim}.BW/{hid}.Data.png',bbox_inches='tight',pad_inches=.1)
        plt.close(f)



        f,ax = plt.subplots(1,1,figsize=(12,8))
        ax.tick_params(labelsize=15)
        ax.set_xticks(np.arange(0,405,30))
        ax.set_yticks(np.arange(0,225,30))
        ax.set_xlim([-15,375])
        ax.set_ylim([-15,195])
        ax.set_xlabel(r'$\phi$-rotation [$^o$]',fontsize=25)
        ax.set_ylabel(r'$\theta$-rotation [$^o$]',fontsize=25)

        plot,bad_p1,bad_d1 = FillNaN(plot)
        plot_padded = PadData(plot)
        plot_padded,bad_p2,bad_d2 = FillNaN(plot_padded)
        if len(bad_d1)>0:
            for i in np.arange(len(bad_p1)):
                ax.add_patch(Rec((bad_p1[i][1]-5,bad_p1[i][0]-5),10,10,color='w'))
        if len(bad_d2)>0:
            for i in np.arange(len(bad_p2)):
                ax.add_patch(Rec((bad_p2[i][1]-5,bad_p2[i][0]-5),10,10,color='w'))

        p = ax.pcolor(Xedge_fill,Yedge_fill,plot_padded,cmap='viridis',edgecolors='k',norm=norm)
        c = f.colorbar(p,ax=ax,pad=0.01,aspect=15)
        c.set_label('b/a',fontsize=25)
        c.ax.tick_params(labelsize=15)
        c.ax.plot([0,1],[np.amin(plot),np.amin(plot)],c='w') 
        c.ax.plot([0,1],[np.amax(plot),np.amax(plot)],c='w') 

        f.savefig(f'../../Images/Interpolation/{sim}.BW/{hid}.Filled.png',bbox_inches='tight',pad_inches=.1)
        plt.close(f)


        #https://scipython.com/blog/non-linear-least-squares-fitting-of-a-two-dimensional-data/
        if True in np.isnan(Z3): Z3,bad_p,bad_d = FillNaN(Z3)
        Zpad = PadData(Z3)
        inter = RGI(points=(Ygrid_inter[:,0],Xgrid_inter[0]),values=Zpad,method='cubic')
        
        # Plot the 3D figure 
        # f = plt.figure()
        # ax = f.gca(projection='3d')
        # ax.set_xlabel(r'$\phi$-rotation [$^o$]',fontsize=15)
        # ax.set_ylabel(r'$\theta$-rotation [$^o$]',fontsize=15)
        # ax.set_zlabel('b/a',fontsize=15)
        # ax.set_xticks(np.arange(0,375,30))
        # ax.set_yticks(np.arange(0,195,30))
        # ax.set_xlim([-15,345])
        # ax.set_ylim([-15,165])
        # ax.set_zlim(0,1)

        # Xi,Yi = np.meshgrid(np.arange(0,360.5,.5),np.arange(0,180.5,.5))
        # Zi = np.zeros(Xi.shape)
        # for x in np.arange(Xi.shape[0]):
        #     for y in np.arange(Xi.shape[1]):
        #         Zi[x][y] = inter((Yi[:,0][x],Xi[0][y]))
        # ax.plot_surface(Xi, Yi, Zi, cmap='viridis',norm=norm)
        # for x in np.arange(len(xang)):
        #     for y in np.arange(len(yang)):
        #         z = Shapes[str(hid)][f'x{x*30:03d}y{y*30:03d}']['b/a']
        #         ax.scatter(y*30,x*30,z,c='k')#c=z,cmap='viridis',norm=norm)

        # plt.show()