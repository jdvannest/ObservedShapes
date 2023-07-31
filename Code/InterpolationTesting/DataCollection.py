import argparse,os,pickle,pymp,pynbody,sys,time,warnings
import numpy as np
import matplotlib.pylab as plt
from pynbody.plot.sph import image
from scipy.optimize import curve_fit
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
def sersic(r, mueff, reff, n):
    return mueff + 2.5*(0.868*n-0.142)*((r/reff)**(1./n) - 1)
warnings.filterwarnings("ignore")

# Storm-3 
# X: 90 - 150
# Y: 150 - 210


SimInfo = pickle.load(open(f'SimulationInfo.BW.pickle','rb'))
simpath = SimInfo['storm']['path']
halos = SimInfo['storm']['halos']
xinit,yinit = 90,150
xfinal,yfinal = 150,210
dx,dy = 10,10 #Angular resolution of rotations

print(f'Loading storm.BW-3')
tstart = time.time()
sim = pynbody.load(simpath)
sim.physical_units()
h = sim.halos()
halo = h[3]
pynbody.analysis.angmom.faceon(halo)
Rhalf = pynbody.analysis.luminosity.half_light_r(halo)
width = 6*Rhalf
ImageSpace = pynbody.filt.Sphere(width*np.sqrt(2)*1.01)
halo.rotate_x(xinit)
halo.rotate_y(yinit)
ImageData,SBData = {},{}
myprint(f'storm.BW-3 loaded.',clear=True)

print(f'Generating images: 0.00%')
xrot,prog = 0,0
while (xinit + xrot*dx)<xfinal:
    yrot = 0
    while (yinit + yrot*dy)<yfinal:
        key = f'x{(xinit + xrot*dx):03d}y{(yinit + yrot*dy):03d}'
        SBData[key] = {}
        SBData[key]['Rhalf'] = Rhalf
        try:
            prof = pynbody.analysis.profile.Profile(halo.s,type='lin',min=.25,max=5*Rhalf,ndim=2,nbins=int((5*Rhalf)/0.1))
            SBData[key]['sb,v'] = prof['sb,v']
            SBData[key]['v_lum_den'] = prof['v_lum_den']
            SBData[key]['rbins'] = prof['rbins']
        except:
            SBData[key]['sb,v'] = np.NaN
            SBData[key]['v_lum_den'] = np.NaN
            SBData[key]['rbins'] = np.NaN
        if not type(SBData[key]['sb,v']) is float:
            try:
                vband = prof['sb,v']
                smooth = np.nanmean(np.pad(vband.astype(float),(0,3-vband.size%3),mode='constant',constant_values=np.nan).reshape(-1,3),axis=1)
                x = np.arange(len(smooth))*0.3 + 0.15
                x[0] = .05
                if True in np.isnan(smooth):
                    x = np.delete(x,np.where(np.isnan(smooth)==True))
                    y = np.delete(smooth,np.where(np.isnan(smooth)==True))
                else: y = smooth
                r0 = x[int(len(x)/2)]
                m0 = np.mean(y[:3])
                par,ign = curve_fit(sersic,x,y,p0=(m0,r0,1),bounds=([10,0,0.5],[40,100,16.5]))
                SBData[key]['Reff'] = par[1]
            except:
                SBData[key]['Reff'] = np.NaN
        else:
            SBData[key]['Reff'] = np.NaN
        f = plt.figure(frameon=False)
        f.set_size_inches(10,10)
        ax = plt.Axes(f, [0., 0., 1., 1.])
        ax.set_axis_off()
        f.add_axes(ax)
        im = image(sim[ImageSpace].s,qty='v_lum_den',width=width,subplot=ax,units='kpc^-2',resolution=1000,show_cbar=False)
        f.savefig(f'Images/storm.3.{key}.png')
        plt.close()
        #Store data
        ImageData[key] = im
        #Progress to next orientation
        prog+=1
        myprint(f'Generating images: {round(prog/((xfinal-xinit)/dx * (yfinal-yinit)/dy)*100,2)}%',clear=True)
        halo.rotate_y(dy)
        yrot+=1
    
    halo.rotate_y(yinit-yfinal)
    halo.rotate_x(dx)
    xrot+=1

tstop = time.time()
print(f'storm.bw-3 imaged in {round((tstop-tstart)/60,2)} minutes.',clear=True)
pickle.dump(SBData,open('Storm.3.Profiles.pickle','wb'))
pickle.dump(ImageData,open('Storm.3.Images.pickle','wb'))