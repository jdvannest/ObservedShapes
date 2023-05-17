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

parser = argparse.ArgumentParser(description='Collect images of all resolved halos from a given simulation. Images will be generated across all orientations.')
parser.add_argument('-f','--feedback',choices=['BW','SB'],default='BW',help='Feedback Model')
parser.add_argument('-s','--simulation',choices=['cptmarvel','elektra','storm','rogue'],required=True,help='Simulation to analyze')
parser.add_argument('-i','--image',action='store_true',help='Generate images in addition to profile data')
parser.add_argument('-n','--numproc',type=int,required=True,help='Number of processors to use')
parser.add_argument('-o','--overwrite',action='store_true',help='Overwrite existing images')
parser.add_argument('-v','--verbose',action='store_true',help='Print halo IDs being analyzed')
args = parser.parse_args()

SimInfo = pickle.load(open(f'SimulationInfo.{args.feedback}.pickle','rb'))
simpath = SimInfo[args.simulation]['path']
halos = SimInfo[args.simulation]['halos']
dx,dy = 30,30 #Angular resolution of rotations
#Check if all halos in sim have been completed
if f'{halos[-1]}.x{180-dx:03d}.y{360-dy}.png' in os.listdir(f'../Images/{args.simulation}.{args.feedback}/{halos[-1]}/') and not args.overwrite:
    print(f'{args.simulation} completed.')
    sys.exit(0)

print(f'Loading {args.simulation}')
tstart = time.time()
sim = pynbody.load(simpath)
sim.physical_units()
h = sim.halos()
if args.image: ImageData = pymp.shared.dict()
SBData = pymp.shared.dict()
myprint(f'{args.simulation} loaded.',clear=True)

prog=pymp.shared.array((1,),dtype=int)
print(f'\tGenerating images: {round(prog[0]/len(halos)*100,2)}%')
with pymp.Parallel(args.numproc) as pl:
    for i in pl.xrange(len(halos)):
        t_start_current = time.time()
        if args.verbose: print(f'\tAnalyzing {halos[i]}...')
        hid = halos[i]
        halo = h[hid]
        pynbody.analysis.angmom.faceon(halo)
        Rhalf = pynbody.analysis.luminosity.half_light_r(halo)
        width = 6*Rhalf
        ImageSpace = pynbody.filt.Sphere(width*np.sqrt(2)*1.01)
        current_image = {}
        current_sb = {}
        xrotation = 0
        while xrotation*dx<180:
            yrotation = 0
            while yrotation*dy<360:
                #Find V-band SB at Reff
                current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}'] = {}
                current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['Rhalf'] = Rhalf
                try:
                    prof = pynbody.analysis.profile.Profile(halo.s,type='lin',min=.25,max=5*Rhalf,ndim=2,nbins=int((5*Rhalf)/0.1))
                    current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['sb,v'] = prof['sb,v']
                    current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['v_lum_den'] = prof['v_lum_den']
                    current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['rbins'] = prof['rbins']
                except:
                    current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['sb,v'] = np.NaN
                    current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['v_lum_den'] = np.NaN
                    current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['rbins'] = np.NaN
                if not type(current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['sb,v']) is float:
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
                        current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['Reff'] = par[1]
                    except:
                        current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['Reff'] = np.NaN
                else:
                    current_sb[f'x{xrotation*dx:03d}y{yrotation*dy:03d}']['Reff'] = np.NaN
                if args.image:
                    #Generate V-band SB image
                    f = plt.figure(frameon=False)
                    f.set_size_inches(10,10)
                    ax = plt.Axes(f, [0., 0., 1., 1.])
                    ax.set_axis_off()
                    f.add_axes(ax)
                    im = image(sim[ImageSpace].s,qty='v_lum_den',width=width,subplot=ax,units='kpc^-2',resolution=1000,show_cbar=False)
                    f.savefig(f'../Images/{args.simulation}.{args.feedback}/{hid}/{hid}.x{xrotation*dx:03d}.y{yrotation*dy:03d}.png')
                    plt.close()
                    #Store data
                    current_image[f'x{xrotation*dx:03d}y{yrotation*dy:03d}'] = im
                #Progress to next orientation
                prog[0]+=1
                if not args.verbose: myprint(f'\tGenerating images: {round(prog[0]/(len(halos)*72)*100,2)}%',clear=True)
                halo.rotate_y(dy)
                yrotation+=1
            halo.rotate_x(dx)
            xrotation+=1
        if args.image: ImageData[str(hid)] = current_image
        SBData[str(hid)] = current_sb
        t_end_current = time.time()
        if args.verbose: print(f'\t\t{hid} done in {round((t_end_current-t_start_current)/60,2)} minutes.')

ImageFile = pickle.load(open(f'../Data/{args.simulation}.{args.feedback}.Images.pickle','rb'))
SBFile = pickle.load(open(f'../Data/{args.simulation}.{args.feedback}.Profiles.pickle','rb'))
for halo in halos:
    if args.image: ImageFile[str(halo)] = ImageData[str(halo)]
    SBFile[str(halo)] = SBData[str(halo)]
if args.image: pickle.dump(ImageFile,open(f'../Data/{args.simulation}.{args.feedback}.Images.pickle','wb'))
pickle.dump(SBFile,open(f'../Data/{args.simulation}.{args.feedback}.Profiles.pickle','wb'))
tstop = time.time()
myprint(f'\t{args.simulation} imaged in {round((tstop-tstart)/60,2)} minutes.',clear=True)
