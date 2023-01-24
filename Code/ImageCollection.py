import argparse,pickle,os,pynbody,sys,time
import numpy as np
import matplotlib.pylab as plt
from pynbody.plot.sph import image
from multiprocessing import Manager,Pool
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)

parser = argparse.ArgumentParser(description='Collect images of all resolved halos from a given simulation. Images will be generated across all orientations.')
parser.add_argument('-f','--feedback',choices=['BW','SB'],default='BW',help='Feedback Model')
parser.add_argument('-s','--simulation',choices=['cptmarvel','elektra','storm','rogue'],required=True,help='Simulation to analyze')
parser.add_argument('-n','--numproc',type=int,required=True,help='Number of processors to use')
parser.add_argument('-o','--overwrite',action='store_true',help='Overwrite existing images')
args = parser.parse_args()

SimInfo = pickle.load(open(f'SimulationInfo.{args.feedback}.pickle','rb'))
simpath = SimInfo[args.simulation]['path']
halos = SimInfo[args.simulation]['path']
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
myprint(f'{args.simulation} loaded.',clear=True)

with Manager() as manager:
    prog = manager.list([0])
    print(f'\tGenerating images: {round(prog[0]/len(halos)*100,2)}%')
    def GenerateImages(hid):
        #Check if halo has been complete
        if f'{hid}.x{180-dx:03d}.y{360-dy}.png' in os.listdir(f'../Images/{args.simulation}.{args.feedback}/{hid}/') and not args.overwrite:
            prog[0]+=1
            myprint(f'\tGenerating images: {round(prog/len(halos)*100,2)}%',clear=True)
            return
        halo = h[hid]
        pynbody.analys.angmom.faceon(halo)
        Rhalf = pynbody.analysis.luminosity.half_light_r(halo)
        width = 3*Rhalf
        ImageSpace = pynbody.filt.Sphere(width*np.sqrt(2)*1.01)
        xrotation,yrotation = 0,0

        while xrotation*dx<180:
            while yrotation*dy<360:
                if f'{hid}.x{xrotation*dx:03d}.y{yrotation*dy:03d}.png' not in os.listdir(f'../Images/{args.simulation}.{args.feedback}/{hid}/') or args.overwrite:
                    f = plt.figure(frameon=False)
                    f.set_size_inches(10,10)
                    ax = plt.Axes(f, [0., 0., 1., 1.])
                    ax.set_axis_off()
                    f.add_axes(ax)
                    im = image(sim[ImageSpace].s,qty='v_lum_den',width=width,subplot=ax,units='kpc^-2',resolution=1000,show_cbar=False)
                    f.savefig(f'../Images/{args.simulation}.{args.feedback}/{hid}/{hid}.x{xrotation*dx:03d}.y{yrotation*dy:03d}.png')
                    plt.close()
                halo.rotate_y(dy)
                yrotation+=1
            halo.rotate_x(dx)
            xrotation+=1
        prog[0]+=1
        myprint(f'\tGenerating images: {round(prog[0]/len(halos)*100,2)}%',clear=True)
    p = Pool(args.numproc)
    p.map(GenerateImages,halos)

tstop = time.time()
print(f'{args.simulation} completed in {round((tstop-tstart)/60,2)} minutes.')