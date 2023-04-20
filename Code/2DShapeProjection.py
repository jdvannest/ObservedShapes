import argparse,pickle,pymp,sys
import numpy as np
from math import pi
from numpy import degrees
from numpy import sin,cos
from matplotlib.patches import Ellipse
import matplotlib.pylab as plt
from ShapeFunctions.Functions import Ellipsoid,Project,kpc2pix,pix2kpc,Project_OLD

def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)

parser = argparse.ArgumentParser(description='Collect images of all resolved halos from a given simulation. Images will be generated across all orientations.')
parser.add_argument('-f','--feedback',choices=['BW','SB'],default='BW',help='Feedback Model')
parser.add_argument('-s','--simulation',choices=['cptmarvel','elektra','storm','rogue'],required=True,help='Simulation to analyze')
parser.add_argument('-n','--numproc',type=int,required=True,help='Number of processors to use')
parser.add_argument('-p','--plot',action='store_true')
#parser.add_argument('-o','--overwrite',action='store_true',help='Overwrite existing images')
args = parser.parse_args()

SimInfo = pickle.load(open(f'SimulationInfo.{args.feedback}.pickle','rb'))
halos = SimInfo[args.simulation]['halos']
Images = pickle.load(open(f'../Data/{args.simulation}.{args.feedback}.Images.pickle','rb'))
Profiles = pickle.load(open(f'../Data/{args.simulation}.{args.feedback}.Profiles.pickle','rb'))
Shapes = pickle.load(open(f'../Data/{args.simulation}.{args.feedback}.3DShapes.pickle','rb'))

ShapeData = pymp.shared.dict()
prog=pymp.shared.array((1,),dtype=int)
print(f'Projecting Ellipsoids for {args.simulation}: 0.000%')
with pymp.Parallel(args.numproc) as pl:
    for i in pl.xrange(len(halos)):
        hid = halos[i]
        current_shape = {}
        for xrot in np.arange(0,180,30):
            for yrot in np.arange(0,360,30):
                #Get 3D Ellipsoid Projection
                rbins = Shapes[str(hid)]['rbins']
                Rhalf = Profiles[str(hid)][f'x000y000']['Rhalf']
                ind_eff = np.argmin(abs(rbins-Rhalf))
                a,ba,ca,Es = [Shapes[str(hid)]['a'][ind_eff],Shapes[str(hid)]['ba'][ind_eff],
                              Shapes[str(hid)]['ca'][ind_eff],Shapes[str(hid)]['Es'][ind_eff]]
                x,y,z = Ellipsoid(a,ba,ca,Es,xrot,yrot)
                x = -x
                ap,bp,cen,phi = Project_OLD(x,y,z)
                atrue,btrue = max([ap,bp]),min([ap,bp])
                current_shape[f'x{xrot:03d}y{yrot:03d}'] = {}
                current_shape[f'x{xrot:03d}y{yrot:03d}']['b/a'] = btrue/atrue
                current_shape[f'x{xrot:03d}y{yrot:03d}']['a'] = atrue
                current_shape[f'x{xrot:03d}y{yrot:03d}']['b'] = btrue

                if args.plot:
                    #Get Isophote
                    Rhalf = Profiles[str(hid)][f'x{xrot:03d}y{yrot:03d}']['Rhalf']
                    rbins = Profiles[str(hid)][f'x{xrot:03d}y{yrot:03d}']['rbins']
                    ind_eff = np.argmin(abs(rbins-Rhalf))
                    v = Profiles[str(hid)][f'x{xrot:03d}y{yrot:03d}']['v_lum_den'][ind_eff]
                    im = Images[str(hid)][f'x{xrot:03d}y{yrot:03d}']
                    #im = np.flip(im,0)
                    iso,tolerance = [[[],[]],.01]
                    if not np.isnan(v):
                        while len(iso[0])==0 and tolerance<.1:
                            iso = np.where((im>v*(1-tolerance)) & (im<v*(1+tolerance)))
                            tolerance+=.01
                    
                    f,ax = plt.subplots(1,1)
                    LogImage = plt.imread(f'../Images/{args.simulation}.{args.feedback}/{hid}/{hid}.x{xrot:03d}.y{yrot:03d}.png')
                    ax.imshow(LogImage)
                    ax.set_xlim([0,1e3])
                    ax.set_ylim([1e3,0])
                    ax.scatter(500,500,c='k',marker='+')
                    ax.scatter(iso[1],iso[0],c='r',marker='.',s=1)
                    x,y = kpc2pix(x,6*Rhalf)+500,kpc2pix(y,6*Rhalf)+500
                    cen = [cen[0]+500,cen[1]+500]
                    ap,bp = kpc2pix(ap,6*Rhalf),kpc2pix(bp,6*Rhalf)
                    ax.add_patch(Ellipse(cen,2*ap,2*bp,angle=degrees(phi),facecolor='None',edgecolor='orange',linewidth=2))
                    plt.plot([-ap*cos(phi)+cen[0],ap*cos(phi)+cen[0]],[-ap*sin(phi)+cen[1],ap*sin(phi)+cen[1]],color='orange')
                    plt.plot([-bp*cos(phi+pi/2)+cen[0],bp*cos(phi+pi/2)+cen[0]],[-bp*sin(phi+pi/2)+cen[1],bp*sin(phi+pi/2)+cen[1]],color='orange')
                    ax.set_title(f'Projected b/a: {round(btrue/atrue,3)}',fontsize=15)
                    f.savefig(f'../Images/{args.simulation}.{args.feedback}/{hid}/{hid}.x{xrot:03d}.y{yrot:03d}.Intrinsic.png',bbox_inches='tight',pad_inches=.1)
                
                prog[0]+=1
                myprint(f'Projecting Ellipsoids for {args.simulation}: {round(prog[0]/(len(halos)*72)*100,3)}%',clear=True)
        ShapeData[str(hid)] = current_shape

#Pymp shared dictionaries 'corrupt' when saving, so transfer to standerd dict objects
ShapeFile = {}
for hid in halos:
    ShapeFile[str(hid)] = ShapeData[str(hid)]
pickle.dump(ShapeFile,open(f'../Data/{args.simulation}.{args.feedback}.ProjectedData.pickle','wb'))