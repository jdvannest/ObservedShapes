from pynbody import array
from pynbody.analysis import profile
import numpy as np
import logging
logger = logging.getLogger('pynbody.analysis.halo')
import numpy as np 
from numpy import matmul as X
from numpy.linalg import eig, inv
from numpy import sin,cos
from math import pi

#Projection Functions
def cart2pol(x, y):
    #Cartesian to Polar coordinate conversion
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    #Polar to Cartesian coordinate conversion
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def Rx(theta):
    #X-axis rotation matrix
    a = np.radians(theta)
    return( [[1,0,0],[0,cos(a),-sin(a)],[0,sin(a),cos(a)]] )

def Ry(theta):
    #Y-axis rotation matrix
    a = np.radians(theta)
    return( [[cos(a),0,sin(a)],[0,1,0],[-sin(a),0,cos(a)]] )

def EllipseFit(x,y):
    #Fit an ellipse to 2D scatter data
    ex = np.array(x)
    ey = np.array(y)
    ex = ex[:,np.newaxis]
    ey = ey[:,np.newaxis]
    D =  np.hstack((ex*ex, ex*ey, ey*ey, ex, ey, np.ones_like(ex)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    return( V[:,n] )

def EllipseAxes(ellipse):
    #Return a and b axes lengths for fit ellipse
    b,c,d,f,g,a = ellipse[1]/2,ellipse[2],ellipse[3]/2,ellipse[4]/2,ellipse[5],ellipse[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1,res2])

def Ellipsoid(b,c,xrot,yrot):
    #Generate 3D scatter points for an ellipsoid with given b/a and c/a axis ratios
    #and x,y rotations
    resolution = 50
    phi = np.linspace(0,2*np.pi,resolution)
    theta = np.linspace(0,np.pi,resolution)
    x,y,z=[[],[],[]]
    for i in np.arange(resolution):
        for j in np.arange(resolution):
            x.append(sin(theta[i])*cos(phi[j]))
            y.append(b*sin(theta[i])*sin(phi[j]))
            z.append(c*cos(theta[i]))
    for i in np.arange(len(x)):
        vec = [x[i],y[i],z[i]]
        rot1 = X(vec,Rx(xrot))
        rot2 = X(rot1,Ry(yrot))
        x[i],y[i],z[i] = rot2
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    return(x,y,z)

def Ellipticity(x,y,z):
    #Return the b/a axis ratios of a flattened 3D rotated ellipsoid
    phi_bins = np.linspace(-pi,pi,50)
    r,phi = cart2pol(np.array(x),np.array(y))
    r_edge,phi_edge = [[],[]]
    for i in np.arange(len(phi_bins)-1):
        binr,binphi = [[],[]]
        for j in np.arange(len(r)):
            if phi_bins[i]<phi[j]<phi_bins[i+1]:
                binr.append(r[j])
                binphi.append(phi[j])
        if len(binr)>0:
            r_edge.append(max(binr))
            phi_edge.append(binphi[binr.index(max(binr))])
    x_edge,y_edge=pol2cart(r_edge,phi_edge)
    Efit = EllipseFit(x_edge,y_edge)
    a,b = EllipseAxes(Efit)
    return(b/a)

def myround(x, base=5):
    return base * round(x/base)

def RandomRotation(res=15):
    phi = np.degrees(np.random.uniform(0,np.pi*2))
    y = int(myround(phi,base=res))
    theta = np.degrees(np.arccos(np.random.uniform(-1,1)))
    x = int(myround(theta,base=res))
    if x>165:
        x = x%180
        y = (180-y)
        if y<0:
            y = y+360
    if y == 360:
        y = 0
    return(x,y)


#Pynbody Shape functions for stellar morphology
def halo_shape_stellar(sim, N=100, rin=None, rout=None, bins='equal', ret_pos=False):
    #-----------------------------FUNCTIONS-----------------------------
    # Define an ellipsoid shell with lengths a,b,c and orientation E:
    def Ellipsoid(r, a,b,c, E):
        x,y,z = np.dot(np.transpose(E),[r[:,0],r[:,1],r[:,2]])
        return (x/a)**2 + (y/b)**2 + (z/c)**2

    # Define moment of inertia tensor:
    MoI = lambda r,m: np.array([[np.sum(m*r[:,i]*r[:,j]) for j in range(3)]\
                               for i in range(3)])

    # Splits 'r' array into N groups containing equal numbers of particles.
    # An array is returned with the radial bins that contain these groups.
    sn = lambda r,N: np.append([r[i*int(len(r)/N):(1+i)*int(len(r)/N)][0]\
                               for i in range(N)],r[-1])

    # Retrieves alignment angle:
    almnt = lambda E: np.arccos(np.dot(np.dot(E,[1.,0.,0.]),[1.,0.,0.]))
    #-----------------------------FUNCTIONS-----------------------------

    if (rout == None): rout = sim.s['r'].max()
    if (rin == None): rin = rout/1E3

    posr = np.array(sim.s['r'])[np.where(sim.s['r'] < rout)]
    pos = np.array(sim.s['pos'])[np.where(sim.s['r'] < rout)]
    mass = np.array(sim.s['mass'])[np.where(sim.s['r'] < rout)]

    rx = [[1.,0.,0.],[0.,0.,-1.],[0.,1.,0.]]
    ry = [[0.,0.,1.],[0.,1.,0.],[-1.,0.,0.]]
    rz = [[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]]

    # Define bins:
    if (bins == 'equal'): # Each bin contains equal number of particles
        mid = sn(np.sort(posr[np.where((posr >= rin) & (posr <= rout))]),N*2)
        rbin = mid[1:N*2+1:2]
        mid = mid[0:N*2+1:2]

    elif (bins == 'log'): # Bins are logarithmically spaced
        mid = profile.Profile(sim.s, type='log', ndim=3, rmin=rin, rmax=rout, nbins=N+1)['rbins']
        rbin = np.sqrt(mid[0:N]*mid[1:N+1])

    elif (bins == 'lin'): # Bins are linearly spaced
        mid = profile.Profile(sim.s, type='lin', ndim=3, rmin=rin, rmax=rout, nbins=N+1)['rbins']
        rbin = 0.5*(mid[0:N]+mid[1:N+1])

    # Define b/a and c/a ratios and angle arrays:
    ba,ca,angle,aout,n,n_i,pos_out,pos_in = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),[],[]
    Es = [0]*N

    # Begin loop through radii:
    for i in range(0,N):

        # Initialise convergence criterion:
        tol = .1
        count = 0

        # Define initial spherical shell:
        a=b=c = rbin[i]
        E = np.identity(3)
        L1,L2 = rbin[i]-mid[i],mid[i+1]-rbin[i]

        # Begin iterative procedure to fit data to shell:
        while True:
            count += 1

            # Collect all particle positions and masses within shell:
            r = pos[np.where((posr < a+L2) & (posr > c-L1*c/a))]
            inner = Ellipsoid(r, a-L1,b-L1*b/a,c-L1*c/a, E)
            outer = Ellipsoid(r, a+L2,b+L2*b/a,c+L2*c/a, E)
            r = r[np.where((inner > 1.) & (outer < 1.))]
            m = mass[np.where((inner > 1.) & (outer < 1.))]
            num = len(r)
            if (count == 1):
                num_i = len(r)
                pos_in.append(r)

            # End iterations if there is no data in range:
            if (len(r) == 0):
                ba[i],ca[i],angle[i],Es[i],n_i[i] = b/a,c/a,almnt(E),E,num_i
                logger.info('No data in range after %i iterations' %count)
                break

            # Calculate shape tensor & diagonalise:
            D = list(np.linalg.eig(MoI(r,m)/np.sum(m)))

            # Purge complex numbers:
            if isinstance(D[1][0,0],complex):
                D[0] = D[0].real ; D[1] = D[1].real
                logger.info('Complex numbers in D removed...')

            # Compute ratios a,b,c from moment of intertia principles:
            anew,bnew,cnew = np.sqrt(abs(D[0])*3.0)

            # The rotation matrix must be reoriented:
            E = D[1]
            if ((bnew > anew) & (anew >= cnew)): E = np.dot(E,rz)
            if ((cnew > anew) & (anew >= bnew)): E = np.dot(np.dot(E,ry),rx)
            if ((bnew > cnew) & (cnew >= anew)): E = np.dot(np.dot(E,rz),rx)
            if ((anew > cnew) & (cnew >= bnew)): E = np.dot(E,rx)
            if ((cnew > bnew) & (bnew >= anew)): E = np.dot(E,ry)
            cnew,bnew,anew = np.sort(np.sqrt(abs(D[0])*3.0))

            # Keep a as semi-major axis and distort b,c by b/a and c/a:
            div = rbin[i]/anew
            anew *= div
            bnew *= div
            cnew *= div

            # Convergence criterion:
            if (np.abs(b/a-bnew/anew) < tol) & (np.abs(c/a-cnew/anew) < tol):
                if (almnt(-E) < almnt(E)): E = -E
                aout[i],ba[i],ca[i],angle[i],Es[i],n[i],n_i[i] = anew,bnew/anew,cnew/anew,almnt(E),E,num,num_i
                pos_out.append(r)
                break

            # Increase tolerance if convergence has stagnated:
            elif (count%10 == 0): tol *= 5.

            # Reset a,b,c for the next iteration:
            a,b,c = anew,bnew,cnew
    if ret_pos:
        return [array.SimArray(rbin,sim.d['pos'].units),aout,ba,ca,angle,Es,n,n_i,np.array(pos_out),np.array(pos_in)]
    else:
        return [array.SimArray(rbin,sim.d['pos'].units),aout,ba,ca,angle,Es,n,n_i]