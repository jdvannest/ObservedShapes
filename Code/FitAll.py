import argparse,os,pickle
config = pickle.load(open('Config.pickle','rb'))

parser = argparse.ArgumentParser(description='Fit Isophotes to all generated images.')
parser.add_argument('-n','--numproc',type=int,required=True,help='Number of processors to use')
parser.add_argument('-o','--overwrite',action='store_true',help='Overwrite existing data')
args = parser.parse_args()

loop = True
while loop:
    method = input('Fit Isophotes, Project 3D Ellipsoids, or Run MCMC (I/3/M): ')
    if method in ['I','3','M']:
        loop = False
        if method=='I':
            script = 'IsophoteFitting.py'
        elif method=='3':
            script = '2DShapeProjection.py -p'
        else:
            script = 'MCMC.py'

over = '-o' if args.overwrite else ''

for feedback in ['BW','SB']:
    sims = pickle.load(open(f'SimulationInfo.{feedback}.pickle','rb'))
    for s in sims:
        os.system(f"{config['python_path']} {script} -f {feedback} -s {s} -n {args.numproc} {over}")
        #os.system(f"/usr/local/anaconda/bin/python {script} -f {feedback} -s {s} -n {args.numproc}")
