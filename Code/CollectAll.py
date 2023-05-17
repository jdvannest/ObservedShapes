import argparse,os,pickle
config = pickle.load(open('Config.pickle','rb'))

parser = argparse.ArgumentParser(description='Collect data from all simulations')
parser.add_argument('-n','--numproc',type=int,required=True,help='Number of processors to use')
parser.add_argument('-v','--verbose',action='store_true',help='Print halo IDs being analyzed')
parser.add_argument('-o','--overwrite',action='store_true',help='Overwrite existing images')
args = parser.parse_args()

overwrite = '-o' if args.overwrite else ''
verbose = '-v' if args.verbose else ''

loop = True
while loop:
    type = input('Collect Images or Shapes (I/S): ')
    if type in ['I','S']:
        loop = False 
if type=='I':
    loop = True
    while loop:
        im = input('Generate images in addition to Profiles? (y/n): ')
        if im in ['y','n']: loop = False
    gen_im = '-i' if im=='y' else ''

for feedback in ['BW','SB']:
    sims = pickle.load(open(f'SimulationInfo.{feedback}.pickle','rb'))
    for s in sims:
        if type=='I':
            os.system(f"{config['python_path']} ImageCollection.py -f {feedback} -s {s} {gen_im} -n {args.numproc} {verbose} {overwrite}")
        elif type=='S':
            os.system(f"{config['python_path']} 3DShapeCollection.py -f {feedback} -s {s} -n {args.numproc} {verbose}")
