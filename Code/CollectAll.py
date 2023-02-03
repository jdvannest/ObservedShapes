import argparse,os,pickle
config = pickle.load(open('Config.pickle','rb'))

parser = argparse.ArgumentParser(description='Collect images of all resolved halos from a given simulation. Images will be generated across all orientations.')
parser.add_argument('-n','--numproc',type=int,required=True,help='Number of processors to use')
parser.add_argument('-o','--overwrite',action='store_true',help='Overwrite existing images')
args = parser.parse_args()

overwrite = '-o' if args.overwrite else ''

for feedback in ['BW','SB']:
    sims = pickle.load(open(f'SimulationInfo.{feedback}.pickle','rb'))
    for s in sims:
        os.system(f"{config['python_path']} ImageCollection.py -f {feedback} -s {s} -n {args.numproc} {overwrite}")