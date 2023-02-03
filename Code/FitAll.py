import argparse,os,pickle
config = pickle.load(open('Config.pickle','rb'))

parser = argparse.ArgumentParser(description='Fit Isophotes to all generated images.')
parser.add_argument('-n','--numproc',type=int,required=True,help='Number of processors to use')
args = parser.parse_args()

for feedback in ['BW','SB']:
    sims = pickle.load(open(f'SimulationInfo.{feedback}.pickle','rb'))
    for s in sims:
        os.system(f"{config['python_path']} IsophoteFitting.py -f {feedback} -s {s} -n {args.numproc}")