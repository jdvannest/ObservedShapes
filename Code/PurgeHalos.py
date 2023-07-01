import argparse,os,pickle,sys

parser = argparse.ArgumentParser(description='Collect images of all resolved halos from a given simulation. Images will be generated across all orientations.')
parser.add_argument('-f','--feedback',choices=['BW','SB'],default='BW',help='Feedback Model')
parser.add_argument('-s','--simulation',choices=['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329'],required=True,help='Simulation to analyze')
args = parser.parse_args()

SimInfo = pickle.load(open(f'SimulationInfo.{args.feedback}.pickle','rb'))
datadir = os.listdir('../Data')
for f in datadir:
    if f.split('.')[0]==args.simulation:
        if f.split('.')[1]==args.feedback:
            s = pickle.load(open(f'../Data/{f}','rb'))
            delete = []
            for halo in s:
                if int(halo) not in SimInfo[args.simulation]['halos']: 
                    delete.append(halo)
            if len(delete)==0: 
                print(f'{args.simulation}.{args.feedback} clean')
                sys.exit(0)
            for halo in delete:
                del s[halo]
                os.system(f'rm -rf ../Images/{args.simulation}.{args.feedback}/{halo}')
            pickle.dump(s,open(f'../Data/{f}','wb'))