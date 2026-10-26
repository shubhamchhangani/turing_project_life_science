"""Generate independent pulse/chase measurements from injected fluxes and pools."""
import sys,csv,json,time
from pathlib import Path
import numpy as np
import dynamic_model as D
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'metabolic-flux-analysis'))
import true_fluxes as T
ROOT=Path(__file__).resolve().parents[2]/'tasks/life-sciences/biology/metabolic-flux-analysis-v2'
POOLS={'G6P':20,'F6P':16,'T3P':7,'P5P':15,'E4P':8,'S7P':18,'PG3':10,'Pyr':18,'Lac':30,'Ala':12,'Ser':14,'AcCoA':5,'OAA':8,'Cit':45,'AKG':22,'Suc':12,'Fum':10,'Mal':30,'Asp':20,'Glu':35}
UNKNOWN=['OAA','Cit','AKG','Mal']
TIMES=[.02,.05,.1,.2,.4,.8,1.6,2.5,3.8,4.5,6,8,12]

def write(folder,name,fields,rows):
    with (folder/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def generate(seed,folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(seed)
    schedule=[]
    for exp,chase in [('glc12',4.0),('gln5',6.5)]:
        for start in [0,chase]:
            row={'experiment':exp,'start_min':start,'co2_13c':.0107}
            for i in range(1,7): row[f'glucose_C{i}']=.99 if exp=='glc12' and i<=2 and start==0 else .0107
            for i in range(1,6): row[f'glutamine_C{i}']=.99 if exp=='gln5' and start==0 else .0107
            schedule.append(row)
    write(folder,'labeling_schedule.csv',list(schedule[0]),schedule)
    write(folder,'pool_sizes.csv',['metabolite','pool_size'],[{'metabolite':m,'pool_size':'unknown' if m in UNKNOWN else POOLS[m]} for m in D.M.BALANCED])
    meas=[]
    for exp in ['glc12','gln5']:
        for t in TIMES:
            for frag,(_,carbs,_) in D.M.FRAGMENTS.items():
                for k in range(len(carbs)+3):meas.append({'experiment':exp,'time_min':t,'fragment':frag,'isotopologue':f'M{k}','fraction':0,'sd':.002})
    start=time.monotonic()
    pred=D.predictions(T.NET,T.EXCH,POOLS,schedule,meas,rtol=2e-10,atol=1e-13)
    print('simulation sec',time.monotonic()-start,flush=True)
    for r,y in zip(meas,pred):r['fraction']=f'{y+rng.normal(0,.002):.7f}'
    write(folder,'mid_measurements.csv',list(meas[0]),meas)
    rates=[]
    for rid in ['GLCup','GLNup','LACout','ALAout','SERout','ASPout','GLUout','P5Pout','ACCOAout']:
        sd=0 if rid=='GLCup' else T.NET[rid]*.03
        rates.append({'reaction_id':rid,'rate':f'{T.NET[rid]+rng.normal(0,sd):.6f}','sd':f'{sd:.6f}'})
    write(folder,'extracellular_rates.csv',list(rates[0]),rates)
    print('measurements',len(meas),'truth SSR',sum(((float(r['fraction'])-y)/.002)**2 for r,y in zip(meas,pred)),flush=True)
    return pred

if __name__=='__main__':generate(int(sys.argv[1]) if len(sys.argv)>1 else 20261023,sys.argv[2] if len(sys.argv)>2 else ROOT/'environment/data')
