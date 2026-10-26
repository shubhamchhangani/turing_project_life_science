"""Author prototype: infer two metabolic phenotypes from bulk pulse-chase MIDs.

Generation and fitting have separate entry points. The fitter reads no truth
file. This is experimental author work, not a submission-ready task package.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('NUMBA_CACHE_DIR','/private/tmp/mfa-mixture-numba')
import argparse, csv, json, shutil, time
from collections import OrderedDict
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import dynamic_model as D
import forward_model as M
import network_model as N

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]/'tasks/life-sciences/biology/metabolic-flux-analysis-v2/environment/data'
POPS=('A','B')

def read_rows(folder,name):
    with (folder/name).open() as f:return list(csv.DictReader(f,delimiter='\t' if name.endswith('.tsv') else ','))

def write_rows(folder,name,fields,rows):
    with (folder/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

class Model:
    def __init__(self,folder):
        self.folder=Path(folder)
        rows=read_rows(self.folder,'network.tsv')
        self.rxns=[r['reaction_id'] for r in rows]
        self.rev={r['reaction_id']:r['reversible']=='yes' for r in rows}
        self.flux=N.FluxModel(self.rxns,self.rev)
        self.weights={r['population']:float(r['protein_fraction']) for r in read_rows(self.folder,'populations.csv')}
        self.known={p:{} for p in POPS};self.unknown={p:[] for p in POPS}
        for r in read_rows(self.folder,'pool_sizes.csv'):
            if r['pool_size']=='unknown':self.unknown[r['population']].append(r['metabolite'])
            else:self.known[r['population']][r['metabolite']]=float(r['pool_size'])
        self.schedule=read_rows(self.folder,'labeling_schedule.csv')
        self.meas=read_rows(self.folder,'mid_measurements.csv') if (self.folder/'mid_measurements.csv').exists() else []
        self.rates=read_rows(self.folder,'extracellular_rates.csv') if (self.folder/'extracellular_rates.csv').exists() else []
        self.width=len(self.flux.free)+len(self.flux.exch)+len(self.unknown['A'])
        self._simulation_cache=OrderedDict()

    def simulate_population(self,pop,exp,net,exchange,pools,schedule,ts,rtol):
        # Exact parameter keys preserve finite-difference derivatives. Cache
        # repeated simulations of the unchanged population during fitting.
        key=(pop,exp,tuple(net.items()),tuple(exchange.items()),tuple(pools.items()),rtol)
        if key in self._simulation_cache:
            self._simulation_cache.move_to_end(key)
            return self._simulation_cache[key]
        states=D.simulate(net,exchange,pools,schedule,ts,rtol=rtol,atol=1e-11)
        self._simulation_cache[key]=states
        if len(self._simulation_cache)>16:self._simulation_cache.popitem(last=False)
        return states

    def decode(self,parameters):
        out={}
        nf=len(self.flux.free);ne=len(self.flux.exch)
        for j,pop in enumerate(POPS):
            p=parameters[j*self.width:(j+1)*self.width]
            net=self.flux.net(p[:nf])
            exchange=dict(zip(self.flux.exch,p[nf:nf+ne]))
            pools={**self.known[pop],**dict(zip(self.unknown[pop],p[nf+ne:]))}
            out[pop]=(net,exchange,pools)
        return out

    def predict(self,decoded,rtol=3e-8):
        predicted={}
        for exp in sorted({r['experiment'] for r in self.meas}):
            ts=sorted({float(r['time_min']) for r in self.meas if r['experiment']==exp})
            schedule=[r for r in self.schedule if r['experiment']==exp]
            sims={}
            for pop in POPS:
                net,exchange,pools=decoded[pop]
                physical={rid:max(value,0) if not self.rev[rid] else value for rid,value in net.items()}
                sims[pop]=self.simulate_population(pop,exp,physical,exchange,pools,schedule,ts,rtol)
            for t in ts:
                mixture={}
                for met in self.known['A'].keys() | set(self.unknown['A']):
                    amounts=[self.weights[p]*decoded[p][2][met] for p in POPS]
                    mixture[met]=sum(a*sims[p][t][met] for p,a in zip(POPS,amounts))/sum(amounts)
                # Pool molecules mix before natural-isotope convolution and
                # peak truncation. Averaging separately normalized raw MIDs
                # would produce the wrong experimental observation model.
                for frag in M.FRAGMENTS:
                    predicted[exp,t,frag]=M.raw_mid(mixture,frag)
        return np.array([predicted[r['experiment'],float(r['time_min']),r['fragment']][int(r['isotopologue'][1:])] for r in self.meas])

    def write_outputs(self,decoded,out):
        out.mkdir(parents=True,exist_ok=True)
        fluxes=[];pools=[]
        for pop,(net,exchange,amounts) in decoded.items():
            f,b=self.flux.directional(net,[exchange[r] for r in self.flux.exch])
            for rid in self.rxns:
                fluxes.append(dict(population=pop,reaction_id=rid,net_flux=f'{net[rid]:.8f}',forward_flux=f'{f[rid]:.8f}',reverse_flux=f'{b.get(rid,0):.8f}'))
            for met in self.unknown[pop]:pools.append(dict(population=pop,metabolite=met,pool_size=f'{amounts[met]:.8f}'))
        write_rows(out,'fluxes.csv',['population','reaction_id','net_flux','forward_flux','reverse_flux'],fluxes)
        write_rows(out,'pools.csv',['population','metabolite','pool_size'],pools)

def generate(folder,seed):
    import runpy
    folder.mkdir(parents=True,exist_ok=True)
    for f in ['network.tsv','fragments.csv','isotope_abundances.csv']:shutil.copy2(BASE/f,folder/f)
    write_rows(folder,'populations.csv',['population','protein_fraction'],[dict(population='A',protein_fraction=.37),dict(population='B',protein_fraction=.63)])
    rows=read_rows(BASE,'pool_sizes.csv')
    b_known=dict(G6P=35,F6P=10,T3P=12,P5P=9,E4P=15,S7P=11,PG3=18,Pyr=32,Lac=18,Ala=21,Ser=9,AcCoA=9,Suc=21,Fum=6,Asp=34,Glu=18)
    pool_rows=[dict(population=p,metabolite=r['metabolite'],pool_size=(r['pool_size'] if p=='A' or r['pool_size']=='unknown' else b_known[r['metabolite']])) for p in POPS for r in rows]
    write_rows(folder,'pool_sizes.csv',['population','metabolite','pool_size'],pool_rows)
    schedule=read_rows(BASE,'labeling_schedule.csv')
    template=schedule[0]
    natural=.0107
    for exp,switches in [('glcU',[(0,True,False),(3,False,False)]),('dual',[(0,True,True),(3,False,True),(7.5,False,False)])]:
        for start,glc,gln in switches:
            row={k:(exp if k=='experiment' else start if k=='start_min' else natural) for k in template}
            for i in range(1,7):row[f'glucose_C{i}']=.99 if glc and (exp=='glcU' or i==1) else natural
            for i in range(1,6):row[f'glutamine_C{i}']=.98 if gln else natural
            schedule.append(row)
    write_rows(folder,'labeling_schedule.csv',list(template),schedule)
    m=Model(folder)
    truth_a=runpy.run_path(str(HERE.parent/'metabolic-flux-analysis/true_fluxes.py'))
    a=(truth_a['NET'],truth_a['EXCH'],{**m.known['A'],'OAA':8,'Cit':45,'AKG':22,'Mal':30})
    bnet=m.flux.net([36,7,8,10,8,10,8,55,65,22])
    b=(bnet,dict(zip(m.flux.exch,[25,35,12,40,24,30,90,20])),{**m.known['B'],'OAA':20,'Cit':20,'AKG':42,'Mal':14})
    decoded={'A':a,'B':b}
    assert all(v>=-1e-8 for p in POPS for r,v in decoded[p][0].items() if not m.rev[r])
    raw=read_rows(BASE,'mid_measurements.csv')
    prototype=[r for r in raw if r['experiment']=='glc12']
    m.meas=[{**r,'experiment':exp,'sd':.001} for exp in ['glc12','gln5','glcU','dual'] for r in prototype]
    y=m.predict(decoded);rng=np.random.default_rng(seed)
    for r,value in zip(m.meas,y):r['fraction']=float(value+rng.normal(0,.001))
    write_rows(folder,'mid_measurements.csv',list(m.meas[0]),m.meas)
    rate_rows=[]
    for r in read_rows(BASE,'extracellular_rates.csv'):
        rid=r['reaction_id'];value=sum(m.weights[p]*decoded[p][0][rid] for p in POPS)
        sd=0 if rid=='GLCup' else .03*abs(value)
        rate_rows.append(dict(reaction_id=rid,rate=value+rng.normal(0,sd),sd=sd))
    write_rows(folder,'extracellular_rates.csv',['reaction_id','rate','sd'],rate_rows)
    expected={p:{'net':d[0],'exchange':d[1],'pools':d[2]} for p,d in decoded.items()}
    (folder.parent/'injected_truth.json').write_text(json.dumps(expected,indent=2)+'\n')
    print('Generated',len(y),'MID observations with seed',seed,flush=True)

def fit(folder,out,starts,seed):
    m=Model(folder);rng=np.random.default_rng(seed)
    y=np.array([float(r['fraction']) for r in m.meas]);sd=np.array([float(r['sd']) for r in m.meas])
    rates={r['reaction_id']:(float(r['rate']),float(r['sd'])) for r in m.rates}
    count=0;begun=time.monotonic()
    def residual(p):
        nonlocal count
        count+=1
        decoded=m.decode(p)
        pred=m.predict(decoded)
        extra=[(sum(m.weights[pop]*decoded[pop][0][r] for pop in POPS)-v)/s for r,(v,s) in rates.items() if s>0]
        penalty=[1e3*min(v,0) for pop in POPS for r,v in decoded[pop][0].items() if not m.rev[r]]
        residuals=np.r_[(pred-y)/sd,extra,penalty]
        if count % 500 == 0:
            print('evaluations',count,'SSR',float(residuals@residuals),'elapsed_sec',time.monotonic()-begun,flush=True)
        return residuals
    nf=len(m.flux.free);ne=len(m.flux.exch);nu=len(m.unknown['A'])
    lo=np.tile(np.r_[np.full(nf,1e-5),np.zeros(ne),np.full(nu,.5)],2)
    hi=np.tile(np.r_[np.full(nf,300),np.full(ne,500),np.full(nu,200)],2)
    best=None
    for start in range(starts):
        points=[]
        for pop in POPS:
            z=[rates[r][0]*rng.uniform(.7,1.3) if r in rates else rng.uniform(10,55) for r in m.flux.free]
            points.extend(np.r_[z,rng.uniform(5,100,ne),rng.uniform(8,60,nu)])
        p=np.clip(points,lo+1e-4,hi-1e-4)
        opt=least_squares(residual,p,bounds=(lo,hi),x_scale='jac',diff_step=1e-4,max_nfev=200,ftol=2e-7,xtol=2e-7,gtol=1e-5)
        print('start',start,'SSR',2*opt.cost,'nfev',opt.nfev,'elapsed_sec',time.monotonic()-begun,flush=True)
        if best is None or opt.cost<best.cost:
            best=opt;m.write_outputs(m.decode(opt.x),out)
            np.save(out/'parameters.npy',opt.x)
            sv=np.linalg.svd(opt.jac,compute_uv=False)
            np.save(out/'jacobian.npy',opt.jac)
            report={'ssr':float(2*opt.cost),'residual_evaluations':count,'elapsed_sec':time.monotonic()-begun,'smallest_jacobian_singular_value':float(sv[-1]),'condition_number':float(sv[0]/sv[-1])}
            (out/'fit_report.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['generate','fit'])
    p.add_argument('--data',type=Path,default=HERE/'draft/data')
    p.add_argument('--out',type=Path,default=HERE/'draft/reference')
    p.add_argument('--seed',type=int,default=20261024)
    p.add_argument('--starts',type=int,default=2)
    a=p.parse_args()
    if a.action=='generate':generate(a.data,a.seed)
    else:fit(a.data,a.out,a.starts,a.seed)
