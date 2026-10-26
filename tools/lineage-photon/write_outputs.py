"""Write the public output contract from a computed result; no expected values."""
import csv,json
from pathlib import Path
from reference import solve,states

def write(model,result,out):
    out.mkdir(parents=True,exist_ok=True)
    rows={
        'evidence.csv':('metric,value',[('log_evidence',result['log_evidence'])]),
        'division_posterior.csv':('state,probability',[(f'{s}:{m}:{p}',v) for (s,m,p),v in zip(states(model),result['parent'])]),
        'daughter_rna.csv':('counts,probability',[(f'{a}:{b}',v) for a,row in enumerate(result['daughter_rna']) for b,v in enumerate(row)]),
        'transcription_counts.csv':('count,probability',[(str(k),v) for k,v in enumerate(result['counts'])]+[('overflow',result['count_overflow'])]),
    }
    for name,(header,records) in rows.items():
        with (out/name).open('w',newline='') as f:
            w=csv.writer(f,lineterminator='\n');w.writerow(header.split(','))
            w.writerows((key,format(value,'.16g')) for key,value in records)

if __name__=='__main__':
    m=json.loads(Path('/root/data/model.json').read_text());r=solve(m);write(m,r,Path('/root/results'))
    print(json.dumps({'reference_seconds':r['seconds']}),flush=True)
