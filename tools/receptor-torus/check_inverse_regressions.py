"""Exercise complete artifact contracts and independently check moment coverage."""
import csv
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import numpy as np

AUTHOR=Path(__file__).resolve().parent
ROOT=AUTHOR.parents[1]
TASK=ROOT/'tasks/life-sciences/biology/inverse-receptor-cooperativity'
sys.path.insert(0,str(TASK/'tests'))
import test_outputs as checks

fitted=np.load(AUTHOR/'draft/field-inverse/newton-fitted.npy')


def encode(rows):
    stream=io.StringIO(newline='')
    csv.writer(stream,lineterminator='\n').writerows(rows)
    return stream.getvalue()


base={'couplings.csv':encode([['x','y','direction','K']]+[[*key,f'{value:.10f}'] for key,value in zip(checks.EDGES,fitted[:-3])]),
      'biases.csv':encode([['receptor_id','h']]+[[key,f'{value:.10f}'] for key,value in zip(checks.SITES,fitted[-3:])])}
cases=[]


def case(name,accept,change):
    cases.append((name,accept,change))


def cell(filename,row,col,value):
    def change(files):
        data=list(csv.reader(io.StringIO(files[filename])))
        data[row][col]=value;files[filename]=encode(data)
    return change


case('genuine_input_only_fit',True,lambda f:None)
case('reordered_rows',True,lambda f:f.update({k:encode([r[0]]+r[1:][::-1]) for k in f for r in [list(csv.reader(io.StringIO(f[k]))) ]}))
case('CRLF',True,lambda f:f.update({k:v.replace('\n','\r\n') for k,v in f.items()}))
def scientific(files):
    for key in files:
        rows=list(csv.reader(io.StringIO(files[key])))
        for row in rows[1:]:row[-1]=f'{float(row[-1]):.12e}'
        files[key]=encode(rows)
case('scientific_notation',True,scientific)
for label,value in [('leading_space',' 0.1'),('trailing_space','0.1 '),('underscore','0_1'),('unicode','٠.1'),('nan','NaN'),('infinity','inf'),('overflow','1e999'),('empty',''),('below_bound','-.350000001'),('above_bound','.900000001')]:
    case(label,False,cell('couplings.csv',1,3,value))
case('bias_whitespace',False,cell('biases.csv',1,1,' .43'))
case('bias_bound',False,cell('biases.csv',1,1,'.800001'))
case('bias_wrong_sign',False,cell('biases.csv',2,1,'.31'))
case('six_hundredth_coupling_wrong',False,cell('couplings.csv',600,3,str(fitted[599]+.01)))
case('leading_zero_coordinate',False,cell('couplings.csv',1,0,'00'))
case('signed_coordinate',False,cell('couplings.csv',1,0,'+0'))
case('out_of_range_coordinate',False,cell('couplings.csv',1,0,'17'))
case('unknown_direction',False,cell('couplings.csv',1,2,'W'))
case('duplicate_edge',False,cell('couplings.csv',2,2,'E'))
case('duplicate_site',False,cell('biases.csv',2,0,'F1'))
case('unknown_site',False,cell('biases.csv',2,0,'F4'))
case('quoted_header',False,lambda f:f.update({'couplings.csv':f['couplings.csv'].replace('x,y,direction,K','"x",y,direction,K',1)}))
case('missing_row',False,lambda f:f.update({'couplings.csv':encode(list(csv.reader(io.StringIO(f['couplings.csv'])))[:-1])}))
case('extra_row',False,lambda f:f.update({'biases.csv':f['biases.csv']+'F3,0.24\n'}))
case('blank_row',False,lambda f:f.update({'couplings.csv':f['couplings.csv']+'\n'}))
case('extra_column',False,lambda f:f.update({'couplings.csv':encode([r+['x'] for r in csv.reader(io.StringIO(f['couplings.csv']))])}))
case('oversize',False,lambda f:f.update({'couplings.csv':'x'*1048577}))
case('invalid_utf8',False,lambda f:f.update({'couplings.csv':b'\xff'}))
case('bare_CR',False,lambda f:f.update({'couplings.csv':f['couplings.csv'].replace('\n','\r')}))
for label in ['missing','fifo','symlink']:
    case(label,False,lambda f:None)
results=[]
cache={}
for name,should_accept,change in cases:
    with tempfile.TemporaryDirectory(prefix='inverse-receptor-check-') as temporary:
        folder=Path(temporary);files=base.copy();change(files)
        for key,content in files.items():
            (folder/key).write_bytes(content if isinstance(content,bytes) else content.encode())
        if name in {'missing','fifo','symlink'}:(folder/'couplings.csv').unlink()
        if name=='fifo':os.mkfifo(folder/'couplings.csv')
        if name=='symlink':(folder/'couplings.csv').symlink_to(folder/'biases.csv')
        start=time.monotonic();reason=''
        try:
            parameters=checks.contract(folder)
            checks.test_01_complete_contract(parameters)
            checks.test_02_every_coupling(parameters)
            checks.test_03_every_bias(parameters)
            key=parameters.tobytes()
            if key not in cache:
                cache[key]=checks.calculate_moments(parameters)
            assert max(abs(cache[key]-checks.observed_moments()))<=2e-6
            accepted=True
        except (AssertionError,ValueError,OSError,UnicodeError,csv.Error) as exc:
            accepted=False;reason=f'{type(exc).__name__}: {exc}'
        results.append({'name':name,'accepted':accepted,'expected':should_accept,'passed':accepted==should_accept,'reason':reason,'seconds':time.monotonic()-start})

# All coordinates remain inside the parameter tolerance, but their collective
# change must be rejected by the complete equilibrium-observable requirement.
shifted=fitted.copy();shifted[:-3]+=4e-5
checks.test_02_every_coupling(shifted);checks.test_03_every_bias(shifted)
moment_error=float(max(abs(checks.calculate_moments(shifted)-checks.observed_moments())))
results.append({'name':'collective_error_within_parameter_tolerance','expected':False,'accepted':moment_error<=2e-6,'passed':moment_error>2e-6,'moment_error':moment_error})
report={'all_passed':all(r['passed'] for r in results),'valid_cases':sum(r['expected'] for r in results),'invalid_cases':sum(not r['expected'] for r in results),'results':results}
(AUTHOR/'draft/field-inverse/regressions.json').write_text(json.dumps(report,indent=2)+'\n')
print({k:v for k,v in report.items() if k!='results'},flush=True)
raise SystemExit(0 if report['all_passed'] else 1)
