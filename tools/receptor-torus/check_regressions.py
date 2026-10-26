"""Behavioral checks using independently calculated thermodynamics, never reference-derived targets."""
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
AUTHOR=Path(__file__).resolve().parent
TASK=ROOT/'tasks/life-sciences/biology/receptor-cluster-thermodynamics'
sys.path.insert(0,str(TASK/'tests'))
import test_outputs as checks

independent=json.loads((AUTHOR/'draft/independent.json').read_text())
expected=(np.array(independent['joint']),independent['log_partition'])
source=AUTHOR/'draft/reference'
base={name:(source/name).read_text() for name in ['joint.csv','thermodynamics.csv']}


def rows(text):
    return list(csv.reader(io.StringIO(text)))


def encode(parsed):
    out=io.StringIO(newline='')
    csv.writer(out,lineterminator='\n').writerows(parsed)
    return out.getvalue()


cases=[]
def case(name,accepted,mutate):
    cases.append((name,accepted,mutate))


def field(name,row,column,value):
    def mutate(files):
        data=rows(files[name]);data[row][column]=value
        files[name]=encode(data)
    return mutate


case('genuine_reference',True,lambda f:None)
case('arbitrary_row_order',True,lambda f:f.update({'joint.csv':encode(rows(f['joint.csv'])[:1]+rows(f['joint.csv'])[1:][::-1])}))
case('CRLF',True,lambda f:f.update({k:v.replace('\n','\r\n') for k,v in f.items()}))
def scientific(files):
    for name in files:
        data=rows(files[name])
        for row in data[1:]:row[1]=f'{float(row[1]):.12e}'
        files[name]=encode(data)
case('scientific_notation',True,scientific)
for label,value in [('leading_space',' 0.1'),('trailing_space','0.1 '),('underscore','0_1'),('unicode_digit','٠.1'),('nan','NaN'),('infinity','inf'),('overflow','1e999'),('empty',''),('negative','-0.000001'),('above_one','1.000001')]:
    case(label,False,field('joint.csv',1,1,value))
case('padded_partition',False,field('thermodynamics.csv',1,1,' 326.091954464613'))
case('wrong_quantity',False,field('thermodynamics.csv',1,0,'free_energy'))
case('bad_partition',False,field('thermodynamics.csv',1,1,str(expected[1]+np.log(2))))
case('pattern_lost_zeroes',False,field('joint.csv',1,0,'0'))
case('duplicate_pattern',False,field('joint.csv',2,0,'000000'))
case('quoted_header',False,lambda f:f.update({'joint.csv':f['joint.csv'].replace('pattern,probability','"pattern",probability',1)}))
case('wrong_header',False,lambda f:f.update({'joint.csv':f['joint.csv'].replace('pattern,probability','pattern,value',1)}))
case('extra_column',False,lambda f:f.update({'joint.csv':encode([r+['extra'] for r in rows(f['joint.csv'])])}))
case('missing_row',False,lambda f:f.update({'joint.csv':encode(rows(f['joint.csv'])[:-1])}))
case('extra_row',False,lambda f:f.update({'joint.csv':f['joint.csv']+'000000,0.1\n'}))
case('blank_row',False,lambda f:f.update({'joint.csv':f['joint.csv']+'\n'}))
case('summary_extra_row',False,lambda f:f.update({'thermodynamics.csv':f['thermodynamics.csv']+'log_partition,0\n'}))
def uniform(files):
    data=rows(files['joint.csv'])
    for r in data[1:]:r[1]='0.015625'
    files['joint.csv']=encode(data)
case('independent_receptors',False,uniform)
def sixth_moment(files):
    data=rows(files['joint.csv'])
    epsilon=min(expected[0])/4
    assert epsilon>1e-7
    for r in data[1:]:
        r[1]=f'{float(r[1])+epsilon*(-1)**r[0].count("1"):.12f}'
    files['joint.csv']=encode(data)
case('sixth_order_only_error',False,sixth_moment)
case('oversized_artifact',False,lambda f:f.update({'joint.csv':'x'*1048577}))
case('invalid_utf8',False,lambda f:f.update({'joint.csv':b'\xff'}))
for name in ['missing_file','fifo','symlink']:
    case(name,False,lambda f:None)

results=[]
for name,should_pass,mutate in cases:
    with tempfile.TemporaryDirectory(prefix='receptor-regression-') as temporary:
        folder=Path(temporary)
        files=base.copy();mutate(files)
        for filename,contents in files.items():
            (folder/filename).write_bytes(contents if isinstance(contents,bytes) else contents.encode())
        if name in {'missing_file','fifo','symlink'}:
            (folder/'joint.csv').unlink()
        if name=='fifo':os.mkfifo(folder/'joint.csv')
        if name=='symlink':(folder/'joint.csv').symlink_to(source/'joint.csv')
        start=time.monotonic()
        reason=''
        try:
            submitted=checks.contract(folder)
            for function in [checks.test_01_complete_contract,checks.test_02_normalization,checks.test_03_inversion_symmetry,checks.test_04_all_probe_marginals]:
                function(submitted)
            checks.test_05_full_joint_distribution(submitted,expected)
            checks.test_06_absolute_partition_function(submitted,expected)
            accepted=True
        except (AssertionError,ValueError,OSError,UnicodeError,csv.Error) as exc:
            accepted=False;reason=f'{type(exc).__name__}: {exc}'
        duration=time.monotonic()-start
        results.append({'case':name,'expected_accept':should_pass,'accepted':accepted,'passed':accepted==should_pass,'seconds':duration,'reason':reason})
report={'all_passed':all(r['passed'] for r in results),'valid_cases':sum(c[1] for c in cases),'invalid_cases':sum(not c[1] for c in cases),'results':results}
(AUTHOR/'draft/regressions.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='results'},indent=2))
raise SystemExit(0 if report['all_passed'] else 1)
