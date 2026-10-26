"""Author-only reproduction of the scientific fixture; requires NumPy and g++."""
import json,struct,subprocess,tempfile
from pathlib import Path
from emissions import tables

if __name__=='__main__':
 here=Path(__file__).parent;data=here.parent/'data';model=json.loads((data/'model.json').read_text());table=tables(model)
 with tempfile.TemporaryDirectory(prefix='independent-partitions-') as tmp:
  d=Path(tmp)
  for source,binary in [('orbit_graph.cpp','graph'),('partition_source.cpp','source')]:subprocess.run(['g++','-O3','-std=c++17',str(here/source),'-o',str(d/binary)],check=True)
  with (d/'emissions.bin').open('wb') as f:f.write(struct.pack('<iii',table.shape[0],2,table.shape[-1]));f.write(table.astype('<f8').tobytes())
  subprocess.run([str(d/'graph'),str(model['span']-1),str(d/'graph.bin')],check=True)
  subprocess.run([str(d/'source'),str(d/'graph.bin'),str(d/'emissions.bin'),str(d/'expectations.json')],check=True)
  generated=json.loads((d/'expectations.json').read_text());expected=json.loads((data/'expectations.json').read_text())
  delta=max(abs(generated['log_evidence']-expected['log_evidence']),max(abs(a-b) for ar,br in zip(generated['read_expectations'],expected['read_expectations']) for a,b in zip(ar,br)))
  assert len(generated['read_expectations'])==len(model['reads']) and delta<1e-9
  print('Independent fixture reproduced; maximum difference:',delta)
