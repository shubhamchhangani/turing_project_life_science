"""Exact input-only marginalization over equivalence classes of phase histories."""
import csv,json,struct,subprocess,tempfile
from pathlib import Path
import numpy as np
from emissions import tables


def main():
 model=json.loads(Path('/root/data/model.json').read_text());tab=tables(model)
 with tempfile.TemporaryDirectory(prefix='phase-posterior-') as temp:
  d=Path(temp)
  for source,binary in [('orbit_graph.cpp','graph'),('infer.cpp','infer')]:
   subprocess.run(['g++','-O3','-std=c++17',str(Path(__file__).parent/source),'-o',str(d/binary)],check=True)
  with (d/'emissions.bin').open('wb') as f:f.write(struct.pack('<iii',tab.shape[0],2,tab.shape[-1]));f.write(tab.astype('<f8').tobytes())
  subprocess.run([str(d/'graph'),str(model['span']-1),str(d/'graph.bin')],check=True)
  subprocess.run([str(d/'infer'),str(d/'graph.bin'),str(d/'emissions.bin'),str(d/'checkpoints.bin'),str(d/'posterior.bin')],check=True)
  answer=np.fromfile(d/'posterior.bin','<f8');out=Path('/root/results');out.mkdir(parents=True,exist_ok=True)
  with (out/'evidence.csv').open('w',newline='') as f:
   writer=csv.writer(f);writer.writerow(['metric','value']);writer.writerow(['log_evidence',f'{answer[0]:.12f}'])
  with (out/'read_posteriors.csv').open('w',newline='') as f:
   writer=csv.writer(f);writer.writerow(['read_id','p_error_free','expected_errors'])
   for r,a in zip(model['reads'],answer[1:].reshape(-1,2)):writer.writerow([r['id'],f'{a[0]:.12f}',f'{a[1]:.12f}'])

if __name__=='__main__':main()
