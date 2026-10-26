"""Screen common-line baseline; truth is used only after estimation for diagnostics."""
import argparse
import json
import os
from pathlib import Path
import time

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from aspire.abinitio import CLSyncVoting
from aspire.source import ArrayImageSource, Simulation
from aspire.utils import Rotation
from aspire.volume import Volume


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--states', type=int, default=1)
    parser.add_argument('--n', type=int, default=128)
    args = parser.parse_args()
    rng = np.random.default_rng(20261025)
    size = 33
    grid = np.stack(np.meshgrid(*([np.arange(size) - size // 2] * 3), indexing='ij'))
    centers = rng.normal(size=(18, 3)) * 4
    densities = []
    for state in range(args.states):
        positions = centers.copy()
        if state:
            positions[:7] += [3, -2, 1]
        density = sum(np.exp(-np.sum((grid - p[:, None, None, None]) ** 2, axis=0) / 3.2)
                      for p in positions)
        densities.append(density)
    sim = Simulation(vols=Volume(np.asarray(densities, dtype=np.float32)),
                     n=args.n, offsets=0, amplitudes=1, seed=418)
    images = sim.images[:].asnumpy().copy()
    # The estimator receives only pixel arrays, never the simulation or its rotations.
    source = ArrayImageSource(images, pixel_size=1.0)
    start = time.monotonic()
    estimate = CLSyncVoting(source, n_theta=180, max_shift=0, mask=False)
    estimate.estimate_rotations()
    rotations = estimate.rotations
    result = {'states': args.states, 'images': args.n,
              'seconds': time.monotonic() - start,
              'rotation_mse_after_global_registration': float(Rotation(rotations).mse(Rotation(sim.rotations))),
              'note': 'Clean synthetic screening only; no release acceptance threshold.'}
    target = Path(__file__).parent / 'draft' / f'screen-{args.states}.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
