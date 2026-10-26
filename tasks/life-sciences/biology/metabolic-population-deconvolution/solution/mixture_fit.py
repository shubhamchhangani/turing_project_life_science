"""Fit only the supplied measurements from independent starting points."""
import sys,time,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from mixture_model import Model, POPS
def fit(folder, out, starts, seed):
    m = Model(folder)
    rng = np.random.default_rng(seed)
    y = np.array([float(r['fraction']) for r in m.meas])
    sd = np.array([float(r['sd']) for r in m.meas])
    rates = {r['reaction_id']: (float(r['rate']), float(r['sd'])) for r in m.rates}
    count = 0
    begun = time.monotonic()

    def residual(p):
        nonlocal count
        count += 1
        decoded = m.decode(p)
        pred = m.predict(decoded)
        extra = [(sum((m.weights[pop] * decoded[pop][0][r] for pop in POPS)) - v) / s for r, (v, s) in rates.items() if s > 0]
        penalty = [1000.0 * min(v, 0) for pop in POPS for r, v in decoded[pop][0].items() if not m.rev[r]]
        residuals = np.r_[(pred - y) / sd, extra, penalty]
        if count % 500 == 0:
            print('evaluations', count, 'SSR', float(residuals @ residuals), 'elapsed_sec', time.monotonic() - begun, flush=True)
        return residuals
    nf = len(m.flux.free)
    ne = len(m.flux.exch)
    nu = len(m.unknown['A'])
    lo = np.tile(np.r_[np.full(nf, 1e-05), np.zeros(ne), np.full(nu, 0.5)], 2)
    hi = np.tile(np.r_[np.full(nf, 300), np.full(ne, 500), np.full(nu, 200)], 2)
    best = None
    for start in range(starts):
        points = []
        for pop in POPS:
            z = [rates[r][0] * rng.uniform(0.7, 1.3) if r in rates else rng.uniform(10, 55) for r in m.flux.free]
            points.extend(np.r_[z, rng.uniform(5, 100, ne), rng.uniform(8, 60, nu)])
        p = np.clip(points, lo + 0.0001, hi - 0.0001)
        opt = least_squares(residual, p, bounds=(lo, hi), x_scale='jac', diff_step=0.0001, max_nfev=200, ftol=2e-07, xtol=2e-07, gtol=1e-05)
        print('start', start, 'SSR', 2 * opt.cost, 'nfev', opt.nfev, 'elapsed_sec', time.monotonic() - begun, flush=True)
        if best is None or opt.cost < best.cost:
            best = opt
            m.write_outputs(m.decode(opt.x), out)
            np.save(out / 'parameters.npy', opt.x)
            sv = np.linalg.svd(opt.jac, compute_uv=False)
            np.save(out / 'jacobian.npy', opt.jac)
            report = {'ssr': float(2 * opt.cost), 'residual_evaluations': count, 'elapsed_sec': time.monotonic() - begun, 'smallest_jacobian_singular_value': float(sv[-1]), 'condition_number': float(sv[0] / sv[-1])}
            (out / 'fit_report.json').write_text(json.dumps(report, indent=2) + '\n')
if __name__ == '__main__':
    fit(Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
