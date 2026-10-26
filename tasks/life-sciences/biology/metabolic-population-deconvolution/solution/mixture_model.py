"""Transient population dynamics and molecule-weighted bulk observations.

This module contains model equations only, with no fitted or injected answers.
"""
import csv
from collections import OrderedDict
from pathlib import Path
import numpy as np
import dynamic_model as D
import forward_model as M
import network_model as N
POPS = ('A', 'B')
def read_rows(folder, name):
    with (folder / name).open() as f:
        return list(csv.DictReader(f, delimiter='\t' if name.endswith('.tsv') else ','))

def write_rows(folder, name, fields, rows):
    with (folder / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

class Model:

    def __init__(self, folder):
        self.folder = Path(folder)
        rows = read_rows(self.folder, 'network.tsv')
        self.rxns = [r['reaction_id'] for r in rows]
        self.rev = {r['reaction_id']: r['reversible'] == 'yes' for r in rows}
        self.flux = N.FluxModel(self.rxns, self.rev)
        self.weights = {r['population']: float(r['protein_fraction']) for r in read_rows(self.folder, 'populations.csv')}
        self.known = {p: {} for p in POPS}
        self.unknown = {p: [] for p in POPS}
        for r in read_rows(self.folder, 'pool_sizes.csv'):
            if r['pool_size'] == 'unknown':
                self.unknown[r['population']].append(r['metabolite'])
            else:
                self.known[r['population']][r['metabolite']] = float(r['pool_size'])
        self.schedule = read_rows(self.folder, 'labeling_schedule.csv')
        self.meas = read_rows(self.folder, 'mid_measurements.csv') if (self.folder / 'mid_measurements.csv').exists() else []
        self.rates = read_rows(self.folder, 'extracellular_rates.csv') if (self.folder / 'extracellular_rates.csv').exists() else []
        self.width = len(self.flux.free) + len(self.flux.exch) + len(self.unknown['A'])
        self._simulation_cache = OrderedDict()

    def simulate_population(self, pop, exp, net, exchange, pools, schedule, ts, rtol):
        key = (pop, exp, tuple(net.items()), tuple(exchange.items()), tuple(pools.items()), rtol)
        if key in self._simulation_cache:
            self._simulation_cache.move_to_end(key)
            return self._simulation_cache[key]
        states = D.simulate(net, exchange, pools, schedule, ts, rtol=rtol, atol=1e-11)
        self._simulation_cache[key] = states
        if len(self._simulation_cache) > 16:
            self._simulation_cache.popitem(last=False)
        return states

    def decode(self, parameters):
        out = {}
        nf = len(self.flux.free)
        ne = len(self.flux.exch)
        for j, pop in enumerate(POPS):
            p = parameters[j * self.width:(j + 1) * self.width]
            net = self.flux.net(p[:nf])
            exchange = dict(zip(self.flux.exch, p[nf:nf + ne]))
            pools = {**self.known[pop], **dict(zip(self.unknown[pop], p[nf + ne:]))}
            out[pop] = (net, exchange, pools)
        return out

    def predict(self, decoded, rtol=3e-08):
        predicted = {}
        for exp in sorted({r['experiment'] for r in self.meas}):
            ts = sorted({float(r['time_min']) for r in self.meas if r['experiment'] == exp})
            schedule = [r for r in self.schedule if r['experiment'] == exp]
            sims = {}
            for pop in POPS:
                net, exchange, pools = decoded[pop]
                physical = {rid: max(value, 0) if not self.rev[rid] else value for rid, value in net.items()}
                sims[pop] = self.simulate_population(pop, exp, physical, exchange, pools, schedule, ts, rtol)
            for t in ts:
                mixture = {}
                for met in self.known['A'].keys() | set(self.unknown['A']):
                    amounts = [self.weights[p] * decoded[p][2][met] for p in POPS]
                    mixture[met] = sum((a * sims[p][t][met] for p, a in zip(POPS, amounts))) / sum(amounts)
                for frag in M.FRAGMENTS:
                    predicted[exp, t, frag] = M.raw_mid(mixture, frag)
        return np.array([predicted[r['experiment'], float(r['time_min']), r['fragment']][int(r['isotopologue'][1:])] for r in self.meas])

    def write_outputs(self, decoded, out):
        out.mkdir(parents=True, exist_ok=True)
        fluxes = []
        pools = []
        for pop, (net, exchange, amounts) in decoded.items():
            f, b = self.flux.directional(net, [exchange[r] for r in self.flux.exch])
            for rid in self.rxns:
                fluxes.append(dict(population=pop, reaction_id=rid, net_flux=f'{net[rid]:.8f}', forward_flux=f'{f[rid]:.8f}', reverse_flux=f'{b.get(rid, 0):.8f}'))
            for met in self.unknown[pop]:
                pools.append(dict(population=pop, metabolite=met, pool_size=f'{amounts[met]:.8f}'))
        write_rows(out, 'fluxes.csv', ['population', 'reaction_id', 'net_flux', 'forward_flux', 'reverse_flux'], fluxes)
        write_rows(out, 'pools.csv', ['population', 'metabolite', 'pool_size'], pools)
