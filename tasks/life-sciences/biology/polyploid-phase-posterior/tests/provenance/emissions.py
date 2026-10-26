import numpy as np

def tables(model):
    span = model['span']
    n = model['variants']
    readlist = model['reads']
    tab = np.empty((n - span + 1, 2, 3, 2 ** span), np.float64)
    for i, r in enumerate(readlist):
        end = r['end']
        for pattern in range(2 ** span):
            likelihood = 1.0
            errors = 0
            for call in r['calls']:
                bit = pattern >> end - call['variant'] & 1
                wrong = bit != call['allele']
                errors += wrong
                likelihood *= call['error_probability'] if wrong else 1 - call['error_probability']
            tab[i // 2, i % 2, :, pattern] = [likelihood, likelihood * (errors == 0), likelihood * errors]
    return tab
