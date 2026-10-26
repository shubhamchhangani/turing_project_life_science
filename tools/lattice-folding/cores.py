"""Sequence-independent compact-core construction, with no stored folds."""
import itertools
from functools import lru_cache

MOVES = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]


def edges(core):
    return sum(tuple(p[k]+d[k] for k in range(3)) in core for p in core for d in MOVES)//2


@lru_cache(maxsize=100000)
def canonical(points):
    images=[]
    for axes in itertools.permutations(range(3)):
        for signs in itertools.product([-1,1],repeat=3):
            transformed=[tuple(p[axes[k]]*signs[k] for k in range(3)) for p in points]
            lo=[min(p[k] for p in transformed) for k in range(3)]
            images.append(tuple(sorted(tuple(p[k]-lo[k] for k in range(3)) for p in transformed)))
    return min(images)


def compact_cores(n, side=4):
    """Peel surface sites, retaining best-contact cores at each layer.

    This is a constructive search, not a claim to exhaust all lattice cores.
    Any threading found in one of these cores is subsequently energy-checked.
    """
    cores={tuple(itertools.product(range(side),repeat=3))}
    history=[]
    for count in range(side**3,n,-1):
        candidates=set()
        best=-1
        for points in sorted(cores):
            core=set(points)
            score=edges(core)
            degrees={p:sum(tuple(p[k]+d[k] for k in range(3)) in core for d in MOVES) for p in core}
            for p,degree in sorted(degrees.items()):
                remaining=score-degree
                if remaining>=best:
                    c=canonical(tuple(sorted(core-{p})))
                    if remaining>best:
                        candidates.clear()
                        best=remaining
                    candidates.add(c)
        cores=candidates
        history.append({'H':count-1,'edges':best,'cores':len(cores)})
    return sorted(cores),history


if __name__=='__main__':
    import json
    import time
    start=time.monotonic()
    c,h=compact_cores(55)
    print(json.dumps(h,indent=2))
    print('seconds',time.monotonic()-start)
