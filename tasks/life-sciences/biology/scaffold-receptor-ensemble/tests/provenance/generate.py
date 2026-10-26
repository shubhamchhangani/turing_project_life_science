"""Original competing receptor sheet with effective scaffold links and local biases."""
import json
from pathlib import Path
import numpy as np


def model(side=32):
    rng=np.random.default_rng(2026102603)
    xy=[[x,y] for y in range(side) for x in range(side)];edges=[]
    for y in range(side):
        for x in range(side):
            a=y*side+x
            for dx,dy in [(1,0),(0,1),(1,1)]:
                if x+dx<side and y+dy<side and (dx!=dy or (x+2*y)%3==0):
                    coupling=round(float(rng.uniform(.3,.7)),7)
                    if rng.random()<.28:coupling=-coupling
                    edges.append([a,(y+dy)*side+x+dx,coupling])
    positions=[(side//2+dx,side//2+dy) for dx,dy in
               [(-5,-4),(-1,-5),(4,-4),(-5,1),(0,0),(5,1),
                (-3,-3),(2,-3),(-3,0),(3,0),(-4,4),(-1,3),(2,4),(5,4),
                (-5,-1),(0,-2),(4,-1),(-2,5),(1,1),(3,6),
                (-6,3),(6,-3),(-3,-6),(2,-6),(-6,-5),(6,6)]]
    sites=[y*side+x for x,y in positions]
    pairs=[(a,b) for a in sites for b in sites if a<b and sum(abs(u-v) for u,v in zip(xy[a],xy[b]))>=4]
    rng.shuffle(pairs)
    for a,b in pairs[:36]:edges.append([a,b,round(float(rng.uniform(.2,.55)*rng.choice([-1,1])),7)])
    fields=[[s,round(float(rng.uniform(-.4,.4)),7)] for s in sites]
    return {'coordinates':xy,'edges':edges,'fields':fields,'ligand_sites':sites[:6],
            'ligand_spins':[1,-1,1,-1,1,-1],'reporter_sites':sites[6:14]}


if __name__=='__main__':
    p=Path(__file__).parent/'draft/model.json';p.write_text(json.dumps(model(),separators=(',',':'))+'\n')
