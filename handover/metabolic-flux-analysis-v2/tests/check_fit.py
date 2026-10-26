"""Evaluate the reported dynamic map against trusted time-course measurements."""
import json,math,sys
from pathlib import Path
import numpy as np
import dynamic_model as D


def weighted_ssr(artifact):
    table=artifact["fluxes"]
    net={r:v["net"] for r,v in table.items()}
    exchange={r:min(v["forward"],v["reverse"]) for r,v in table.items()}
    net={r:max(v,0) if not D.M.REACTIONS[r][2] else v for r,v in net.items()}
    known,unknown,schedule,meas,rates=D.read_data(Path(__file__).with_name("data"))
    pools={**known,**artifact["pools"]}
    pred=D.predictions(net,exchange,pools,schedule,meas,rtol=2e-8,atol=1e-11)
    residuals=[(p-float(r["fraction"]))/float(r["sd"]) for p,r in zip(pred,meas)]
    residuals += [(table[r["reaction_id"]]["net"]-float(r["rate"]))/float(r["sd"]) for r in rates if float(r["sd"])>0]
    assert all(math.isfinite(float(v)) for v in residuals)
    return math.fsum(float(v)**2 for v in residuals)

if __name__=='__main__': print(weighted_ssr(json.load(sys.stdin)))
