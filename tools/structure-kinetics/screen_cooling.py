"""Check coupled forward physics; this is not an inverse solution or model trial."""
import json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

TREF=300.
ATTEMPT=1e6
ENERGIES=np.array([0.,.25,.55,.85])
BARRIERS=np.array([[0.,3.4,4.6,5.2],[3.4,0.,3.8,4.8],[4.6,3.8,0.,3.6],[5.2,4.8,3.6,0.]])

def generator(temperature):
    # Column probability convention; symmetric transition-state barriers imply
    # detailed balance with pi proportional to exp(-E*TREF/temperature).
    q=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if i!=j:q[j,i]=ATTEMPT*np.exp((ENERGIES[i]-BARRIERS[i,j])*TREF/temperature)
    q[np.diag_indices(4)]=-q.sum(axis=0)
    return q

def stationary(t):
    p=np.exp(-ENERGIES*TREF/t);return p/p.sum()

def temperature(time,duration):
    return 300.-220.*time/duration

def midpoint(duration,steps):
    p=stationary(300.)
    for j in range(steps):p=expm(generator(temperature((j+.5)*duration/steps,duration))*duration/steps)@p
    return p

if __name__=='__main__':
    started=time.monotonic();records=[]
    for duration in [4e-5,6e-4]:
        adaptive=solve_ivp(lambda t,p:generator(temperature(t,duration))@p,(0.,duration),stationary(300.),rtol=2e-11,atol=2e-13,method='DOP853').y[:,-1]
        coarse=midpoint(duration,512);fine=midpoint(duration,1024)
        assert np.max(np.abs(adaptive-fine))<3e-7
        assert np.max(np.abs(fine-adaptive))<np.max(np.abs(coarse-adaptive))
        records.append(dict(cooling_seconds=duration,frozen_populations=adaptive.tolist(),
            midpoint_1024_max_error=float(np.max(np.abs(adaptive-fine))),
            total_variation_from_room_temperature=float(np.abs(adaptive-stationary(300.)).sum()/2),
            total_variation_from_final_temperature=float(np.abs(adaptive-stationary(80.)).sum()/2)))
    for t in [80.,150.,300.]:
        q=generator(t);pi=stationary(t)
        assert np.max(np.abs(q@pi))<1e-10
        assert np.max(np.abs(q*pi[None,:]-(q*pi[None,:]).T))<1e-10
    result=dict(status='forward_coupling_screen_passed_not_inverse_validation',records=records,
                reference_temperature_population=stationary(300.).tolist(),
                final_temperature_equilibrium=stationary(80.).tolist(),seconds=time.monotonic()-started)
    out=Path(__file__).parent/'draft';out.mkdir(exist_ok=True)
    (out/'cooling-screen.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
