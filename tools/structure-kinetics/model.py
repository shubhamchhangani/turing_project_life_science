"""Public, parameterized forward physics; no planted parameters or inverse hints."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

def coordinates(config, angles):
    core=np.asarray(config['core']);domain=np.asarray(config['mobile_domain'])
    hinge=np.asarray(config['hinge']);result=[]
    for azimuth,elevation in angles:
        rotation=Rotation.from_euler('zy',[azimuth,elevation]).as_matrix()
        result.append(np.vstack([core,domain@rotation.T+hinge]))
    return np.asarray(result)

def kinetics(config,energies,barriers,temperature):
    e=np.r_[0.,energies];q=np.zeros((4,4));b=iter(barriers)
    for i in range(4):
        for j in range(i+1,4):
            barrier=next(b)
            q[j,i]=config['attempt_rate']*np.exp((e[i]-barrier)*config['reference_temperature']/temperature)
            q[i,j]=config['attempt_rate']*np.exp((e[j]-barrier)*config['reference_temperature']/temperature)
    q[np.diag_indices(4)]=-q.sum(axis=0)
    pi=np.exp(-e*config['reference_temperature']/temperature);pi/=pi.sum()
    return q,pi

def frozen(config,energies,barriers,duration):
    _,initial=kinetics(config,energies,barriers,config['reference_temperature'])
    def rhs(t,p):
        T=config['reference_temperature']+(config['freeze_temperature']-config['reference_temperature'])*t/duration
        return kinetics(config,energies,barriers,T)[0]@p
    return solve_ivp(rhs,(0,duration),initial,method='DOP853',rtol=2e-10,atol=1e-12).y[:,-1]

def projection(config,xyz,rotation,shift,defocus):
    n=config['image_size'];step=config['pixel_size'];sigma=config['gaussian_sigma']
    xy=(xyz@rotation.T)[:,:2]+np.asarray(shift)
    axis=(np.arange(n)-(n-1)/2)*step
    xx,yy=np.meshgrid(axis,axis,indexing='xy')
    weights=np.asarray(config['bead_weights'])
    image=np.sum(weights[:,None,None]*np.exp(-((xx[None]-xy[:,0,None,None])**2+(yy[None]-xy[:,1,None,None])**2)/(2*sigma*sigma)),axis=0)/(2*np.pi*sigma*sigma)
    freq=np.fft.fftfreq(n,d=step);fx,fy=np.meshgrid(freq,freq,indexing='xy');r2=fx*fx+fy*fy
    transfer=np.cos(np.pi*defocus*r2)*np.exp(-config['transfer_decay']*r2)
    return np.fft.ifft2(np.fft.fft2(image)*transfer).real

def photon_model(config,angles,energies,barriers,temperature,pair):
    q,pi=kinetics(config,energies,barriers,temperature)
    off=config['dye_off'];on=config['dye_on']
    dye=np.array([[-on,off],[on,-off]]) # state 0 dark, 1 bright
    Q=np.kron(q,np.eye(2))+np.kron(np.eye(4),dye)
    initial=np.kron(pi,[off/(on+off),on/(on+off)])
    xyz=coordinates(config,angles);a,b=config['label_pairs'][pair]
    distances=np.linalg.norm(xyz[:,a]-xyz[:,b],axis=1)
    efficiencies=1/(1+(distances/config['forster_distance'])**6)
    rates=np.empty((3,8));bg=config['photon_background'];brightness=config['photon_brightness']
    for k,E in enumerate(efficiencies):
        for bright in [0,1]:
            effective=E*bright;i=2*k+bright
            rates[0,i]=bg+brightness*(1-effective)
            rates[1,i]=bg+config['acceptor_detection']*brightness*effective
            rates[2,i]=bg+config['direct_acceptor_brightness']*bright
    return Q,initial,rates
