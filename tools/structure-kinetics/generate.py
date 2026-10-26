"""Generate author screening data; public observations and private truth are separate."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from model import coordinates,frozen,projection,photon_model

def config():
    return dict(schema='author-screen-1',states=4,
        core=[[-17,-9,-4],[-13,11,6],[-4,-12,13],[3,8,-13],[7,-3,3],[-10,0,-15]],
        mobile_domain=[[25,0,0],[27,6,3],[22,-5,6]],hinge=[3,1,0],
        bead_weights=[1.5,1.2,1.1,1.4,.9,1.3,1.,.8,.7],
        gaussian_sigma=2.7,image_size=32,pixel_size=2.5,transfer_decay=9.,image_noise_sd=.00035,
        reference_temperature=300.,freeze_temperature=80.,attempt_rate=1e6,
        cooling_seconds=[4e-5,6e-4],label_pairs=[[0,6],[1,6],[2,6],[3,7]],
        forster_distance=35.,dye_on=3500.,dye_off=1200.,photon_background=700.,
        photon_brightness=140000.,acceptor_detection=.82,direct_acceptor_brightness=45000.,
        temperatures=[285.,315.],trace_seconds=.008,
        geometry_contract='azimuth followed by elevation: scipy intrinsic/extrinsic convention explicitly z then y extrinsic; xyz=R_y(elevation) R_z(azimuth) local + hinge',
        units=dict(length='angstrom',time='second',temperature='kelvin',angle='radian'),
        parameter_bounds=dict(azimuth=[-2.7,2.7],elevation=[-1.2,1.2],energies=[-.5,1.5],barriers=[2.,6.5]))

def simulate(rng,Q,pi,rates,duration):
    state=int(rng.choice(len(pi),p=pi));t=0.;times=[];channels=[]
    changes=0
    while True:
        transitions=Q[:,state].copy();transitions[state]=0
        hazards=np.r_[transitions,rates[:,state]];total=hazards.sum()
        t+=rng.exponential(1/total)
        if t>=duration:break
        choice=int(rng.choice(len(hazards),p=hazards/total))
        if choice<len(pi):state=choice;changes+=1
        else:times.append(t);channels.append(choice-len(pi))
    return np.array(times),np.array(channels,dtype=np.int8),changes

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--images',type=int,default=128);parser.add_argument('--replicates',type=int,default=3)
    parser.add_argument('--duration',type=float,default=.008);parser.add_argument('--output',type=Path,default=Path(__file__).parent/'draft');a=parser.parse_args()
    assert a.images>0 and a.replicates>0 and a.duration>0
    started=time.monotonic();root=a.output;public=root/'public';private=root/'private'
    for p in [public,private]:p.mkdir(parents=True,exist_ok=True)
    if (public/'manifest.json').exists():raise SystemExit('Screen already exists; use a new directory for a new dataset')
    c=config();c['trace_seconds']=a.duration;rng=np.random.default_rng(2026102607)
    # Never import these parameters in the inverse solver.
    angles=np.array([[-1.2,-.55],[-.3,.7],[.85,-.35],[1.8,.65]])
    energies=np.array([.25,.55,.85]);barriers=np.array([3.4,4.6,5.2,3.8,4.8,3.6])
    xyz=coordinates(c,angles);weights=[frozen(c,energies,barriers,d) for d in c['cooling_seconds']]
    rotations=Rotation.random(a.images,random_state=rng).as_matrix();shifts=rng.uniform(-2,2,(a.images,2));defocus=rng.uniform(35,90,a.images)
    protocol=np.arange(a.images)%2;classes=np.array([rng.choice(4,p=weights[g]) for g in protocol])
    images=np.array([projection(c,xyz[k],r,s,d)+rng.normal(0,c['image_noise_sd'],(c['image_size'],)*2) for k,r,s,d in zip(classes,rotations,shifts,defocus)],dtype=np.float32)
    np.savez_compressed(public/'images.npz',images=images,defocus=defocus,protocol=protocol)
    metadata=[];all_t=[];all_ch=[];all_ids=[];event_counts=[]
    for T in c['temperatures']:
        for pair in range(len(c['label_pairs'])):
            Q,pi,rates=photon_model(c,angles,energies,barriers,T,pair)
            for rep in range(a.replicates):
                tid=len(metadata);ts,cs,n=simulate(rng,Q,pi,rates,c['trace_seconds'])
                metadata.append(dict(trace_id=tid,temperature=T,label_pair=pair,duration=c['trace_seconds']))
                all_t.extend(ts);all_ch.extend(cs);all_ids.extend([tid]*len(ts));event_counts.append(n)
    np.savez_compressed(public/'photons.npz',time=np.array(all_t),channel=np.array(all_ch,dtype=np.int8),trace_id=np.array(all_ids,dtype=np.int16))
    (public/'config.json').write_text(json.dumps(c,indent=2)+'\n');(public/'traces.json').write_text(json.dumps(metadata,indent=2)+'\n')
    np.savez_compressed(private/'truth.npz',angles=angles,energies=energies,barriers=barriers,rotations=rotations,shifts=shifts,classes=classes,frozen_populations=weights)
    manifest=dict(status='prototype_dataset_not_a_task_package',images=a.images,traces=len(metadata),photons=len(all_t),
        seconds=time.monotonic()-started,files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(public.iterdir()) if p.is_file()})
    (public/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
