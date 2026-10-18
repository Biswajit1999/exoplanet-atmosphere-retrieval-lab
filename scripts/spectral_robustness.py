"""Independent conditional WLS sensitivity tests for WASP-39b benchmark bins.

Gaussian feature templates are not opacities; coefficients are not abundances.
This is an unconstrained linear fit, unlike the browser's nonnegative grid.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def templates(wavelength, width_scale=1.0):
    w=np.asarray(wavelength,dtype=float)
    if not np.all(np.isfinite(w)) or np.any(w <= 0) or width_scale <= 0:
        raise ValueError('Positive finite wavelength and width_scale required')
    def g(c, width):
        return np.exp(-0.5*((w-c)/(width*width_scale))**2)
    return (0.62*g(1.15,.08)+g(1.40,.12)+0.42*g(1.90,.14)+0.48*g(2.70,.18),
            g(4.30,.16),g(4.05,.07))


def design_matrix(points, include_slope=True, width_scale=1.0):
    if len(points) < 10:
        raise ValueError('At least ten bins required')
    wavelength=np.asarray([p['wavelength_um'] for p in points],float)
    groups=sorted({p['instrument'] for p in points})
    h2o,co2,so2=templates(wavelength,width_scale)
    columns=[np.array([float(p['instrument']==g) for p in points]) for g in groups]
    names=['offset:'+g for g in groups]
    if include_slope:
        columns.append(np.log(wavelength/2.2));names.append('slope')
    columns.extend([h2o,co2,so2]);names.extend(['H2O_template','CO2_template','SO2_nuisance'])
    return np.column_stack(columns),names


def weighted_fit(points, error_floor=0.0, include_slope=True, width_scale=1.0):
    """Unconstrained WLS with a numerical rank check and formal covariance."""
    if not np.isfinite(error_floor) or error_floor < 0:
        raise ValueError('error_floor must be nonnegative and finite')
    X,names=design_matrix(points,include_slope,width_scale)
    y=np.asarray([p['transit_depth_ppm'] for p in points],float)
    sig=np.asarray([p['uncertainty_ppm'] for p in points],float)
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(sig)) or np.any(sig <= 0):
        raise ValueError('Finite depths and positive uncertainties required')
    sig=np.hypot(sig,error_floor)
    A=X/sig[:,None]
    beta,_,rank,s=np.linalg.lstsq(A,y/sig,rcond=None)
    if rank != len(names) or s[-1] <= s[0]*1e-10:
        raise ValueError('Singular or poorly identified spectral design')
    residual=y-X@beta
    cov=np.linalg.inv(A.T@A)
    dof=len(y)-len(names)
    return {'n_bins':len(y),'dof':dof,'chi2':float(np.sum((residual/sig)**2)),
            'reduced_chi2':float(np.sum((residual/sig)**2)/dof),
            'condition_number':float(s[0]/s[-1]),
            'coefficients':dict(zip(names,map(float,beta))),
            'formal_1sigma':dict(zip(names,map(float,np.sqrt(np.diag(cov))))),
            'warning':'Unconstrained conditional template fit, not molecular abundance retrieval'}


def floor_sensitivity(points, floors=(0,50,100,200,300)):
    """Quantify how extra unverified white noise changes fit and uncertainty."""
    return [{'floor_ppm':float(f),**weighted_fit(points,error_floor=f)} for f in floors]


def leave_one_instrument_out(points):
    """Measure influence of each product on conditional template coefficients."""
    groups=sorted({p['instrument'] for p in points})
    baseline=weighted_fit(points)
    out=[]
    for group in groups:
        kept=[p for p in points if p['instrument']!=group]
        if len(kept) < 10:
            continue
        try:
            fit=weighted_fit(kept)
        except ValueError:
            out.append({'excluded':group,'identifiable':False})
            continue
        out.append({'excluded':group,'identifiable':True,'n_bins':len(kept),
                    'delta_H2O_ppm':fit['coefficients']['H2O_template']-baseline['coefficients']['H2O_template'],
                    'delta_CO2_ppm':fit['coefficients']['CO2_template']-baseline['coefficients']['CO2_template'],
                    'reduced_chi2':fit['reduced_chi2']})
    return {'baseline':baseline,'omissions':out}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('benchmark',type=Path)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    data=json.loads(args.benchmark.read_text())
    report={'source':data.get('dataset'), 'covariance_assumption':'independent bins, optional scalar floor',
            'floor_sensitivity':floor_sensitivity(data['points']),
            'instrument_influence':leave_one_instrument_out(data['points'])}
    text=json.dumps(report,indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+'\n')
    else:
        print(text)

if __name__=='__main__':
    main()
