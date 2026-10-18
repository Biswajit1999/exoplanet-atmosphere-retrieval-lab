import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from spectral_robustness import templates,design_matrix,weighted_fit,floor_sensitivity,leave_one_instrument_out


def synthetic_points(noise=0.0,seed=12):
    rng=np.random.default_rng(seed)
    wave=np.linspace(.7,5.1,180)
    h,c,s=templates(wave)
    out=[]
    for i,w in enumerate(wave):
        group=['A','B','C'][i%3]
        offset={'A':20000,'B':20500,'C':19800}[group]
        y=offset+80*np.log(w/2.2)+250*h[i]+500*c[i]-70*s[i]+rng.normal(0,noise)
        out.append({'wavelength_um':float(w),'transit_depth_ppm':float(y),
                    'uncertainty_ppm':30.0,'instrument':group})
    return out


def test_template_peak():
    _,co2,_=templates([4.30])
    assert np.isclose(co2[0],1.0)


def test_fit_recovers_injected_templates():
    r=weighted_fit(synthetic_points())
    assert abs(r['coefficients']['H2O_template']-250)<1e-7
    assert abs(r['coefficients']['CO2_template']-500)<1e-7
    assert abs(r['coefficients']['SO2_nuisance']+70)<1e-7
    assert r['chi2']<1e-15


def test_slope_and_offsets():
    r=weighted_fit(synthetic_points())
    assert abs(r['coefficients']['slope']-80)<1e-7
    assert abs(r['coefficients']['offset:A']-20000)<1e-7


def test_error_floor_increases_formal_uncertainty():
    p=synthetic_points(noise=12)
    a=weighted_fit(p,error_floor=0)
    b=weighted_fit(p,error_floor=100)
    assert b['formal_1sigma']['H2O_template']>a['formal_1sigma']['H2O_template']
    assert b['chi2']<a['chi2']


def test_floor_sweep():
    f=floor_sensitivity(synthetic_points(noise=20),floors=(0,50,100))
    assert len(f)==3
    assert f[0]['chi2']>=f[1]['chi2']>=f[2]['chi2']


def test_omission_recovery_noiseless():
    out=leave_one_instrument_out(synthetic_points())
    assert len(out['omissions'])==3
    assert all(x['identifiable'] for x in out['omissions'])
    assert all(abs(x['delta_CO2_ppm'])<1e-6 for x in out['omissions'])


def test_bad_sigma():
    p=synthetic_points();p[0]['uncertainty_ppm']=0
    with pytest.raises(ValueError):weighted_fit(p)


def test_singular_design():
    p=synthetic_points()
    for q in p:q['wavelength_um']=2.0
    with pytest.raises(ValueError):weighted_fit(p)


def test_negative_floor_rejected():
    with pytest.raises(ValueError):weighted_fit(synthetic_points(),error_floor=-1)
