import numpy as np
import pytest
from hea_md.io import Curve
from hea_md.mechanics import AnalysisConfig,analyze,elastic_fit,centered_mean,engineering_stress_from_pressure

def bilinear():
 x=np.linspace(0,.1,1001)
 return Curve(x,np.where(x<=.04,130*x,5.2+10*(x-.04)))

@pytest.mark.parametrize('mean',[1,11,21,31])
@pytest.mark.parametrize('confirm',[1,3,5])
def test_analytical_offset(mean,confirm):
 r=analyze(bilinear(),AnalysisConfig(smooth_points=mean,confirm_points=confirm))
 assert r['elastic_fit']['slope_GPa']==pytest.approx(130)
 p=r['proof_stress']['point'];exact=.04+130*.002/(130-10)
 assert p['strain']==pytest.approx(exact,abs=1e-12)
 assert p['stress_GPa']==pytest.approx(5.2+10*(exact-.04),abs=1e-12)

def test_linear_has_no_proof():
 x=np.linspace(0,.1,1001);r=analyze(Curve(x,130*x-.1))
 assert r['proof_stress']['point'] is None
 assert r['peak_common_window']['at_boundary']
 assert r['elastic_fit']['intercept_GPa']==pytest.approx(-.1)

@pytest.mark.parametrize('spike',[.4,-.4,1.,-1.])
def test_isolated_excursions(spike):
 c=bilinear();c.stress[330]+=spike
 p=analyze(c)['proof_stress']['point']
 assert p['strain']>.04

def test_directed_crossing_ignores_upward_spike():
 c=bilinear();c.stress[330]+=.4
 r=analyze(c,AnalysisConfig(smooth_points=1,confirm_points=1))
 assert r['proof_stress']['point']['strain']>.04

def test_persistence_rejects_downward_spike():
 c=bilinear();c.stress[330]-=.4
 r=analyze(c,AnalysisConfig(smooth_points=1,confirm_points=3))
 assert r['proof_stress']['point']['strain']>.04

def test_multiple_loading_peaks_are_separated():
 x=np.linspace(0,.2,2001);y=np.where(x<.08,130*x,1+100*(x-.08))
 r=analyze(Curve(x,y))
 assert r['peak_full_record']['stress_GPa']>r['peak_common_window']['stress_GPa']
 assert r['peak_full_record']['strain']>.1

@pytest.mark.parametrize('kwargs',[{'fit_max':0},{'fit_min':-.1},{'fit_max':.2},{'offset':0},{'analysis_max':1},
 {'offset':float('nan')},{'smooth_points':2},{'smooth_points':True},{'confirm_points':0},{'confirm_points':1.2}])
def test_invalid_config(kwargs):
 with pytest.raises(ValueError):AnalysisConfig(**kwargs).validate()

def test_short_input_not_extrapolated():
 with pytest.raises(ValueError):analyze(Curve(np.linspace(0,.05,501),np.linspace(0,5,501)))

def test_fit_coverage_and_sign():
 with pytest.raises(ValueError):elastic_fit(bilinear(),0,.2)
 with pytest.raises(ValueError):elastic_fit(Curve(np.array([0,.01,.02]),np.array([0,-1,-2])),0,.02)

def test_irregular_grid_rejects_sample_mean():
 c=Curve(np.array([0,.01,.021,.03,.04,.05]),np.arange(6.))
 with pytest.raises(ValueError):centered_mean(c,3)

def test_mean_has_no_padded_edges():
 c=bilinear();x,y=centered_mean(c,21)
 assert len(x)==981 and x[0]==pytest.approx(.001)
 assert y[0]==pytest.approx(np.mean(c.stress[:21]))

def test_fit_matches_lstsq():
 c=bilinear();x=c.strain[:201];y=c.stress[:201]
 e,b=np.linalg.lstsq(np.column_stack([x,np.ones(len(x))]),y,rcond=None)[0]
 r=elastic_fit(c,0,.02)
 assert r['slope_GPa']==pytest.approx(e)
 assert r['intercept_GPa']==pytest.approx(b,abs=1e-12)

@pytest.mark.parametrize('pressure,area,initial,expected',[(10000,1,1,1),(20000,1.1,1,2.2),(-10000,2,1,-2)])
def test_engineering_stress(pressure,area,initial,expected):
 assert engineering_stress_from_pressure(pressure,area,initial)==pytest.approx(expected)

@pytest.mark.parametrize('params',[(1,0,1),(1,1,0),(float('inf'),1,1)])
def test_bad_areas(params):
 with pytest.raises(ValueError):engineering_stress_from_pressure(*params)
