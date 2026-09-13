"""Transformation properties, empty outcomes, configuration coupling and provenance."""
import importlib.util
import json
from pathlib import Path
import shutil
from decimal import Decimal, localcontext
import numpy as np
import pytest

from hea_md.io import Curve
from hea_md.mechanics import AnalysisConfig, analyze
from hea_md.campaign import calculate, configuration
from hea_md.plotting import offset_series


def bilinear():
    x=np.linspace(0,.1,1001)
    return Curve(x,np.where(x<=.04,130*x,5.2+10*(x-.04)))


@pytest.mark.parametrize('factor',[.1,.5,2.,10.])
def test_stress_scaling_preserves_strains(factor):
    c=bilinear();a=analyze(c);b=analyze(Curve(c.strain,c.stress*factor))
    assert b['elastic_fit']['slope_GPa']==pytest.approx(factor*a['elastic_fit']['slope_GPa'])
    for field in ['proof_stress']:
        assert b[field]['point']['strain']==pytest.approx(a[field]['point']['strain'])
        assert b[field]['point']['stress_GPa']==pytest.approx(factor*a[field]['point']['stress_GPa'])
    assert b['peak_common_window']['strain']==a['peak_common_window']['strain']


@pytest.mark.parametrize('shift',[-3.,-.1,1.,10.])
def test_stress_translation_preserves_slope_and_crossing(shift):
    c=bilinear();a=analyze(c);b=analyze(Curve(c.strain,c.stress+shift))
    assert b['elastic_fit']['slope_GPa']==pytest.approx(a['elastic_fit']['slope_GPa'])
    assert b['elastic_fit']['intercept_GPa']==pytest.approx(a['elastic_fit']['intercept_GPa']+shift)
    assert b['proof_stress']['point']['strain']==pytest.approx(a['proof_stress']['point']['strain'])


def test_common_window_independent_of_later_samples():
    c=bilinear();x=np.linspace(.1001,.2,1000)
    b=analyze(Curve(np.r_[c.strain,x],np.r_[c.stress,np.full(1000,300.)]))
    a=analyze(c)
    for key in ['elastic_fit','proof_stress','peak_common_window','peak_averaged_common_window']:
        assert a[key]==b[key]
    assert b['peak_full_record']['stress_GPa']==300.


@pytest.mark.parametrize('mean,hi,offset',[(1,.01,.002),(11,.015,.004),(31,.02,.003)])
def test_figure_geometry_follows_config(mean,hi,offset):
    c=bilinear();cfg=AnalysisConfig(fit_max=hi,offset=offset,smooth_points=mean)
    r=analyze(c,cfg);p=offset_series(c,r)
    np.testing.assert_allclose(p['elastic']-p['offset'],r['elastic_fit']['slope_GPa']*offset)
    assert len(p['x'])==1001-mean+1
    assert p['x'][0]==pytest.approx((mean//2)*.0001)
    assert p['point']==r['proof_stress']['point']


@pytest.mark.parametrize('key,value',[
    ('sensitivity_fit_max',[]),('sensitivity_fit_max',[.01,.01]),('sensitivity_fit_max',[.11]),
    ('sensitivity_fit_max',[True]),('sensitivity_smooth_points',[2]),('sensitivity_confirm_points',[0]),
    ('sensitivity_confirm_points',[1.1]),('sensitivity_smooth_points',[float('nan')])])
def test_sensitivity_grid_rejects_invalid_settings(root,tmp_path,key,value):
    (tmp_path/'configs').mkdir()
    obj=json.loads((root/'configs/analysis.json').read_text());obj[key]=value
    (tmp_path/'configs/analysis.json').write_text(json.dumps(obj))
    with pytest.raises(ValueError):configuration(tmp_path)


def test_all_linear_records_produce_undefined_ranges(root,monkeypatch):
    import hea_md.campaign as module
    real=module.load_curve
    def linear(path):
        curve=real(path)
        return Curve(curve.strain,130*curve.strain-.1,curve.source,curve.sha256,curve.origin_adjustments)
    monkeypatch.setattr(module,'load_curve',linear)
    obj=calculate(root)
    assert all(row['proof_GPa'] is None for row in obj['sensitivity'])
    assert all(row['proof_min_GPa'] is None and row['proof_max_GPa'] is None for row in obj['method_sensitivity_ranges'])



