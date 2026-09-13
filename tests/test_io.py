import numpy as np
import pytest
from hea_md.io import Curve,load_curve,write_json,sha256_file

@pytest.mark.parametrize('text',[
 'strain,stress\n0,-.1\n.01,1.2\n.02,2.5\n',
 '# Fix print output\n0,-.1\n.01,1.2\n.02,2.5\n',
 '0 -.1\n.01 1.2\n.02 2.5\n'])
def test_read_retains_stress(tmp_path,text):
 p=tmp_path/'c.csv';p.write_text(text);c=load_curve(p)
 assert c.stress[0]==-.1
 assert c.sha256==sha256_file(p)
 assert len(c.strain)==3

@pytest.mark.parametrize('text',[
 '0,1\n.01,nan\n.02,3','0,1\n.01,inf\n.02,3','0,1\n.01,2,3\n.02,4',
 '0,1\n.01,2\n.01,3','0,1\n.02,2\n.01,3','0,1\n-0.01,2\n.02,3',
 '0,1\n.01,2','strange,header\n0,1\n.01,2\n.02,3',
 'strain,stress\nstrain,stress\n0,1\n.01,2\n.02,3',
 '0,1\nstrain,stress\n.01,2\n.02,3','0,1\n.5,2\n1,3'])
def test_bad_inputs_rejected(tmp_path,text):
 p=tmp_path/'c.csv';p.write_text(text)
 with pytest.raises(ValueError):load_curve(p)

@pytest.mark.parametrize('origin',[-8e-16,7e-16,0.])
def test_roundoff_is_logged(tmp_path,origin):
 p=tmp_path/'c.csv';p.write_text(f'{origin},-.1\n.01,1\n.02,2\n')
 before=p.read_bytes();c=load_curve(p)
 assert c.strain[0]==0
 assert bool(c.origin_adjustments)==(origin!=0)
 assert p.read_bytes()==before

@pytest.mark.parametrize('unit,factor',[('GPa',1),('MPa',1000),('Pa',1e9),('bar',1e4)])
def test_units(tmp_path,unit,factor):
 p=tmp_path/'c.csv';p.write_text(f'0,0\n1,{factor}\n2,{2*factor}\n')
 c=load_curve(p,strain_unit='percent',stress_unit=unit)
 np.testing.assert_allclose(c.strain,[0,.01,.02]);np.testing.assert_allclose(c.stress,[0,1,2])

@pytest.mark.parametrize('bad',[(np.array([0,.1]),np.array([0,1])),(np.array([0,.1,.2]),np.array([0,1])),(np.array([[0,.1,.2]]),np.array([[0,1,2]]))])
def test_curve_shape(bad):
 with pytest.raises(ValueError):Curve(*bad).validate()

def test_nonfinite_json(tmp_path):
 with pytest.raises(ValueError):write_json(tmp_path/'c.json',{'value':float('nan')})
