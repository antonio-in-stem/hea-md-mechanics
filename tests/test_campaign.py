from pathlib import Path
import json
import shutil
import numpy as np
import pytest
from hea_md.campaign import load_campaign,calculate,export,stage,resolve_inside
from hea_md.cli import compare_results,main
from hea_md.io import sha256_file,load_curve

EXPECTED={
 'eq':(136.59949909291953,4.330156661888059,5.25312875597326),
 'zr24':(130.39953024687136,4.240413937815951,5.00605727384936),
 'zr28':(125.59144700580251,4.196665734040464,4.96721503734773),
 'hf16':(137.75647980505246,4.601311564826266,5.64690051427443),
 'hf12':(138.98147498367473,4.716089899122435,5.80298557484108),
 'ta16':(130.18380002101208,3.92043651130245,4.69088472347828),
 'ta12':(122.337191784548,3.7512178082963,4.17166611306709),
 'nb24':(145.1934742465064,4.718626790270606,5.77899793087597),
 'nb28':(153.96821256601586,4.924986902333556,6.3858758238348)}

@pytest.mark.parametrize('cid',list(EXPECTED))
def test_numerical_regression(campaign_result,cid):
 r=next(r for r in campaign_result['summary'] if r['case_id']==cid)
 np.testing.assert_allclose([r['slope_GPa'],r['proof_stress_GPa'],r['peak_common_GPa']],EXPECTED[cid],rtol=1e-10,atol=1e-11)

@pytest.mark.parametrize('cid',list(EXPECTED))
def test_independent_lstsq(root,cid):
 c=next(c for c in load_campaign(root)['cases'] if c['id']==cid)
 a=np.loadtxt(root/'data/raw'/c['curve_file'],delimiter=',',comments='#')
 mask=(a[:,0]>=-1e-12)&(a[:,0]<=.02+1e-12);x,y=a[mask].T
 slope=np.linalg.lstsq(np.column_stack([x,np.ones(len(x))]),y,rcond=None)[0][0]
 assert slope==pytest.approx(EXPECTED[cid][0],rel=1e-11)

def test_stored_analysis(root,campaign_result):
 compare_results(json.loads((root/'results/analysis.json').read_text()),campaign_result)

def test_coverage_and_sensitivity(campaign_result):
 assert len(campaign_result['cases'])==9
 assert sum(r['input']['samples'] for r in campaign_result['cases'])==16009
 assert len(campaign_result['sensitivity'])==540
 assert all(r['proof_stress']['status']=='found' for r in campaign_result['cases'])
 assert all(r['elastic_fit']['n_fit']==201 for r in campaign_result['cases'])
 assert sum(r['proof_GPa'] is None for r in campaign_result['sensitivity'])==3

def test_zirconium_identity(root):
 m={c['id']:c for c in load_campaign(root)['cases']}
 assert m['zr28']['curve_file']=='strain_stress_1.csv'
 assert m['zr24']['curve_file']=='strain_stress_2.csv'

def clone(root,tmp_path):
 for d in ['configs','simulation','data']:shutil.copytree(root/d,tmp_path/d)
 return tmp_path

@pytest.mark.parametrize('change',['missing','modified','duplicate_id','wrong_composition'])
def test_case_identity_fail_closed(root,tmp_path,change):
 p=clone(root,tmp_path);c=p/'configs/campaign.json';m=json.loads(c.read_text())
 if change=='missing':(p/'data/raw/strain_stress_3.csv').unlink()
 elif change=='modified':(p/'data/raw/strain_stress_3.csv').write_text('0,0\n.01,1\n.02,2\n')
 elif change=='duplicate_id':m['cases'][1]['id']=m['cases'][0]['id'];c.write_text(json.dumps(m))
 else:m['cases'][0]['nominal_atomic_percent']['Hf']=21;c.write_text(json.dumps(m))
 with pytest.raises(ValueError):load_campaign(p)

def test_path_escape(tmp_path):
 with pytest.raises(ValueError):resolve_inside(tmp_path,'../outside')

def test_export_no_overwrite(root,tmp_path):
 out=tmp_path/'out';out.mkdir();(out/'marker').write_text('x')
 with pytest.raises(FileExistsError):export(root,out)

def test_output_not_raw(root):
 with pytest.raises(ValueError):export(root,root/'data/raw')

def test_stage_requires_potentials(root,tmp_path):
 with pytest.raises(FileNotFoundError):stage(root,'eq',tmp_path/'run',tmp_path/'potential')

def test_stage_no_overwrite(root,tmp_path):
 p=tmp_path/'run';p.mkdir()
 with pytest.raises(FileExistsError):stage(root,'eq',p,tmp_path)

def test_unknown_case(root,tmp_path):
 with pytest.raises(ValueError):stage(root,'something',tmp_path/'run',tmp_path)

def test_stage_records_but_does_not_run(root,tmp_path):
 pot=tmp_path/'pot';pot.mkdir()
 # Format-only fixture, not a usable interatomic model.
 (pot/'library.meam').write_text("'Hf' 'Nb' 'Ta' 'Ti' 'Zr'\n")
 (pot/'HfNbTaTiZr.meam').write_text('\n'.join(f'nn2({i},{j})=1' for i in range(1,6) for j in range(i,6)))
 out=tmp_path/'run';r=stage(root,'eq',out,pot)
 assert r['lammps_executed'] is False
 assert (out/'in.compression.lmp').is_file()
 assert len(list((out/'outputs').iterdir()))==0
 assert sha256_file(pot/'library.meam')==r['potential_sha256']['library.meam']

def test_equal_config_paired_trends(campaign_result):
 rs=campaign_result['paired_sensitivity']
 nb=next(r for r in rs if r['case_id']=='nb28' and r['metric']=='slope_GPa')
 assert nb['minimum_change_percent']>12
 hf=next(r for r in rs if r['case_id']=='hf16' and r['metric']=='slope_GPa')
 assert hf['minimum_change_percent']<0<hf['maximum_change_percent']

def test_comparison_rejects_false_results():
 with pytest.raises(ValueError):compare_results({'x':1.},{'x':1.1})
 with pytest.raises(ValueError):compare_results({'x':[]},{'x':[1]})
