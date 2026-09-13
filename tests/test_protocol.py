import re
import pytest
from hea_md.campaign import load_campaign

@pytest.mark.parametrize('cid',['eq','zr24','zr28','hf16','hf12','ta16','ta12','nb24','nb28'])
def test_declared_protocol(root,cid):
 c=next(c for c in load_campaign(root)['cases'] if c['id']==cid)
 text='\n'.join(line.split('#',1)[0].strip() for line in (root/c['input_file']).read_text().splitlines())
 assert 'units metal' in text and 'boundary p p p' in text
 assert 'timestep 0.001' in text
 assert 'orient z 0 0 1' in text
 assert 'variable strainrate equal -1.0e9' in text
 assert 'v_strainrate/1.0e12' in text
 assert 'z erate ${strainrate1} units box remap x' in text
 assert '((0.0001*pzz)*(lx*ly))/(v_ly0*v_lx0)' in text
 assert 'x 0.0 0.0 $(1000*dt) y 0.0 0.0 $(1000*dt)' in text
 runs=[int(x) for x in re.findall(r'^run (\d+)',text,re.M)]
 assert runs[:6]==[30000,300000,30000,300000,300000,50000]
 assert sum(runs[:6])*.001==1010
 assert runs[-1]==(100000 if cid in ['eq','zr28'] else 200000)
 selections=[float(x) for x in re.findall(r'^set type [1-4] type/fraction [2-5] ([\d.]+) 1234567',text,re.M)]
 assert len(selections)==4
 f2,f3,f4,f5=selections
 for k,v in zip(['Hf','Nb','Ta','Ti','Zr'],[1-f2,f2-f3,f3-f4,f4-f5,f5]):
  assert v*100==pytest.approx(c['nominal_atomic_percent'][k])

def test_no_z_barostat_during_deformation(root):
 c=load_campaign(root)['cases'][0];text=(root/c['input_file']).read_text()
 line=next(line for line in text.splitlines() if 'fix' in line and 'npt temp 300.0' in line)
 assert ' z ' not in line and ' iso ' not in line
