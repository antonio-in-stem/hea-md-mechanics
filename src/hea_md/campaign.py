"""Explicit composition-to-file mappings and deterministic campaign analysis."""
from __future__ import annotations
import csv
import itertools
import json
import math
from pathlib import Path
import shutil
import re
from .io import load_curve,sha256_file,write_json
from .mechanics import AnalysisConfig,analyze

ELEMENTS={"Hf","Nb","Ta","Ti","Zr"}


def resolve_inside(root:Path,relative:str)->Path:
    root=Path(root).resolve()
    path=(root/relative).resolve()
    if not path.is_relative_to(root): raise ValueError("Path escapes the project directory.")
    return path


def load_campaign(root:Path)->dict:
    root=Path(root).resolve()
    obj=json.loads((root/'configs/campaign.json').read_text(encoding='utf-8'))
    cases=obj.get('cases',[])
    required_units = {'strain_unit': 'fraction', 'stress_unit': 'GPa',
                      'strain_convention': 'engineering compression positive',
                      'stress_convention': 'engineering compression positive'}
    for key, value in required_units.items():
        if obj.get(key) != value:
            raise ValueError(f'Unsupported campaign unit or convention: {key}.')
    if obj.get('schema_version')!=1 or len(cases)!=9:raise ValueError("Expected nine campaign cases.")
    for key in ['id','curve_file','input_file']:
        if len({c[key] for c in cases})!=len(cases):raise ValueError(f"Duplicate {key}.")
    for c in cases:
        if not re.fullmatch(r'[a-z][a-z0-9]*',c['id']):raise ValueError("Invalid case ID.")
        if Path(c['curve_file']).name!=c['curve_file']:raise ValueError("Curve must be a local basename.")
        comp=c['nominal_atomic_percent']
        if set(comp)!=ELEMENTS or not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) and 0<=v<=100 for v in comp.values()) or abs(sum(comp.values())-100)>1e-9:
            raise ValueError("Invalid nominal composition.")
        for rel,key in [(f"data/raw/{c['curve_file']}",'sha256'),(c['input_file'],'input_sha256')]:
            p=resolve_inside(root,rel)
            if not p.is_file() or sha256_file(p)!=c[key]:raise ValueError(f"Missing or modified input: {rel}.")
    if sum(c['id']=='eq' for c in cases)!=1:raise ValueError("One equiatomic reference is required.")
    return obj


def configuration(root:Path)->tuple[AnalysisConfig,dict]:
    a=json.loads((Path(root)/'configs/analysis.json').read_text(encoding='utf-8'))
    cfg=AnalysisConfig(**{k:a[k] for k in AnalysisConfig.__dataclass_fields__})
    cfg.validate()
    grids = ['sensitivity_fit_max', 'sensitivity_smooth_points', 'sensitivity_confirm_points']
    allowed = set(AnalysisConfig.__dataclass_fields__) | set(grids)
    if set(a) != allowed:
        raise ValueError("Unexpected or missing analysis settings.")
    for key in grids:
        values = a[key]
        if not isinstance(values, list) or not values:
            raise ValueError(f"{key} must be a nonempty list.")
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in values):
            raise ValueError(f"{key} requires finite numbers.")
        if len(set(values)) != len(values):
            raise ValueError(f"{key} contains duplicate settings.")
    for hi, mean, confirm in itertools.product(*(a[key] for key in grids)):
        AnalysisConfig(cfg.fit_min, hi, cfg.analysis_max, cfg.offset, mean, confirm).validate()
    return cfg,a


def write_csv(path:Path,rows:list[dict])->None:
    if not rows:raise ValueError("No rows to write.")
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)


def calculate(root:Path)->dict:
    campaign=load_campaign(root)
    cfg,settings=configuration(root)
    case_results=[]; sensitivity=[]
    for c in campaign['cases']:
        curve=load_curve(Path(root)/'data/raw'/c['curve_file'])
        if len(curve.strain)!=c['samples'] or abs(curve.strain[-1]-c['maximum_strain'])>1e-12:
            raise ValueError(f"Unexpected coverage for {c['id']}.")
        result=analyze(curve,cfg);result['case_id']=c['id'];case_results.append(result)
        for hi,mean,n in itertools.product(settings['sensitivity_fit_max'],settings['sensitivity_smooth_points'],settings['sensitivity_confirm_points']):
            alt=AnalysisConfig(cfg.fit_min,hi,cfg.analysis_max,cfg.offset,mean,n)
            out=analyze(curve,alt);point=out['proof_stress']['point']
            sensitivity.append({'case_id':c['id'],'fit_max':hi,'smooth_points':mean,'confirm_points':n,
               'slope_GPa':out['elastic_fit']['slope_GPa'],'proof_status':out['proof_stress']['status'],
               'proof_GPa':None if point is None else point['stress_GPa'],
               'proof_strain':None if point is None else point['strain'],
               'averaged_peak_GPa':out['peak_averaged_common_window']['stress_GPa']})
    reference=next(r for r in case_results if r['case_id']=='eq')
    rows=[]
    for c,r in zip(campaign['cases'],case_results):
        fit=r['elastic_fit'];p=r['proof_stress']['point'];rp=reference['proof_stress']['point']
        rows.append({'case_id':c['id'],'label':c['label'],'samples':r['input']['samples'],
          'slope_GPa':fit['slope_GPa'],'intercept_GPa':fit['intercept_GPa'],'fit_R2':fit['r_squared'],'fit_RMSE_GPa':fit['rmse_GPa'],
          'proof_stress_GPa':None if p is None else p['stress_GPa'],'proof_strain':None if p is None else p['strain'],
          'peak_common_GPa':r['peak_common_window']['stress_GPa'],'peak_common_strain':r['peak_common_window']['strain'],
          'peak_full_GPa':r['peak_full_record']['stress_GPa'],'peak_full_strain':r['peak_full_record']['strain'],
          'maximum_recorded_strain':r['input']['strain_max'],
          'slope_change_percent':100*(fit['slope_GPa']/reference['elastic_fit']['slope_GPa']-1),
          'proof_change_percent':None if p is None or rp is None else 100*(p['stress_GPa']/rp['stress_GPa']-1),
          'peak_change_percent':100*(r['peak_common_window']['stress_GPa']/reference['peak_common_window']['stress_GPa']-1)})
    spans=[]
    for c in campaign['cases']:
        sub=[s for s in sensitivity if s['case_id']==c['id']]
        valid_proofs = [s['proof_GPa'] for s in sub if s['proof_GPa'] is not None]
        spans.append({'case_id':c['id'], 'slope_min_GPa':min(s['slope_GPa'] for s in sub),'slope_max_GPa':max(s['slope_GPa'] for s in sub),
                      'proof_min_GPa':min(valid_proofs) if valid_proofs else None,
                      'proof_max_GPa':max(valid_proofs) if valid_proofs else None,
                      'averaged_peak_min_GPa':min(s['averaged_peak_GPa'] for s in sub),
                      'averaged_peak_max_GPa':max(s['averaged_peak_GPa'] for s in sub),
                      'configurations':len(sub),'missing_proof_results':sum(s['proof_GPa'] is None for s in sub)})
    lookup={(r['fit_max'],r['smooth_points'],r['confirm_points']):r for r in sensitivity if r['case_id']=='eq'}
    paired=[]
    for c in campaign['cases']:
        sub=[r for r in sensitivity if r['case_id']==c['id']]
        for metric in ['slope_GPa','proof_GPa','averaged_peak_GPa']:
            differences=[]
            for r in sub:
                ref=lookup[r['fit_max'],r['smooth_points'],r['confirm_points']][metric]
                if r[metric] is not None and ref is not None:
                    differences.append(100*(r[metric]/ref-1))
            paired.append({'case_id':c['id'],'metric':metric,
                'minimum_change_percent':min(differences) if differences else None,
                'maximum_change_percent':max(differences) if differences else None,
                'positive_configurations':sum(x>0 for x in differences),
                'negative_configurations':sum(x<0 for x in differences),'valid_configurations':len(differences)})
    return {'schema_version':1,'analysis_settings':settings,'reference_case':'eq','cases':case_results,
            'summary':rows,'sensitivity':sensitivity,'method_sensitivity_ranges':spans,'paired_sensitivity':paired,
            'interpretation':'Method-sensitivity ranges are not confidence intervals; each composition has one loading record.'}


def export(root:Path,output:Path,*,overwrite:bool=False)->dict:
    root,output=Path(root).resolve(),Path(output).resolve()
    for protected in ['data','simulation','src','tests','configs']:
        if output.is_relative_to(root/protected):raise ValueError("Output may not overwrite project inputs or source.")
    if output==root:raise ValueError("Use a results directory, not the project root.")
    if output.exists() and any(output.iterdir()) and not overwrite:raise FileExistsError("Output exists. Select another directory or use --overwrite.")
    value=calculate(root)
    output.mkdir(parents=True,exist_ok=True)
    write_json(output/'analysis.json',value)
    for name,key in [('summary','summary'),('sensitivity','sensitivity'),('method-sensitivity','method_sensitivity_ranges'),('paired-sensitivity','paired_sensitivity')]:
        write_csv(output/f'{name}.csv',value[key])
    return value


def stage(root:Path,case_id:str,destination:Path,potential_dir:Path)->dict:
    root,destination,potential_dir=Path(root).resolve(),Path(destination).resolve(),Path(potential_dir).resolve()
    campaign=load_campaign(root)
    cases=[c for c in campaign['cases'] if c['id']==case_id]
    if not cases:raise ValueError("Unknown case ID.")
    if destination.exists():raise FileExistsError("Choose an unused run directory.")
    for folder in ['data','simulation','src','tests','configs','results','manuscript']:
        if destination.is_relative_to(root/folder):raise ValueError("Runs must not overwrite project files.")
    files=[potential_dir/n for n in ['library.meam','HfNbTaTiZr.meam']]
    if any(not p.is_file() or p.stat().st_size==0 for p in files):raise FileNotFoundError("Supply library.meam and HfNbTaTiZr.meam.")
    if not ELEMENTS.issubset(set(re.findall(r"['\"]([A-Za-z]+)['\"]",'\n'.join(line.split('#', 1)[0] for line in files[0].read_text(encoding='utf-8').splitlines())))):
        raise ValueError("The MEAM library must contain Hf, Nb, Ta, Ti and Zr.")
    params='\n'.join(line.split('#',1)[0] for line in files[1].read_text(encoding='utf-8').splitlines())
    nn={(min(int(i),int(j)),max(int(i),int(j))):float(v) for i,j,v in re.findall(r'nn2\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([\d.eE+\-]+)',params)}
    if not all(nn.get((i,j))==1. for i in range(1,6) for j in range(i,6)):
        raise ValueError("This 2NN protocol requires explicit nn2=1 for all 15 element pairs.")
    c=cases[0];destination.mkdir(parents=True)
    for d in ['potential','outputs']:(destination/d).mkdir()
    shutil.copyfile(root/c['input_file'],destination/'in.compression.lmp')
    for f in files:shutil.copyfile(f,destination/'potential'/f.name)
    record={'case_id':case_id,'input_sha256':c['input_sha256'],'potential_sha256':{p.name:sha256_file(p) for p in files},
            'lammps_executed':False,'command':['lmp','-in','in.compression.lmp','-log','log.lammps'],
            'element_order':['Hf','Nb','Ta','Ti','Zr']}
    write_json(destination/'run-inputs.json',record)
    return record
