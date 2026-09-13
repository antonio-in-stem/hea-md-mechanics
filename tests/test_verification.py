"""Export-corruption tests exercise the verification boundary."""
from pathlib import Path
import csv
import json
import shutil
import pytest

from hea_md.verification import verify, compare_results, EXPORTS


@pytest.mark.parametrize('filename',list(EXPORTS))
@pytest.mark.parametrize('fault',['number','nan','label','missing_row','extra_row','header','extra_column','empty_number'])
def test_corrupted_export_fails(root,tmp_path,filename,fault):
    output=tmp_path/'results'
    shutil.copytree(root/'results',output)
    path=output/filename
    with path.open(newline='',encoding='utf-8') as stream:
        rows=list(csv.reader(stream))
    numeric=next(i for i,name in enumerate(rows[0]) if name.endswith('GPa') or name=='minimum_change_percent')
    if fault=='number':rows[1][numeric]='999'
    elif fault=='nan':rows[1][numeric]='NaN'
    elif fault=='label':rows[1][0]='not_a_case'
    elif fault=='missing_row':rows.pop()
    elif fault=='extra_row':rows.append(rows[-1])
    elif fault=='header':rows[0][numeric]='incorrect_column'
    elif fault=='extra_column':rows[1].append('1')
    elif fault=='empty_number':rows[1][numeric]=''
    with path.open('w',newline='',encoding='utf-8') as stream:csv.writer(stream).writerows(rows)
    with pytest.raises(ValueError):verify(root,output/'analysis.json')


def test_undefined_offset_cannot_be_filled_with_zero(root,tmp_path):
    output=tmp_path/'results';shutil.copytree(root/'results',output)
    path=output/'sensitivity.csv'
    with path.open(newline='',encoding='utf-8') as stream:
        reader=csv.DictReader(stream);fields=reader.fieldnames;rows=list(reader)
    row=next(r for r in rows if r['proof_GPa']=='');row['proof_GPa']='0'
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fields);writer.writeheader();writer.writerows(rows)
    with pytest.raises(ValueError,match='empty field'):verify(root,output/'analysis.json')


@pytest.mark.parametrize('left,right',[(1,True),(True,1),(1.,float('nan')),(float('inf'),float('inf')),
                                     (None,0),({'a':1},{'a':1,'b':2}),([1,2],[2,1])])
def test_strict_json_types(left,right):
    with pytest.raises(ValueError):compare_results(left,right)


def test_roundoff_tolerance_is_allowed():
    compare_results({'v':136.59949909291953},{'v':136.59949909291956})
