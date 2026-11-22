from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
c=json.loads((R/'results/certificate_d71.json').read_text('utf8'))
m=json.loads((R/'publication/zenodo_metadata_blueprint.json').read_text('utf8'))
assert c['status']=='PASS' and c['dimension']==m['scientific_dimension']==71
assert c['vertices']==m['norm_one_vertices']==76
assert c['pairs']==m['exact_labelled_squared_pairs']==2850
assert c['affine_rank']==c['centered_gram_rank']==71
assert c['eigenvalue_38_multiplicity']==m['model_eigenvalue_38_multiplicity']==65
assert c['eigenvalue_75_over_2_multiplicity']==m['model_eigenvalue_75_over_2_multiplicity']==4
assert c['parent_tag_commit']=='0c56da5cd75613453226fc8ed252dacc01bab22b'
assert c['parent_doi']==m['scientific_parent_doi']=='10.5281/zenodo.24433371'
assert m['related_identifiers']==[{'identifier':c['parent_doi'],'relation':'isDerivedFrom','resource_type':'publication-preprint','scheme':'doi'}]
assert len(m['keywords'])==9 and len(set(m['keywords']))==9
assert 'affine Euclidean dimension seventy-one' in m['keywords'] and 'sixty-fivefold ideal eigenvalue' in m['keywords']
assert not set(m['keywords']).intersection({'affine Euclidean dimension seventy','sixty-fourfold ideal eigenvalue','affine Euclidean dimension seventy-two','sixty-sixfold ideal eigenvalue'})
assert len(m['references'])==2 and c['parent_doi'] in m['references'][0] and '10.5281/zenodo.24407952' in m['references'][1]
assert all(x in (R/'publication/POLYBRUTUS_71D_NOTE.md').read_text('utf8') for x in ('## References',c['parent_doi'],'10.5281/zenodo.24407952'))
print('D71 METADATA 9 KEYWORDS / 2 REFERENCES / DIRECT D70 DOI PASS')
