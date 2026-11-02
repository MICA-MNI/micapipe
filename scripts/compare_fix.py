#!/usr/bin/env python3
"""Compare paired FIX classifications and cleaned NIfTI timeseries."""
import argparse,ast,json
from pathlib import Path
import nibabel as nib
import numpy as np

def labels(path):
    text=Path(path).read_text().strip().splitlines()[-1]
    result=ast.literal_eval(text)
    return {int(x) for x in result}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('legacy',type=Path);p.add_argument('python',type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    old_labels=list(a.legacy.glob('fix4melview*thr20.txt'))
    new_labels=list(a.python.glob('fix4melview*thr20.txt'))
    if len(old_labels)!=1 or len(new_labels)!=1:
        raise SystemExit('Expected one threshold-20 label file per backend')
    old,new=labels(old_labels[0]),labels(new_labels[0])
    n=nib.load(a.legacy/'filtered_func_data.ica/melodic_IC.nii.gz').shape[3]
    if any(x<1 or x>n for x in old|new):
        raise SystemExit('ICA labels outside the component range')
    agreement=1-len(old^new)/n
    jaccard=len(old&new)/len(old|new) if old|new else 1.0
    x=nib.load(a.legacy/'filtered_func_data_clean.nii.gz')
    y=nib.load(a.python/'filtered_func_data_clean.nii.gz')
    mask=nib.load(a.legacy/'mask.nii.gz').get_fdata()>0
    if x.shape!=y.shape or not np.allclose(x.affine,y.affine,atol=1e-5):
        raise SystemExit('Cleaned volume geometry differs')
    # Compare fluctuations after voxelwise demeaning, so static tissue intensity
    # cannot conceal differences in denoised functional signal.
    u=x.get_fdata(dtype=np.float32)[mask].astype(np.float64)
    v=y.get_fdata(dtype=np.float32)[mask].astype(np.float64)
    if not (np.isfinite(u).all() and np.isfinite(v).all()):
        raise SystemExit('Non-finite cleaned data')
    raw_pearson=float(np.corrcoef(u.ravel(),v.ravel())[0,1])
    u-=u.mean(axis=1,keepdims=True);v-=v.mean(axis=1,keepdims=True)
    norm_u=float(np.sum(u*u));norm_v=float(np.sum(v*v))
    if norm_u<=0 or norm_v<=0:
        raise SystemExit('Cleaned data have no temporal variance')
    pearson=float(np.sum(u*v)/np.sqrt(norm_u*norm_v))
    nrmse=float(np.sqrt(np.sum((u-v)**2)/norm_u))
    passed=agreement>=.95 and jaccard>=.90 and pearson>=.999 and nrmse<=.01
    result=dict(components=n,legacy_noise=sorted(old),python_noise=sorted(new),
                disagreement=sorted(old^new),label_agreement=agreement,noise_jaccard=jaccard,
                cleaned_pearson=pearson,cleaned_nrmse=nrmse,raw_pearson=raw_pearson,
                provisional_engineering_gate_passed=passed,
                thresholds=dict(label_agreement=.95,noise_jaccard=.90,cleaned_pearson=.999,cleaned_nrmse=.01))
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
