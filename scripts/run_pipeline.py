"""Launch database-dependent research stages in a caller-owned work directory."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

REPO=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['hetionet-prepare','hetionet-kge','hetionet-rgcn','hetionet-summary','primekg-train','primekg-summary','mimic-extract','mimic-severity','mimic-assemble'])
    p.add_argument('--workdir',type=Path,required=True)
    p.add_argument('--seeds',type=int,nargs='+',default=[1,2,3])
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--max-epochs',type=int)
    a=p.parse_args();a.workdir=a.workdir.resolve()
    if a.threads<1:p.error('threads must be positive')
    if "'" in str(a.workdir):p.error('Work directory cannot contain a single quote in the SQL extraction stages')
    env=dict(os.environ,KG_AUDIT_WORKDIR=str(a.workdir),PYTHONHASHSEED='0',OPENBLAS_NUM_THREADS='1',
             OMP_NUM_THREADS=str(a.threads),THREADS=str(a.threads))
    env['PYTHONPATH']=str(REPO/'src')+os.pathsep+env.get('PYTHONPATH','')
    names={'hetionet-prepare':'hetionet_prepare.py','hetionet-kge':'hetionet_train_kge.py',
           'hetionet-rgcn':'hetionet_train_rgcn.py','hetionet-summary':'hetionet_summarize.py',
           'primekg-train':'primekg_train.py','primekg-summary':'primekg_summarize.py',
           'mimic-extract':'mimic_extract.py','mimic-severity':'mimic_severity.py','mimic-assemble':'mimic_assemble.py'}
    if a.stage in {'hetionet-kge','hetionet-rgcn'}:
        for seed in a.seeds:
            current=dict(env,SEED=str(seed),SUFFIX=f'_seed{seed}',MODELS='complex,distmult,rotate',
                         SPLITS='random,coldstart,scaffold',EPOCHS=str(a.max_epochs or 20))
            subprocess.run([sys.executable,str(REPO/'pipelines'/names[a.stage])],env=current,check=True)
    else:
        env['SEEDS']=' '.join(map(str,a.seeds))
        if a.max_epochs is not None:env['MAX_EPOCHS']=str(a.max_epochs)
        subprocess.run([sys.executable,str(REPO/'pipelines'/names[a.stage])],env=env,check=True)


if __name__=='__main__':main()
