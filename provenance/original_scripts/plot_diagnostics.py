"""Export chemical-character and available occupation-comparison figures."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'validation/figures';out.mkdir(exist_ok=True)
    rows=list(csv.DictReader((ROOT/'orbitals/orbital_character_shells.csv').open()))[78:98]
    groups=[(0,'Fe','3d','Fe atom 0: 3d'),(5,'Fe','3d','Fe atom 5: 3d'),
            (1,'S','3p','Bridge S atom 1: 3p'),(4,'S','3p','Bridge S atom 4: 3p'),
            *[(a,'S','3p',f'Terminal S atom {a}: 3p') for a in [2,3,6,7]]]
    weights=np.array([[float(row[f'atom{a}_{el}_{shell}_meta_lowdin_weight']) for row in rows] for a,el,shell,_ in groups])
    weights=np.vstack([weights,1-weights.sum(axis=0)])
    fig,ax=plt.subplots(figsize=(11,4.7),constrained_layout=True)
    im=ax.imshow(weights,vmin=0,vmax=1,aspect='auto',cmap='viridis')
    ax.set_yticks(range(9),[x[3] for x in groups]+['Other AO character'])
    ax.set_xticks(range(20));ax.set_xlabel('Initial localized CAS20 orbital (zero-based)')
    ax.set_title('Chemical character before correlated natural-orbital rotation')
    for x in [1.5,6.5,12.5,17.5]:ax.axvline(x,color='white',lw=1)
    fig.colorbar(im,ax=ax,label='Meta-Löwdin weight')
    for suffix in ['png','pdf']:fig.savefig(out/f'initial_CAS20_character.{suffix}',dpi=180)
    plt.close(fig)
    for stage in ['provisional','production']:
        path=ROOT/'validation'/f'occupation_comparison_{stage}.json'
        if not path.exists():continue
        data=json.loads(path.read_text());ref=np.array(data['reference']['occupations']);ours=np.array(data['reconstructed']['occupations'])
        fig,(ax,bx)=plt.subplots(2,1,figsize=(8,6),sharex=True,constrained_layout=True,gridspec_kw={'height_ratios':[3,1]})
        x=np.arange(1,21)
        ax.plot(x,ref,'o-',label='Historical FCIDUMP');ax.plot(x,ours,'s--',label=f'Reconstructed ({stage})')
        ax.axhline(1.97,color='grey',ls=':',lw=1,label='1.97 screening value')
        ax.set_ylabel('Spin-summed natural occupation');ax.legend(fontsize=8);ax.set_ylim(0,2.07)
        ax.set_title('Singlet CAS(30e,20o) DMRG density spectra')
        bx.axhline(0,color='grey',lw=1);bx.bar(x,ours-ref)
        bx.set_ylabel('Difference');bx.set_xlabel('Descending occupation rank');bx.set_xticks(x)
        for suffix in ['png','pdf']:fig.savefig(out/f'occupation_comparison_{stage}.{suffix}',dpi=180)
        plt.close(fig)
    print('Diagnostic figures saved in',out)

if __name__=='__main__':main()
