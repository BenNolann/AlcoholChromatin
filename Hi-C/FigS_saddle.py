import os,sys,pickle,collections
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"]="none"
PROJ="/Zulu/jordan/alcoholATAC"; L=f"{PROJ}/bensalcohol_out/loops"
os.chdir(PROJ)
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"svg"); os.makedirs(OUT,exist_ok=True)
sys.path.insert(0,L)
RED="#C44E52"; BLU="#4C72B0"; GRN="#238b45"; GRY="#999999"; ORA="#E1812C"
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.svg",bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    print(f"  -> svg/{name}.svg")

"""Supplementary - saddle plots on the CRUSH axis and the raw GC axis (50-200 kb band).
Input: loops/saddle_newcrush.pkl  (upstream/saddle2.py, built on the crush-rs 5 kb tracks)"""
S=pickle.load(open(f"{L}/saddle_newcrush.pkl","rb"))
def quad(m,k=4):
    return dict(BB=np.nanmean(m[:k,:k]),AA=np.nanmean(m[-k:,-k:]),
                AB=(np.nanmean(m[:k,-k:])+np.nanmean(m[-k:,:k]))/2)
for axis,tag in (("crush","CRUSH"),("gc","GC")):
    fig,ax=plt.subplots(1,4,figsize=(14.5,3.9))
    mats={c:S[c][(axis,0)] for c in ("ctrl","etoh","with")}
    vmin=np.nanpercentile(np.log2([m for m in mats.values()]),1)
    vmax=np.nanpercentile(np.log2([m for m in mats.values()]),99)
    for i,(c,lab) in enumerate([("ctrl","Control"),("etoh","EtOH"),("with","Withdrawal")]):
        m=mats[c]; q=quad(m); st=(q["AA"]+q["BB"])/(2*q["AB"])
        ax[i].imshow(np.log2(m),cmap="RdBu_r",vmin=vmin,vmax=vmax,origin="lower")
        ax[i].set_title(f"{lab}\nstrength={st:.3f}",fontsize=10)
        ax[i].set_xticks([]); ax[i].set_yticks([])
        if i==0: ax[i].set_ylabel(f"B $\\rightarrow$ A quantile ({tag})",fontsize=9)
    d=np.log2(mats["etoh"]/mats["ctrl"]); mm=np.nanpercentile(np.abs(d),98)
    im=ax[3].imshow(d,cmap="bwr",vmin=-mm,vmax=mm,origin="lower")
    ax[3].set_title("log2 EtOH / Ctl",fontsize=10); ax[3].set_xticks([]); ax[3].set_yticks([])
    plt.colorbar(im,ax=ax[3],fraction=.046)
    qc,qe=quad(mats["ctrl"]),quad(mats["etoh"])
    fig.suptitle(f"Saddle, {tag} axis, 50-200 kb   |   AA {100*(qe['AA']/qc['AA']-1):+.1f}%   "
                 f"BB {100*(qe['BB']/qc['BB']-1):+.1f}%   AB {100*(qe['AB']/qc['AB']-1):+.1f}%",fontsize=11,y=1.02)
    plt.tight_layout(); save(fig,f"FigS_saddle_{axis}")
