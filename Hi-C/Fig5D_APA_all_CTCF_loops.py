import os,sys,pickle,collections
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
# Ensure PDF renders text as TrueType (actual text boxes)
matplotlib.rcParams['pdf.fonttype'] = 42


PROJ="/Zulu/jordan/alcoholATAC"; L=f"{PROJ}/bensalcohol_out/loops"
os.chdir(PROJ)
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
OUT=os.path.join(HERE,"svg"); os.makedirs(OUT,exist_ok=True)
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from bullseye import bullseye
BLU="#5a72ad"; RED="#ab5054"; PRP="#7d2c88"; GRY="#999999"; ORA="#E1812C"
W=10; FL=2000; NB=80
COND=[("ctrl","Control"),("etoh","EtOH"),("with","Withdrawal")]
COL={"Control":BLU,"EtOH":RED,"Withdrawal":PRP}
STYLE={"Control":dict(lw=1.9,ls="-"),"EtOH":dict(lw=1.9,ls="-"),"Withdrawal":dict(lw=1.9,ls=(0,(4,2.2)))}
BR=LinearSegmentedColormap.from_list("bright_red",[(1,1,1),(1,0,0)])
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.svg",bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    fig.savefig(f"{OUT}/{name}.pdf")

    print(f"  -> svg/{name}.svg")
def apa_row(fig,gs,key,APA,label):
    MATS={c:APA[c][key][0]/max(APA[c][key][1],1) for c,_ in COND}
    vmax=max(np.nanmax(m) for m in MATS.values()); vmin=min(np.nanmin(m) for m in MATS.values())
    for i,(c,lab) in enumerate(COND):
        ax=fig.add_subplot(gs[0,2*i:2*i+2],projection="polar")
        m=MATS[c]; sc=m[W-1:W+2,W-1:W+2].mean()/max(m[:5,-5:].mean(),1e-9)
        th,r,C,_=bullseye(m); ax.grid(False)
        ax.pcolormesh(th,r,C,cmap=BR,vmin=vmin,vmax=vmax,rasterized=True)
        print(vmin)
        print(vmax)
        ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_theta_offset(np.pi/2)
        ax.spines["polar"].set_visible(False)
        ax.set_title(f"{lab}\nAPA={sc:.2f}",fontsize=10,pad=12)
        if i==0: ax.text(-0.18,0.5,f"{label}\nn={APA['ctrl'][key][1]:,}",transform=ax.transAxes,
                         rotation=90,va="center",ha="center",fontsize=9)
def prof_row(fig,gs,PROF):
    x=(np.arange(NB)+.5)*(2*FL/NB)-FL
    for j,prot in enumerate(("CTCF","RAD21")):
        ax=fig.add_subplot(gs[1,3*j:3*j+3])
        for lab in ("Control","EtOH","Withdrawal"):
            ax.plot(x/1000,PROF[(prot,lab)],color=COL[lab],label=lab,**STYLE[lab])
        ax.axvline(0,color="grey",lw=.7,ls=":")
        ax.set_title(prot,fontsize=10.5); ax.set_xlabel("distance from CTCF motif (kb)",fontsize=8.5)
        if j==0: ax.set_ylabel("mean BPM",fontsize=9); ax.legend(fontsize=8,frameon=False)
        ax.set_xlim(-2,2); ax.tick_params(labelsize=8)

"""Fig 5D - SIPMeta bullseye APA of all CTCF loops + CTCF/RAD21 profiles at their anchors.
Inputs: loops/apa_loops.pkl (upstream/panel_apa.py), loops/panelb_profiles.pkl (upstream/panel_b_profiles.py)"""
APA=pickle.load(open(f"{L}/apa_loops.pkl","rb"))
PROF=pickle.load(open(f"{L}/panelb_profiles.pkl","rb"))["prof"]
fig=plt.figure(figsize=(10.5,7.8))
gs=GridSpec(2,6,figure=fig,height_ratios=[1.30,1.0],hspace=.36,wspace=.85)
apa_row(fig,gs,"all CTCF loops",APA,"all CTCF loops")
prof_row(fig,gs,PROF)
fig.suptitle("CTCF loops are unchanged; occupancy at their anchors is not",fontsize=12,y=0.98)

save(fig,"Fig5D_APA_all_CTCF_loops")
for p in ("CTCF","RAD21"):
    c=NB//2; k={l:PROF[(p,l)][c-4:c+4].mean() for l in ("Control","EtOH","Withdrawal")}
    print(f"     {p:6s} ctrl={k['Control']:.3f}  etoh={k['EtOH']:.3f} ({100*(k['EtOH']/k['Control']-1):+.1f}%)"
          f"  with={k['Withdrawal']:.3f} ({100*(k['Withdrawal']/k['Control']-1):+.1f}%)")

