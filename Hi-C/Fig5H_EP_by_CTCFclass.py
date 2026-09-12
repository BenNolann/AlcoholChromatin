import os,sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
# Ensure PDF renders text as TrueType (actual text boxes)
matplotlib.rcParams['pdf.fonttype'] = 42

PROJ="/Zulu/jordan/alcoholATAC"; L=f"{PROJ}/bensalcohol_out/loops"
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"svg")
os.makedirs(OUT,exist_ok=True)
BLU="#5a72ad"; RED="#ab5054"; PRP="#7d2c88"; GRY="#999999"; ORA="#E1812C"
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.svg",bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    fig.savefig(f"{OUT}/{name}.pdf")

"""Fig 5H - E-P contact change by CTCF/ZF5 class at the anchors.
Writes BOTH the MA-normalised and the un-normalised version; the MA one matches 5G/5I.
Input: loops/ep_unbiased_MA.tsv"""
D=pd.read_csv(f"{L}/ep_unbiased_MA.tsv",sep="\t")
ORD=["GGG","GAT","otherCTCF","noCTCF"]; LAB=["GGG\nanchor","GAT\nanchor","other CTCF\nanchor","no CTCF"]
def build(ce,cw,name,note):
    fig,ax=plt.subplots(figsize=(6.6,4.7)); w=.36
    for i,(col,c,lab) in enumerate([(ce,RED,"EtOH / Ctl"),(cw,PRP,"Withdraw / Ctl")]):
        mu=[D.loc[D["ctcfcls"]==k,col].mean() for k in ORD]
        er=[1.96*D.loc[D["ctcfcls"]==k,col].sem() for k in ORD]
        ax.bar(np.arange(4)+(i-.5)*w,mu,w,yerr=er,capsize=4,color=c,label=lab,
               edgecolor="k",lw=.6,error_kw=dict(lw=1.1))
    ax.axhline(0,color="k",lw=.9)
    yl=ax.get_ylim(); ax.set_ylim(yl[0],yl[1]*1.30)
    for j,k in enumerate(ORD):
        ax.annotate(f"n={int((D['ctcfcls']==k).sum()):,}",xy=(j,ax.get_ylim()[1]),xytext=(0,-9),
                    textcoords="offset points",ha="center",va="top",fontsize=7.5,color="#444444")
    ax.set_xticks(range(4)); ax.set_xticklabels(LAB,fontsize=9)
    ax.set_ylabel(f"E-P contact log2FC\n{note}",fontsize=9)
    ax.set_title("E-P gain is largest at GGG CTCF anchors,\nbut vanishes at matched GC (z=+0.17, p=0.87)",fontsize=9.5)
    ax.legend(fontsize=9,frameon=False,loc="lower left")
    save(fig,name)
    for k in ORD:
        s=D[D["ctcfcls"]==k]
        print(f"     {k:10s} n={len(s):7,}  EtOH={s[ce].mean():+.4f}+-{1.96*s[ce].sem():.4f}"
              f"   With={s[cw].mean():+.4f}+-{1.96*s[cw].sem():.4f}")
print("  MA-normalised (matches 5G/5I):")
build("L2E_ma","L2W_ma","Fig5H_EP_by_CTCFclass_MA","MA-normalised, (o+1)/(e+1)")
print("  un-normalised:")
build("L2E","L2W","Fig5H_EP_by_CTCFclass","(o+1)/(e+1) distance-normalised")
