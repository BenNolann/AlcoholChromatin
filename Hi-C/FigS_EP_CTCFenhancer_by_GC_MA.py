import os,sys,pickle,collections
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
# Ensure PDF renders text as TrueType (actual text boxes)
matplotlib.rcParams['pdf.fonttype'] = 42

PROJ="/Zulu/jordan/alcoholATAC"; L=f"{PROJ}/bensalcohol_out/loops"
os.chdir(PROJ)
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"svg"); os.makedirs(OUT,exist_ok=True)
sys.path.insert(0,L)
BLU="#5a72ad"; RED="#ab5054"; PRP="#7d2c88"; GRY="#999999"; ORA="#E1812C"
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.svg",bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    fig.savefig(f"{OUT}/{name}.pdf")

"""Supplementary - E-P contacts with a CTCF site at the ENHANCER, by GC, MA-normalised.
Plotted on the same y-axis as Fig 5G for direct comparison.
Input: loops/ep_unbiased_MA.tsv"""
D=pd.read_csv(f"{L}/ep_unbiased_MA.tsv",sep="\t")
S=D[D["ctcf_at_enh"]].copy(); x=np.arange(1,6)
def rng(T):
    q=np.nanquantile(T["gc"],[i/5 for i in range(1,5)]); b=np.digitize(T["gc"],q)
    lo=min(min(T.loc[b==k,c].mean()-1.96*T.loc[b==k,c].sem() for k in range(5)) for c in ("L2E_ma","L2W_ma"))
    hi=max(max(T.loc[b==k,c].mean()+1.96*T.loc[b==k,c].sem() for k in range(5)) for c in ("L2E_ma","L2W_ma"))
    return lo,hi
l1=rng(D); l2=rng(S); pad=.08*(max(l1[1],l2[1])-min(l1[0],l2[0]))
YL=(min(l1[0],l2[0])-pad,max(l1[1],l2[1])+pad)
q=np.nanquantile(S["gc"],[i/5 for i in range(1,5)]); b=np.digitize(S["gc"],q)
fig,ax=plt.subplots(figsize=(5.8,4.7))
for col,c,lab in [("L2E_ma",RED,"EtOH / Ctl"),("L2W_ma",PRP,"Withdraw / Ctl")]:
    m=[S.loc[b==k,col].mean() for k in range(5)]; e=[1.96*S.loc[b==k,col].sem() for k in range(5)]
    ax.errorbar(x,m,yerr=e,fmt="-o",color=c,capsize=3,lw=2,label=lab)
ax.axhline(0,color="k",lw=.9); ax.set_xticks(x); ax.set_ylim(*YL)
ax.set_xlabel("GC quintile of the E-P pair",fontsize=9)
ax.set_ylabel("E-P contact log2FC\nMA-normalised, (o+1)/(e+1)",fontsize=9)
se=S[b==4]["L2E_ma"].mean()-S[b==0]["L2E_ma"].mean(); sw=S[b==4]["L2W_ma"].mean()-S[b==0]["L2W_ma"].mean()
ax.set_title(f"E-P contacts with CTCF at the enhancer\nn={len(S):,}   Q5-Q1: {se:+.3f} EtOH, {sw:+.3f} withdrawal ({100*sw/se:.0f}% retained)",fontsize=9.5)
ax.legend(fontsize=9,frameon=False)
save(fig,"FigS_EP_CTCFenhancer_by_GC_MA")
