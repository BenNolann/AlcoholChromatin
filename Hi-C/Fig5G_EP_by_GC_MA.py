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


"""Fig 5G - all enhancer-promoter contacts by GC, MA-normalised.
Input: loops/ep_unbiased_MA.tsv (upstream/ep_unbiased.py -> upstream/ma_normalise.py)"""
D=pd.read_csv(f"{L}/ep_unbiased_MA.tsv",sep="\t")
x=np.arange(1,6)
q=np.nanquantile(D["gc"],[i/5 for i in range(1,5)]); b=np.digitize(D["gc"],q)
fig,ax=plt.subplots(figsize=(5.8,4.7))
for col,c,lab in [("L2E_ma",RED,"EtOH / Ctl"),("L2W_ma",PRP,"Withdraw / Ctl")]:
    m=[D.loc[b==k,col].mean() for k in range(5)]
    e=[1.96*D.loc[b==k,col].sem() for k in range(5)]
    ax.errorbar(x,m,yerr=e,fmt="-o",color=c,capsize=3,lw=2,label=lab)
ax.axhline(0,color="k",lw=.9); ax.set_xticks(x)
ax.set_xlabel("GC quintile of the E-P pair",fontsize=9)
ax.set_ylabel("E-P contact log2FC\nMA-normalised, (o+1)/(e+1)",fontsize=9)
se=D[b==4]["L2E_ma"].mean()-D[b==0]["L2E_ma"].mean(); sw=D[b==4]["L2W_ma"].mean()-D[b==0]["L2W_ma"].mean()
ax.set_title(f"All enhancer-promoter contacts\nn={len(D):,}   Q5-Q1: {se:+.3f} EtOH, {sw:+.3f} withdrawal ({100*sw/se:.0f}% retained)",fontsize=9.5)
ax.legend(fontsize=9,frameon=False)
save(fig,"Fig5G_EP_by_GC_MA")
