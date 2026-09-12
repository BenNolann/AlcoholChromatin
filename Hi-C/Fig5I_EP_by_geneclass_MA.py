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

"""Fig 5I - E-P contact change by gene trajectory class (MA-normalised) + guanine-gradient retention.
Inputs: loops/ep_unbiased_MA.tsv, loops/genes_classed.tsv"""
D=pd.read_csv(f"{L}/ep_unbiased_MA.tsv",sep="\t")
D=D[D["genes"].notna()&(D["genes"]!=".")].copy()
D["g"]=D["genes"].str.split(","); X=D.explode("g"); X["gene"]=X["g"].str.strip()
G=pd.read_csv(f"{L}/genes_classed.tsv",sep="\t")[["gene","class"]]
X=X.merge(G,on="gene",how="inner")
X["grp"]=np.where(X["class"]=="Lasting","Lasting",
         np.where(X["class"]=="Acute","Acute",
         np.where(X["class"]=="notDE","Unchanged",None)))
S=X[X["grp"].notna()].copy()
ORD=["Lasting","Acute","Unchanged"]; COL={"Lasting":RED,"Acute":ORA,"Unchanged":GRY}
fig,ax=plt.subplots(1,2,figsize=(11.6,4.8)); w=.34
for i,(col,hatch) in enumerate([("L2E_ma",""),("L2W_ma","//")]):
    mu=[S.loc[S["grp"]==g,col].mean() for g in ORD]; er=[1.96*S.loc[S["grp"]==g,col].sem() for g in ORD]
    ax[0].bar(np.arange(3)+(i-.5)*w,mu,w,yerr=er,capsize=4,hatch=hatch,
              color=[COL[g] for g in ORD],edgecolor="k",lw=.7,error_kw=dict(lw=1.1))
ax[0].axhline(0,color="k",lw=.9); ax[0].set_xticks(range(3))
ax[0].set_xticklabels([f"{g}\nn={int((S['grp']==g).sum()):,} pairs" for g in ORD],fontsize=9)
ax[0].set_ylabel("E-P contact log2FC\nMA-normalised, (o+1)/(e+1)",fontsize=9.5)
ax[0].set_title("Contacts at persistent genes stay elevated;\ncontacts elsewhere return to baseline",fontsize=10.5)
h=[plt.Rectangle((0,0),1,1,fc="white",ec="k",hatch=x) for x in ("","//")]
ax[0].legend(h,["EtOH / Ctl","Withdrawal / Ctl"],fontsize=9,frameon=False,loc="upper right")
res={}
for g in ORD:
    s=S[S["grp"]==g]
    q=np.nanquantile(s["gc"],[i/5 for i in range(1,5)]); b=np.digitize(s["gc"],q)
    res[g]=[(s[b==4][c].mean()-s[b==0][c].mean(),
             1.96*np.sqrt(s[b==4][c].sem()**2+s[b==0][c].sem()**2)) for c in ("L2E_ma","L2W_ma")]
for i,hatch in enumerate(["","//"]):
    mu=[res[g][i][0] for g in ORD]; er=[res[g][i][1] for g in ORD]
    ax[1].bar(np.arange(3)+(i-.5)*w,mu,w,yerr=er,capsize=4,hatch=hatch,
              color=[COL[g] for g in ORD],edgecolor="k",lw=.7,error_kw=dict(lw=1.1))
for j,g in enumerate(ORD):
    top=res[g][1][0]+res[g][1][1]
    ax[1].annotate(f"{100*res[g][1][0]/res[g][0][0]:.0f}% retained",xy=(j+.5*w,top),xytext=(0,6),
                   textcoords="offset points",ha="center",fontsize=8.5,fontweight="bold")
ax[1].axhline(0,color="k",lw=.9); ax[1].set_xticks(range(3)); ax[1].set_xticklabels(ORD,fontsize=9.5)
ax[1].set_ylabel("guanine gradient of E-P contacts\n(Q5 - Q1 of GC)",fontsize=9.5)
ax[1].set_title("The guanine gradient is induced equally\nbut persists only at lasting genes",fontsize=10.5)
ax[1].set_ylim(0,0.098)
plt.tight_layout(); save(fig,"Fig5I_EP_by_geneclass_MA")
