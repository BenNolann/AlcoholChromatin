"""Fig 5F - GGG CTCF sites are depleted from convergent loop anchors.
Enrichment vs a permutation null matched on motif score x control CTCF occupancy.
Inputs: loops/ctcf_anchor_annot.bed, humanZF/fc_table.tsv, loops/loop_table.tsv,
        loops/conv_anchor_zf5.tsv"""
import os,sys,collections
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"]="none"
PROJ="/Zulu/jordan/alcoholATAC"; L=f"{PROJ}/bensalcohol_out/loops"; os.chdir(PROJ)
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"svg"); os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(21); NPERM=20000
GRN="#238b45"; RED="#C44E52"; GRY="#999999"
ut={}
for l in open("bensalcohol_out/humanZF/fc_table.tsv"):
    if l.startswith("sid"): continue
    p=l.rstrip("\n").split("\t"); ut[p[0]]=float(p[3])
S=[]
for l in open(f"{L}/ctcf_anchor_annot.bed"):
    p=l.rstrip("\n").split("\t")
    if p[3] in ut: S.append((p[0],int(p[1]),p[3],float(p[4]),p[5],p[6],ut[p[3]]))
BIN=25000; idx=collections.defaultdict(list)
for i,s in enumerate(S): idx[(s[0],s[1]//BIN)].append(i)
H=None;R=[]
for l in open(f"{L}/loop_table.tsv"):
    p=l.rstrip("\n").split("\t")
    if H is None: H={k:i for i,k in enumerate(p)}; continue
    R.append(p)
C=[l.rstrip("\n").split("\t") for l in open(f"{L}/conv_anchor_zf5.tsv")][1:]
isC=np.zeros(len(S),bool)
def hits(ch,a,b):
    o=[]
    for k in range(a//BIN,b//BIN+1):
        for i in idx.get((ch,k),[]):
            if a<=S[i][1]<b: o.append(i)
    return o
for r,cc in zip(R,C):
    if cc[0]!="none":
        for i in hits(r[H['chr1']],int(r[H['x1']]),int(r[H['x2']])):
            if S[i][4]=='-': isC[i]=True
    if cc[1]!="none":
        for i in hits(r[H['chr2']],int(r[H['y1']]),int(r[H['y2']])):
            if S[i][4]=='+': isC[i]=True
zf5=np.array([s[5] for s in S]); sc=np.array([s[3] for s in S]); u=np.array([s[6] for s in S])
def qb(v,n): return np.digitize(v,np.quantile(v,[i/n for i in range(1,n)]))
cell=qb(sc,12)*12+qb(u,12)
def enr(feat,target):
    f=np.asarray(feat,bool); obs=int((f&target).sum()); null=np.zeros(NPERM)
    for c in np.unique(cell):
        m=cell==c; K=int((m&target).sum()); n=int((m&f).sum()); N=int(m.sum())
        if K==0 or n==0: continue
        null+=rng.hypergeometric(K,N-K,n,size=NPERM)
    return obs/null.mean(),obs/np.percentile(null,97.5),obs/np.percentile(null,2.5),(obs-null.mean())/null.std(ddof=1)
fig,ax=plt.subplots(figsize=(5.2,4.6)); COL={"GGG":GRN,"GAT":RED,"other":GRY}
for i,k in enumerate(("GGG","GAT","other")):
    r_,lo,hi,zz=enr(zf5==k,isC)
    ax.bar(i,np.log2(r_),color=COL[k],width=.6)
    ax.errorbar(i,np.log2(r_),yerr=[[np.log2(r_)-np.log2(lo)],[np.log2(hi)-np.log2(r_)]],
                color='k',capsize=4,lw=1.2)
    ax.text(i,0.185,f"{r_:.3f}\nz={zz:+.1f}",ha='center',va='top',fontsize=9)
    print(f"  {k:6s} obs/exp={r_:.3f}  z={zz:+.2f}  n={int((zf5==k).sum())}")
ax.axhline(0,color='k',lw=1); ax.set_xticks(range(3)); ax.set_xticklabels(("GGG","GAT","other"))
ax.set_ylim(-0.19,0.21)
ax.set_ylabel("log2 enrichment as convergent loop anchor",fontsize=9)
ax.set_title("GGG CTCF sites are depleted from loop anchors\n(matched on motif score $\\times$ occupancy)",fontsize=10)
fig.savefig(f"{OUT}/Fig5F_anchor_depletion.svg",bbox_inches='tight')
fig.savefig(f"{OUT}/Fig5F_anchor_depletion.png",dpi=200,bbox_inches='tight')
print("  -> svg/Fig5F_anchor_depletion.svg")
