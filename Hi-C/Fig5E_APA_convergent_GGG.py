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
        print(vmin)
        print(vmax)
        th,r,C,_=bullseye(m); ax.grid(False)
        ax.pcolormesh(th,r,C,cmap=BR,vmin=vmin,vmax=vmax,rasterized=True)
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

"""Fig 5E - SIPMeta bullseye APA of convergent GGG/GGG loops + CTCF/RAD21 profiles
at the convergent GGG anchors.  Same construction and colour scale rules as Fig 5D.
Inputs: loops/apa_loops.pkl, loops/ctcf_anchor_annot.bed, loops/loop_table.tsv,
        loops/conv_anchor_zf5.tsv, chipseq/{CTCF,RAD21}_{UT,E6,E6R6}*.bw"""
import pyBigWig
APA=pickle.load(open(f"{L}/apa_loops.pkl","rb"))
KEY="convergent GGG/GGG"
S=[]
for l in open(f"{L}/ctcf_anchor_annot.bed"):
    p=l.rstrip("\n").split("\t"); S.append((p[0],int(p[1]),float(p[4]),p[5],p[6]))
BIN=25000; idx=collections.defaultdict(list)
for i,s in enumerate(S): idx[(s[0],s[1]//BIN)].append(i)
H=None;R=[]
for l in open(f"{L}/loop_table.tsv"):
    p=l.rstrip("\n").split("\t")
    if H is None: H={k:i for i,k in enumerate(p)}; continue
    R.append(p)
C=[l.rstrip("\n").split("\t") for l in open(f"{L}/conv_anchor_zf5.tsv")][1:]
def bstrand(ch,a,b,strand):
    bi=None
    for k in range(a//BIN,b//BIN+1):
        for i in idx.get((ch,k),[]):
            if a<=S[i][1]<b and S[i][3]==strand and (bi is None or S[i][2]>S[bi][2]): bi=i
    return bi
anch=set()
for r,cc in zip(R,C):
    if cc[0]=="GGG" and cc[1]=="GGG":
        i=bstrand(r[H['chr1']],int(r[H['x1']]),int(r[H['x2']]),'-')
        j=bstrand(r[H['chr2']],int(r[H['y1']]),int(r[H['y2']]),'+')
        if i is not None: anch.add(i)
        if j is not None: anch.add(j)
anch=sorted(anch); print(f"     convergent GGG/GGG anchor sites: {len(anch)}")
BWF={"CTCF":[("Control","chipseq/CTCF_UT_10bp_BPM.bw"),("EtOH","chipseq/CTCF_E6_10bp_MABPM.bw"),
             ("Withdrawal","chipseq/CTCF_E6R6_10bp_MABPM.bw")],
     "RAD21":[("Control","chipseq/RAD21_UT_10bp_BPM.bw"),("EtOH","chipseq/RAD21_E6_10bp_MABPM.bw"),
              ("Withdrawal","chipseq/RAD21_E6R6_10bp_MABPM.bw")]}
PROF={}
for prot,items in BWF.items():
    for lab,f in items:
        bw=pyBigWig.open(f); ok=set(bw.chroms()); acc=np.zeros(NB); n=0
        for i in anch:
            ch,pos=S[i][0],S[i][1]
            key=ch if ch in ok else ch.replace("chr","")
            if key not in ok: continue
            a,b=pos-FL,pos+FL
            if a<0 or b>=bw.chroms()[key]: continue
            try: raw=np.array(bw.values(key,a,b),dtype=np.float32)
            except RuntimeError: continue
            if raw.size!=2*FL: continue
            acc+=np.nan_to_num(raw).reshape(NB,(2*FL)//NB).mean(axis=1); n+=1
        PROF[(prot,lab)]=acc/max(n,1)
fig=plt.figure(figsize=(10.5,7.8))
gs=GridSpec(2,6,figure=fig,height_ratios=[1.30,1.0],hspace=.36,wspace=.85)
apa_row(fig,gs,KEY,APA,"convergent GGG/GGG")
prof_row(fig,gs,PROF)
fig.suptitle("Convergent GGG/GGG loops behave like the bulk",fontsize=12,y=0.98)
save(fig,"Fig5E_APA_convergent_GGG")
for p in ("CTCF","RAD21"):
    c=NB//2; k={l:PROF[(p,l)][c-4:c+4].mean() for l in ("Control","EtOH","Withdrawal")}
    print(f"     {p:6s} ctrl={k['Control']:.3f}  etoh={k['EtOH']:.3f} ({100*(k['EtOH']/k['Control']-1):+.1f}%)"
          f"  with={k['Withdrawal']:.3f} ({100*(k['Withdrawal']/k['Control']-1):+.1f}%)")
