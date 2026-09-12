import numpy as np, pandas as pd, subprocess, collections, pickle
OUT="bensalcohol_out/loops"; RES=10000; W=10; OFFMAX=200
HIC={"ctrl":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/ctrl_merged_nodups_q30_592m.hic",
     "etoh":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/etoh_merged_nodups_q30_592m.hic",
     "with":"/Zulu/bnolan/Projects/Personal/Ethanol/GEO/HiC/T2_mega_q30.hic"}
# ---- loop sets ----
H=None;R=[]
for l in open(f"{OUT}/loop_table.tsv"):
    p=l.rstrip("\n").split("\t")
    if H is None: H={k:i for i,k in enumerate(p)}; continue
    R.append(p)
C=[l.rstrip("\n").split("\t") for l in open(f"{OUT}/conv_anchor_zf5.tsv")][1:]
cza=np.array([c[0] for c in C]); czb=np.array([c[1] for c in C])
SETS={"all CTCF loops":np.ones(len(R),bool),
      "convergent GGG/GGG":(cza=="GGG")&(czb=="GGG"),
      "convergent GAT/GAT":(cza=="GAT")&(czb=="GAT"),
      "convergent any-CTCF":(cza!="none")&(czb!="none")}
for k,v in SETS.items(): print(f"{k}: n={int(v.sum())}")
rows=[]
for i,r in enumerate(R):
    rows.append((r[H['chr1']].replace("chr",""),(int(r[H['x1']])+int(r[H['x2']]))//2,
                 (int(r[H['y1']])+int(r[H['y2']]))//2))
rows=pd.DataFrame(rows,columns=["chrn","m1","m2"])
for k,v in SETS.items(): rows[k]=v
def band(cond,ch):
    cmd=(f"straw NONE {HIC[cond]} {ch} {ch} BP {RES} 2>/dev/null | "
         f"awk '{{d=$2-$1; if(d<0)d=-d; if(d<={(OFFMAX+1)*RES}) print $1,$2,$3}}'")
    p=subprocess.run(["bash","-c",cmd],capture_output=True,text=True)
    if not p.stdout.strip(): return None,0
    df=pd.read_csv(pd.io.common.StringIO(p.stdout),sep=" ",header=None,names=["i","j","v"])
    bi=df["i"].to_numpy()//RES; bj=df["j"].to_numpy()//RES
    lo=np.minimum(bi,bj); hi=np.maximum(bi,bj); off=hi-lo; nb=int(hi.max())+1
    M=np.zeros((nb,OFFMAX+1),dtype=np.float64); ok=off<=OFFMAX
    np.add.at(M,(lo[ok],off[ok]),df["v"].to_numpy()[ok]); return M,nb
APA={c:{k:[np.zeros((2*W+1,2*W+1)),0] for k in SETS} for c in HIC}
for ch in [str(c) for c in range(1,23)]+["X"]:
    sub=rows[rows["chrn"]==ch]
    if not len(sub): continue
    for cond in HIC:
        M,nb=band(cond,ch)
        if M is None: continue
        exp=np.array([M[:nb-o,o].sum()/max(nb-o,1) for o in range(OFFMAX+1)]); exp[exp<=0]=np.nan
        for _,rr in sub.iterrows():
            i,j=rr["m1"]//RES,rr["m2"]//RES
            if i>j: i,j=j,i
            if j-i>OFFMAX-W or i<W or j+W>=nb: continue
            win=np.full((2*W+1,2*W+1),np.nan)
            for a in range(-W,W+1):
                ii=i+a
                for b in range(-W,W+1):
                    jj=j+b
                    if ii<0 or jj>=nb or jj<ii: continue
                    o=jj-ii
                    if o<=OFFMAX and np.isfinite(exp[o]): win[a+W,b+W]=M[ii,o]/exp[o]
            g=np.nan_to_num(win)
            for k in SETS:
                if rr[k]: APA[cond][k][0]+=g; APA[cond][k][1]+=1
    print(f"  chr{ch} done",flush=True)
pickle.dump({c:{k:(v[0],v[1]) for k,v in APA[c].items()} for c in APA},open(f"{OUT}/apa_loops.pkl","wb"))
print("\nAPA scores (centre 3x3 / corner):")
for k in SETS:
    line=f"  {k:22s}"
    for c in ("ctrl","etoh","with"):
        m=APA[c][k][0]/max(APA[c][k][1],1)
        s=m[W-1:W+2,W-1:W+2].mean()/max(m[:5,-5:].mean(),1e-9)
        line+=f"{c}={s:.3f}  "
    line+=f"n={APA['ctrl'][k][1]}"
    print(line)
