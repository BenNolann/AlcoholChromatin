import numpy as np, pandas as pd, subprocess, pyBigWig, pickle
OUT="bensalcohol_out/loops"; RES=10000; NQ=20
BANDS=[(50000,200000),(200000,500000),(500000,2000000)]
OFFMAX=200
HIC={"ctrl":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/ctrl_merged_nodups_q30_592m.hic",
     "etoh":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/etoh_merged_nodups_q30_592m.hic",
     "with":"/Zulu/bnolan/Projects/Personal/Ethanol/GEO/HiC/T2_mega_q30.hic"}
CR={c:pyBigWig.open(f)
    for c,f in [("ctrl","/tmp/claude-1000/-Zulu-jordan-alcoholATAC/0804ae8d-3e68-44fe-ac8d-cac57f0fb3df/scratchpad/crush/ctrl/ctrlmergedCrush_5000.bw"),("etoh","/tmp/claude-1000/-Zulu-jordan-alcoholATAC/0804ae8d-3e68-44fe-ac8d-cac57f0fb3df/scratchpad/crush/etoh/etohmergedCrush_5000.bw"),("with","/tmp/claude-1000/-Zulu-jordan-alcoholATAC/0804ae8d-3e68-44fe-ac8d-cac57f0fb3df/scratchpad/crush/with/withmergedCrush_5000.bw")]}
gcbw=pyBigWig.open("chipseq/hg38.gc5Base.bw"); gok=set(gcbw.chroms())
CHR=[str(c) for c in range(1,23)]+["X"]

def track(bw,ch,nb,prefix_ok):
    key=("chr"+ch) if ("chr"+ch) in prefix_ok else ch
    if key not in prefix_ok: return None
    L=bw.chroms()[key]
    n=min(nb,L//RES)
    v=np.array(bw.stats(key,0,n*RES,nBins=n,type="mean"),dtype=float)
    out=np.full(nb,np.nan); out[:n]=v; return out

def band(cond,ch):
    cmd=(f"straw NONE {HIC[cond]} {ch} {ch} BP {RES} 2>/dev/null | "
         f"awk '{{d=$2-$1; if(d<0)d=-d; if(d<={(OFFMAX+1)*RES}) print $1,$2,$3}}'")
    p=subprocess.run(["bash","-c",cmd],capture_output=True,text=True)
    if not p.stdout.strip(): return None,0
    df=pd.read_csv(pd.io.common.StringIO(p.stdout),sep=" ",header=None,names=["i","j","v"])
    bi=df["i"].to_numpy()//RES; bj=df["j"].to_numpy()//RES
    lo=np.minimum(bi,bj); hi=np.maximum(bi,bj); off=hi-lo; nb=int(hi.max())+1
    M=np.zeros((nb,OFFMAX+1),dtype=np.float64)
    ok=off<=OFFMAX
    np.add.at(M,(lo[ok],off[ok]),df["v"].to_numpy()[ok])
    return M,nb

# ---- pass 1: collect the axis values genome-wide to set global quantiles ----
AX={"crush":{}, "gc":{}}
NB={}
for ch in CHR:
    M,nb=band("ctrl",ch)
    if M is None: continue
    NB[ch]=nb
    AX["crush"][ch]=track(CR["ctrl"],ch,nb,set(CR["ctrl"].chroms()))
    AX["gc"][ch]=track(gcbw,ch,nb,gok)
for axis in ("crush","gc"):
    allv=np.concatenate([v for v in AX[axis].values() if v is not None])
    q=np.nanquantile(allv,[i/NQ for i in range(1,NQ)])
    AX[axis]["_q"]=q
    print(f"{axis} axis quantile edges ready (n={np.isfinite(allv).sum()})")

# ---- pass 2: saddle per condition, per band, per axis ----
res={}
for cond in HIC:
    S={(a,b):np.zeros((NQ,NQ)) for a in ("crush","gc") for b in range(len(BANDS))}
    N={(a,b):np.zeros((NQ,NQ)) for a in ("crush","gc") for b in range(len(BANDS))}
    for ch in CHR:
        M,nb=band(cond,ch)
        if M is None: continue
        exp=np.array([M[:nb-o,o].sum()/max(nb-o,1) for o in range(OFFMAX+1)]); exp[exp<=0]=np.nan
        for axis in ("crush","gc"):
            a=AX[axis].get(ch)
            if a is None: continue
            a=a[:nb] if len(a)>=nb else np.pad(a,(0,nb-len(a)),constant_values=np.nan)
            qi=np.digitize(a,AX[axis]["_q"]); qi=np.where(np.isfinite(a),qi,-1)
            for bidx,(lo_d,hi_d) in enumerate(BANDS):
                for o in range(max(1,lo_d//RES),min(OFFMAX,hi_d//RES)+1):
                    if not np.isfinite(exp[o]): continue
                    v=M[:nb-o,o]/exp[o]
                    q1=qi[:nb-o]; q2=qi[o:nb]
                    m=(q1>=0)&(q2>=0)
                    if not m.any(): continue
                    np.add.at(S[(axis,bidx)],(q1[m],q2[m]),v[m])
                    np.add.at(N[(axis,bidx)],(q1[m],q2[m]),1)
    res[cond]={k:(S[k]/np.maximum(N[k],1)) for k in S}
    print(f"  {cond} saddle done",flush=True)
pickle.dump(res,open(f"{OUT}/saddle_newcrush.pkl","wb"))

def strength(m,k=4):
    BB=np.nanmean(m[:k,:k]); AA=np.nanmean(m[-k:,-k:])
    AB=(np.nanmean(m[:k,-k:])+np.nanmean(m[-k:,:k]))/2
    return (AA+BB)/(2*AB), AA, BB, AB
print("\n=== COMPARTMENTALISATION STRENGTH  (AA+BB)/2AB ===")
for axis in ("crush","gc"):
    print(f"\n  axis = {axis.upper()}")
    print(f"    {'band':16s}{'ctrl':>9s}{'etoh':>9s}{'with':>9s}{'EtOH/Ctl':>11s}")
    for bidx,(lo_d,hi_d) in enumerate(BANDS):
        s={c:strength(res[c][(axis,bidx)])[0] for c in HIC}
        print(f"    {lo_d//1000}-{hi_d//1000}kb{'':6s}{s['ctrl']:9.4f}{s['etoh']:9.4f}{s['with']:9.4f}{s['etoh']/s['ctrl']:11.4f}")
