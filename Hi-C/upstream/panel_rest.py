import numpy as np, pandas as pd, subprocess, pickle, pyBigWig
OUT="bensalcohol_out/loops"; RES=10000
HIC={"ctrl":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/ctrl_merged_nodups_q30_592m.hic",
     "etoh":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/etoh_merged_nodups_q30_592m.hic",
     "with":"/Zulu/bnolan/Projects/Personal/Ethanol/GEO/HiC/T2_mega_q30.hic"}
P={}
# ---------- (a) example region matrix, 250 kb bins, chr1 ----------
EX=("1",0,100_000_000); EXRES=250000
mats={}
for c in HIC:
    cmd=f"straw NONE {HIC[c]} {EX[0]}:{EX[1]}:{EX[2]} {EX[0]}:{EX[1]}:{EX[2]} BP {EXRES} 2>/dev/null"
    p=subprocess.run(["bash","-c",cmd],capture_output=True,text=True)
    df=pd.read_csv(pd.io.common.StringIO(p.stdout),sep="\t",header=None,names=["i","j","v"])
    n=(EX[2]-EX[1])//EXRES+1
    M=np.zeros((n,n))
    bi=(df["i"].to_numpy()-EX[1])//EXRES; bj=(df["j"].to_numpy()-EX[1])//EXRES
    ok=(bi>=0)&(bj>=0)&(bi<n)&(bj<n)
    np.add.at(M,(bi[ok],bj[ok]),df["v"].to_numpy()[ok]); M=M+M.T-np.diag(np.diag(M))
    mats[c]=M
P["example"]={"mats":mats,"region":EX,"res":EXRES}
print("example region done",flush=True)
# ---------- (b) contact decay curves ----------
dec={c:None for c in HIC}
for c in HIC:
    tot=np.zeros(3000); cnt=np.zeros(3000)
    for ch in ["1","2","3"]:
        cmd=(f"straw NONE {HIC[c]} {ch} {ch} BP {RES} 2>/dev/null | awk '{{d=$2-$1; if(d<0)d=-d; print int(d/{RES}), $3}}'")
        p=subprocess.run(["bash","-c",cmd],capture_output=True,text=True)
        a=np.loadtxt(pd.io.common.StringIO(p.stdout))
        if a.ndim==1: a=a[None,:]
        o=a[:,0].astype(int); v=a[:,1]
        m=o<3000
        np.add.at(tot,o[m],v[m]); np.add.at(cnt,o[m],1)
    dec[c]=tot/np.maximum(cnt,1)
    print(f"  decay {c} done",flush=True)
P["decay"]=dec
# ---------- (f) compartment change by GC within A ----------
cr=pyBigWig.open("bensalcohol_out/crush_rs/ctrlmergedCrush_5000.bw")
gcbw=pyBigWig.open("chipseq/hg38.gc5Base.bw"); gok=set(gcbw.chroms())
OFF=200
res_f={}
for c in HIC: res_f[c]={}
CH=[str(x) for x in range(1,23)]
acc={c:{} for c in HIC}
for ch in CH:
    key="chr"+ch if "chr"+ch in set(cr.chroms()) else ch
    if key not in set(cr.chroms()): continue
    L=cr.chroms()[key]; nb=L//RES
    A=np.array(cr.stats(key,0,nb*RES,nBins=nb,type="mean"),dtype=float)
    gk="chr"+ch if "chr"+ch in gok else ch
    G=np.array(gcbw.stats(gk,0,min(nb*RES,gcbw.chroms()[gk]),nBins=min(nb,gcbw.chroms()[gk]//RES),type="mean"),dtype=float)
    n=min(len(A),len(G)); A=A[:n]; G=G[:n]
    isA=A>np.nanmedian(A)
    for c in HIC:
        cmd=(f"straw NONE {HIC[c]} {ch} {ch} BP {RES} 2>/dev/null | "
             f"awk '{{d=$2-$1; if(d<0)d=-d; if(d>=50000 && d<=500000) print $1,$2,$3}}'")
        p=subprocess.run(["bash","-c",cmd],capture_output=True,text=True)
        if not p.stdout.strip(): continue
        df=pd.read_csv(pd.io.common.StringIO(p.stdout),sep=" ",header=None,names=["i","j","v"])
        bi=df["i"].to_numpy()//RES; bj=df["j"].to_numpy()//RES
        lo=np.minimum(bi,bj); hi=np.maximum(bi,bj); off=hi-lo
        M=np.zeros((n,OFF+1)); ok=(off<=OFF)&(hi<n)
        np.add.at(M,(lo[ok],off[ok]),df["v"].to_numpy()[ok])
        exp=np.array([M[:n-o,o].sum()/max(n-o,1) for o in range(OFF+1)]); exp[exp<=0]=np.nan
        for o in range(5,51):
            if not np.isfinite(exp[o]): continue
            v=M[:n-o,o]/exp[o]
            a1=isA[:n-o]; a2=isA[o:n]; g1=G[:n-o]; g2=G[o:n]
            m=a1&a2&np.isfinite(g1)&np.isfinite(g2)
            if not m.any(): continue
            gm=(g1[m]+g2[m])/2
            acc[c].setdefault("gc",[]).append(gm); acc[c].setdefault("oe",[]).append(v[m])
    print(f"  compartment chr{ch} done",flush=True)
CG={}
for c in HIC:
    gcv=np.concatenate(acc[c]["gc"]); oev=np.concatenate(acc[c]["oe"]); CG[c]=(gcv,oev)
q=np.nanquantile(CG["ctrl"][0],[i/5 for i in range(1,5)])
P["AA_by_gc"]={c:[np.nanmean(CG[c][1][np.digitize(CG[c][0],q)==k]) for k in range(5)] for c in HIC}
P["AA_by_gc_edges"]=q
print("AA by GC quintile:",{c:np.round(P['AA_by_gc'][c],4) for c in HIC})
pickle.dump(P,open(f"{OUT}/panels.pkl","wb"))
print(f"wrote {OUT}/panels.pkl")
