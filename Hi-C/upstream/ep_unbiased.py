"""E-P pairs defined ONLY by annotation (enhancer bin x promoter bin, 20kb-2Mb).
No FitHiC significance selection anywhere -> no ascertainment bias."""
import numpy as np, pandas as pd, subprocess, collections, pyBigWig
RES=10000; LO,HI=20000,2000000; OFFMAX=HI//RES
OUT="bensalcohol_out/loops"
HIC={"ctrl":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/ctrl_merged_nodups_q30_592m.hic",
     "etoh":"/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/etoh_merged_nodups_q30_592m.hic",
     "with":"/Zulu/bnolan/Projects/Personal/Ethanol/GEO/HiC/T2_mega_q30.hic"}
# ---- annotation ----
P=collections.defaultdict(set); E=collections.defaultdict(set); genes={}
for l in open(f"{OUT}/ep_bins.tsv"):
    if l.startswith("chr\t"): continue
    p=l.rstrip("\n").split("\t"); ch=p[0].replace("chr",""); b=int(p[1])
    if p[2]=="P": P[ch].add(b); genes[(ch,b)]=p[3]
    else: E[ch].add(b)
print(f"promoter bins={sum(len(v) for v in P.values())}  enhancer bins={sum(len(v) for v in E.values())}")
# ---- ZF5 per bin ----
best={}
for l in open(f"{OUT}/ctcf_anchor_annot.bed"):
    p=l.rstrip("\n").split("\t"); k=(p[0].replace("chr",""),int(p[1])//RES); sc=float(p[4])
    if k not in best or sc>best[k][0]: best[k]=(sc,p[6])
gcbw=pyBigWig.open("chipseq/hg38.gc5Base.bw"); gok=set(gcbw.chroms())
rows=[]
for ch in sorted(P, key=lambda c:(len(c),c)):
    if ch in ("Y","M") or ch not in [str(i) for i in range(1,23)]+["X"]: continue
    pb=sorted(P[ch]); eb=E[ch]
    pairs=[(min(a,b),max(a,b)) for a in pb for b in range(a-OFFMAX,a+OFFMAX+1)
           if b in eb and LO<=abs(b-a)*RES<=HI]
    pairs=sorted(set(pairs))
    if not pairs: continue
    M={}
    for cond in HIC:
        cmd=(f"straw NONE {HIC[cond]} {ch} {ch} BP {RES} 2>/dev/null | "
             f"awk '{{d=$2-$1; if(d<0)d=-d; if(d<={OFFMAX*RES}) print $1,$2,$3}}'")
        q=subprocess.run(["bash","-c",cmd],capture_output=True,text=True)
        if not q.stdout.strip(): break
        df=pd.read_csv(pd.io.common.StringIO(q.stdout),sep=" ",header=None,names=["i","j","v"])
        bi=df["i"].to_numpy()//RES; bj=df["j"].to_numpy()//RES
        lo=np.minimum(bi,bj); hi=np.maximum(bi,bj); off=hi-lo; nb=int(hi.max())+1
        X=np.zeros((nb,OFFMAX+1)); ok=off<=OFFMAX
        np.add.at(X,(lo[ok],off[ok]),df["v"].to_numpy()[ok])
        exp=np.array([X[:nb-o,o].sum()/max(nb-o,1) for o in range(OFFMAX+1)])
        M[cond]=(X,exp,nb)
    if len(M)<3: continue
    gk="chr"+ch if "chr"+ch in gok else ch
    L=gcbw.chroms()[gk]; nbg=L//RES
    G=np.array(gcbw.stats(gk,0,nbg*RES,nBins=nbg,type="mean"),dtype=float)
    for a,b in pairs:
        o=b-a
        if any(a>=M[c][2] for c in M): continue
        vals={}
        for c in M:
            X,exp,nb=M[c]; vals[c]=(X[a,o]+1)/(exp[o]+1)
        z1=best.get((ch,a),(0,"none"))[1]; z2=best.get((ch,b),(0,"none"))[1]
        gcv=(G[a]+G[b])/2 if a<len(G) and b<len(G) else np.nan
        pb_=a if a in P[ch] else b
        rows.append((("chr"+ch),a*RES,b*RES,o*RES,vals["ctrl"],vals["etoh"],vals["with"],
                     z1,z2,gcv,genes.get((ch,pb_),".")))
    print(f"  chr{ch}: {len(pairs)} E-P pairs",flush=True)
D=pd.DataFrame(rows,columns=["chr","m1","m2","dist","dn_ctrl","dn_etoh","dn_with","z1","z2","gc","genes"])
D["L2E"]=np.log2(D["dn_etoh"]/D["dn_ctrl"]); D["L2W"]=np.log2(D["dn_with"]/D["dn_ctrl"])
anyG=(D["z1"]=="GGG")|(D["z2"]=="GGG"); anyA=(D["z1"]=="GAT")|(D["z2"]=="GAT")
D["ctcfcls"]=np.select([anyG&~anyA,anyA&~anyG,anyG&anyA,(D["z1"]!="none")|(D["z2"]!="none")],
                       ["GGG","GAT","GGG+GAT","otherCTCF"],default="noCTCF")
D.to_csv(f"{OUT}/ep_unbiased.tsv",sep="\t",index=False)
print(f"\nTOTAL E-P pairs (annotation-defined, no significance filter): {len(D):,}")
print(f"  mean L2E={D['L2E'].mean():+.4f}  median={D['L2E'].median():+.4f}")
print(f"  mean L2W={D['L2W'].mean():+.4f}  median={D['L2W'].median():+.4f}")
print(f"wrote {OUT}/ep_unbiased.tsv")
