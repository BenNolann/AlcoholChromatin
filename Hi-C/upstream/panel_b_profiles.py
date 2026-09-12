import numpy as np, collections, pyBigWig, pickle
OUT="bensalcohol_out/loops"; FL=2000; NB=80   # +/-2 kb, 50 bp bins
S=[]
for l in open(f"{OUT}/ctcf_anchor_annot.bed"):
    p=l.rstrip("\n").split("\t"); S.append((p[0],int(p[1]),float(p[4]),p[5]))
BIN=25000; idx=collections.defaultdict(list)
for i,s in enumerate(S): idx[(s[0],s[1]//BIN)].append(i)
H=None;R=[]
for l in open(f"{OUT}/loop_table.tsv"):
    p=l.rstrip("\n").split("\t")
    if H is None: H={k:i for i,k in enumerate(p)}; continue
    R.append(p)
def best(ch,a,b):
    bi=None
    for k in range(a//BIN,b//BIN+1):
        for i in idx.get((ch,k),[]):
            if a<=S[i][1]<b and (bi is None or S[i][2]>S[bi][2]): bi=i
    return bi
anch=set()
for r in R:
    for c,a,b in ((r[H['chr1']],int(r[H['x1']]),int(r[H['x2']])),(r[H['chr2']],int(r[H['y1']]),int(r[H['y2']]))):
        i=best(c,a,b)
        if i is not None: anch.add(i)
anch=sorted(anch)
print(f"CTCF loop-anchor sites for profiles: {len(anch)}",flush=True)
BW={
 "CTCF":  [("Control","chipseq/CTCF_UT_10bp_BPM.bw"),("EtOH","chipseq/CTCF_E6_10bp_MABPM.bw"),("Withdrawal","chipseq/CTCF_E6R6_10bp_MABPM.bw")],
 "RAD21": [("Control","chipseq/RAD21_UT_10bp_BPM.bw"),("EtOH","chipseq/RAD21_E6_10bp_MABPM.bw"),("Withdrawal","chipseq/RAD21_E6R6_10bp_MABPM.bw")],
 "RNAPIIs2p":[("Control","chipseq/RNAPII_UT_10bp_BPM.bw"),("EtOH","chipseq/RNAPII_E6_10bp_MABPM.bw"),("Withdrawal","chipseq/RNAPII_E6R6_10bp_MABPM.bw")],
}
PROF={}
for prot,items in BW.items():
    for lab,f in items:
        bw=pyBigWig.open(f); ok=set(bw.chroms())
        acc=np.zeros(NB); n=0
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
        print(f"  {prot} {lab}: n={n}",flush=True)
pickle.dump({"prof":PROF,"flank":FL,"nb":NB,"n":len(anch)},open(f"{OUT}/panelb_profiles.pkl","wb"))
print(f"wrote {OUT}/panelb_profiles.pkl")
