import numpy as np, pandas as pd, collections
OUT="bensalcohol_out/loops"; RES=10000
E="/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/RNA/merged_expression_results.tsv"
H=None;G=[]
for l in open(E):
    p=l.rstrip("\n").split("\t")
    if H is None: H={k:i for i,k in enumerate(p)}; continue
    try:
        cT=float(p[H['cTPM']]); ch=p[H['chr']]; sd=p[H['strand']].strip()
        tss=int(p[H['start']]) if sd=='+' else int(p[H['stop']])
        le=float(p[H['log2FoldChange.x']]); lw=float(p[H['log2FoldChange.y']])
        pe=p[H['padj.x']]; pw=p[H['padj.y']]; g=p[H['symbol']].strip()
    except: continue
    if not ch.startswith("chr") or cT<=1: continue
    pe=float(pe) if pe not in ("NA","") else 1.0; pw=float(pw) if pw not in ("NA","") else 1.0
    G.append((ch,tss,g,le,lw,pe,pw,cT))
print(f"expressed genes: {len(G)}")
# gene trajectory class, manuscript rule
def gcls(le,lw,pe,pw):
    if not((pe<0.05)or(pw<0.05)): return "notDE"
    ce=abs(le)>0.585; cw=abs(lw)>0.585
    if ce and not cw: return "Acute"
    if ce and cw:     return "Lasting"
    if cw and not ce: return "Delayed"
    return "NotAffected"
prom={}   # (chr,bin) -> list of gene records
for ch,tss,g,le,lw,pe,pw,cT in G:
    for b in {(tss-2000)//RES,(tss+2000)//RES}:
        prom.setdefault((ch,b),[]).append((g,le,lw,gcls(le,lw,pe,pw)))
print(f"promoter bins: {len(prom)}")
# ---- enhancers: union H3K27ac, excluding promoter-proximal ----
D="/Zulu/bnolan/Projects/Personal/Ethanol/Assays/CUTNTAG/nf_k27ac_k27me3_k9me3_9-3/macs3"
tssidx=collections.defaultdict(list)
for ch,tss,*_ in G: tssidx[(ch,tss//RES)].append(tss)
def near_tss(ch,p,w=2000):
    for b in range((p-w)//RES,(p+w)//RES+1):
        for t in tssidx.get((ch,b),[]):
            if abs(t-p)<=w: return True
    return False
enh=set(); npk=0
for s in ("k27acctrl_A","k27aceth_A","k27acrec_A"):
    for l in open(f"{D}/{s}/{s}_peaks.narrowPeak"):
        p=l.rstrip("\n").split("\t")
        if len(p)<10 or not p[0].startswith("chr"): continue
        summit=int(p[1])+int(p[9]); npk+=1
        if near_tss(p[0],summit): continue
        enh.add((p[0],summit//RES))
print(f"H3K27ac peaks read: {npk}   distal enhancer bins: {len(enh)}")
enh_only=enh-set(prom.keys())
print(f"enhancer bins after removing promoter bins: {len(enh_only)}")
with open(f"{OUT}/ep_bins.tsv","w") as o:
    o.write("chr\tbin\ttype\tgenes\n")
    for (ch,b),recs in prom.items():
        o.write(f"{ch}\t{b}\tP\t{','.join(r[0] for r in recs)}\n")
    for ch,b in enh_only: o.write(f"{ch}\t{b}\tE\t.\n")
pd.DataFrame([(ch,tss,g,le,lw,gcls(le,lw,pe,pw)) for ch,tss,g,le,lw,pe,pw,cT in G],
             columns=["chr","tss","gene","l2FC_etoh","l2FC_with","class"]).to_csv(f"{OUT}/genes_classed.tsv",sep="\t",index=False)
import collections as C
print("gene classes:",dict(C.Counter(gcls(le,lw,pe,pw) for ch,tss,g,le,lw,pe,pw,cT in G)))
print(f"wrote {OUT}/ep_bins.tsv and genes_classed.tsv")
