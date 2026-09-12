import numpy as np, collections, os
OUT="bensalcohol_out/loops"
LOOPS="/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/Hi-C/mergednodups/refinedloops/filtered_loops_with_oe.bedpe"

# ---- CTCF anchor annotation, indexed by (chr, 25kb bin) ----
BIN=25000
idx=collections.defaultdict(list)
for l in open(f"{OUT}/ctcf_anchor_annot.bed"):
    p=l.rstrip("\n").split("\t")
    ch=p[0]; pos=int(p[1]); sc=float(p[4]); st=p[5]; zf5=p[6]; l2=float(p[8])
    idx[(ch,pos//BIN)].append((sc,zf5,st,l2,pos))
print("annotated CTCF sites:",sum(len(v) for v in idx.values()))

def anchor_call(ch,s,e):
    """best CTCF site (max motif score) whose summit falls in [s,e)."""
    best=None
    for b in range(s//BIN, e//BIN + 1):
        for rec in idx.get((ch,b),[]):
            if s<=rec[4]<e and (best is None or rec[0]>best[0]): best=rec
    return best

rows=[];miss=0
for l in open(LOOPS):
    p=l.rstrip("\n").split("\t")
    if p[0]=="0": continue                     # header
    c1="chr"+p[0]; x1,x2=int(p[1]),int(p[2])
    c2="chr"+p[3]; y1,y2=int(p[4]),int(p[5])
    oe_c,oe_e,oe_w=float(p[9]),float(p[10]),float(p[11])
    if not all(np.isfinite([oe_c,oe_e,oe_w])) or min(oe_c,oe_e,oe_w)<=0: miss+=1; continue
    a=anchor_call(c1,x1,x2); b=anchor_call(c2,y1,y2)
    za = a[1] if a else "none"; zb = b[1] if b else "none"
    sa = a[0] if a else np.nan; sb = b[0] if b else np.nan
    ora= a[2] if a else "."   ; orb= b[2] if b else "."
    dist=((y1+y2)//2)-((x1+x2)//2)
    rows.append([c1,x1,x2,c2,y1,y2,dist,oe_c,oe_e,oe_w,za,zb,sa,sb,ora,orb])
print(f"loops kept={len(rows)}  dropped(nonpositive OE)={miss}")

d=np.array([r[6] for r in rows],float)
oc=np.array([r[7] for r in rows]); oe=np.array([r[8] for r in rows]); ow=np.array([r[9] for r in rows])
l2e=np.log2(oe/oc); l2w=np.log2(ow/oc)
def dec(v):
    q=np.quantile(v,[i/10 for i in range(1,10)]); return np.digitize(v,q)+1
de,dw=dec(l2e),dec(l2w)

with open(f"{OUT}/loop_table.tsv","w") as o:
    o.write("chr1\tx1\tx2\tchr2\ty1\ty2\tdist\toe_ctrl\toe_etoh\toe_with\t"
            "zf5_a\tzf5_b\tscore_a\tscore_b\tstrand_a\tstrand_b\tlog2FC_etoh\tlog2FC_with\tdec_etoh\tdec_with\n")
    for i,r in enumerate(rows):
        o.write("\t".join(map(str,r))+f"\t{l2e[i]:.5f}\t{l2w[i]:.5f}\t{de[i]}\t{dw[i]}\n")

print(f"\nlog2FC EtOH/Ctl : median={np.median(l2e):+.4f}  up={int((l2e>0).sum())} down={int((l2e<0).sum())}")
print(f"log2FC With/Ctl : median={np.median(l2w):+.4f}  up={int((l2w>0).sum())} down={int((l2w<0).sum())}")
print("decile median log2FC(EtOH):"," ".join(f"{np.median(l2e[de==k]):+.2f}" for k in range(1,11)))
print("loop distance: median=%.0f kb  range %.0f-%.0f kb"%(np.median(d)/1e3,d.min()/1e3,d.max()/1e3))
ca=collections.Counter(r[10] for r in rows); cb=collections.Counter(r[11] for r in rows)
print("anchor A ZF5:",dict(ca)); print("anchor B ZF5:",dict(cb))
both=collections.Counter(tuple(sorted((r[10],r[11]))) for r in rows)
print("\nanchor-pair classes (top):")
for k,v in both.most_common(10): print(f"   {k}: {v}")
