import numpy as np, collections
OUT="bensalcohol_out/loops"; BIN=25000
# per-anchor best-scoring motif SEPARATELY for each strand
idx=collections.defaultdict(list)
for l in open(f"{OUT}/ctcf_anchor_annot.bed"):
    p=l.rstrip("\n").split("\t")
    idx[(p[0],int(p[1])//BIN)].append((float(p[4]),p[5],p[6],int(p[1])))  # score,strand,zf5,pos
def best_by_strand(ch,s,e):
    b={'+':None,'-':None}
    for k in range(s//BIN,e//BIN+1):
        for r in idx.get((ch,k),[]):
            if s<=r[3]<e:
                if b[r[1]] is None or r[0]>b[r[1]][0]: b[r[1]]=r
    return b

H=None;R=[]
for l in open(f"{OUT}/loop_table.tsv"):
    p=l.rstrip("\n").split("\t")
    if H is None: H={k:i for i,k in enumerate(p)}; continue
    R.append(p)
print(f"loops={len(R)}")

A=[];B=[]
for r in R:
    A.append(best_by_strand(r[H['chr1']],int(r[H['x1']]),int(r[H['x2']])))
    B.append(best_by_strand(r[H['chr2']],int(r[H['y1']]),int(r[H['y2']])))

# Which pairing is the convergent one?  Compare loop anchors against DISTANCE-MATCHED
# random anchor pairs drawn from the same anchor pool (shuffle which right anchor goes
# with which left anchor, preserving each anchor's own motif content).
rng=np.random.default_rng(3)
def counts(Alist,Blist):
    c=collections.Counter()
    for a,b in zip(Alist,Blist):
        for sa in '+-':
            for sb in '+-':
                if a[sa] and b[sb]: c[sa+sb]+=1
    return c
obs=counts(A,B)
perm=collections.Counter()
NSH=20
for _ in range(NSH):
    o=rng.permutation(len(B))
    c=counts(A,[B[i] for i in o])
    for k,v in c.items(): perm[k]+=v/NSH
print("\nanchor motif-strand pairing (left,right):   observed   shuffled   obs/exp")
for k in ("++","+-","-+","--"):
    print(f"  {k}:  {obs[k]:7d}   {perm[k]:9.0f}   {obs[k]/perm[k]:.3f}")
# Convergent = the dominant pairing.  The shuffle preserves each anchor's own strand
# content, so obs/exp has no power here; the 12.6x count asymmetry between the two
# anti-parallel pairings is itself the convergence signal.  HOMER's ctcf.motif is the
# reverse complement of the canonical CTCF orientation, so '-'/'+' here = canonical
# convergent (each anchor motif pointing INTO the loop).
conv = max(("+-","-+"), key=lambda k: obs[k])
print(f"\n=> convergent configuration in this annotation = left '{conv[0]}' / right '{conv[1]}'")

# assign ZF5 at the CONVERGENT sites only
za=[];zb=[];sa=[];sb=[]
for a,b in zip(A,B):
    ra=a[conv[0]]; rb=b[conv[1]]
    za.append(ra[2] if ra else "none"); zb.append(rb[2] if rb else "none")
    sa.append(ra[0] if ra else np.nan); sb.append(rb[0] if rb else np.nan)
za=np.array(za);zb=np.array(zb)
both=(za!="none")&(zb!="none")
print(f"\nloops with a convergent CTCF motif at BOTH anchors: {int(both.sum())}")
cc=collections.Counter(tuple(sorted((x,y))) for x,y in zip(za[both],zb[both]))
print("convergent-pair ZF5 classes:")
for k,v in cc.most_common(): print(f"   {k}: {v}")
with open(f"{OUT}/conv_anchor_zf5.tsv","w") as o:
    o.write("zf5_conv_a\tzf5_conv_b\tscore_a\tscore_b\n")
    for i in range(len(R)):
        o.write(f"{za[i]}\t{zb[i]}\t{sa[i]}\t{sb[i]}\n")
print(f"\nwrote {OUT}/conv_anchor_zf5.tsv")
