import re, os, numpy as np
# CTCF peak -> ZF5 triplet + motif score + coords.  Same conventions as zf5_2dmatched.py:
# HOMER ctcf.motif is the reverse complement of the manuscript ZF numbering, so
# ZF5 GGG = CCC at [8:11] of the oriented 20-mer, GAT = ATC.
H="bensalcohol_out/humanZF"; OUT="bensalcohol_out/loops"; os.makedirs(OUT,exist_ok=True)
COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return "".join(COMP[b] for b in reversed(s))
PAT=re.compile(r'(-?\d+)\(([ACGTN]+),([+-]),[-0-9.]+\)'); L=20

coord={}
for l in open(f"{H}/peaks_pm100.bed"):
    p=l.rstrip("\n").split("\t")
    coord[p[3]]=(p[0],(int(p[1])+int(p[2]))//2)      # sid -> (chr, summit)
score={}
for l in open(f"{H}/annot_ctcf.txt"):
    if l.startswith("PeakID"): continue
    p=l.rstrip("\n").split("\t")
    try: score[p[0]]=float(p[-1])
    except: pass
fc={}
for l in open(f"{H}/fc_table.tsv"):
    if l.startswith("sid"): continue
    p=l.rstrip("\n").split("\t"); fc[p[0]]=(float(p[5]),p[6])   # log2fc(E6/UT), decile cat

seq={};strand={}
for l in open(f"{H}/annot_seq.txt"):
    if l.startswith("PeakID"): continue
    p=l.rstrip("\n").split("\t"); ins=PAT.findall(p[-1])
    if not ins: continue
    b=min(ins,key=lambda m:abs(int(m[0])-100))
    s=b[1] if b[2]=='+' else rc(b[1])
    if len(s)==L and 'N' not in s:
        seq[p[0]]=s; strand[p[0]]=b[2]

n=0
with open(f"{OUT}/ctcf_anchor_annot.bed","w") as o:
    for sid,s in seq.items():
        if sid not in coord or sid not in score: continue
        ch,su=coord[sid]; tri=s[8:11]
        zf5 = "GGG" if tri=="CCC" else ("GAT" if tri=="ATC" else "other")
        l2,cat = fc.get(sid,(float('nan'),"NA"))
        o.write(f"{ch}\t{su}\t{su+1}\t{sid}\t{score[sid]:.4f}\t{strand[sid]}\t{zf5}\t{tri}\t{l2:.4f}\t{cat}\n")
        n+=1
print(f"wrote {OUT}/ctcf_anchor_annot.bed  n={n}")
import collections
c=collections.Counter(("GGG" if seq[s][8:11]=="CCC" else "GAT" if seq[s][8:11]=="ATC" else "other")
                      for s in seq if s in coord and s in score)
print("ZF5 class counts:",dict(c))
