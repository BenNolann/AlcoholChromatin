"""Assign which anchor of each E-P pair is the enhancer, and record its ZF5 class.
ep_unbiased.tsv -> ep_unbiased_EP.tsv   (adds z_enh, z_prom, ctcf_at_enh)"""
import numpy as np, pandas as pd
L="/Zulu/jordan/alcoholATAC/bensalcohol_out/loops"; RES=10000
D=pd.read_csv(f"{L}/ep_unbiased.tsv",sep="\t")
Eb=set()
for l in open(f"{L}/ep_bins.tsv"):
    if l.startswith("chr\t"): continue
    p=l.rstrip("\n").split("\t")
    if p[2]=="E": Eb.add((p[0],int(p[1])))
k1=list(zip(D["chr"],D["m1"]//RES)); k2=list(zip(D["chr"],D["m2"]//RES))
m1E=np.array([k in Eb for k in k1]); m2E=np.array([k in Eb for k in k2])
ok=m1E^m2E
D=D[ok].copy(); sel=m1E[ok]
D["z_enh"]=np.where(sel,D["z1"],D["z2"]); D["z_prom"]=np.where(sel,D["z2"],D["z1"])
D["ctcf_at_enh"]=D["z_enh"]!="none"
D.to_csv(f"{L}/ep_unbiased_EP.tsv",sep="\t",index=False)
print(f"pairs={len(D):,}  CTCF at enhancer={int(D['ctcf_at_enh'].sum()):,}")
