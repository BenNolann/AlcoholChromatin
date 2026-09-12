"""MA-normalise the E-P contact log-ratios: remove the intensity-dependent (M vs A) trend
by subtracting a running median of M across A. Applied to both comparisons for symmetry."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
OUT="bensalcohol_out/loops"
D=pd.read_csv(f"{OUT}/ep_unbiased_EP.tsv",sep="\t")
for c in ("ctrl","etoh","with"): D=D[D[f"dn_{c}"]>0]
D=D.reset_index(drop=True)
def ma_norm(dn_num,dn_den,nbin=200):
    M=np.log2(dn_num/dn_den); A=0.5*np.log2(dn_num*dn_den)
    edges=np.nanquantile(A,np.linspace(0,1,nbin+1)); edges[0]-=1e-9; edges[-1]+=1e-9
    idx=np.digitize(A,edges[1:-1])
    med=np.full(nbin,np.nan); ctr=np.full(nbin,np.nan)
    for k in range(nbin):
        m=idx==k
        if m.sum()>20: med[k]=np.median(M[m]); ctr[k]=np.median(A[m])
    ok=np.isfinite(med)&np.isfinite(ctr)
    trend=np.interp(A,ctr[ok],med[ok])
    return M-trend, M, A
D["L2E_ma"],M_e,A_e=ma_norm(D["dn_etoh"].values,D["dn_ctrl"].values)
D["L2W_ma"],M_w,A_w=ma_norm(D["dn_with"].values,D["dn_ctrl"].values)
D["A"]=A_e
print(f"n={len(D):,}")
print(f"  before MA:  mean M etoh={M_e.mean():+.4f}   with={M_w.mean():+.4f}")
print(f"  after  MA:  mean M etoh={D['L2E_ma'].mean():+.4f}   with={D['L2W_ma'].mean():+.4f}")
qa=np.nanquantile(A_e,[i/10 for i in range(1,10)]); b=np.digitize(A_e,qa)
print("\n  residual A-trend after MA (should be ~flat):")
print(f"    EtOH  D1={D.loc[b==0,'L2E_ma'].mean():+.4f}  D10={D.loc[b==9,'L2E_ma'].mean():+.4f}")
print(f"    With  D1={D.loc[b==0,'L2W_ma'].mean():+.4f}  D10={D.loc[b==9,'L2W_ma'].mean():+.4f}")
# ---- key results recomputed ----
print("\n=== GC gradient, MA-normalised ===")
q=np.nanquantile(D["gc"],[i/5 for i in range(1,5)]); bg=np.digitize(D["gc"],q)
print(f"  {'':<5}{'EtOH':>10}{'+-95%':>9}{'With':>10}{'+-95%':>9}")
for k in range(5):
    s=D[bg==k]
    print(f"  Q{k+1:<4}{s['L2E_ma'].mean():+10.4f}{1.96*s['L2E_ma'].sem():9.4f}"
          f"{s['L2W_ma'].mean():+10.4f}{1.96*s['L2W_ma'].sem():9.4f}")
se=D[bg==4]["L2E_ma"].mean()-D[bg==0]["L2E_ma"].mean(); sw=D[bg==4]["L2W_ma"].mean()-D[bg==0]["L2W_ma"].mean()
print(f"  Q5-Q1: EtOH={se:+.4f}  With={sw:+.4f}  -> {100*sw/se:.0f}% retained")
print("\n=== by gene class, MA-normalised ===")
D["gene"]=D["genes"].str.split(","); X=D.explode("gene"); X["gene"]=X["gene"].str.strip()
G=pd.read_csv(f"{OUT}/genes_classed.tsv",sep="\t")[["gene","class"]]
X=X.merge(G,on="gene",how="inner")
print(f"  {'class':11s}{'n pairs':>10}{'EtOH':>10}{'+-95%':>9}{'With':>10}{'+-95%':>9}")
for c in ("Lasting","Delayed","NotAffected","Acute","notDE"):
    s=X[X["class"]==c]
    if len(s)<500: continue
    print(f"  {c:11s}{len(s):10,}{s['L2E_ma'].mean():+10.4f}{1.96*s['L2E_ma'].sem():9.4f}"
          f"{s['L2W_ma'].mean():+10.4f}{1.96*s['L2W_ma'].sem():9.4f}")
print("\n  GC gradient within class (MA-normalised):")
print(f"  {'class':11s}{'Q5-Q1 EtOH':>13}{'Q5-Q1 With':>13}{'retained':>10}")
for c in ("Lasting","Acute","notDE"):
    s=X[X["class"]==c]
    qq=np.nanquantile(s["gc"],[i/5 for i in range(1,5)]); bb=np.digitize(s["gc"],qq)
    e=s[bb==4]["L2E_ma"].mean()-s[bb==0]["L2E_ma"].mean()
    w=s[bb==4]["L2W_ma"].mean()-s[bb==0]["L2W_ma"].mean()
    print(f"  {c:11s}{e:+13.4f}{w:+13.4f}{100*w/e:9.0f}%")
D.to_csv(f"{OUT}/ep_unbiased_MA.tsv",sep="\t",index=False)
print(f"\nwrote {OUT}/ep_unbiased_MA.tsv (L2E_ma / L2W_ma)")
