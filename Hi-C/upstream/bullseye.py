"""Faithful port of SIPMeta's bullseye transform (Rowley lab bullseye.py).
Groups matrix cells by Manhattan distance from the centre and lays each ring
uniformly around the polar angle, which removes the Manhattan-distance artefact
that makes square APA plots look diamond-shaped."""
import numpy as np, math
from bisect import bisect_left
def _closest(lst,x):
    pos=bisect_left(lst,x)
    if pos==0: return 0
    if pos==len(lst): return -1
    return pos if lst[pos]-x < x-lst[pos-1] else pos-1
def bullseye(mat,uniform=True):
    mat=np.asarray(mat,float); n=len(mat)
    if n%2==0: raise ValueError("bullseye needs an odd-sized matrix (central bin)")
    c=(n-1)/2
    rings=math.floor(math.sqrt(((n//2)**2)/2))
    maxrings=math.ceil(math.sqrt(rings**2+rings**2))
    slices=(n-1)*100
    md=np.zeros((n,n)); at=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            md[i,j]=abs(j-int(c))+abs(i-c)
            at[i,j]=np.arctan2(j-int(c), i-c)
    m_at=[[] for _ in range(n)]; m_sc=[[] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for d in range(n):
                if md[b,d]==a: m_sc[a].append(mat[b,d]); m_at[a].append(at[b,d])
    rad=np.linspace(0,maxrings,maxrings+1); azm=np.linspace(0,2*np.pi,slices)
    r,th=np.meshgrid(rad,azm)
    C=np.zeros([slices,int(maxrings)+1])
    for ring in range(int(maxrings)+1):
        if ring>=len(m_at) or not m_at[ring]:
            C[:,ring]=C[:,ring-1] if ring else 0; continue
        thetas=[(x+np.pi)*((slices/2)/np.pi) for x in m_at[ring]]
        scores=list(m_sc[ring])
        thetas.append(0); scores.append(scores[0])
        sc_sorted=[x for _,x in sorted(zip(thetas,scores))]
        th_sorted=sorted(thetas)
        if uniform and ring:
            sps=slices/(len(scores)-1)
            th_sorted=[sps*t for t in range(len(th_sorted))]
        for i in range(slices):
            C[i,ring]=sc_sorted[_closest(th_sorted,i)]
    return th,r,C,rings
