import os,sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

# Ensure PDF renders text as TrueType (actual text boxes)
matplotlib.rcParams['pdf.fonttype'] = 42

PROJ="/Zulu/jordan/alcoholATAC"; L=f"{PROJ}/bensalcohol_out/loops"
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"svg")
os.makedirs(OUT,exist_ok=True)
BLU="#5a72ad"; RED="#ab5054"; PRP="#7d2c88"; GRY="#999999"; ORA="#E1812C"
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.svg",bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    print(f"  -> svg/{name}.svg")

"""Fig 5B - A-compartment self-association shifts toward G-rich.
Input: loops/panels.pkl  (key 'AA_by_gc'), produced by upstream/panel_rest.py"""
import pickle
P=pickle.load(open(f"{L}/panels.pkl","rb"))
x=np.arange(1,6)
fig,ax=plt.subplots(figsize=(5.6,4.6))
for c,lab,col in [("ctrl","Ctl",BLU),("etoh","EtOH",RED),("with","Withdraw",PRP)]:
    ax.plot(x,P["AA_by_gc"][c],"-o",color=col,label=lab,lw=2)
d=100*(np.array(P["AA_by_gc"]["etoh"])/np.array(P["AA_by_gc"]["ctrl"])-1)
yl=ax.get_ylim(); ax.set_ylim(yl[0],yl[1]+0.10*(yl[1]-yl[0]))
for xi,di in zip(x,d):
    ax.annotate(f"{di:+.1f}%",xy=(xi,P["AA_by_gc"]["etoh"][xi-1]),xytext=(0,10),
                textcoords="offset points",ha="center",fontsize=9,
                color=RED if di>0 else BLU,fontweight="bold")
ax.set_xticks(x); ax.set_xlabel("GC quintile of A-compartment bins",fontsize=9)
ax.set_ylabel("mean O/E, A-A contacts (50-500 kb)",fontsize=9)
ax.set_title("A-compartment self-association shifts toward G-rich",fontsize=10)
ax.legend(fontsize=9)
#save(fig,"Fig5B_compartment_by_GC")
# Save as pdf
plt.savefig("Fig5B_compartment_by_GC.pdf")

