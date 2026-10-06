"""Parse PrimeKG edges.csv from stdin into compact int32 arrays (CSV never stored)."""
import csv, os, sys, numpy as np
OUT = "work/data/primekg"
rel2i, X, Y, R, disp = {}, [], [], [], {}
rd = csv.reader(sys.stdin)
hdr = next(rd); print(hdr, flush=True)
ir, idr, ix, iy = hdr.index("relation"), hdr.index("display_relation"), hdr.index("x_index"), hdr.index("y_index")
for n, row in enumerate(rd):
    k = rel2i.setdefault(row[ir], len(rel2i)); disp[row[ir]] = row[idr]
    X.append(int(row[ix])); Y.append(int(row[iy])); R.append(k)
    if n % 1000000 == 0: print(n, flush=True)
rels = sorted(rel2i, key=rel2i.get)
np.savez(f"{OUT}/edges.tmp.npz", x=np.array(X, np.int32), y=np.array(Y, np.int32), r=np.array(R, np.int16),
         rels=np.array(rels), disp=np.array([disp[k] for k in rels]))
os.replace(f"{OUT}/edges.tmp.npz", f"{OUT}/edges.npz")
print("done", n + 1, flush=True)
