#!/usr/bin/env python
"""Compare numeric tokens, citation keys, placeholders and figure/table references between two versions of each file."""
import collections, re, sys
from pathlib import Path
old_dir, new_dir = Path(sys.argv[1]), Path(sys.argv[2])
NUM = re.compile(r"[−\-+]?\d[\d,]*(?:\.\d+)?%?")
CIT = re.compile(r"\{\{[^}]+\}\}")
REF = re.compile(r"(?:Figs?\.\s*\d+[A-F–\-]*(?:\s*and\s*\d+)?|Table\s+S?\d+[a-z]?|Section[s]?\s+S\d+(?:–S\d+)?|PLACEHOLDER|PKG_[A-Z_]+|TABLE\d_PLACEHOLDER)")
ok = True
for f in sorted(old_dir.glob("*.md")):
    a, b = f.read_text(), (new_dir / f.name).read_text()
    for name, rx in [("numbers", NUM), ("citations", CIT), ("refs", REF)]:
        ca, cb = collections.Counter(rx.findall(a)), collections.Counter(rx.findall(b))
        if ca != cb:
            ok = False
            print(f.name, name, "removed:", dict(ca - cb), "added:", dict(cb - ca))
    ha = [l for l in a.splitlines() if l.startswith("#")]; hb = [l for l in b.splitlines() if l.startswith("#")]
    if ha != hb:
        ok = False; print(f.name, "headings changed")
print("ALL PRESERVED" if ok else "DIFFERENCES FOUND")
