#!/usr/bin/env python
"""Additional file 1, Section S22: placebo negatives and testing-intensity policies (from make_tables_v5.py output)."""
import re, sys


def fix_minus(t):
    """Typographic minus for negative numbers; a rounded zero carries no sign."""
    t = re.sub(r"(?<![\w.])[-−]0\.000\b", "0.000", t)
    return re.sub(r"(?<![\w.])-(?=\d)", "−", t)


T = sys.argv[1] if len(sys.argv) > 1 else "work/tables_v5"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/manuscript/supp_S22.md"
text = ["## S22. Placebo negatives and the testing-intensity policy", "",
        "Every restricting policy can be written as masking all stopped pairs and negating a set S of them. For each S, a degree-matched placebo "
        "rewired S by double-edge swaps (30 attempted swaps per pair), so that every compound and disease received as many negatives as under the real "
        "policy, and a uniform placebo drew the same number of pairs at random. Placebo pairs were restricted to unlabelled pairs without a stopped trial "
        "that were not palliative or off-label pairs; fitting positives of the partition were excluded, held-out treatments were not. One placebo draw was "
        "made per partition and set (random seed derived from the partition seed). No original pair remained after rewiring, except that in 6 of the 60 "
        "degree-matched placebo sets of the Hetionet held-out tasks one pair could not be swapped out and was dropped. The testing-intensity policies negated a scientific stop only if the "
        "pair had at most one, two or four registered trials of any status that started before 2015, and masked all other stopped pairs; two was the "
        "pre-specified primary threshold. Scripts: pipelines/17_placebo_intensity_hetionet.py and pipelines/18_placebo_intensity_primekg.py.", "",
        "**Table S22a.** Change in per-disease AP relative to masking when the negated set of a policy (stopped pairs) or a placebo set of the same size is "
        "added on top of masking. 95% paired disease-cluster bootstrap intervals (1,000 draws). Hetionet: means over five partitions (MF: and three "
        "initialisation seeds); PrimeKG: three partitions.", "", fix_minus(open(f"{T}/S22a_placebo.md").read()), "",
        "**Table S22b.** Testing-intensity policies for all scorers: per-disease AP and change relative to no write-back and to masking.", "",
        fix_minus(open(f"{T}/S22b_intensity.md").read()), ""]
open(OUT, "w").write("\n".join(text))
print("written", OUT)
