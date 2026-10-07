## S22. Placebo negatives and the testing-intensity policy

Every restricting policy can be written as masking all stopped pairs and negating a set S of them. For each S, a degree-matched placebo rewired S by double-edge swaps (30 attempted swaps per pair), so that every compound and disease received as many negatives as under the real policy, and a uniform placebo drew the same number of pairs at random. Placebo pairs were restricted to unlabelled pairs without a stopped trial that were not palliative or off-label pairs; fitting positives of the partition were excluded, held-out treatments were not. One placebo draw was made per partition and set (random seed derived from the partition seed). No original pair remained after rewiring, except that in 6 of the 60 degree-matched placebo sets of the Hetionet held-out tasks one pair could not be swapped out and was dropped. The testing-intensity policies negated a scientific stop only if the pair had at most one, two or four registered trials of any status that started before 2015, and masked all other stopped pairs; two was the pre-specified primary threshold. Scripts: pipelines/17_placebo_intensity_hetionet.py and pipelines/18_placebo_intensity_primekg.py.

**Table S22a.** Change in per-disease AP relative to masking when the negated set of a policy (stopped pairs) or a placebo set of the same size is added on top of masking. 95% paired disease-cluster bootstrap intervals (1,000 draws). Hetionet: means over five partitions (MF: and three initialisation seeds); PrimeKG: three partitions.

| Graph | Outcome | Scorer | Negated set | Negatives added | Stopped pairs: Δ vs masking (95% CI) | Placebo, same drug and disease counts: Δ vs masking (95% CI) | Placebo, uniform: Δ vs masking (95% CI) |
|---|---|---|---|---|---|---|---|
| Hetionet | Random edge | Graph head | Flat negative | 1752 | −0.119 (−0.152, −0.087) | −0.002 (−0.010, +0.006) | −0.001 (−0.005, +0.003) |
| Hetionet | Random edge | Graph head | Typed | 372 | −0.042 (−0.065, −0.021) | +0.001 (−0.003, +0.005) | +0.000 (−0.001, +0.001) |
| Hetionet | Random edge | Graph head | Typed + scope | 202 | −0.027 (−0.048, −0.007) | −0.003 (−0.009, +0.001) | 0.000 (−0.002, +0.001) |
| Hetionet | Random edge | Graph head | Typed + role | 191 | −0.019 (−0.035, −0.006) | −0.001 (−0.006, +0.003) | −0.001 (−0.002, +0.001) |
| Hetionet | Random edge | Graph head | Typed + role + scope | 97 | −0.013 (−0.026, −0.002) | 0.000 (−0.003, +0.002) | 0.000 (−0.001, +0.000) |
| Hetionet | Random edge | Graph head | Typed + ≤2 trials | 77 | −0.005 (−0.015, +0.001) | −0.001 (−0.002, +0.000) | −0.002 (−0.004, 0.000) |
| Hetionet | Random edge | Label-only MF | Flat negative | 1752 | −0.096 (−0.123, −0.069) | −0.012 (−0.023, +0.000) | +0.001 (−0.004, +0.007) |
| Hetionet | Random edge | Label-only MF | Typed | 372 | −0.038 (−0.057, −0.022) | −0.009 (−0.016, −0.003) | +0.002 (−0.004, +0.008) |
| Hetionet | Random edge | Label-only MF | Typed + scope | 202 | −0.031 (−0.046, −0.018) | +0.002 (−0.003, +0.008) | 0.000 (−0.001, +0.001) |
| Hetionet | Random edge | Label-only MF | Typed + role | 191 | −0.016 (−0.027, −0.006) | −0.006 (−0.012, −0.001) | 0.000 (−0.003, +0.002) |
| Hetionet | Random edge | Label-only MF | Typed + role + scope | 97 | −0.008 (−0.014, −0.003) | +0.002 (−0.002, +0.007) | −0.001 (−0.004, +0.000) |
| Hetionet | Random edge | Label-only MF | Typed + ≤2 trials | 77 | +0.002 (−0.003, +0.007) | −0.001 (−0.004, +0.002) | −0.001 (−0.003, +0.000) |
| Hetionet | Compound disjoint | Graph head | Flat negative | 1763 | −0.129 (−0.153, −0.103) | +0.013 (+0.002, +0.025) | +0.009 (+0.002, +0.017) |
| Hetionet | Compound disjoint | Graph head | Typed | 377 | −0.063 (−0.088, −0.039) | +0.001 (−0.007, +0.009) | −0.004 (−0.007, −0.001) |
| Hetionet | Compound disjoint | Graph head | Typed + scope | 203 | −0.033 (−0.046, −0.019) | +0.006 (+0.001, +0.013) | +0.001 (−0.004, +0.006) |
| Hetionet | Compound disjoint | Graph head | Typed + role | 194 | −0.033 (−0.049, −0.018) | −0.005 (−0.012, +0.001) | +0.000 (−0.001, +0.002) |
| Hetionet | Compound disjoint | Graph head | Typed + role + scope | 99 | −0.015 (−0.026, −0.005) | +0.000 (−0.002, +0.003) | 0.000 (−0.002, +0.000) |
| Hetionet | Compound disjoint | Graph head | Typed + ≤2 trials | 76 | +0.003 (0.000, +0.007) | −0.001 (−0.003, +0.000) | −0.001 (−0.003, +0.001) |
| Hetionet | Compound disjoint | Label-only MF | Flat negative | 1763 | −0.026 (−0.037, −0.015) | −0.020 (−0.031, −0.009) | +0.001 (−0.010, +0.011) |
| Hetionet | Compound disjoint | Label-only MF | Typed | 377 | −0.019 (−0.027, −0.011) | −0.016 (−0.025, −0.007) | +0.002 (−0.002, +0.006) |
| Hetionet | Compound disjoint | Label-only MF | Typed + scope | 203 | −0.011 (−0.017, −0.006) | −0.009 (−0.015, −0.003) | +0.000 (−0.004, +0.005) |
| Hetionet | Compound disjoint | Label-only MF | Typed + role | 194 | −0.016 (−0.024, −0.009) | −0.017 (−0.024, −0.010) | +0.001 (−0.003, +0.005) |
| Hetionet | Compound disjoint | Label-only MF | Typed + role + scope | 99 | −0.007 (−0.012, −0.002) | −0.005 (−0.010, +0.000) | +0.003 (+0.001, +0.005) |
| Hetionet | Compound disjoint | Label-only MF | Typed + ≤2 trials | 76 | −0.005 (−0.009, −0.002) | −0.005 (−0.009, −0.001) | +0.001 (−0.001, +0.003) |
| Hetionet | External approvals | Graph head | Flat negative | 1649 | −0.019 (−0.043, +0.003) | −0.001 (−0.003, 0.000) | −0.001 (−0.002, +0.000) |
| Hetionet | External approvals | Graph head | Typed | 340 | −0.008 (−0.024, +0.001) | −0.001 (−0.002, 0.000) | −0.001 (−0.002, +0.000) |
| Hetionet | External approvals | Graph head | Typed + scope | 183 | −0.005 (−0.021, +0.003) | −0.001 (−0.002, 0.000) | +0.000 (0.000, +0.000) |
| Hetionet | External approvals | Graph head | Typed + role | 176 | −0.001 (−0.004, +0.001) | 0.000 (−0.001, +0.001) | 0.000 (−0.001, +0.000) |
| Hetionet | External approvals | Graph head | Typed + role + scope | 89 | +0.001 (0.000, +0.002) | 0.000 (0.000, 0.000) | +0.001 (0.000, +0.003) |
| Hetionet | External approvals | Graph head | Typed + ≤2 trials | 76 | 0.000 (−0.001, +0.000) | +0.000 (−0.001, +0.001) | 0.000 (0.000, +0.000) |
| Hetionet | External approvals | Label-only MF | Flat negative | 1649 | −0.010 (−0.016, −0.004) | −0.003 (−0.005, +0.000) | −0.001 (−0.002, +0.001) |
| Hetionet | External approvals | Label-only MF | Typed | 340 | −0.003 (−0.005, +0.000) | −0.001 (−0.004, +0.002) | −0.002 (−0.003, −0.001) |
| Hetionet | External approvals | Label-only MF | Typed + scope | 183 | +0.000 (−0.002, +0.003) | +0.000 (−0.001, +0.002) | −0.001 (−0.002, +0.000) |
| Hetionet | External approvals | Label-only MF | Typed + role | 176 | −0.002 (−0.004, +0.000) | −0.001 (−0.002, +0.000) | 0.000 (0.000, +0.000) |
| Hetionet | External approvals | Label-only MF | Typed + role + scope | 89 | +0.001 (−0.001, +0.003) | 0.000 (−0.001, +0.000) | 0.000 (−0.001, +0.000) |
| Hetionet | External approvals | Label-only MF | Typed + ≤2 trials | 76 | 0.000 (−0.001, +0.000) | +0.000 (0.000, +0.001) | +0.001 (0.000, +0.002) |
| PrimeKG | Random edge | Graph head | Flat negative | 3665 | −0.006 (−0.011, −0.002) | −0.002 (−0.006, +0.001) | −0.001 (−0.004, +0.003) |
| PrimeKG | Random edge | Graph head | Typed | 686 | −0.003 (−0.007, +0.000) | −0.001 (−0.004, +0.003) | +0.001 (−0.001, +0.003) |
| PrimeKG | Random edge | Graph head | Typed + role + scope | 273 | +0.000 (−0.002, +0.002) | −0.001 (−0.003, +0.001) | −0.012 (−0.018, −0.008) |
| PrimeKG | Random edge | Graph head | Typed + ≤2 trials | 204 | −0.002 (−0.004, +0.001) | +0.001 (−0.001, +0.004) | 0.000 (−0.002, +0.001) |
| PrimeKG | Random edge | Label-only MF | Flat negative | 3665 | −0.013 (−0.017, −0.008) | −0.001 (−0.003, +0.001) | +0.001 (−0.001, +0.003) |
| PrimeKG | Random edge | Label-only MF | Typed | 686 | −0.005 (−0.008, −0.003) | +0.000 (−0.001, +0.002) | 0.000 (−0.001, +0.001) |
| PrimeKG | Random edge | Label-only MF | Typed + role + scope | 273 | −0.001 (−0.003, +0.000) | +0.001 (0.000, +0.002) | 0.000 (−0.002, +0.001) |
| PrimeKG | Random edge | Label-only MF | Typed + ≤2 trials | 204 | −0.001 (−0.003, +0.001) | +0.000 (−0.001, +0.001) | +0.000 (−0.001, +0.001) |
| PrimeKG | Compound disjoint | Graph head | Flat negative | 3665 | −0.047 (−0.055, −0.039) | −0.036 (−0.043, −0.029) | +0.009 (+0.005, +0.014) |
| PrimeKG | Compound disjoint | Graph head | Typed | 681 | −0.020 (−0.024, −0.016) | −0.018 (−0.022, −0.014) | −0.019 (−0.023, −0.015) |
| PrimeKG | Compound disjoint | Graph head | Typed + role + scope | 272 | −0.012 (−0.017, −0.007) | +0.005 (+0.002, +0.008) | −0.004 (−0.005, −0.002) |
| PrimeKG | Compound disjoint | Graph head | Typed + ≤2 trials | 204 | −0.012 (−0.014, −0.009) | −0.002 (−0.003, +0.000) | −0.004 (−0.006, −0.003) |
| PrimeKG | Compound disjoint | Label-only MF | Flat negative | 3665 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) | +0.006 (+0.006, +0.006) |
| PrimeKG | Compound disjoint | Label-only MF | Typed | 681 | 0.000 (0.000, +0.000) | 0.000 (0.000, +0.000) | +0.001 (+0.001, +0.001) |
| PrimeKG | Compound disjoint | Label-only MF | Typed + role + scope | 272 | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Compound disjoint | Label-only MF | Typed + ≤2 trials | 204 | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | External approvals | Graph head | Flat negative | 3455 | −0.002 (−0.005, +0.003) | +0.000 (−0.003, +0.004) | 0.000 (−0.006, +0.006) |
| PrimeKG | External approvals | Graph head | Typed | 631 | −0.001 (−0.002, 0.000) | +0.001 (−0.001, +0.006) | +0.005 (+0.001, +0.010) |
| PrimeKG | External approvals | Graph head | Typed + role + scope | 255 | +0.003 (+0.001, +0.004) | −0.001 (−0.004, +0.004) | −0.006 (−0.012, +0.001) |
| PrimeKG | External approvals | Graph head | Typed + ≤2 trials | 200 | 0.000 (−0.001, 0.000) | −0.003 (−0.009, +0.002) | 0.000 (0.000, +0.000) |
| PrimeKG | External approvals | Label-only MF | Flat negative | 3455 | −0.032 (−0.047, −0.021) | −0.002 (−0.008, +0.004) | +0.001 (−0.004, +0.005) |
| PrimeKG | External approvals | Label-only MF | Typed | 631 | −0.008 (−0.013, −0.003) | −0.001 (−0.005, +0.002) | −0.001 (−0.004, +0.001) |
| PrimeKG | External approvals | Label-only MF | Typed + role + scope | 255 | −0.006 (−0.013, −0.002) | −0.001 (−0.003, +0.001) | +0.000 (−0.003, +0.004) |
| PrimeKG | External approvals | Label-only MF | Typed + ≤2 trials | 200 | −0.002 (−0.006, +0.001) | +0.000 (−0.002, +0.002) | 0.000 (−0.001, +0.001) |

**Table S22b.** Testing-intensity policies for all scorers: per-disease AP and change relative to no write-back and to masking.

| Graph | Policy | Negated pairs (recorded or approved) | Outcome | Scorer | Per-disease AP | Δ vs no write-back (95% CI) | Δ vs masking (95% CI) |
|---|---|---|---|---|---|---|---|
| Hetionet | Typed + ≤1 trial | 58 (2%) | Random edge | Degree reference | 0.132 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | Random edge | Label-only MF | 0.337 | +0.019 (+0.009, +0.029) | −0.001 (−0.005, +0.002) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | Random edge | Graph head | 0.355 | +0.019 (−0.002, +0.040) | −0.006 (−0.016, +0.001) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | Random edge | Hybrid | 0.357 | +0.018 (−0.003, +0.038) | −0.005 (−0.014, +0.000) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | Compound disjoint | Degree reference | 0.036 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | Compound disjoint | Label-only MF | 0.084 | −0.002 (−0.006, +0.001) | −0.003 (−0.006, +0.000) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | Compound disjoint | Graph head | 0.367 | +0.004 (−0.002, +0.010) | 0.000 (−0.001, +0.001) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | Compound disjoint | Hybrid | 0.367 | +0.004 (−0.002, +0.010) | 0.000 (−0.001, +0.001) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | External approvals | Degree reference | 0.020 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | External approvals | Label-only MF | 0.040 | +0.004 (+0.001, +0.006) | 0.000 (−0.001, +0.000) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | External approvals | Graph head | 0.057 | +0.008 (−0.008, +0.022) | +0.000 (0.000, +0.000) |
| Hetionet | Typed + ≤1 trial | 58 (2%) | External approvals | Hybrid | 0.057 | +0.011 (+0.001, +0.024) | 0.000 (0.000, +0.000) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Random edge | Degree reference | 0.132 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Random edge | Label-only MF | 0.340 | +0.023 (+0.012, +0.033) | +0.002 (−0.003, +0.007) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Random edge | Graph head | 0.356 | +0.019 (−0.002, +0.040) | −0.005 (−0.015, +0.001) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Random edge | Hybrid | 0.357 | +0.019 (−0.002, +0.039) | −0.004 (−0.015, +0.002) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Compound disjoint | Degree reference | 0.036 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Compound disjoint | Label-only MF | 0.082 | −0.005 (−0.008, −0.001) | −0.005 (−0.009, −0.002) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Compound disjoint | Graph head | 0.369 | +0.007 (0.000, +0.015) | +0.003 (0.000, +0.007) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | Compound disjoint | Hybrid | 0.369 | +0.007 (−0.001, +0.015) | +0.003 (0.000, +0.007) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | External approvals | Degree reference | 0.020 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | External approvals | Label-only MF | 0.040 | +0.004 (+0.001, +0.006) | 0.000 (−0.001, +0.000) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | External approvals | Graph head | 0.056 | +0.007 (−0.008, +0.021) | 0.000 (−0.001, +0.000) |
| Hetionet | Typed + ≤2 trials | 77 (3%) | External approvals | Hybrid | 0.057 | +0.011 (+0.001, +0.024) | 0.000 (0.000, +0.000) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Random edge | Degree reference | 0.132 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Random edge | Label-only MF | 0.339 | +0.022 (+0.008, +0.035) | +0.001 (−0.009, +0.011) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Random edge | Graph head | 0.356 | +0.020 (−0.002, +0.042) | −0.004 (−0.015, +0.002) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Random edge | Hybrid | 0.357 | +0.018 (−0.002, +0.038) | −0.004 (−0.015, +0.002) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Compound disjoint | Degree reference | 0.036 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Compound disjoint | Label-only MF | 0.078 | −0.008 (−0.012, −0.004) | −0.008 (−0.013, −0.005) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Compound disjoint | Graph head | 0.369 | +0.007 (−0.002, +0.015) | +0.002 (−0.003, +0.008) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | Compound disjoint | Hybrid | 0.369 | +0.007 (−0.002, +0.015) | +0.002 (−0.003, +0.008) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | External approvals | Degree reference | 0.020 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | External approvals | Label-only MF | 0.040 | +0.004 (+0.001, +0.007) | 0.000 (−0.001, +0.001) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | External approvals | Graph head | 0.056 | +0.007 (−0.009, +0.021) | −0.001 (−0.002, +0.000) |
| Hetionet | Typed + ≤4 trials | 126 (5%) | External approvals | Hybrid | 0.057 | +0.010 (+0.001, +0.023) | −0.001 (−0.003, +0.000) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | Random edge | Degree reference | 0.114 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | Random edge | Label-only MF | 0.620 | 0.000 (−0.002, +0.002) | −0.001 (−0.003, +0.000) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | Random edge | Graph head | 0.122 | −0.007 (−0.012, −0.003) | −0.010 (−0.015, −0.005) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | Compound disjoint | Degree reference | 0.007 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | Compound disjoint | Label-only MF | 0.007 | +0.000 (0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | Compound disjoint | Graph head | 0.147 | +0.003 (−0.001, +0.007) | −0.004 (−0.006, −0.002) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | External approvals | Degree reference | 0.063 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | External approvals | Label-only MF | 0.154 | +0.007 (−0.001, +0.017) | 0.000 (−0.002, +0.001) |
| PrimeKG | Typed + ≤1 trial | 203 (5%) | External approvals | Graph head | 0.045 | −0.005 (−0.011, +0.001) | −0.005 (−0.011, +0.000) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | Random edge | Degree reference | 0.114 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | Random edge | Label-only MF | 0.621 | +0.001 (−0.002, +0.003) | −0.001 (−0.003, +0.001) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | Random edge | Graph head | 0.130 | +0.001 (−0.002, +0.004) | −0.002 (−0.004, +0.001) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | Compound disjoint | Degree reference | 0.007 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | Compound disjoint | Label-only MF | 0.007 | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | Compound disjoint | Graph head | 0.139 | −0.005 (−0.009, 0.000) | −0.012 (−0.014, −0.009) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | External approvals | Degree reference | 0.063 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | External approvals | Label-only MF | 0.152 | +0.006 (−0.003, +0.015) | −0.002 (−0.006, +0.001) |
| PrimeKG | Typed + ≤2 trials | 310 (7%) | External approvals | Graph head | 0.050 | 0.000 (−0.001, +0.001) | 0.000 (−0.001, 0.000) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | Random edge | Degree reference | 0.114 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | Random edge | Label-only MF | 0.622 | +0.002 (−0.001, +0.004) | +0.000 (−0.002, +0.002) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | Random edge | Graph head | 0.128 | −0.002 (−0.005, +0.001) | −0.004 (−0.007, −0.001) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | Compound disjoint | Degree reference | 0.007 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | Compound disjoint | Label-only MF | 0.007 | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | Compound disjoint | Graph head | 0.149 | +0.005 (+0.002, +0.009) | −0.002 (−0.004, +0.000) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | External approvals | Degree reference | 0.063 | +0.000 (+0.000, +0.000) | +0.000 (+0.000, +0.000) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | External approvals | Label-only MF | 0.152 | +0.005 (−0.003, +0.015) | −0.002 (−0.007, +0.001) |
| PrimeKG | Typed + ≤4 trials | 430 (8%) | External approvals | Graph head | 0.052 | +0.001 (−0.001, +0.006) | +0.001 (−0.001, +0.005) |
