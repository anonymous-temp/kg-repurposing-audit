<!-- TABLE4 -->
**Table 4** Write-back policies in Hetionet under the graph head

| Policy | Negated pairs (recorded or approved) | RE held-out negated / masked | RE per-disease AP (Δ; 95% CI) | CD held-out negated / masked | CD per-disease AP (Δ; 95% CI) | E2 per-disease AP (Δ; 95% CI) | E3 AUROC |
|---|---|---|---|---|---|---|---|
| No write-back | 0 | 0.0 / 0.0 | 0.336 | 0.0 / 0.0 | 0.362 | 0.049 | 0.658 |
| Flat negative | 2,025 (30%) | 73.0 / 0.0 | 0.242 (−0.095; −0.125, −0.069) | 75.4 / 0.0 | 0.238 (−0.124; −0.147, −0.100) | 0.037 (−0.012; −0.029, +0.001) | 0.664 |
| Typed | 459 (39%) | 23.6 / 49.4 | 0.318 (−0.018; −0.043, +0.005) | 24.2 / 51.2 | 0.304 (−0.058; −0.083, −0.036) | 0.049 (0.000; −0.022, +0.017) | 0.665 |
| Typed + scope | 252 (39%) | 14.4 / 58.6 | 0.333 (−0.003; −0.027, +0.021) | 14.0 / 61.4 | 0.334 (−0.028; −0.042, −0.014) | 0.051 (+0.002; −0.021, +0.019) | 0.658 |
| Typed + role | 234 (36%) | 10.4 / 62.6 | 0.341 (+0.005; −0.019, +0.026) | 13.8 / 61.6 | 0.334 (−0.028; −0.044, −0.013) | 0.055 (+0.006; −0.008, +0.020) | 0.655 |
| Typed + role + scope | 120 (33%) | 5.2 / 67.8 | 0.347 (+0.011; −0.011, +0.032) | 7.8 / 67.6 | 0.352 (−0.010; −0.023, +0.002) | 0.057 (+0.008; −0.006, +0.022) | 0.655 |
| Typed + ≤2 trials | 77 (3%) | 0.6 / 72.4 | 0.356 (+0.019; −0.002, +0.040) | 0.2 / 75.2 | 0.369 (+0.007; 0.000, +0.015) | 0.056 (+0.007; −0.008, +0.021) | 0.655 |
| Mask all | 0 | 0.0 / 73.0 | 0.360 (+0.024; +0.008, +0.043) | 0.0 / 75.4 | 0.367 (+0.004; −0.002, +0.011) | 0.056 (+0.008; −0.008, +0.022) | 0.658 |

Negated pairs: all pairs that a policy writes back as explicit negatives, with the share that are recorded treatments or approved indications. Held-out treatments negated or masked are means per partition (151 random-edge test positives; 140–162 compound-disjoint test positives). Per-disease AP is followed by its difference from no write-back and the 95% paired disease-cluster bootstrap interval. Typed + ≤2 trials negates efficacy or safety stops only for pairs with at most two registered trials before 2015. E3 AUROC: external approved indications ranked above later scientific failures. RE, random-edge task; CD, compound-disjoint task; E2, external approved indications.

<!-- TABLE5 -->
**Table 5** Write-back policies in PrimeKG

| Policy | Negated pairs (recorded or approved) | RE held-out negated / masked | RE per-disease AP (Δ; 95% CI) | CD held-out negated / masked | CD per-disease AP (Δ; 95% CI) | E2 per-disease AP (Δ; 95% CI) | E3 AUROC |
|---|---|---|---|---|---|---|---|
| No write-back | 0 | 0.0 / 0.0 | 0.620 | 0.0 / 0.0 | 0.144 | 0.147 | 0.651 |
| Flat negative | 5,547 (22%) | 152.0 / 0.0 | 0.609 (−0.011; −0.016, −0.007) | 148.0 / 0.0 | 0.104 (−0.040; −0.047, −0.032) | 0.122 (−0.025; −0.037, −0.015) | 0.658 |
| Typed | 1,010 (28%) | 40.0 / 112.0 | 0.616 (−0.004; −0.007, −0.001) | 35.7 / 112.3 | 0.131 (−0.013; −0.018, −0.008) | 0.146 (0.000; −0.008, +0.009) | 0.662 |
| Typed + scope | 845 (30%) | 37.7 / 114.3 | 0.616 (−0.004; −0.007, −0.002) | 33.7 / 114.3 | 0.126 (−0.018; −0.024, −0.012) | 0.147 (+0.000; −0.009, +0.010) | 0.661 |
| Typed + role | 503 (23%) | 13.7 / 138.3 | 0.620 (+0.000; −0.002, +0.002) | 13.0 / 135.0 | 0.150 (+0.006; +0.001, +0.010) | 0.146 (−0.001; −0.010, +0.010) | 0.659 |
| Typed + role + scope | 408 (25%) | 13.3 / 138.7 | 0.620 (0.000; −0.002, +0.002) | 13.0 / 135.0 | 0.139 (−0.005; −0.011, +0.001) | 0.148 (+0.001; −0.008, +0.011) | 0.657 |
| Typed + ≤2 trials | 310 (7%) | 2.7 / 149.3 | 0.621 (+0.001; −0.002, +0.003) | 2.0 / 146.0 | 0.139 (−0.005; −0.009, 0.000) | 0.152 (+0.006; −0.003, +0.015) | 0.649 |
| Mask all | 0 | 0.0 / 152.0 | 0.621 (+0.001; 0.000, +0.003) | 0.0 / 148.0 | 0.151 (+0.007; +0.003, +0.011) | 0.154 (+0.007; 0.000, +0.017) | 0.647 |

Columns as in Table 4. Per-disease AP refers to the strongest scorer in each outcome set: label-only matrix factorisation in the random-edge task and for external approved indications (E2, E3), and the graph head in the compound-disjoint task. Held-out treatments are means over three partitions (1,878 random-edge test positives; 1,643–2,003 compound-disjoint test positives).
