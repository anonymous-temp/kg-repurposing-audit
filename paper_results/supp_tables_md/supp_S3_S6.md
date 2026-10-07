## S3. Hyper-parameter selection

Hyper-parameters were selected on mean validation AP over the five partitions of each task, without write-back. Grid: MF dimension 16 or 32; weight decay 10⁻⁵, 10⁻⁴ or 10⁻³; embedding checkpoint after 5, 10, 15 or 20 pretraining epochs (graph head); 200, 400, 800 or 1,600 training epochs. The hybrid used the selected MF dimension and checkpoint and re-selected weight decay and epochs. The full grid (840 configuration–epoch–partition evaluations) is in paper_results/selection_grid.tsv.

Table S3. Selected configurations and their mean validation AP.

| Task | Scorer | MF dimension | Embedding checkpoint (epochs) | Weight decay | Training epochs | Validation AP |
|---|---|---|---|---|---|---|
| Random edge | Label-only MF | 32 | — | 1e-04 | 1600 | 0.104 |
| Random edge | Graph head | — | 15 | 1e-05 | 1600 | 0.211 |
| Random edge | Hybrid | 32 | 15 | 1e-04 | 1600 | 0.205 |
| Compound disjoint | Label-only MF | 32 | — | 1e-05 | 200 | 0.118 |
| Compound disjoint | Graph head | — | 15 | 1e-03 | 1600 | 0.340 |
| Compound disjoint | Hybrid | 32 | 15 | 1e-03 | 1600 | 0.340 |

For the random-edge graph head with the selected checkpoint and weight decay, mean validation AP was 0.166, 0.184, 0.194 and 0.211 at 200, 400, 800 and 1,600 epochs. The largest budget was selected, and performance may increase with longer training. The graph head is initialised at zero and is deterministic; MF and the hybrid were fitted with three initialisation seeds.

## S4. Full results for held-out treatments and external approved indications

Values are means (SD over partitions) for E1 and means over initialisation seeds for E2. 'Sci X / other Y' gives the action applied to pairs with a scientific stop and to other stopped pairs. 'Typed + scope' negates only same-concept scientific stops and masks all other stopped pairs. Negative weight is the weight of explicit negatives relative to unlabelled pairs.

Table S4a. Random-edge task (E1).

| Policy | Negative weight | Scorer | Pooled AP | Per-disease AP | MRR | Hits@10 | AUROC | Held-out negated | Held-out masked |
|---|---|---|---|---|---|---|---|---|---|
| No write-back | 10.0 | Degree reference | 0.042 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.688 (0.008) | +0.000 | +0.000 |
| No write-back | 10.0 | Disease degree | 0.054 (0.003) | 0.027 (0.002) | 0.191 (0.013) | 0.404 (0.034) | 0.725 (0.006) | +0.000 | +0.000 |
| No write-back | 10.0 | Graph head | 0.257 (0.025) | 0.336 (0.019) | 0.427 (0.026) | 0.662 (0.022) | 0.847 (0.018) | +0.000 | +0.000 |
| No write-back | 10.0 | Hybrid | 0.253 (0.027) | 0.339 (0.026) | 0.437 (0.031) | 0.680 (0.023) | 0.849 (0.018) | +0.000 | +0.000 |
| No write-back | 10.0 | MF | 0.178 (0.019) | 0.317 (0.014) | 0.387 (0.028) | 0.691 (0.023) | 0.851 (0.005) | +0.000 | +0.000 |
| Sci ignore / other mask | 10.0 | Degree reference | 0.042 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.687 (0.009) | +0.000 | +49.400 |
| Sci ignore / other mask | 10.0 | Graph head | 0.251 (0.046) | 0.342 (0.019) | 0.445 (0.037) | 0.687 (0.025) | 0.849 (0.019) | +0.000 | +49.400 |
| Sci ignore / other mask | 10.0 | Hybrid | 0.249 (0.040) | 0.349 (0.023) | 0.443 (0.040) | 0.681 (0.025) | 0.852 (0.020) | +0.000 | +49.400 |
| Sci ignore / other mask | 10.0 | MF | 0.185 (0.023) | 0.329 (0.017) | 0.396 (0.021) | 0.708 (0.032) | 0.856 (0.006) | +0.000 | +49.400 |
| Sci ignore / other negate | 10.0 | Degree reference | 0.043 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.692 (0.009) | +49.400 | +0.000 |
| Sci ignore / other negate | 10.0 | Graph head | 0.176 (0.015) | 0.274 (0.026) | 0.378 (0.022) | 0.652 (0.024) | 0.837 (0.016) | +49.400 | +0.000 |
| Sci ignore / other negate | 10.0 | Hybrid | 0.169 (0.013) | 0.274 (0.022) | 0.381 (0.020) | 0.651 (0.013) | 0.839 (0.016) | +49.400 | +0.000 |
| Sci ignore / other negate | 10.0 | MF | 0.133 (0.013) | 0.261 (0.018) | 0.306 (0.017) | 0.614 (0.027) | 0.819 (0.009) | +49.400 | +0.000 |
| Sci mask / other ignore | 10.0 | Degree reference | 0.042 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.687 (0.009) | +0.000 | +23.600 |
| Sci mask / other ignore | 10.0 | Graph head | 0.278 (0.026) | 0.349 (0.017) | 0.446 (0.016) | 0.689 (0.020) | 0.849 (0.018) | +0.000 | +23.600 |
| Sci mask / other ignore | 10.0 | Hybrid | 0.273 (0.026) | 0.354 (0.017) | 0.459 (0.028) | 0.688 (0.018) | 0.852 (0.018) | +0.000 | +23.600 |
| Sci mask / other ignore | 10.0 | MF | 0.187 (0.019) | 0.327 (0.013) | 0.399 (0.019) | 0.693 (0.023) | 0.852 (0.006) | +0.000 | +23.600 |
| Mask all | 10.0 | Degree reference | 0.042 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.687 (0.009) | +0.000 | +73.000 |
| Mask all | 10.0 | Graph head | 0.264 (0.039) | 0.360 (0.024) | 0.463 (0.030) | 0.698 (0.031) | 0.853 (0.019) | +0.000 | +73.000 |
| Mask all | 10.0 | Hybrid | 0.265 (0.038) | 0.361 (0.025) | 0.469 (0.034) | 0.689 (0.023) | 0.856 (0.019) | +0.000 | +73.000 |
| Mask all | 10.0 | MF | 0.190 (0.022) | 0.338 (0.016) | 0.406 (0.016) | 0.704 (0.038) | 0.857 (0.006) | +0.000 | +73.000 |
| Sci mask / other negate | 10.0 | Degree reference | 0.043 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.692 (0.009) | +49.400 | +23.600 |
| Sci mask / other negate | 10.0 | Graph head | 0.202 (0.022) | 0.290 (0.024) | 0.397 (0.018) | 0.654 (0.026) | 0.839 (0.016) | +49.400 | +23.600 |
| Sci mask / other negate | 10.0 | Hybrid | 0.194 (0.018) | 0.288 (0.021) | 0.400 (0.016) | 0.660 (0.012) | 0.840 (0.016) | +49.400 | +23.600 |
| Sci mask / other negate | 10.0 | MF | 0.143 (0.018) | 0.268 (0.019) | 0.312 (0.015) | 0.618 (0.032) | 0.822 (0.009) | +49.400 | +23.600 |
| Sci negate / other ignore | 10.0 | Degree reference | 0.043 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.689 (0.009) | +23.600 | +0.000 |
| Sci negate / other ignore | 10.0 | Graph head | 0.232 (0.032) | 0.304 (0.016) | 0.402 (0.030) | 0.646 (0.031) | 0.841 (0.018) | +23.600 | +0.000 |
| Sci negate / other ignore | 10.0 | Hybrid | 0.226 (0.033) | 0.304 (0.025) | 0.404 (0.031) | 0.654 (0.032) | 0.843 (0.018) | +23.600 | +0.000 |
| Sci negate / other ignore | 10.0 | MF | 0.152 (0.017) | 0.288 (0.013) | 0.351 (0.030) | 0.655 (0.027) | 0.833 (0.007) | +23.600 | +0.000 |
| Typed | 3.0 | Degree reference | 0.042 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.687 (0.009) | +23.600 | +49.400 |
| Typed | 3.0 | Graph head | 0.238 (0.046) | 0.327 (0.022) | 0.431 (0.038) | 0.668 (0.030) | 0.847 (0.019) | +23.600 | +49.400 |
| Typed | 3.0 | Hybrid | 0.237 (0.044) | 0.332 (0.022) | 0.428 (0.041) | 0.669 (0.028) | 0.849 (0.019) | +23.600 | +49.400 |
| Typed | 3.0 | MF | 0.177 (0.023) | 0.314 (0.015) | 0.384 (0.024) | 0.702 (0.030) | 0.853 (0.006) | +23.600 | +49.400 |
| Typed | 10.0 | Degree reference | 0.043 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.689 (0.009) | +23.600 | +49.400 |
| Typed | 10.0 | Graph head | 0.227 (0.047) | 0.318 (0.024) | 0.418 (0.038) | 0.654 (0.022) | 0.843 (0.019) | +23.600 | +49.400 |
| Typed | 10.0 | Hybrid | 0.222 (0.043) | 0.319 (0.024) | 0.412 (0.040) | 0.657 (0.030) | 0.845 (0.019) | +23.600 | +49.400 |
| Typed | 10.0 | MF | 0.162 (0.023) | 0.299 (0.016) | 0.362 (0.032) | 0.669 (0.034) | 0.840 (0.005) | +23.600 | +49.400 |
| Typed | 30.0 | Degree reference | 0.043 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.691 (0.009) | +23.600 | +49.400 |
| Typed | 30.0 | Graph head | 0.221 (0.050) | 0.314 (0.026) | 0.405 (0.041) | 0.637 (0.025) | 0.839 (0.020) | +23.600 | +49.400 |
| Typed | 30.0 | Hybrid | 0.216 (0.045) | 0.313 (0.025) | 0.402 (0.041) | 0.642 (0.030) | 0.841 (0.020) | +23.600 | +49.400 |
| Typed | 30.0 | MF | 0.151 (0.023) | 0.290 (0.009) | 0.345 (0.037) | 0.624 (0.028) | 0.820 (0.009) | +23.600 | +49.400 |
| Flat negative | 3.0 | Degree reference | 0.043 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.689 (0.009) | +73.000 | +0.000 |
| Flat negative | 3.0 | Graph head | 0.207 (0.023) | 0.282 (0.015) | 0.394 (0.031) | 0.648 (0.028) | 0.841 (0.016) | +73.000 | +0.000 |
| Flat negative | 3.0 | Hybrid | 0.205 (0.023) | 0.288 (0.020) | 0.396 (0.025) | 0.658 (0.023) | 0.843 (0.017) | +73.000 | +0.000 |
| Flat negative | 3.0 | MF | 0.151 (0.011) | 0.285 (0.015) | 0.350 (0.021) | 0.656 (0.023) | 0.837 (0.007) | +73.000 | +0.000 |
| Flat negative | 10.0 | Degree reference | 0.043 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.693 (0.009) | +73.000 | +0.000 |
| Flat negative | 10.0 | Graph head | 0.154 (0.020) | 0.242 (0.023) | 0.352 (0.031) | 0.625 (0.035) | 0.833 (0.016) | +73.000 | +0.000 |
| Flat negative | 10.0 | Hybrid | 0.145 (0.018) | 0.241 (0.017) | 0.346 (0.032) | 0.623 (0.023) | 0.833 (0.016) | +73.000 | +0.000 |
| Flat negative | 10.0 | MF | 0.112 (0.005) | 0.241 (0.022) | 0.278 (0.027) | 0.585 (0.046) | 0.803 (0.012) | +73.000 | +0.000 |
| Flat negative | 30.0 | Degree reference | 0.044 (0.007) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.701 (0.008) | +73.000 | +0.000 |
| Flat negative | 30.0 | Graph head | 0.130 (0.018) | 0.221 (0.022) | 0.320 (0.034) | 0.596 (0.044) | 0.824 (0.016) | +73.000 | +0.000 |
| Flat negative | 30.0 | Hybrid | 0.119 (0.014) | 0.216 (0.017) | 0.319 (0.034) | 0.597 (0.032) | 0.823 (0.016) | +73.000 | +0.000 |
| Flat negative | 30.0 | MF | 0.083 (0.004) | 0.200 (0.022) | 0.225 (0.030) | 0.460 (0.038) | 0.742 (0.018) | +73.000 | +0.000 |
| Typed + scope | 10.0 | Degree reference | 0.042 (0.006) | 0.132 (0.017) | 0.191 (0.013) | 0.404 (0.034) | 0.688 (0.008) | +14.400 | +58.600 |
| Typed + scope | 10.0 | Graph head | 0.234 (0.038) | 0.333 (0.028) | 0.429 (0.031) | 0.662 (0.019) | 0.847 (0.018) | +14.400 | +58.600 |
| Typed + scope | 10.0 | Hybrid | 0.234 (0.037) | 0.336 (0.022) | 0.429 (0.034) | 0.663 (0.024) | 0.850 (0.018) | +14.400 | +58.600 |
| Typed + scope | 10.0 | MF | 0.170 (0.020) | 0.307 (0.021) | 0.377 (0.024) | 0.669 (0.041) | 0.845 (0.005) | +14.400 | +58.600 |

Table S4b. Compound-disjoint task (E1).

| Policy | Negative weight | Scorer | Pooled AP | Per-disease AP | MRR | Hits@10 | AUROC | Held-out negated | Held-out masked |
|---|---|---|---|---|---|---|---|---|---|
| No write-back | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +0.000 | +0.000 |
| No write-back | 10.0 | Disease degree | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +0.000 | +0.000 |
| No write-back | 10.0 | Graph head | 0.276 (0.028) | 0.362 (0.075) | 0.479 (0.039) | 0.690 (0.054) | 0.851 (0.021) | +0.000 | +0.000 |
| No write-back | 10.0 | Hybrid | 0.276 (0.028) | 0.362 (0.075) | 0.479 (0.039) | 0.690 (0.054) | 0.851 (0.021) | +0.000 | +0.000 |
| No write-back | 10.0 | MF | 0.075 (0.007) | 0.086 (0.007) | 0.182 (0.021) | 0.416 (0.021) | 0.709 (0.025) | +0.000 | +0.000 |
| Sci ignore / other mask | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +0.000 | +51.200 |
| Sci ignore / other mask | 10.0 | Graph head | 0.273 (0.025) | 0.366 (0.068) | 0.478 (0.036) | 0.698 (0.050) | 0.854 (0.020) | +0.000 | +51.200 |
| Sci ignore / other mask | 10.0 | Hybrid | 0.273 (0.025) | 0.366 (0.068) | 0.478 (0.036) | 0.698 (0.050) | 0.854 (0.020) | +0.000 | +51.200 |
| Sci ignore / other mask | 10.0 | MF | 0.075 (0.006) | 0.086 (0.007) | 0.182 (0.021) | 0.417 (0.018) | 0.709 (0.025) | +0.000 | +51.200 |
| Sci ignore / other negate | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +51.200 | +0.000 |
| Sci ignore / other negate | 10.0 | Graph head | 0.186 (0.038) | 0.301 (0.060) | 0.380 (0.018) | 0.630 (0.060) | 0.819 (0.025) | +51.200 | +0.000 |
| Sci ignore / other negate | 10.0 | Hybrid | 0.186 (0.038) | 0.301 (0.060) | 0.380 (0.018) | 0.630 (0.060) | 0.819 (0.025) | +51.200 | +0.000 |
| Sci ignore / other negate | 10.0 | MF | 0.030 (0.003) | 0.060 (0.010) | 0.091 (0.009) | 0.203 (0.025) | 0.372 (0.040) | +51.200 | +0.000 |
| Sci mask / other ignore | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +0.000 | +24.200 |
| Sci mask / other ignore | 10.0 | Graph head | 0.279 (0.026) | 0.366 (0.075) | 0.483 (0.038) | 0.696 (0.054) | 0.852 (0.021) | +0.000 | +24.200 |
| Sci mask / other ignore | 10.0 | Hybrid | 0.279 (0.026) | 0.366 (0.075) | 0.483 (0.038) | 0.696 (0.054) | 0.852 (0.021) | +0.000 | +24.200 |
| Sci mask / other ignore | 10.0 | MF | 0.075 (0.007) | 0.087 (0.007) | 0.181 (0.021) | 0.416 (0.021) | 0.709 (0.025) | +0.000 | +24.200 |
| Mask all | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +0.000 | +75.400 |
| Mask all | 10.0 | Graph head | 0.275 (0.026) | 0.367 (0.070) | 0.480 (0.040) | 0.698 (0.042) | 0.856 (0.019) | +0.000 | +75.400 |
| Mask all | 10.0 | Hybrid | 0.275 (0.026) | 0.367 (0.070) | 0.480 (0.040) | 0.698 (0.042) | 0.856 (0.019) | +0.000 | +75.400 |
| Mask all | 10.0 | MF | 0.075 (0.006) | 0.087 (0.008) | 0.182 (0.020) | 0.416 (0.018) | 0.709 (0.025) | +0.000 | +75.400 |
| Sci mask / other negate | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +51.200 | +24.200 |
| Sci mask / other negate | 10.0 | Graph head | 0.188 (0.037) | 0.303 (0.058) | 0.379 (0.018) | 0.634 (0.050) | 0.820 (0.025) | +51.200 | +24.200 |
| Sci mask / other negate | 10.0 | Hybrid | 0.188 (0.037) | 0.303 (0.058) | 0.379 (0.018) | 0.634 (0.050) | 0.820 (0.025) | +51.200 | +24.200 |
| Sci mask / other negate | 10.0 | MF | 0.030 (0.003) | 0.060 (0.010) | 0.092 (0.009) | 0.202 (0.026) | 0.372 (0.040) | +51.200 | +24.200 |
| Sci negate / other ignore | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +24.200 | +0.000 |
| Sci negate / other ignore | 10.0 | Graph head | 0.192 (0.059) | 0.304 (0.081) | 0.402 (0.045) | 0.645 (0.042) | 0.822 (0.030) | +24.200 | +0.000 |
| Sci negate / other ignore | 10.0 | Hybrid | 0.189 (0.061) | 0.299 (0.084) | 0.397 (0.047) | 0.644 (0.047) | 0.822 (0.030) | +24.200 | +0.000 |
| Sci negate / other ignore | 10.0 | MF | 0.049 (0.003) | 0.068 (0.010) | 0.130 (0.013) | 0.317 (0.028) | 0.513 (0.058) | +24.200 | +0.000 |
| Typed | 3.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +24.200 | +51.200 |
| Typed | 3.0 | Graph head | 0.216 (0.046) | 0.322 (0.081) | 0.431 (0.035) | 0.678 (0.050) | 0.840 (0.023) | +24.200 | +51.200 |
| Typed | 3.0 | Hybrid | 0.217 (0.046) | 0.322 (0.081) | 0.431 (0.035) | 0.678 (0.050) | 0.840 (0.023) | +24.200 | +51.200 |
| Typed | 3.0 | MF | 0.050 (0.004) | 0.068 (0.009) | 0.132 (0.012) | 0.316 (0.026) | 0.542 (0.056) | +24.200 | +51.200 |
| Typed | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +24.200 | +51.200 |
| Typed | 10.0 | Graph head | 0.193 (0.054) | 0.304 (0.082) | 0.408 (0.040) | 0.648 (0.047) | 0.824 (0.030) | +24.200 | +51.200 |
| Typed | 10.0 | Hybrid | 0.193 (0.054) | 0.302 (0.083) | 0.407 (0.040) | 0.643 (0.046) | 0.824 (0.030) | +24.200 | +51.200 |
| Typed | 10.0 | MF | 0.048 (0.003) | 0.068 (0.009) | 0.130 (0.012) | 0.316 (0.026) | 0.512 (0.059) | +24.200 | +51.200 |
| Typed | 30.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +24.200 | +51.200 |
| Typed | 30.0 | Graph head | 0.182 (0.059) | 0.291 (0.081) | 0.391 (0.052) | 0.629 (0.049) | 0.816 (0.032) | +24.200 | +51.200 |
| Typed | 30.0 | Hybrid | 0.182 (0.060) | 0.292 (0.081) | 0.392 (0.053) | 0.625 (0.048) | 0.815 (0.029) | +24.200 | +51.200 |
| Typed | 30.0 | MF | 0.050 (0.003) | 0.067 (0.008) | 0.135 (0.009) | 0.311 (0.028) | 0.500 (0.058) | +24.200 | +51.200 |
| Flat negative | 3.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +75.400 | +0.000 |
| Flat negative | 3.0 | Graph head | 0.171 (0.058) | 0.292 (0.089) | 0.377 (0.043) | 0.640 (0.055) | 0.803 (0.041) | +75.400 | +0.000 |
| Flat negative | 3.0 | Hybrid | 0.180 (0.051) | 0.293 (0.081) | 0.388 (0.038) | 0.633 (0.054) | 0.808 (0.034) | +75.400 | +0.000 |
| Flat negative | 3.0 | MF | 0.030 (0.003) | 0.059 (0.009) | 0.082 (0.007) | 0.159 (0.023) | 0.387 (0.047) | +75.400 | +0.000 |
| Flat negative | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +75.400 | +0.000 |
| Flat negative | 10.0 | Graph head | 0.126 (0.022) | 0.238 (0.044) | 0.330 (0.040) | 0.604 (0.033) | 0.787 (0.034) | +75.400 | +0.000 |
| Flat negative | 10.0 | Hybrid | 0.126 (0.021) | 0.239 (0.043) | 0.334 (0.044) | 0.595 (0.032) | 0.784 (0.029) | +75.400 | +0.000 |
| Flat negative | 10.0 | MF | 0.030 (0.002) | 0.060 (0.013) | 0.078 (0.006) | 0.151 (0.018) | 0.339 (0.041) | +75.400 | +0.000 |
| Flat negative | 30.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +75.400 | +0.000 |
| Flat negative | 30.0 | Graph head | 0.111 (0.024) | 0.224 (0.046) | 0.292 (0.034) | 0.551 (0.042) | 0.773 (0.030) | +75.400 | +0.000 |
| Flat negative | 30.0 | Hybrid | 0.111 (0.025) | 0.225 (0.049) | 0.294 (0.032) | 0.549 (0.049) | 0.770 (0.029) | +75.400 | +0.000 |
| Flat negative | 30.0 | MF | 0.029 (0.002) | 0.059 (0.013) | 0.077 (0.002) | 0.151 (0.012) | 0.321 (0.041) | +75.400 | +0.000 |
| Typed + scope | 10.0 | Degree reference | 0.080 (0.009) | 0.036 (0.001) | 0.204 (0.023) | 0.442 (0.033) | 0.741 (0.018) | +14.000 | +61.400 |
| Typed + scope | 10.0 | Graph head | 0.239 (0.048) | 0.334 (0.074) | 0.446 (0.046) | 0.674 (0.055) | 0.845 (0.022) | +14.000 | +61.400 |
| Typed + scope | 10.0 | Hybrid | 0.239 (0.048) | 0.334 (0.074) | 0.446 (0.046) | 0.674 (0.055) | 0.845 (0.022) | +14.000 | +61.400 |
| Typed + scope | 10.0 | MF | 0.058 (0.009) | 0.075 (0.003) | 0.153 (0.024) | 0.344 (0.036) | 0.566 (0.062) | +14.000 | +61.400 |

Table S4c. Approved indications absent from Hetionet (E2) and later scientific failures (E3).

| Policy | Negative weight | Scorer | Pooled AP | Per-disease AP | AUROC | Median percentile, scientific stratum | Median percentile, other-stop stratum | Median percentile, no stop | E3 AUROC |
|---|---|---|---|---|---|---|---|---|---|
| No write-back | 10.0 | Degree reference | +0.006 | +0.020 | 0.692 | 0.878 | 0.874 | +0.790 | +0.501 |
| No write-back | 10.0 | Disease degree | +0.009 | +0.004 | 0.695 | 0.794 | 0.779 | +0.838 | +0.553 |
| No write-back | 10.0 | Graph head | +0.025 | +0.049 | 0.666 | 0.879 | 0.778 | +0.756 | +0.658 |
| No write-back | 10.0 | Hybrid | +0.025 | +0.046 | 0.674 | 0.83 | 0.797 | +0.785 | +0.662 |
| No write-back | 10.0 | MF | +0.014 | +0.037 | 0.734 | 0.941 | 0.914 | +0.863 | +0.560 |
| Sci ignore / other mask | 10.0 | Degree reference | +0.006 | +0.020 | 0.692 | 0.878 | 0.874 | +0.790 | +0.501 |
| Sci ignore / other mask | 10.0 | Graph head | +0.030 | +0.051 | 0.677 | 0.845 | 0.833 | +0.765 | +0.662 |
| Sci ignore / other mask | 10.0 | Hybrid | +0.029 | +0.056 | 0.683 | 0.814 | 0.85 | +0.784 | +0.666 |
| Sci ignore / other mask | 10.0 | MF | +0.015 | +0.040 | 0.744 | 0.945 | 0.934 | +0.864 | +0.559 |
| Sci ignore / other negate | 10.0 | Degree reference | +0.006 | +0.020 | 0.695 | 0.867 | 0.867 | +0.790 | +0.505 |
| Sci ignore / other negate | 10.0 | Graph head | +0.016 | +0.042 | 0.654 | 0.866 | 0.767 | +0.745 | +0.657 |
| Sci ignore / other negate | 10.0 | Hybrid | +0.015 | +0.040 | 0.662 | 0.846 | 0.801 | +0.755 | +0.662 |
| Sci ignore / other negate | 10.0 | MF | +0.010 | +0.032 | 0.671 | 0.889 | 0.823 | +0.863 | +0.560 |
| Sci mask / other ignore | 10.0 | Degree reference | +0.006 | +0.020 | 0.692 | 0.878 | 0.874 | +0.790 | +0.501 |
| Sci mask / other ignore | 10.0 | Graph head | +0.028 | +0.047 | 0.666 | 0.894 | 0.786 | +0.757 | +0.653 |
| Sci mask / other ignore | 10.0 | Hybrid | +0.027 | +0.046 | 0.675 | 0.872 | 0.808 | +0.780 | +0.661 |
| Sci mask / other ignore | 10.0 | MF | +0.014 | +0.039 | 0.743 | 0.949 | 0.919 | +0.859 | +0.555 |
| Mask all | 10.0 | Degree reference | +0.006 | +0.020 | 0.692 | 0.878 | 0.874 | +0.790 | +0.501 |
| Mask all | 10.0 | Graph head | +0.031 | +0.056 | 0.678 | 0.867 | 0.836 | +0.760 | +0.658 |
| Mask all | 10.0 | Hybrid | +0.030 | +0.057 | 0.684 | 0.828 | 0.852 | +0.790 | +0.665 |
| Mask all | 10.0 | MF | +0.015 | +0.041 | 0.748 | 0.955 | 0.928 | +0.865 | +0.557 |
| Sci mask / other negate | 10.0 | Degree reference | +0.006 | +0.020 | 0.695 | 0.867 | 0.867 | +0.790 | +0.505 |
| Sci mask / other negate | 10.0 | Graph head | +0.019 | +0.038 | 0.655 | 0.876 | 0.776 | +0.744 | +0.653 |
| Sci mask / other negate | 10.0 | Hybrid | +0.017 | +0.039 | 0.662 | 0.879 | 0.789 | +0.760 | +0.659 |
| Sci mask / other negate | 10.0 | MF | +0.010 | +0.033 | 0.672 | 0.896 | 0.824 | +0.863 | +0.558 |
| Sci negate / other ignore | 10.0 | Degree reference | +0.006 | +0.020 | 0.693 | 0.872 | 0.872 | +0.790 | +0.502 |
| Sci negate / other ignore | 10.0 | Graph head | +0.023 | +0.042 | 0.662 | 0.827 | 0.75 | +0.761 | +0.662 |
| Sci negate / other ignore | 10.0 | Hybrid | +0.023 | +0.039 | 0.67 | 0.777 | 0.793 | +0.781 | +0.665 |
| Sci negate / other ignore | 10.0 | MF | +0.013 | +0.035 | 0.686 | 0.864 | 0.904 | +0.859 | +0.558 |
| Typed | 3.0 | Degree reference | +0.006 | +0.020 | 0.692 | 0.878 | 0.874 | +0.790 | +0.501 |
| Typed | 3.0 | Graph head | +0.028 | +0.049 | 0.675 | 0.823 | 0.835 | +0.763 | +0.663 |
| Typed | 3.0 | Hybrid | +0.027 | +0.050 | 0.682 | 0.763 | 0.843 | +0.784 | +0.667 |
| Typed | 3.0 | MF | +0.014 | +0.039 | 0.732 | 0.928 | 0.933 | +0.864 | +0.563 |
| Typed | 10.0 | Degree reference | +0.006 | +0.020 | 0.693 | 0.878 | 0.874 | +0.790 | +0.502 |
| Typed | 10.0 | Graph head | +0.027 | +0.049 | 0.673 | 0.82 | 0.831 | +0.766 | +0.665 |
| Typed | 10.0 | Hybrid | +0.025 | +0.050 | 0.68 | 0.751 | 0.845 | +0.786 | +0.668 |
| Typed | 10.0 | MF | +0.014 | +0.038 | 0.693 | 0.874 | 0.906 | +0.860 | +0.561 |
| Typed | 30.0 | Degree reference | +0.006 | +0.020 | 0.694 | 0.867 | 0.867 | +0.790 | +0.505 |
| Typed | 30.0 | Graph head | +0.026 | +0.049 | 0.671 | 0.811 | 0.833 | +0.765 | +0.666 |
| Typed | 30.0 | Hybrid | +0.024 | +0.049 | 0.677 | 0.74 | 0.839 | +0.782 | +0.668 |
| Typed | 30.0 | MF | +0.014 | +0.039 | 0.673 | 0.749 | 0.906 | +0.860 | +0.551 |
| Flat negative | 3.0 | Degree reference | +0.006 | +0.020 | 0.693 | 0.872 | 0.872 | +0.790 | +0.502 |
| Flat negative | 3.0 | Graph head | +0.019 | +0.039 | 0.658 | 0.833 | 0.768 | +0.755 | +0.660 |
| Flat negative | 3.0 | Hybrid | +0.019 | +0.036 | 0.666 | 0.795 | 0.777 | +0.768 | +0.665 |
| Flat negative | 3.0 | MF | +0.012 | +0.035 | 0.684 | 0.893 | 0.893 | +0.863 | +0.558 |
| Flat negative | 10.0 | Degree reference | +0.006 | +0.020 | 0.695 | 0.867 | 0.867 | +0.790 | +0.505 |
| Flat negative | 10.0 | Graph head | +0.015 | +0.037 | 0.652 | 0.814 | 0.753 | +0.748 | +0.664 |
| Flat negative | 10.0 | Hybrid | +0.014 | +0.033 | 0.659 | 0.789 | 0.771 | +0.757 | +0.666 |
| Flat negative | 10.0 | MF | +0.010 | +0.031 | 0.662 | 0.828 | 0.774 | +0.855 | +0.560 |
| Flat negative | 30.0 | Degree reference | +0.006 | +0.020 | 0.698 | 0.861 | 0.871 | +0.800 | +0.507 |
| Flat negative | 30.0 | Graph head | +0.013 | +0.033 | 0.644 | 0.791 | 0.708 | +0.741 | +0.664 |
| Flat negative | 30.0 | Hybrid | +0.012 | +0.031 | 0.649 | 0.765 | 0.742 | +0.748 | +0.665 |
| Flat negative | 30.0 | MF | +0.008 | +0.026 | 0.621 | 0.661 | 0.664 | +0.853 | +0.541 |
| Typed + scope | 10.0 | Degree reference | +0.006 | +0.020 | 0.692 | 0.878 | 0.874 | +0.790 | +0.501 |
| Typed + scope | 10.0 | Graph head | +0.030 | +0.051 | 0.677 | 0.864 | 0.836 | +0.760 | +0.658 |
| Typed + scope | 10.0 | Hybrid | +0.030 | +0.052 | 0.683 | 0.801 | 0.843 | +0.788 | +0.665 |
| Typed + scope | 10.0 | MF | +0.015 | +0.041 | 0.706 | 0.942 | 0.926 | +0.860 | +0.559 |

## S5. Policy and scorer contrasts with bootstrap intervals

Table S5. Differences from no write-back for every scorer (pooled and per-disease AP) with 95% paired disease-cluster bootstrap intervals (1,000 draws).

| Set | Contrast | Pooled AP difference | Per-disease AP difference | Pooled AP 95% CI | Per-disease AP 95% CI |
|---|---|---|---|---|---|
| E1 random edge | Flat negative vs no write-back (degree reference) | +0.000 | +0.000 | (−0.000, +0.002) | (+0.000, +0.000) |
| E1 random edge | Flat negative vs no write-back (MF) | -0.066 | -0.076 | (−0.100, −0.041) | (−0.101, −0.054) |
| E1 random edge | Flat negative vs no write-back (graph head) | -0.102 | -0.095 | (−0.134, −0.061) | (−0.125, −0.069) |
| E1 random edge | Flat negative vs no write-back (hybrid) | -0.108 | -0.098 | (−0.142, −0.065) | (−0.129, −0.071) |
| E1 random edge | Mask all vs no write-back (degree reference) | -0.000 | +0.000 | (−0.000, +0.000) | (+0.000, +0.000) |
| E1 random edge | Mask all vs no write-back (MF) | +0.012 | +0.020 | (+0.002, +0.021) | (+0.010, +0.030) |
| E1 random edge | Mask all vs no write-back (graph head) | +0.007 | +0.024 | (−0.019, +0.040) | (+0.008, +0.043) |
| E1 random edge | Mask all vs no write-back (hybrid) | +0.012 | +0.023 | (−0.011, +0.041) | (+0.004, +0.040) |
| E1 random edge | Typed vs no write-back (degree reference) | +0.000 | +0.000 | (−0.000, +0.001) | (+0.000, +0.000) |
| E1 random edge | Typed vs no write-back (MF) | -0.016 | -0.018 | (−0.037, +0.006) | (−0.036, −0.001) |
| E1 random edge | Typed vs no write-back (graph head) | -0.029 | -0.018 | (−0.059, +0.008) | (−0.043, +0.005) |
| E1 random edge | Typed vs no write-back (hybrid) | -0.031 | -0.019 | (−0.058, +0.002) | (−0.042, +0.003) |
| random | Typed + scope - no_writeback [degree] | +0.000 | +0.000 | (−0.000, +0.000) | (+0.000, +0.000) |
| random | Typed + scope - no_writeback [mf] | -0.008 | -0.011 | (−0.027, +0.008) | (−0.026, +0.003) |
| random | Typed + scope - no_writeback [graph] | -0.022 | -0.003 | (−0.050, +0.014) | (−0.027, +0.021) |
| random | Typed + scope - no_writeback [hybrid] | -0.019 | -0.003 | (−0.044, +0.013) | (−0.025, +0.018) |
| E1 compound disjoint | Flat negative vs no write-back (degree reference) | +0.000 | +0.000 | (+0.000, +0.000) | (+0.000, +0.000) |
| E1 compound disjoint | Flat negative vs no write-back (MF) | -0.045 | -0.026 | (−0.065, −0.023) | (−0.037, −0.015) |
| E1 compound disjoint | Flat negative vs no write-back (graph head) | -0.149 | -0.124 | (−0.188, −0.105) | (−0.147, −0.100) |
| E1 compound disjoint | Flat negative vs no write-back (hybrid) | -0.150 | -0.123 | (−0.190, −0.104) | (−0.147, −0.098) |
| E1 compound disjoint | Mask all vs no write-back (degree reference) | +0.000 | +0.000 | (+0.000, +0.000) | (+0.000, +0.000) |
| E1 compound disjoint | Mask all vs no write-back (MF) | -0.000 | +0.001 | (−0.001, +0.001) | (−0.001, +0.003) |
| E1 compound disjoint | Mask all vs no write-back (graph head) | -0.001 | +0.004 | (−0.013, +0.012) | (−0.002, +0.011) |
| E1 compound disjoint | Mask all vs no write-back (hybrid) | -0.001 | +0.004 | (−0.013, +0.012) | (−0.002, +0.011) |
| E1 compound disjoint | Typed vs no write-back (degree reference) | +0.000 | +0.000 | (+0.000, +0.000) | (+0.000, +0.000) |
| E1 compound disjoint | Typed vs no write-back (MF) | -0.027 | -0.018 | (−0.046, −0.011) | (−0.027, −0.011) |
| E1 compound disjoint | Typed vs no write-back (graph head) | -0.083 | -0.058 | (−0.127, −0.045) | (−0.083, −0.036) |
| E1 compound disjoint | Typed vs no write-back (hybrid) | -0.082 | -0.060 | (−0.125, −0.046) | (−0.085, −0.038) |
| compound | Typed + scope - no_writeback [degree] | +0.000 | +0.000 | (+0.000, +0.000) | (+0.000, +0.000) |
| compound | Typed + scope - no_writeback [mf] | -0.017 | -0.011 | (−0.026, −0.011) | (−0.017, −0.004) |
| compound | Typed + scope - no_writeback [graph] | -0.037 | -0.028 | (−0.060, −0.012) | (−0.042, −0.014) |
| compound | Typed + scope - no_writeback [hybrid] | -0.036 | -0.028 | (−0.060, −0.012) | (−0.042, −0.014) |
| E2 | Flat negative vs no write-back (degree reference) | +0.000 | +0.000 | (−0.000, +0.000) | (+0.000, +0.000) |
| E2 | Flat negative vs no write-back (MF) | -0.004 | -0.006 | (−0.008, −0.001) | (−0.012, −0.001) |
| E2 | Flat negative vs no write-back (graph head) | -0.010 | -0.012 | (−0.019, −0.003) | (−0.029, +0.001) |
| E2 | Flat negative vs no write-back (hybrid) | -0.011 | -0.013 | (−0.020, −0.002) | (−0.031, −0.001) |
| E2 | Mask all vs no write-back (degree reference) | -0.000 | +0.000 | (−0.000, +0.000) | (+0.000, +0.000) |
| E2 | Mask all vs no write-back (MF) | +0.001 | +0.004 | (−0.001, +0.003) | (+0.001, +0.007) |
| E2 | Mask all vs no write-back (graph head) | +0.006 | +0.008 | (−0.003, +0.016) | (−0.008, +0.022) |
| E2 | Mask all vs no write-back (hybrid) | +0.005 | +0.011 | (−0.002, +0.015) | (+0.001, +0.024) |
| E2 | Typed vs no write-back (degree reference) | +0.000 | +0.000 | (−0.000, +0.000) | (+0.000, +0.000) |
| E2 | Typed vs no write-back (MF) | +0.000 | +0.002 | (−0.001, +0.003) | (−0.002, +0.005) |
| E2 | Typed vs no write-back (graph head) | +0.002 | -0.000 | (−0.005, +0.010) | (−0.022, +0.017) |
| E2 | Typed vs no write-back (hybrid) | +0.000 | +0.004 | (−0.006, +0.008) | (−0.016, +0.019) |
| e2 | Typed + scope - no_writeback [degree] | -0.000 | +0.000 | (−0.000, +0.000) | (+0.000, +0.000) |
| e2 | Typed + scope - no_writeback [mf] | +0.001 | +0.004 | (−0.001, +0.003) | (+0.001, +0.008) |
| e2 | Typed + scope - no_writeback [graph] | +0.005 | +0.002 | (−0.003, +0.014) | (−0.021, +0.019) |
| e2 | Typed + scope - no_writeback [hybrid] | +0.004 | +0.006 | (−0.003, +0.013) | (−0.014, +0.023) |

## S6. Stop-reason categories and registry status

Table S6. Hetionet pairs with a stopped phase 1–3 trial, by Open Targets stop-reason category, and the share that are recorded treatments or have approval evidence. A trial can carry several categories.

| Category | Stopped trials | Pairs | Recorded or approved | Share |
|---|---|---|---|---|
| Insufficient_Enrollment | 1,991 | 1,340 | 416 | 31% |
| Business_Administrative | 1,367 | 960 | 340 | 35% |
| No reason given | 774 | 696 | 293 | 42% |
| Negative | 486 | 434 | 171 | 39% |
| Logistics_Resources | 308 | 362 | 141 | 39% |
| Study_Design | 305 | 356 | 158 | 44% |
| Invalid_Reason | 259 | 322 | 153 | 48% |
| Safety_Sideeffects | 229 | 255 | 118 | 46% |
| Another_Study | 205 | 229 | 106 | 46% |
| Study_Staff_Moved | 184 | 245 | 116 | 47% |
| Regulatory | 103 | 121 | 63 | 52% |
| Uncategorised | 100 | 153 | 82 | 54% |
| Covid19 | 70 | 94 | 34 | 36% |
| No_Context | 52 | 45 | 31 | 69% |

By registry status (phase 1–3): terminated 4,386 trials on 2,119 pairs (26% recorded or approved); withdrawn 1,468 trials on 1,181 pairs (34%); suspended 129 trials on 179 pairs (53%).
