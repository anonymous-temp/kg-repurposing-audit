# Reproducibility

- **Inputs.** `paper_results/input_sha256.txt` lists the SHA-256 of every input file used for the manuscript.
- **Determinism.** Partitions use fixed seeds (42, 1, 2, 3, 4); seed 42 reproduces the test sets of release 0.2.0. Embedding pretraining uses seed 1; scorers use initialisation seeds 1–3; the bootstrap uses seed 20261006. CPU PyTorch results can differ in the last digits across library versions and thread counts.
- **Selection.** Scorer hyper-parameters are selected on validation AP without write-back and then fixed for every policy (`paper_results/selection.json`, full grid in `selection_grid.tsv`). Test sets never enter selection.
- **Idempotence.** Each pipeline stage writes its own output and skips existing outputs, so an interrupted run resumes where it stopped.
- **Levels of reproduction.** (1) Tables and figures from `paper_results/` with `pipelines/09_make_tables.py` and `08_make_figures.py`; (2) metrics and intervals from the released candidate-level scores; (3) full rerun from the public inputs with the numbered pipelines.
- **Recorded runtime.** Python 3.12.3, PyTorch 2.14.0 (CPU), NumPy 2.2.6, pandas 2.2.3, scikit-learn 1.6.1, PyArrow 25.0, LinkML 1.11.1; 4-core CPU.
