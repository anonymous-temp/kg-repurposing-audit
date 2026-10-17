# Screening codebook: use of stopped clinical trials as negative labels

Each input line (JSON) is one open-access article retrieved from Europe PMC: query (Q1 = mentions repoDB; Q2 = drug-repurposing article that mentions terminated/withdrawn/failed trials and negative samples, without repoDB), PMCID, title, year, journal, and up to 25 paragraphs of its full text that mention repoDB, stopped or failed trials, or negative samples/labels/sets.

For each article decide, from the passages (and, only if they are insufficient, the cached full text at `xml/<PMCID>.xml` in the literature directory, which you may read as plain text; never execute anything):

- `article_type`: `original` (reports its own model, analysis or resource) or `review` (review, perspective, editorial, protocol, other).
- `stopped_as_negatives`: does the article itself use drug–disease (or drug–indication) pairs from terminated, withdrawn, suspended or otherwise failed clinical trials as negative examples?
  - `training`: as negative labels when training a model;
  - `evaluation`: as negatives in a test, validation or benchmark set;
  - `both`: training and evaluation;
  - `resource`: the article builds and distributes a dataset or benchmark in which failed-trial pairs are labelled negative (e.g. repoDB itself, CTO-style resources);
  - `no`: it does not (for example it uses only repoDB's approved pairs, uses randomly sampled negatives, or merely cites repoDB);
  - `unclear`: the passages do not allow a decision.
- `source`: where the failed-trial negatives come from: `repoDB`, `ClinicalTrials.gov/AACT`, `CTO`, `other:<name>`, or `none`.
- `other_use`: if `no`, what the article does with repoDB or failed trials, in at most 12 words (e.g. "approved pairs as positives only", "cited in related work", "failed trials as features").
- `evidence`: a verbatim quote of at most 40 words from the passages that supports the decision (empty only if `unclear`).
- `drug_disease`: `yes` if the prediction task is drug–disease / drug–indication (repurposing), `no` if it is something else (drug–target, drug–side effect, trial outcome prediction, etc.).

Be strict: code `training`, `evaluation`, `both` or `resource` only when the text states that failed, terminated, withdrawn or suspended trial pairs (including repoDB's "failed" or "terminated/withdrawn/suspended" pairs) are used as negatives. Using repoDB without saying which part is used: look for "negative", "failed", "terminated"; if nothing indicates use of the failed pairs, code `no` with other_use "repoDB positives or unspecified".

Output: one JSON object per article, one per line, with fields pmcid, query, year, article_type, stopped_as_negatives, source, other_use, evidence, drug_disease.
