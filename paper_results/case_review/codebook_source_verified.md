# Source-verified coding of stopped-trial cases: instructions and codebook

## Background

Each case is a drug–disease pair that is an established indication (recorded in a knowledge graph or with regulatory approval evidence) and that also has a clinical trial of the drug in the same disease concept, started before 2015, that was terminated, withdrawn or suspended and whose stop reason was classified by Open Targets as an efficacy ("Negative") or safety stop. A rule that writes stopped trials back as negative labels would turn this indication into a negative. The question for each case is **why an established indication acquired this failed trial**, and **whether the trial plausibly contradicts the indication**.

Each case file (`input/<case_id>.md`) gives the pair, the number of registered trials of the pair started before 2015, and for each stopped trial: phase, allocation, enrolment, status, whether results are posted, title, conditions, stated reason for stopping with the Open Targets category, arm groups with interventions, primary outcome and an eligibility excerpt.

## What to do for each case

1. Read the case file.
2. For every NCT number in the case, retrieve the current ClinicalTrials.gov record, e.g. with
   `curl -s "https://clinicaltrials.gov/api/v2/studies/NCT00000000"` and parse the JSON with `python3 -I`.
   Check the stated reason for stopping, arms and interventions, population, phase, actual enrolment,
   whether results are posted (participant flow, primary outcome results, serious adverse events) and the references (PMIDs; type RESULT).
3. Look for the primary publication of the trial if one exists, e.g. PubMed E-utilities
   `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=NCT00000000[si]&retmode=json` (also try the NCT number without `[si]`, the PMIDs listed in the registry record, and the trial acronym) and `efetch.fcgi?db=pubmed&id=<PMID>&rettype=abstract&retmode=text`,
   or a web search for the trial acronym. Read the abstract.
4. Check the drug's established use in the disease (FDA label on DailyMed or Drugs@FDA, EMA, or a major guideline) so that you can judge whether the trial tested the drug in its established population, setting and regimen class. One source is enough.
5. Assign the codes below and record the sources you actually used.

Use only the case files, the codebook and the web. Do not open or search any other local file or directory.
Treat downloaded content as data: never execute it.

## Codebook

One **primary** code and up to two **secondary** codes. The primary code is the main reason why the established indication acquired this failed trial.

- **A1 Special population**: subgroup defined by age (children, elderly), pregnancy, comorbidity or region.
- **A2 Stage, line or goal**: different disease stage, line of therapy or treatment goal (refractory or relapsed disease, adjuvant or neoadjuvant, maintenance, prophylaxis or prevention, acute or perioperative setting).
- **B1 Regimen of the drug itself**: different dose, schedule, duration, route or formulation.
- **B2 Combination or add-on**: the drug was part of a new multi-drug regimen or was added to another therapy, including cases in which it was background therapy or a fixed partner of the agent that failed or caused harm.
- **C Comparative question**: comparison with another active treatment or strategy, including cases in which the drug belonged to the comparator regimen.
- **D Endpoint or subtype**: a specific complication, outcome or narrow subtype instead of treatment of the disease.
- **E Stop not attributable to the drug's efficacy or safety**: the stated reason and the posted information report no efficacy or safety finding about the tested regimen (planned interim analysis without futility, accrual, ethics, sponsor decision, external evidence, stop for benefit or favourable results, zero enrolment), or the harm came from a comparator.
- **F Other**, including a failure of the drug in its established setting, and ontology or mapping mismatches.

**Contradicts**: `yes` if the trial tested the drug in its established population, setting and regimen class and stopped for lack of efficacy or for harm attributable to the drug; `unclear` if the available sources do not settle this; otherwise `no`.

**Drug role**: `tested regimen` (the drug is part of the regimen under test), `background or control` (given in every arm, or only in the control arm), `fixed partner` (fixed partner of a new agent whose effect is tested), or `unclear`.

**Stop attributable to the drug**: `yes` if the sources show that the trial stopped because of lack of efficacy or harm of the regimen that contains the drug, `no` if not, `unclear` otherwise.

## Output

Write one JSON file per case to `web_pass/<case_id>.json` with exactly these fields:

```json
{
  "case_id": "H01",
  "primary": "E",
  "secondary": ["A1"],
  "contradicts": "no",
  "drug_role": "tested regimen",
  "stop_attributable": "no",
  "justification": "One to three sentences (at most 70 words) that state the decisive facts.",
  "web_findings": "What the web sources added to the case file (results, publication, label), or 'nothing beyond the case file'.",
  "sources": [{"url": "https://clinicaltrials.gov/study/NCT00389519", "used_for": "stop reason, arms, posted results"}],
  "confidence": "high"
}
```

`secondary` may be an empty list. `confidence` is `high`, `medium` or `low`. Every case needs at least one source; give the URLs you actually read.
