## S1. Mapping Open Targets 26.09 to Hetionet v1.0

**Compounds.** Each Open Targets drug record lists DrugBank identifiers among its cross-references. A ChEMBL identifier was assigned to a Hetionet compound when either the record itself or its parent (salt) record carried the compound's DrugBank identifier. This mapped 1,460 of the 1,552 Hetionet compounds to 2,901 ChEMBL identifiers. Unmapped compounds remained in all grids and received no clinical evidence.

**Diseases.** 132 of the 137 Hetionet diseases had a current DOID cross-reference in the Open Targets disease index. The remaining five were resolved as follows; each assignment was checked against the Open Targets term name and its obsolete cross-references.

Table S1. Hetionet diseases without a current DOID cross-reference in Open Targets 26.09.

| Hetionet disease | DOID | Open Targets term used | Basis |
|---|---|---|---|
| multiple sclerosis | DOID:2377 | EFO_0803536 | obsolete DOID cross-reference; in release 26.09 this term is the parent of relapsing–remitting and chronic progressive multiple sclerosis |
| azoospermia | DOID:14227 | HP_0000027 (Azoospermia) | exact name |
| endogenous depression | DOID:1595 | MONDO_0002009 (major depressive disorder) | obsolete DOID cross-reference |
| vascular cancer | DOID:175 | EFO_0003967 (vascular sarcoma) | obsolete DOID cross-reference |
| pleural cancer | DOID:9917 | MONDO_0006294 (pleural cancer) | exact name |

An Open Targets condition term without a direct mapping was assigned to the Hetionet disease of its most specific mapped ancestor (largest number of ancestors). Of the 2,950 condition terms assigned, 142 were direct and 2,827 were assigned through an ancestor; 19 terms reached two Hetionet diseases at equal depth and were assigned to both. A median of eight condition terms was assigned to each Hetionet disease. The depth gap between a trial condition and its Hetionet disease defined the scope match; 46.0% of stopped-trial report–pair rows were same-concept assignments.

**Clinical reports.** After mapping, 56,825 clinical reports involved at least one mapped compound and one mapped condition: 47,259 trial reports from ClinicalTrials.gov, 7,743 drug-label reports (DailyMed), 982 curated-resource reports and 841 regulatory reports. Approval-stage reports from non-trial sources came from DailyMed (7,743 reports), TTD (497), EMA (364), PMDA (275), ATC (203) and FDA (54). Stopped trials comprised 4,850 terminated, 1,693 withdrawn and 135 suspended reports.

## S2. Partitions

Table S2. Partitions of the recorded treatments.

| Task | Seed | Fitting positives | Validation positives | Test positives | Test compounds | Grid pairs | Test compounds without fitting label | Test positives from such compounds |
|---|---|---|---|---|---|---|---|---|
| Random edge | 42 | 544 | 60 | 151 | 116 | 8,676 | 59 | 64 |
| Random edge | 1 | 544 | 60 | 151 | 111 | 8,304 | 57 | 61 |
| Random edge | 2 | 544 | 60 | 151 | 107 | 7,996 | 52 | 58 |
| Random edge | 3 | 544 | 60 | 151 | 121 | 9,038 | 56 | 60 |
| Random edge | 4 | 544 | 60 | 151 | 107 | 7,991 | 46 | 49 |
| Compound disjoint | 42 | 509 | 84 | 162 | 77 | 5,913 | 77 | 162 |
| Compound disjoint | 1 | 555 | 55 | 145 | 77 | 5,920 | 77 | 145 |
| Compound disjoint | 2 | 533 | 63 | 159 | 77 | 5,913 | 77 | 159 |
| Compound disjoint | 3 | 501 | 114 | 140 | 77 | 5,918 | 77 | 140 |
| Compound disjoint | 4 | 564 | 48 | 143 | 77 | 5,924 | 77 | 143 |

Of the Hetionet compounds with a recorded treatment, 67% have only one. Random-edge holdout therefore leaves 46–59 test compounds per partition without any fitting label, and about 40% of random-edge test positives pose a compound cold-start problem. Seed 42 reproduces the test sets of software release v0.2.0.
