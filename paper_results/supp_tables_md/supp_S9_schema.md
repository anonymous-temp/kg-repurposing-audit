## S9. Schema classes, Biolink and PROV-O mappings

The handoff schema (schemas/handoff.linkml.yaml, version 0.3.0) maps classes and slots to Biolink Model 4.4.5 and PROV-O where an equivalent term exists. Terms marked "local" have no equivalent in Biolink 4.4.5; the stop-reason vocabulary follows Razuvayevskaya et al. (2024).

Table S8. Schema elements and mappings.

| Schema element | Mapping | Note |
|---|---|---|
| TreatmentStrategy | biolink:ChemicalOrDrugOrTreatmentToDiseaseOrPhenotypicFeatureAssociation | drug = biolink:subject; indication = biolink:object |
| TreatmentStrategy.population | biolink:population_context_qualifier | |
| TreatmentStrategy.knowledge_level | biolink:knowledge_level (KnowledgeLevelEnum) | 'prediction' for model output |
| TreatmentStrategy.agent_type | biolink:agent_type (AgentTypeEnum) | 'computational_model' for model output |
| TreatmentStrategy.evidence | biolink:has_evidence | |
| CandidatePrediction | prov:Entity | source_sha256 = prov:wasDerivedFrom |
| EvidenceAssertion | biolink:Association | source = biolink:primary_knowledge_source; source_date = prov:generatedAtTime |
| EvidenceAssertion.predicate | biolink:treats, biolink:in_clinical_trials_for, biolink:studied_to_treat, biolink:contraindicated_in | |
| EvidenceAssertion.negated | biolink:negated | |
| EvidenceAssertion.clinical_approval_status | biolink:clinical_approval_status (ClinicalApprovalStatusEnum) | |
| EvidenceAssertion.research_phase | biolink:max_research_phase (ResearchPhaseEnum) | |
| FailureAnnotation | biolink:ClinicalTrial | for event_type = stopped_trial |
| FailureAnnotation.registry_status | biolink:clinical_trial_overall_status | subset TERMINATED, WITHDRAWN, SUSPENDED |
| FailureAnnotation.trial_phase | biolink:clinical_trial_phase | |
| FailureAnnotation.trial_start_date | biolink:clinical_trial_start_date | |
| FailureAnnotation.event_type | local | stopped trial, development discontinuation, market withdrawal, corporate event |
| FailureAnnotation.stop_reason_categories | local (17 classes + Uncategorised) | Razuvayevskaya et al. 2024 |
| FailureAnnotation.failure_type | local | efficacy, safety, operational, design, uninformative, success, not reported |
| FailureAnnotation.scope_match | local | same, narrower or broader concept; unmapped |
| FailureAnnotation.implication | local | retain, negate, qualify, defer (rule in kg_audit.failures) |
| ActionabilityConstraint | local | eight review domains with pass, concern, fail or unknown |

**Validation.** The JSON Schema in schemas/handoff.schema.json was generated from the LinkML source with LinkML 1.11.1 (JsonSchemaGenerator, top class TreatmentStrategy). Unit tests check that the enumerations in the LinkML file equal the constants used by the rule implementation, that every example validates against the JSON Schema, and that the checker reproduces the stored outputs.

**Checker semantics.** validate_record rejects a failure annotation whose stated implication differs from the one implied by its failure type, scope match and recorded-indication flag. An operational stop therefore cannot be recorded as negating a treatment. writeback_action returns store_scoped_counterevidence only for a same-concept efficacy or safety failure with a dated source. Other scientific failures are stored as context qualifications, and all other reasons are deferred. assess_strategy reports scoped counterevidence only when the drug, indication, population, regimen, comparator and endpoint of the failure all match the strategy under review. These comparisons ignore case and repeated whitespace. Evidence dated after the requested review date is never used.
