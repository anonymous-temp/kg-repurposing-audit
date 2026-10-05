"""Declared evidence-handling invariants, independent of clinical classification."""
import copy
import unittest
from kg_audit.evidence import assess_strategy, validate_record, writeback_action, DOMAINS


def complete_record():
    return {'id':'synthetic:strategy1','drug':'synthetic:drug1','indication':'synthetic:disease1',
            'population':'synthetic population','comparator':'synthetic comparator','regimen':'synthetic regimen',
            'endpoint':'synthetic endpoint','required_exposure':'specified for this synthetic fixture',
            'assessment_date':'2026-01-01','evidence':[{'id':'synthetic:source1','source':'https://example.org/synthetic',
            'source_date':'2025-12-01'}],
            'constraints':[{'domain':d,'status':'pass','evidence_ids':['synthetic:source1']} for d in DOMAINS],
            'failures':[]}


class EvidenceTests(unittest.TestCase):
    def test_complete_fixture_can_pass_documentation_gate(self):
        self.assertEqual(assess_strategy(complete_record(),'2026-01-01')['status'],'ready_for_evidence_review')

    def test_missing_every_domain_never_passes(self):
        for domain in DOMAINS:
            record=complete_record();record['constraints']=[c for c in record['constraints'] if c['domain']!=domain]
            with self.subTest(domain=domain):
                self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'incomplete')

    def test_empty_or_unknown_critical_values_never_pass(self):
        for key in ['population','comparator','regimen','endpoint','required_exposure']:
            for value in [None,'','unknown','not reported']:
                record=complete_record();record[key]=value
                self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'incomplete')

    def test_future_evidence_cannot_support_current_readiness(self):
        record=complete_record();record['evidence'][0]['source_date']='2026-02-01'
        self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'incomplete')

    def test_future_failure_does_not_supersede_past_assessment(self):
        record=complete_record();record['failures']=[{'domain':'efficacy','implication':'negate',
            'source_date':'2026-02-01','evidence_ids':['synthetic:source2'],
            'scope':{'drug':'synthetic:drug1','indication':'synthetic:disease1','population':'synthetic population','endpoint':'synthetic endpoint','regimen':'synthetic regimen','comparator':'synthetic comparator'}}]
        record['evidence'].append({'id':'synthetic:source2','source':'https://example.org/future','source_date':'2026-02-01'})
        self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'ready_for_evidence_review')
        self.assertEqual(assess_strategy(record,'2026-03-01')['status'],'scoped_counterevidence')

    def test_missing_timestamp_is_not_silently_current(self):
        record=complete_record();record['evidence'][0].pop('source_date')
        self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'incomplete')

    def test_unknown_exposure_with_explanatory_suffix_stays_unknown(self):
        record=complete_record();record['required_exposure']='unknown for the proposed effect'
        self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'incomplete')

    def test_future_assessment_is_not_applied_retroactively(self):
        record=complete_record()
        self.assertEqual(assess_strategy(record,'2025-12-15')['status'],'not_yet_assessed')

    def test_percentile_interval_may_exclude_observed_point(self):
        record=complete_record();record['evidence'][0]['effect']={'measure':'bootstrap difference','estimate':.1,'lower':.11,'upper':.2}
        self.assertEqual(validate_record(record),[])

    def test_orphan_evidence_is_validation_error(self):
        record=complete_record();record['constraints'][0]['evidence_ids']=['missing']
        self.assertTrue(validate_record(record))

    def test_contradictory_domain_records_do_not_overwrite(self):
        record=complete_record();record['constraints'].append({'domain':DOMAINS[0],'status':'fail','evidence_ids':['synthetic:source1']})
        self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'conflicting_records')

    def test_unscoped_negative_is_not_written_as_global_negative(self):
        self.assertEqual(writeback_action({'domain':'efficacy','implication':'negate'}),'defer')

    def test_registry_status_alone_cannot_make_negative(self):
        for status in ['TERMINATED','WITHDRAWN','SUSPENDED','UNKNOWN','COMPLETED']:
            self.assertEqual(writeback_action({'registry_status':status}),'defer')

    def test_non_efficacy_failure_does_not_negate(self):
        for domain in ['commercial','recruitment','execution','unknown']:
            self.assertEqual(writeback_action({'domain':domain,'implication':'negate'}),'defer')

    def test_scoped_negative_preserves_scope(self):
        f={'domain':'efficacy','implication':'negate','source_date':'2025-12-01','evidence_ids':['synthetic:source1'],
           'scope':{'drug':'d','indication':'i','population':'p','regimen':'r','endpoint':'e','comparator':'c'}}
        self.assertEqual(writeback_action(f),'store_scoped_counterevidence')

    def test_different_regimen_or_comparator_does_not_supersede(self):
        for field in ('regimen','comparator'):
            record=complete_record()
            scope={k:record[k] for k in ('drug','indication','population','regimen','endpoint','comparator')}
            scope[field]='a different context'
            record['failures']=[{'domain':'efficacy','implication':'negate','source_date':'2025-12-01',
                                 'evidence_ids':['synthetic:source1'],'scope':scope}]
            self.assertEqual(assess_strategy(record,'2026-01-01')['status'],'ready_for_evidence_review')


if __name__=='__main__':unittest.main()
