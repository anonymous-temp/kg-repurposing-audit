import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import urlparse,parse_qs
import pandas as pd
import jsonschema
from kg_audit import cli
from kg_audit.bridge import export_records
from kg_audit.registry import fetch_sample,reason_cues


class IOTests(unittest.TestCase):
    def test_cli_demo_and_paired_comparison(self):
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            cli.demo(tmp);frame=pd.read_csv(Path(tmp)/'synthetic_scores.csv')
            self.assertFalse(frame.duplicated(['compound','disease']).any())
            output=Path(tmp)/'comparison.json'
            with patch('sys.argv',['kg-audit','compare','--scores',str(Path(tmp)/'synthetic_scores.csv'),'--model','model','--reference','reference','--bootstrap','100','--output',str(output)]):cli.main()
            self.assertIn('input_sha256',json.loads(output.read_text()))

    def test_cli_evidence_preserves_missingness(self):
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            record=Path(tmp)/'record.json';record.write_text(json.dumps({'id':'synthetic:1','drug':'synthetic:d','indication':'synthetic:i'}))
            output=Path(tmp)/'result.json'
            with patch('sys.argv',['kg-audit','evidence','--record',str(record),'--as-of','2026-01-01','--output',str(output)]):cli.main()
            self.assertEqual(json.loads(output.read_text())['status'],'incomplete')

    def test_bridge_has_per_disease_rank_domain_and_no_inferred_regimen(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'scores.tsv';pd.DataFrame({'compound':['a','b','a'],'disease':['x','x','y'],'s':[.2,.9,.3]}).to_csv(path,sep='\t',index=False)
            out=Path(tmp)/'records.jsonl';report=export_records(path,'s',out,'synthetic','random','uniform',1,'2026-01-01')
            rows=[json.loads(x) for x in out.read_text().splitlines()]
            self.assertEqual([r['prediction']['rank'] for r in rows],[2,1,1])
            self.assertEqual([r['prediction']['candidate_count'] for r in rows],[2,2,1])
            self.assertTrue(all('regimen' not in r for r in rows));self.assertEqual(report['biomedical_fields_inferred'],0)
            schema=json.loads((Path(__file__).resolve().parents[1]/'schemas/handoff.schema.json').read_text())
            for row in rows:jsonschema.validate(row,schema)

    def test_registry_snapshot_is_bounded_cached_and_not_a_negative_label(self):
        def response(request,timeout=60):
            status=parse_qs(urlparse(request.full_url).query)['filter.overallStatus'][0]
            payload={'totalCount':1,'studies':[{'protocolSection':{'identificationModule':{'nctId':'SYNTHETIC_'+status},
                'statusModule':{'overallStatus':status,'whyStopped':'Funding ended','lastUpdatePostDateStruct':{'date':'2025-01-01'}},
                'armsInterventionsModule':{'interventions':[{'type':'DRUG','name':'synthetic drug'}]}}}]}
            return io.BytesIO(json.dumps(payload).encode())
        with tempfile.TemporaryDirectory() as tmp,patch('urllib.request.urlopen',side_effect=response) as call,patch('time.sleep'):
            data=fetch_sample(tmp,'2026-01-01',1);self.assertEqual(call.call_count,3)
            self.assertTrue(all(x['global_negative_edges_added']==0 for x in data['summary'].values()))
            fetch_sample(tmp,'2026-01-01',1);self.assertEqual(call.call_count,3)
            with self.assertRaises(ValueError):fetch_sample(tmp,'2026-01-02',1)

    def test_reason_cues_are_nonexclusive(self):
        self.assertEqual(reason_cues('Safety and recruitment concerns'),['operational','safety'])
        self.assertEqual(reason_cues('Insufficient efficacy'),['efficacy'])
        self.assertEqual(reason_cues(''),[])
