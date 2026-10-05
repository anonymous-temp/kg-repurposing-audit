"""Validate local JSONL exports without printing record contents."""
import argparse
from collections import Counter
import json
from pathlib import Path
import jsonschema
from kg_audit.evidence import assess_strategy


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('records',type=Path);p.add_argument('--as-of',required=True)
    a=p.parse_args();schema=json.loads((Path(__file__).resolve().parents[1]/'schemas/handoff.schema.json').read_text())
    validator=jsonschema.Draft7Validator(schema,format_checker=jsonschema.FormatChecker())
    status=Counter();errors=0;count=0
    for line in a.records.open():
        record=json.loads(line);errors+=sum(1 for _ in validator.iter_errors(record));count+=1
        status[assess_strategy(record,a.as_of)['status']]+=1
    print(json.dumps({'records':count,'structural_errors':errors,'documentation_status':dict(status)},indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
