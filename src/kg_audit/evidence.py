"""Conservative documentation rules, not clinical eligibility decisions.

No trial registry status or workflow failure creates a global negative indication.
Dates constrain the information available at a requested review time.
"""
from datetime import date
import math

DOMAINS = ('mechanism','systemic_exposure','tissue_exposure','regimen',
           'population','comparator_endpoint','safety','evidence_path')
STATUSES = {'pass','concern','fail','unknown'}
UNKNOWN = {'','unknown','not reported','not available','unspecified','n/a','none'}


def known(value):
    text=str(value).strip().lower()
    return value is not None and text not in UNKNOWN and not text.startswith(('unknown ', 'not reported ', 'not available '))


def iso_date(value):
    return date.fromisoformat(value) if isinstance(value,str) and value else None


def validate_record(record):
    """Return structural errors. Missing scientific evidence is not fabricated."""
    errors=[]
    if not isinstance(record,dict):return ['record must be an object']
    for key in ('id','drug','indication'):
        if not known(record.get(key)):errors.append(f'missing {key}')
    for key in ('evidence','constraints','failures'):
        if not isinstance(record.get(key,[]),list):errors.append(f'{key} must be a list')
    if errors:return errors
    evidence=record.get('evidence',[])
    if any(not isinstance(x,dict) for x in evidence):return ['evidence entries must be objects']
    ids=[x.get('id') for x in evidence]
    if None in ids or len(set(ids))!=len(ids):errors.append('evidence identifiers must be present and unique')
    for entry in evidence:
        if not known(entry.get('source')):errors.append('evidence source missing')
        if entry.get('source_date'):
            try:iso_date(entry['source_date'])
            except ValueError:errors.append('invalid evidence source_date')
        effect=entry.get('effect')
        if effect is not None:
            if not isinstance(effect,dict):errors.append('effect must be an object');continue
            values=[effect.get(k) for k in ('estimate','lower','upper')]
            if any(not isinstance(v,(int,float)) or isinstance(v,bool) or not math.isfinite(v) for v in values):
                errors.append('effect and interval must be finite numbers')
            elif values[1]>values[2]:errors.append('effect interval bounds are reversed')
    for group in ('constraints','failures'):
        for item in record.get(group,[]):
            if not isinstance(item,dict):errors.append(f'{group} entry must be an object');continue
            refs=item.get('evidence_ids',[])
            if not isinstance(refs,list) or any(x not in ids for x in refs):errors.append('orphan evidence reference')
            if group=='constraints':
                if item.get('domain') not in DOMAINS:errors.append('unknown constraint domain')
                if item.get('status') not in STATUSES:errors.append('invalid constraint status')
            if item.get('source_date'):
                try:iso_date(item['source_date'])
                except ValueError:errors.append('invalid failure source_date')
    return errors


def writeback_action(failure):
    """Return a non-destructive evidence action; a flat negative is never inferred."""
    implication=failure.get('implication','defer')
    if implication=='retain':
        return 'retain_existing_evidence' if failure.get('evidence_ids') else 'defer'
    if implication=='qualify':
        return 'store_context_qualification' if failure.get('scope') and failure.get('evidence_ids') else 'defer'
    if implication=='negate':
        scope=failure.get('scope') or {}
        scoped=all(known(scope.get(k)) for k in ('drug','indication','population','regimen','comparator','endpoint'))
        if failure.get('domain')=='efficacy' and scoped and failure.get('evidence_ids') and failure.get('source_date'):
            return 'store_scoped_counterevidence'
    return 'defer'


def assess_strategy(record,as_of):
    """Assess completeness and provenance as of an ISO date.

    A 'ready' result means a complete record for evidence review, not an actionable
    prescription, therapeutic efficacy, regulatory approval or expert consensus.
    """
    cutoff=iso_date(as_of)
    if cutoff is None:raise ValueError('as_of must be an ISO date')
    errors=validate_record(record)
    if errors:return {'status':'invalid_record','reasons':errors,'as_of':as_of}
    if record.get('assessment_date'):
        try:assessment_date=iso_date(record['assessment_date'])
        except ValueError:return {'status':'invalid_record','reasons':['invalid assessment_date'],'as_of':as_of}
        if assessment_date>cutoff:return {'status':'not_yet_assessed','reasons':['assessment occurs after requested review date'],'as_of':as_of}
    evidence={x['id']:x for x in record.get('evidence',[])}
    def available(refs):
        return bool(refs) and all(iso_date(evidence[i].get('source_date')) is not None
                                  and iso_date(evidence[i]['source_date'])<=cutoff for i in refs)
    current=[x for x in record.get('constraints',[]) if available(x.get('evidence_ids',[]))]
    for domain in DOMAINS:
        values={x['status'] for x in current if x['domain']==domain}
        if len(values)>1:return {'status':'conflicting_records','reasons':[domain],'as_of':as_of}
    for failure in record.get('failures',[]):
        day=iso_date(failure.get('source_date'))
        if day is not None and day<=cutoff and available(failure.get('evidence_ids',[])):
            action=writeback_action(failure)
            if action=='store_scoped_counterevidence':
                scope=failure['scope']
                if all(scope.get(k)==record.get(k) for k in ('drug','indication','population','regimen','comparator','endpoint')):
                    return {'status':'scoped_counterevidence','reasons':['counterevidence retained with its regimen and population scope'],
                            'as_of':as_of,'scope':scope}
    missing=[k for k in ('population','comparator','regimen','endpoint','required_exposure') if not known(record.get(k))]
    observed={x['domain']:x['status'] for x in current}
    missing += [domain for domain in DOMAINS if domain not in observed or observed[domain]=='unknown']
    if missing:return {'status':'incomplete','reasons':sorted(set(missing)),'as_of':as_of}
    pending=[domain for domain,status in observed.items() if status!='pass']
    if pending:return {'status':'requires_evidence_review','reasons':pending,'as_of':as_of}
    return {'status':'ready_for_evidence_review','reasons':[],'as_of':as_of}
