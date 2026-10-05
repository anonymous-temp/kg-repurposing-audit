"""Fail if the tracked release tree contains databases, credentials or local paths."""
from pathlib import Path
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
DENIED_SUFFIXES={'.duckdb','.db','.sqlite','.sqlite3','.parquet','.feather','.pt','.pth','.npz','.gz','.zip','.docx'}
PATTERNS={
 'private absolute path':re.compile(r'/Users/[A-Za-z0-9_.-]+/|/home/[A-Za-z0-9_.-]+/|[A-Z]:\\Users\\'),
 'private key':re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
 'token':re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{25,}|sk-(?:proj-)?[A-Za-z0-9_-]{30,}|AKIA[A-Z0-9]{16})\b'),
}


def main():
    git=subprocess.run(['git','ls-files','-z'],cwd=ROOT,capture_output=True,check=False)
    files=[ROOT/part.decode() for part in git.stdout.split(b'\0') if part] if git.returncode==0 else []
    if not files:
        files=[p for p in ROOT.rglob('*') if p.is_file() and not any(x in {'.git','.venv','__pycache__','.pytest_cache','.ruff_cache'} for x in p.parts)]
    errors=[]
    for path in files:
        rel=path.relative_to(ROOT)
        if path.is_symlink():errors.append(f'{rel}: symlink');continue
        if path.suffix.lower() in DENIED_SUFFIXES:errors.append(f'{rel}: excluded binary/database type');continue
        if any(x.lower() in {'data','datasets','outputs','checkpoints'} for x in rel.parts):errors.append(f'{rel}: data/output directory')
        if path.stat().st_size>2_000_000:errors.append(f'{rel}: unexpectedly large file')
        if path.name.startswith('.env'):errors.append(f'{rel}: environment file')
        text=path.read_text(errors='replace')
        for label,pattern in PATTERNS.items():
            if pattern.search(text):errors.append(f'{rel}: {label}')
        # The only tabular tracked files are aggregate publication summaries.
        if path.suffix in {'.csv','.tsv'} and rel.parts[0]!='paper_results':errors.append(f'{rel}: unexpected tabular data')
    if errors:
        print('\n'.join(errors));return 1
    print(f'PASS: {len(files)} source/metadata files checked; no excluded database, model, document, credential or local-path artifact found.')
    return 0


if __name__=='__main__':sys.exit(main())
