"""Shared helpers for the supplementary-figure scripts b_*.py: markdown-table parser, number parsing, Wilson interval."""
import math, os, re, sys

MS = os.environ.get("SUPP_MD", "paper_results/supp_tables_md")                     # generated supplementary tables (markdown)
OUT = os.environ.get("FIG_OUT", "work/figures_supp")
SEL = os.environ.get("SELECTION", "paper_results/selection")
CASES = os.environ.get("CASES", "paper_results/case_review/case_taxonomy_verified.tsv")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import figstyle as S  # noqa: E402

os.makedirs(OUT, exist_ok=True)
NUM = re.compile(r"-?\d+(?:\.\d+)?")


def num(s):
    """First number in a cell: handles thousands separators and the Unicode minus."""
    t = s.replace("−", "-").replace(",", "")
    m = NUM.search(t)
    if not m:
        raise ValueError(f"no number in {s!r}")
    return float(m.group())


def nums(s):
    return [float(x) for x in NUM.findall(s.replace("−", "-").replace(",", ""))]


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def read_table(path, caption_regex):
    """Return (header, rows) of the first pipe table following the line that matches caption_regex at its start
    (a leading ** is allowed)."""
    lines = open(path, encoding="utf-8").read().splitlines()
    pat = re.compile(r"^\**" + caption_regex)
    i = next(k for k, l in enumerate(lines) if pat.search(l))
    j = i + 1
    while not lines[j].startswith("|"):
        j += 1
    header = cells(lines[j])
    j += 2
    rows = []
    while j < len(lines) and lines[j].startswith("|"):
        rows.append(cells(lines[j]))
        j += 1
    return header, rows


def wilson(k, n, z=1.959964):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def lines_after(path, caption_regex, n=3):
    """The first n non-empty text lines after a table (used for sentences under a table)."""
    lines = open(path, encoding="utf-8").read().splitlines()
    pat = re.compile(r"^\**" + caption_regex)
    i = next(k for k, l in enumerate(lines) if pat.search(l))
    j = i + 1
    while not lines[j].startswith("|"):
        j += 1
    while j < len(lines) and lines[j].startswith("|"):
        j += 1
    out = [l for l in lines[j:] if l.strip()]
    return out[:n]
