import re
import warnings
def normalize_taxon_name(name):
    name = name.strip("\"'")
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)

def clean_sequence(sequence):
    return re.sub(r'[\(\[\{].*?[\)\]\}]', '?', sequence)

def statenum(records, position):
    return len({r.seq[position].upper() for r in records if r.seq[position] not in '-?'})

def statemax(records, position):
    states = {
        int(c) for r in records 
        if (c := r.seq[position].upper()) not in '-?'
        and c.isdigit()
    }
    return max(states) + 1 if states else 0

def collapse_ranges(positions):
    positions = sorted(positions)
    ranges = []
    start = prev = positions[0]
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            ranges.append((start, prev))
            start = prev = pos
    ranges.append((start, prev))
    return ranges

def parse_taxon_line(line):

    if line[0] in {'"', "'"}:
        match = re.match(r'(["\'])(.*?)\1\s+(.*)', line)
        raw_name, seq_part = match.groups()[1:] if match else (None, None)
    else:
        match = re.match(r'(\S+)\s+(.*)', line)
        raw_name, seq_part = match.groups() if match else (None, None)

    if raw_name is None:
        return None, None

    name = normalize_taxon_name(raw_name)
    seq = clean_sequence(seq_part.replace(' ', '').replace('\t', ''))
    return name, seq

def clean_warn(message, category=UserWarning):
    original_format = warnings.formatwarning
    warnings.formatwarning = lambda msg, cat, *_: f"⚠️  {cat.__name__}: {msg}\n"
    warnings.warn(message, category)
    warnings.formatwarning = original_format
