import re
import warnings
def normalize_taxon_name(name):
    name = name.strip("\"'")
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)

def clean_sequence(sequence):
    return re.sub(r'[\(\[\{].*?[\)\]\}]', '?', sequence)

def clean_seq_keep_poly(sequence):
    clean_seq = clean_sequence(sequence)
    char_index = 0
    poly_dict = {}

    i = 0
    while i < len(sequence):
        if sequence[i] in '{[(':
            closing = {'{': '}', '[': ']', '(': ')'}[sequence[i]]
            j = i + 1
            while j < len(sequence) and sequence[j] != closing:
                j += 1
            if j == len(sequence):
                raise ValueError("Unmatched bracket in sequence")

            states = set(c for c in sequence[i+1:j] if c.isdigit())
            if states:
                poly_dict[char_index] = states
            i = j + 1
        else:
            i += 1
        char_index += 1

    return clean_seq, poly_dict


def statenum(records, position, poly=None):
    states = set()
    for r in records:
        c = r.seq[position].upper()


        if poly and r.id in poly and position in poly[r.id]:
            states.update(poly[r.id][position])
        

        elif c not in '-?':
            states.add(c)

    return len(states)


def statemax(records, position, poly=None):
    states = set()
    for r in records:
        c = r.seq[position].upper()

        if poly and r.id in poly and position in poly[r.id]:
            states.update({s for s in poly[r.id][position] if s.isdigit()})


        elif c not in '-?' and c.isdigit():
            states.add(c)

    return max(map(int, states)) + 1 if states else 0


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
    raw_seq = seq_part.replace(' ', '').replace('\t', '')
    return name, raw_seq

def clean_warn(message, category=UserWarning):
    original_format = warnings.formatwarning
    warnings.formatwarning = lambda msg, cat, *_: f"⚠️  {cat.__name__}: {msg}\n"
    warnings.warn(message, category)
    warnings.formatwarning = original_format
