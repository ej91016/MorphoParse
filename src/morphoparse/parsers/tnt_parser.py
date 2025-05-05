from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from morphoparse.utils import clean_warn, parse_taxon_line
import re

def parse_tnt(file_path):
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    seq_data = {}
    ntax=seq_len=None
    try:
        i = next(i for i, line in enumerate(lines) if line.lower().startswith("xread"))
    except StopIteration:
        raise ValueError("No 'xread' statement found.")
    i += 1

    if i < len(lines):
        line = lines[i]
        if line[0] in {"'", '"'}:
            quote = line[0]
            if line.endswith(quote) and line.count(quote) % 2 == 0:
                i += 1
            else:
                i += 1
                while i < len(lines):
                    if lines[i].endswith(quote) and lines[i].count(quote) % 2 == 1:
                        i += 1
                        break
                    i += 1

    while i < len(lines) and not re.match(r'^\d+\s+\d+$', lines[i]):
        i += 1
    if i < len(lines) and re.match(r'^\d+\s+\d+$', lines[i]):
        seq_len, ntax = map(int, lines[i].split())
        i += 1

    while i < len(lines):
        line = lines[i]
        if line == ';':
            break
        if line.startswith('['):
            i += 1
            continue
        name, seq = parse_taxon_line(line)
        if name:
            seq_data[name] = seq_data.get(name, '') + seq
        i += 1

    if not seq_data:
        raise ValueError("No sequence data found in TNT file, make sure nchar and ntaxa are specified")
    if len(seq_data) != ntax:
        clean_warn(f"Expected {ntax} sequences, found {len(seq_data)}")

    lengths = {len(seq) for seq in seq_data.values()}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent sequence lengths: {lengths}")
    if seq_len not in lengths:
        clean_warn(f"Expected {seq_len} characters, found {list(lengths)[0]}")

    return [SeqRecord(Seq(seq), id=name) for name, seq in seq_data.items()]
