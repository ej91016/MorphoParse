from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from morphoparse.utils import clean_warn, parse_taxon_line
import re

def parse_nexus(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    ntax = nchar = None
    for line in lines:
        if ntax is None:
            ntax_match = re.search(r'ntax\s*=\s*(\d+)', line, re.I)
            if ntax_match:
                ntax = int(ntax_match.group(1))
        if nchar is None:
            nchar_match = re.search(r'nchar\s*=\s*(\d+)', line, re.I)
            if nchar_match:
                nchar = int(nchar_match.group(1))
        if ntax is not None and nchar is not None:
            break

    if ntax is None:
        clean_warn("Could not find ntax in the file.")         
    if nchar is None:
        clean_warn("Could not find nchar in the file.")

    matrix_start = next((i for i, l in enumerate(lines) if l.strip().lower() == 'matrix'), None)
    if matrix_start is None:
        raise ValueError("No MATRIX found")

    seq_data = {}
    for line in lines[matrix_start + 1:]:
        line = line.strip()
        if not line or line.startswith('['):
            continue
        if line == ';':
            break
        name, seq = parse_taxon_line(line)
        if name:
            seq_data[name] = seq_data.get(name, '') + seq

    if ntax is not None and len(seq_data) != ntax:
        clean_warn(f"Expected {ntax} taxa, found {len(seq_data)}")

    lengths = {len(seq) for seq in seq_data.values()}
    if nchar is not None:
        if len(lengths) != 1:
            raise ValueError(f"Inconsistent sequence lengths: {lengths}")
        if nchar not in lengths:
            clean_warn(f"Expected {nchar} characters, found {len(seq)}")

    return [SeqRecord(Seq(seq), id=name) for name, seq in seq_data.items()]