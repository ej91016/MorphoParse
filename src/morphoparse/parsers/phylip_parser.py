from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from morphoparse.utils import clean_warn, parse_taxon_line

def parse_phylip(file_path):
    with open(file_path, 'r') as f:
        ntax, seq_len = map(int, f.readline().strip().split())
        seq_data = {}
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, seq = parse_taxon_line(line)
            if name:
                seq_data[name] = seq_data.get(name, '') + seq

    if len(seq_data) != ntax:
        clean_warn(f"Expected {ntax} sequences, found {len(seq_data)}")
    lengths = {len(seq) for seq in seq_data.values()}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent sequence lengths: {lengths}")
    if seq_len not in lengths:
        clean_warn(f"Expected {seq_len} characters, found {list(lengths)[0]}")
    
    return [SeqRecord(Seq(seq), id=name) for name, seq in seq_data.items()]
