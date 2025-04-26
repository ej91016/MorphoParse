from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from morphoparse.utils import clean_sequence, normalize_taxon_name, parse_taxon_line

def parse_phylip(file_path):
    with open(file_path, 'r') as f:
        num_seqs, seq_len = map(int, f.readline().strip().split())
        seq_data = {}
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, seq = parse_taxon_line(line)
            if name:
                seq_data[name] = seq_data.get(name, '') + seq

    if len(seq_data) != num_seqs:
        raise ValueError(f"Expected {num_seqs} sequences, found {len(seq_data)}")
    lengths = {len(seq) for seq in seq_data.values()}
    if len(lengths) != 1 or seq_len not in lengths:
        raise ValueError("Sequence lengths are inconsistent or don't match header")
    return [SeqRecord(Seq(seq), id=name) for name, seq in seq_data.items()]
