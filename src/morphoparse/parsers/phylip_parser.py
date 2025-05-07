from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from morphoparse.utils import clean_warn, parse_taxon_line, clean_sequence, clean_seq_keep_poly

def parse_phylip(file_path, keep=False):
    with open(file_path, 'r') as f:
        ntax, seq_len = map(int, f.readline().strip().split())
        seq_data = {}
        poly_data = {} if keep else None

        for line in f:
            line = line.strip()
            if not line:
                continue

            name, raw_seq = parse_taxon_line(line)
            if name:
                prior_len = len(seq_data.get(name, ''))

                if keep:
                    clean_seq, poly = clean_seq_keep_poly(raw_seq)
                    if poly:
                        for k, v in poly.items():
                            poly_data.setdefault(name, {})[prior_len + k] = v
                else:
                    clean_seq = clean_sequence(raw_seq)

                seq_data[name] = seq_data.get(name, '') + clean_seq

    if len(seq_data) != ntax:
        clean_warn(f"Expected {ntax} sequences, found {len(seq_data)}")

    lengths = {len(seq) for seq in seq_data.values()}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent sequence lengths: {lengths}")
    if seq_len not in lengths:
        clean_warn(f"Expected {seq_len} characters, found {list(lengths)[0]}")

    seq_records = [SeqRecord(Seq(seq), id=name) for name, seq in seq_data.items()]
    return seq_records, poly_data
