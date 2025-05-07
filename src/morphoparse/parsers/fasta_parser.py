from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from morphoparse.utils import clean_sequence, clean_seq_keep_poly, normalize_taxon_name 

def parse_fasta(file_path, keep=False):
    records, name, seq_lines = [], None, []
    expected_length = first_taxon = None
    poly_data = {} if keep else None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if name is not None:
                    record, poly = _finalize_record(name, seq_lines, expected_length, first_taxon, keep)
                    expected_length = expected_length or len(record.seq)
                    first_taxon = first_taxon or name
                    records.append(record)
                    if keep:
                        poly_data[name] = poly
                name = normalize_taxon_name(line[1:].strip())
                seq_lines = []
            else:
                seq_lines.append(line)

        if name:
            record, poly = _finalize_record(name, seq_lines, expected_length, first_taxon, keep)
            records.append(record)
            if keep:
                poly_data[name] = poly
    return records, poly_data

def _finalize_record(name, seq_lines, expected_length, first_taxon, keep=False):
    full_seq = ''.join(seq_lines)
    if keep:
        cleaned_seq, poly = clean_seq_keep_poly(full_seq)
    else:
        cleaned_seq = clean_sequence(full_seq)
        poly = None

    seq_len = len(cleaned_seq)
    if expected_length is not None and seq_len != expected_length:
        raise ValueError(
            f"Inconsistent sequence lengths: expected {expected_length} (from {first_taxon}), got {seq_len} (from {name})"
        )

    record = SeqRecord(Seq(cleaned_seq), id=name)
    return record, poly

