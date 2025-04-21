from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from morphoparse.utils import clean_sequence, normalize_taxon_name

def parse_fasta(file_path):
    records, name, seq_lines = [], None, []
    expected_length = first_taxon = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if name is not None:
                    record = _finalize_record(name, seq_lines, expected_length, first_taxon)
                    expected_length = expected_length or len(record.seq)
                    first_taxon = first_taxon or name
                    records.append(record)
                name = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)

        if name:
            record = _finalize_record(name, seq_lines, expected_length, first_taxon)
            records.append(record)

    return records

def _finalize_record(name, seq_lines, expected_length, first_taxon):
    full_seq = ''.join(seq_lines)
    cleaned_seq = clean_sequence(full_seq)
    seq_len = len(cleaned_seq)

    if expected_length is not None and seq_len != expected_length:
        raise ValueError(
            f"Inconsistent sequence lengths: expected {expected_length} (from {first_taxon}), got {seq_len} (from {name})"
        )

    return SeqRecord(Seq(cleaned_seq), id=normalize_taxon_name(name))
