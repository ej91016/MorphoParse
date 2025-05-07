from pathlib import Path

def write_fasta(records, output_file, poly=None):
    with open(f"{output_file}.fasta", 'w') as f:
        for record in records:
            seq = list(str(record.seq))

            if poly and record.id in poly:
                for i, states in poly[record.id].items():
                    seq[i] = "{" + "".join(sorted(states)) + "}"

            f.write(f">{record.id}\n{''.join(seq)}\n")

    print(f"Created FASTA file: {Path(output_file).name}.fasta")