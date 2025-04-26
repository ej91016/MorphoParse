from pathlib import Path

def write_fasta(records, output_file):
    with open(f"{output_file}.fasta", 'w') as f:
        for record in records:
            f.write(f">{record.id}\n{record.seq}\n")
    print(f"Created FASTA file: {Path(output_file).name}.fasta")