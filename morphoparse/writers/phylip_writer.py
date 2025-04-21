from pathlib import Path

def write_phylip(records, output_file):
    padding = max(len(r.id) for r in records) + 2
    with open(f"{output_file}.phy", 'w') as f:
        f.write(f"{len(records)} {len(records[0].seq)}\n")
        for r in records:
            f.write(f"{r.id.ljust(padding)}{r.seq}\n")
    print(f"Created PHYLIP file: {Path(output_file).name}.phy")