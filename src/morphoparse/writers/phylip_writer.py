from pathlib import Path

def write_phylip(records, output_file, poly=None):
    padding = max(len(r.id) for r in records) + 2
    num_taxa = len(records)
    seq_len = len(records[0].seq)

    with open(f"{output_file}.phy", 'w') as f:
        f.write(f"{num_taxa} {seq_len}\n")

        for r in records:
            name = r.id.ljust(padding)
            seq = list(str(r.seq))

            if poly and r.id in poly:
                for i, states in poly[r.id].items():
                    seq[i] = "{" + "".join(sorted(states)) + "}"

            f.write(f"{name}{''.join(seq)}\n")

    print(f"Created PHYLIP file: {Path(output_file).name}.phy")