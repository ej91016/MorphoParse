from pathlib import Path

def write_tnt(records, output_file, poly=None):
    nchar = len(records[0].seq)
    ntax = len(records)
    padding = max(len(r.id) for r in records) + 2

    if poly:
        out_seqs = []
        for r in records:
            seq = list(str(r.seq))
            if r.id in poly:
                for i, states in poly[r.id].items():
                    seq[i] = '[' + ''.join(sorted(states)) + ']'
            out_seqs.append(''.join(seq))
    else:
        out_seqs = [str(r.seq) for r in records]

    with open(f"{output_file}.tnt", 'w') as f:
        f.write("xread\n")
        f.write(f"'{'Produced with morphoparse'}'\n")
        f.write(f"{nchar} {ntax}\n")
        for r, seq in zip(records, out_seqs):
            f.write(f"{r.id.ljust(padding)} {seq}\n")
        f.write(";\n")
    print(f"Created TNT file: {Path(output_file).name}.tnt")