from pathlib import Path

def write_nexus(records, output_file, poly=None, gap='-', missing='?'):
    padding = max(len(r.id) for r in records) + 2
    nchar = len(records[0].seq)
    ntax = len(records)

    if poly:
        out_seqs = []
        for r in records:
            seq = list(str(r.seq))
            if r.id in poly:
                for i, states in poly[r.id].items():
                    seq[i] = '{' + ''.join(sorted(states)) + '}'
            out_seqs.append(''.join(seq))
    else:
        out_seqs = [str(r.seq) for r in records]

    symbols = sorted({c for s in out_seqs for c in s if c not in {gap, missing} and not c in '{}'})
    
    with open(f"{output_file}.nex", 'w') as f:
        f.write("#NEXUS\n\nBEGIN DATA;\n")
        f.write(f"  DIMENSIONS NTAX={ntax} NCHAR={nchar};\n")
        f.write(f"  FORMAT DATATYPE=STANDARD GAP={gap} MISSING={missing} SYMBOLS=\"{''.join(symbols)}\";\n")
        f.write("  MATRIX\n")
        for r, seq in zip(records, out_seqs):
            f.write(f"  {r.id.ljust(padding)}{seq}\n")
        f.write("  ;\nEND;\n")
    print(f"Created NEXUS file: {Path(output_file).name}.nex")
