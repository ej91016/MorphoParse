from pathlib import Path

def write_nexus(records, output_file, gap='-', missing='?'):
    padding = max(len(r.id) for r in records) + 2
    nchar = len(records[0].seq)
    ntax = len(records)
    symbols = sorted({c for r in records for c in str(r.seq) if c not in {gap, missing}})
    with open(f"{output_file}.nex", 'w') as f:
        f.write("#NEXUS\n\nBEGIN DATA;\n")
        f.write(f"  DIMENSIONS NTAX={ntax} NCHAR={nchar};\n")
        f.write(f"  FORMAT DATATYPE=STANDARD GAP={gap} MISSING={missing} SYMBOLS=\"{''.join(symbols)}\";\n")
        f.write("  MATRIX\n")
        for r in records:
            f.write(f"  {r.id.ljust(padding)}{r.seq}\n")
        f.write("  ;\nEND;\n")
    print(f"Created NEXUS file: {Path(output_file).name}.nex")