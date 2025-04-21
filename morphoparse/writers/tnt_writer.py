from pathlib import Path

def write_tnt(records, output_file):
    nchar = len(records[0].seq)
    ntax = len(records)
    padding = max(len(r.id) for r in records) + 2
    with open(f"{output_file}.tnt", 'w') as f:
        f.write("xread\n")
        f.write(f"'{'Produced with morphoparse'}'\n")
        f.write(f"{nchar} {ntax}\n")
        for r in records:
            f.write(f"{r.id.ljust(padding)} {r.seq}\n")
        f.write(";\n")
    print(f"Created TNT file: {Path(output_file).name}.tnt")