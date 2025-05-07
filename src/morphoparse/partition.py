from pathlib import Path
from Bio.Seq import Seq
from morphoparse.utils import statenum, statemax, collapse_ranges
from morphoparse.file_utils import write_file
from morphoparse.weight import write_paup_weights, write_tnt_weights
from morphoparse.poly import Poly
import math


def create_partition_file(records, poly, output_file, ver="raxml", asc_correction=False, out_format="phylip",
                          paup=False, tnt=False, max=False):

    state_groups = {}
    for pos in range(len(records[0].seq)):
        r = statemax(records, pos, poly) if max else statenum(records, pos, poly)
        state_groups.setdefault(r, []).append(pos + 1)

    suffix_map = {
        "raxml": "_raxml.models", 
        "raxmlng": "_raxmlng.models", 
        "iqtree": "_iqtree.nex"
    }
    suffix = suffix_map.get(ver, "_raxml.models")
    out_file = f"{Path(output_file)}{suffix}"
    
    partition_content = []
    iqtree_partitions = []
    paup_lines = []
    tnt_lines = []
    subset_data = {}

    for state, positions in sorted(state_groups.items()):
        ranges = collapse_ranges(positions)
        range_strs = [f"{s}-{e}" if s != e else str(s) for s, e in ranges]
        weight = round(10 * math.log(state)) if state != 0 else 0

        if ver != "iqtree":
            model_template = (
                "MULTI{0}_MK+G4+ASC_LEWIS, m{0}= " if ver == "raxmlng" and asc_correction else
                "MULTI{0}_MK+G4, m{0}= " if ver == "raxmlng" else
                "ASC_MULTI, m{0}= " if asc_correction else
                "MULTI, m{0}= "
            )
            partition_content.append(f"{model_template.format(state)}{','.join(range_strs)}")

        elif ver == "iqtree":
            cols = [pos - 1 for pos in positions]
            subset_records = []
            for rec in records:
                new_seq = "".join(rec.seq[i] for i in cols)
                rec_copy = rec[:]
                rec_copy.seq = Seq(new_seq)
                subset_records.append(rec_copy)
            subset_data[state] = subset_records
            out_ext = out_format if out_format == 'fasta' else out_format[:3]
            part_path = Path(f"{output_file}_m{state}").with_suffix('')
            iqtree_partitions.append(f"    charset m{state} = {part_path}.{out_ext}: *;")

        if paup:
            paup_lines.append(f"{weight}:{' '.join(range_strs)}")
        if tnt:
            tnt_lines.append(f"ccode /{weight or 1} {' '.join(str(p-1) for p in positions)};")

    print(f"SSA partition written to {Path(out_file).name}")
    with open(out_file, 'w') as f:
        if ver == "iqtree":
            f.write("#nexus\nbegin sets;\n")
            f.write("\n".join(iqtree_partitions) + "\n")
            f.write("    charpartition mine =\n")
            iq_template = "MK+ASC+G" if asc_correction else "MK+G"
            f.write(",\n".join(f"      {iq_template}: m{state}" for state in sorted(state_groups)) + ";\n")
            f.write("end;\n")
            poly_all = Poly(poly)
            for state, records in subset_data.items():
                subset_indices = [pos - 1 for pos in state_groups[state]]
                subset_poly = poly_all.subset(subset_indices) if poly else None
                write_file(records, f"{output_file}_m{state}", out_format, subset_poly)
        else:
            f.write("\n".join(partition_content) + "\n")

    if paup and paup_lines:
        write_paup_weights(output_file, paup_lines)
    if tnt and tnt_lines:
        write_tnt_weights(output_file, tnt_lines)