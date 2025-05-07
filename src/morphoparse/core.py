from pathlib import Path
from Bio.Seq import Seq
from morphoparse.parsers import fasta_parser, phylip_parser, nexus_parser, tnt_parser
from morphoparse.writers import fasta_writer, phylip_writer, nexus_writer, tnt_writer
from morphoparse.utils import statenum, statemax, collapse_ranges, detect_format
from morphoparse.remap import remap_sparse_states
from morphoparse.weight import write_paup_weights, write_tnt_weights
from morphoparse.poly import Poly
import math
import os

def parse_sequences(input, file_format, keep):
    if input.stat().st_size == 0:
        raise ValueError("Empty file")    
    parsers = {
        'fasta': fasta_parser.parse_fasta,
        'phylip': phylip_parser.parse_phylip,
        'nexus': nexus_parser.parse_nexus,
        'tnt': tnt_parser.parse_tnt,
    }
    return parsers[file_format](input, keep)

def write_file(records, output, out_format, poly):
    writers = {
        'fasta': fasta_writer.write_fasta,
        'phylip': phylip_writer.write_phylip,
        'nexus': nexus_writer.write_nexus,
        'tnt': tnt_writer.write_tnt,
    }
    writers[out_format](records, output, poly)

def create_partition_file(records, poly, output_file, ver="raxml", asc_correction=False, out_format="phylip",
                          paup=False, tnt=False, max=False, ):

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

def run_pipeline(args):
    input = Path(args.input)
    if not input.exists():
        raise FileNotFoundError(f"Input file not found: {input.absolute()}")
    format = detect_format(input,args.format)
    print(f"MorphoPase will parse {input.name} in {format.upper()} format")
    out_format = (args.out_format or format).lower()
    keep_poly=args.keep_poly
    records, poly_data = parse_sequences(input, format, keep_poly)
    output = args.output or input.with_suffix('')
    remove_missing = args.remap or args.remove_missing
    remove_uninformative = args.remap or args.remove_mono
    reorder_states = args.remap or args.reorder
    if remove_missing or remove_uninformative or reorder_states:
        logfile = f"{output}_remap.txt"
        records, poly_data = remap_sparse_states(
            records, poly_data, logfile,
            remove_missing,
            remove_uninformative,
            reorder_states,
            
        )

    write_file(records, f"{output}_mparse", out_format, poly_data)
    create_partition_file(
        records,
        poly_data,
        output,
        ver=args.software,
        asc_correction=args.asc,
        paup=args.paup,
        tnt=args.tnt,
        max=not reorder_states,
        out_format=out_format
    )

#this needs update
def run_morphoparse(
    input_file,
    output_prefix,
    input_format='phylip',
    output_format=None,
    remap=False,
    remove_missing=False,
    remove_mono=False,
    reorder=False,
    asc=False,
    software=False,
    paup=False,
    tnt=False
):
    output_format = (output_format or input_format).lower()
    records = parse_sequences(input_file, input_format)
    remove_missing = remap or remove_missing
    remove_uninformative = remap or remove_mono
    reorder_states = remap or reorder

    if remove_missing or remove_mono or reorder:
        logfile = f"{output_prefix}_remap.log"
        records = remap_sparse_states(
            records,
            logfile,
            remove_missing,
            remove_uninformative,
            reorder_states
        )

    write_file(records, f"{output_prefix}_mparse", output_format)
    create_partition_file(
        records,
        output_prefix,
        ver=software,
        asc_correction=asc,
        paup=paup,
        tnt=tnt,
        max=not (remap or reorder),
        out_format=output_format
    )