"""
Orchestrates high-level logic: parsing, remapping, writing, partitioning.
"""
from pathlib import Path
from morphoparse.parsers import fasta_parser, phylip_parser, nexus_parser, tnt_parser
from morphoparse.writers import fasta_writer, phylip_writer, nexus_writer, tnt_writer
from morphoparse.utils import statenum, statemax, collapse_ranges
from morphoparse.remap import remap_sparse_states
import math
import os

def parse_sequences(file_path, file_format):
    if os.path.getsize(file_path) == 0:
        raise ValueError("Empty file")    
    parsers = {
        'fasta': fasta_parser.parse_fasta,
        'phylip': phylip_parser.parse_phylip,
        'nexus': nexus_parser.parse_nexus,
        'tnt': tnt_parser.parse_tnt,
    }
    return parsers[file_format](file_path)

def write_file(records, output, out_format):
    writers = {
        'fasta': fasta_writer.write_fasta,
        'phylip': phylip_writer.write_phylip,
        'nexus': nexus_writer.write_nexus,
        'tnt': tnt_writer.write_tnt,
    }
    writers[out_format](records, output)

def create_partition_file(records, output_file, ver="raxml", asc_correction=False, paup=False, tnt=False, max=False, out_format="phylip"):
    state_groups = {}
    for pos in range(len(records[0].seq)):
        r = statemax(records, pos) if max else statenum(records, pos)
        state_groups.setdefault(r, []).append(pos + 1)

    suffix_map = {"raxml": "_raxml.models", "raxmlng": "_raxmlng.models", "iqtree": "_iqtree.nex"}
    suffix = suffix_map.get(ver, "_raxml.models")
    out_file = f"{Path(output_file)}{suffix}"
    print(f"SSA partition written to {Path(out_file).name}")

    paup_lines, tnt_lines = [], []

    if ver == "iqtree":
        part_lines = []
        out_ext = out_format if out_format=='fasta' else out_format[:3]

        for state, positions in sorted(state_groups.items()):
            cols = [pos - 1 for pos in positions]
            subset_records = []
            for rec in records:
                new_seq = "".join(rec.seq[i] for i in cols)
                rec_copy = rec[:]
                rec_copy.seq = new_seq
                subset_records.append(rec_copy)

            part_path = Path(f"{output_file}_m{state}").with_suffix('')
            write_file(subset_records, part_path, out_format)
            part_lines.append(f"    charset m{state} = {part_path}.{out_ext}: *;")

        # Write NEXUS partition file
        with open(out_file, "w") as f:
            iq_template=(
                "MK+ASC+G" if asc_correction else "MK+G")
            f.write("#nexus\nbegin sets;\n")
            for line in part_lines:
                f.write(f"{line}\n")
            f.write(f"    charpartition  mine =\n")
            for i, state in enumerate(sorted(state_groups), 1):
                lend = ',' if i<len(state_groups) else ';' 
                f.write(f"      {iq_template}: m{state}{lend}\n")
            f.write("end;\n")

    else:
        # RAxML/RAxML-NG style
        model_template = (
            "MULTI{0}_MK+G4+ASC_LEWIS, m{0}= " if ver == "raxmlng" and asc_correction else
            "MULTI{0}_MK+G4, m{0}= " if ver == "raxmlng" else
            "ASC_MULTI, m{0}= " if asc_correction else
            "MULTI, m{0}= "
        )

        with open(out_file, 'w') as f:
            for state, positions in sorted(state_groups.items()):
                ranges = collapse_ranges(positions)
                range_strs = [f"{s}-{e}" if s != e else str(s) for s, e in ranges]
                f.write(f"{model_template.format(state)}{','.join(range_strs)}\n")

                weight = round(10 * math.log(state)) if state != 0 else 0
                if paup:
                    paup_lines.append(f"{weight}:{' '.join(range_strs)}")
                if tnt:
                    tnt_lines.append(f"ccode /{weight or 1} {' '.join(str(p-1) for p in positions)};")

    if paup and paup_lines:
        paup_file = f"{Path(output_file)}_paup_weights.txt"
        with open(paup_file, 'w') as f:
            f.write("BEGIN PAUP;\n")
            for line in paup_lines:
                weight, rest = line.split(":")
                positions = rest.strip().split()
                prefix = f"  WEIGHTS {weight}:"
                current_line, current_length = [prefix], len(prefix)
                for pos in positions:
                    next_piece = f" {pos}"
                    if current_length + len(next_piece) + 1 > 250:
                        f.write("".join(current_line) + ";\n")
                        current_line, current_length = [prefix, f" {pos}"], len(prefix) + len(f" {pos}")
                    else:
                        current_line.append(f" {pos}")
                        current_length += len(next_piece)
                f.write("".join(current_line) + ";\n")
            f.write("END;\n")
        print(f"SSA weight scheme written to {paup_file}")

    if tnt and tnt_lines:
        tnt_file = f"{Path(output_file)}_tnt_weights.txt"
        with open(tnt_file, 'w') as f:
            for line in tnt_lines:
                f.write(line + "\n")
        print(f"SSA weight scheme written to {tnt_file}")


def run_pipeline(args):
    args.out_format = (args.out_format or args.format).lower()
    records = parse_sequences(args.input, args.format)
    args.output = args.output or Path(args.input).with_suffix('')
    remove_missing = args.remap or args.remove_missing
    remove_uninformative = args.remap or args.remove_mono
    reorder_states = args.remap or args.reorder

    if remove_missing or remove_uninformative or reorder_states:
        logfile = f"{args.output}_remap.txt"
        records = remap_sparse_states(
            records, logfile,
            remove_missing,
            remove_uninformative,
            reorder_states
        )

    write_file(records, f"{args.output}_clean", args.out_format)
    create_partition_file(
        records,
        args.output,
        ver=args.software,
        asc_correction=args.asc,
        paup=args.paup,
        tnt=args.tnt,
        max=not reorder_states,
        out_format=args.out_format
    )

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

    write_file(records, f"{output_prefix}_clean", output_format)
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