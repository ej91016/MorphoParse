from pathlib import Path
from morphoparse.utils import detect_format
from morphoparse.file_utils import parse_sequences, write_file
from morphoparse.remap import remap_sparse_states
from morphoparse.partition import create_partition_file

def run_pipeline(args):
    input = Path(args.input)
    if not input.exists():
        raise FileNotFoundError(f"Input file not found: {input.absolute()}")
    format = detect_format(input,args.format)
    print(f"MorphoPase will parse {input.name} in {format.upper()} format")
    out_format = (args.out_format or format).lower()
    keep_poly=args.poly
    records, poly_data = parse_sequences(input, format, keep_poly)
    output = args.output or input.with_suffix('')
    remove_missing = args.remap or args.remove_missing
    remove_mono = args.remap or args.remove_mono
    reorder = args.remap or args.reorder
    if remove_missing or remove_mono or reorder:
        logfile = f"{output}_remap.txt"
        records, poly_data = remap_sparse_states(
            records,
            poly_data,
            logfile,
            remove_missing,
            remove_mono,
            reorder,
            
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
        max=not reorder,
        out_format=out_format
    )

#this needs update
def run_morphoparse(
    input_file,
    output_prefix,
    input_format=None,
    output_format=None,
    poly=False,
    remap=False,
    remove_missing=False,
    remove_mono=False,
    reorder=False,
    asc=False,
    software=False,
    paup=False,
    tnt=False
):
    input = Path(input_file)
    if not input.exists():
        raise FileNotFoundError(f"Input file not found: {input.absolute()}")
    format = detect_format(input,input_format)
    print(f"MorphoPase will parse {input.name} in {format.upper()} format")
    output_format = (output_format or input_format).lower()
    records, poly_data = parse_sequences(input, format, poly)
    remove_missing = remap or remove_missing
    remove_uninformative = remap or remove_mono
    reorder_states = remap or reorder

    if remove_missing or remove_mono or reorder:
        logfile = f"{output_prefix}_remap.log"
        records = remap_sparse_states(
            records,
            poly_data,
            logfile,
            remove_missing,
            remove_uninformative,
            reorder_states
        )

    write_file(records, f"{output_prefix}_mparse", output_format, poly_data)
    create_partition_file(
        records,
        poly_data,
        output_prefix,
        ver=software,
        asc_correction=asc,
        paup=paup,
        tnt=tnt,
        max=not (remap or reorder),
        out_format=output_format
    )
    return{"records": records, "poly_data": poly_data}