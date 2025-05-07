from morphoparse.parsers import parse_fasta, parse_phylip, parse_nexus, parse_tnt
from morphoparse.writers import write_fasta, write_phylip, write_nexus, write_tnt

def parse_sequences(input, file_format, keep):
    if input.stat().st_size == 0:
        raise ValueError("Empty file")    
    parsers = {
        'fasta': parse_fasta,
        'phylip': parse_phylip,
        'nexus': parse_nexus,
        'tnt': parse_tnt,
    }
    return parsers[file_format](input, keep)

def write_file(records, output, out_format, poly):
    writers = {
        'fasta': write_fasta,
        'phylip': write_phylip,
        'nexus': write_nexus,
        'tnt': write_tnt,
    }
    writers[out_format](records, output, poly)
