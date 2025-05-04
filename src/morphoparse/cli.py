"""
Handles command-line interface parsing.
"""
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="MorphoParse v1.0.0: Parsing, partitioning, and weighting morphological character matrices for phylogenetic analysis.")
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-o", "--output", help="Output prefix")
    parser.add_argument("-f", "--format", choices=['fasta','phylip','nexus','tnt'], default='phylip', 
        help="Input format (default: phylip)")
    parser.add_argument("-g", "--out_format", choices=['fasta','phylip','nexus','tnt'], 
        help="Output format (default: same as input)")

    parser.add_argument("-r", "--remap", action="store_true",
                    help="Remap characters: removes + reorder")
    parser.add_argument("--remove-missing", action="store_true",
                    help="Remove characters with only missing data (?, -)")
    parser.add_argument("--remove-mono", action="store_true",
                    help="Remove characters with only one unambiguous state")
    parser.add_argument("--reorder", action="store_true",
                    help="Renumber states to 0,1,... consistently per site")
    
    parser.add_argument("-a", "--asc", action="store_true", 
        help="Apply ASC correction (recommended)")
    parser.add_argument("-s", "--software", default='raxml', choices=['raxml','raxmlng','iqtree'],
        help="Indicate program to format for (default: RAxML)")

    parser.add_argument("-p", "--paup", action="store_true", 
        help="Generate SSA weight scheme for PAUP*")
    parser.add_argument("-t", "--tnt", action="store_true", 
        help="Generate SSA weight scheme for TNT")
    parser.add_argument("--version", action="version", version="MorphoParse 1.0.0 (2025-04-21)")

    return parser.parse_args()