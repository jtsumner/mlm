import argparse
from Bio import SeqIO

def parse_scaffolds(sample, scaffolds, parsed, min_length, assembler):
    with open(parsed, "w+") as pa:
        with open(scaffolds, "r") as scaf:
            for seq_rec in SeqIO.parse(scaf, 'fasta'):
                if len(seq_rec.seq) >= int(min_length):
                    if assembler == "spades":
                        seq_rec.id = seq_rec.id.split('_')[1] # get contig number
                    new_id = "{}_{}".format(sample, seq_rec.id)
                    fa = '>{}\n{}\n'.format(new_id, seq_rec.seq)
                    pa.write(fa)

# Create an argument parser
parser = argparse.ArgumentParser(description='Parse scaffolds')

# Add command line arguments
parser.add_argument('--sample', required=True, help='Sample name')
parser.add_argument('--scaffolds', required=True, help='Path to the input scaffolds file')
parser.add_argument('--parsed', required=True, help='Path to the output parsed file')
parser.add_argument('--min_length', required=True, help='Passing length of scaffolds')
parser.add_argument('--assembler', default='other', help='Assembler name (default: other)')

# Parse the command line arguments
args = parser.parse_args()

# Call the function with command line arguments
parse_scaffolds(args.sample, args.scaffolds, args.parsed, args.min_length, args.assembler)

## USAGE ##
## python parse_contigs.py --sample SAMPLE_NAME --scaffolds INPUT_FILE --parsed OUTPUT_FILE --min_length INT
