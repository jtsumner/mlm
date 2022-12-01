from Bio import SeqIO

def parse_scaffolds(sample, scaffolds, parsed):
    with open(parsed, "w+") as pa:
        with open(scaffolds, "r") as scaf:
            for seq_rec in SeqIO.parse(scaf, 'fasta'):
                if len(seq_rec.seq) >= 1000:
                    new_id = "{}_{}".format(sample, seq_rec.id)
                    fa = '>{}\n{}\n'.format(new_id, seq_rec.seq)
                    pa.write(fa)

parse_scaffolds(snakemake.wildcards.sample, snakemake.input[0], snakemake.output[0])
