from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir
from os.path import isfile, isdir, join
from os import getcwd


def parse_files(files, parsed_names):
    for i in range(len(files)):
        with open(parsed_names[i], 'w') as fi:
            for seq_rec in SeqIO.parse(files[i], 'fasta'):
                if len(seq_rec.seq) >= 1000:
                    fa = '>{}\n{}\n'.format(seq_rec.id, seq_rec.seq)
                    fi.write(fa)


def make_df(parsed_names, sample_names):
    frames = {}
    for i in range(len(parsed_names)):
        contig_dict = {'Metagenome': [], 'Contig': [], 'Length': [], 'GC_content': []}
        for seq_rec in SeqIO.parse(parsed_names[i], 'fasta'):
            contig_dict['Metagenome'].append(sample_names[i])
            contig_dict['Contig'].append(seq_rec.id)
            contig_dict['Length'].append(len(seq_rec.seq))
            contig_dict['GC_content'].append(GC(seq_rec.seq))
        frames[names[i]] = pd.DataFrame(contig_dict)

    metagens = pd.concat(list(frames.values()))
    return metagens


def boxplot_contig(axes, data, xticklabels, ylabel, title, color='black'):
    axes.boxplot(data)
    axes.set_xticklabels(xticklabels, rotation=90)
    axes.set_ylabel(ylabel)
    axes.set_title(title)


def parse_length(dataframe, names, names_col, threshold, threshold_col):
    length_sub = []
    for name in names:
        dataframe_sub = dataframe[dataframe[names_col] == name]
        parsed_df_sub = dataframe_sub[dataframe_sub[threshold_col] > threshold]
        length_sub.append(parsed_df_sub[threshold_col])
    return length_sub


def make_figure(metagens, names):
    base_lengths = parse_length(metagens, names, 'Metagenome', 0, 'Length')
    base_lengths_g500 = parse_length(metagens, names, 'Metagenome', 500, 'Length')
    base_lengths_g1000 = parse_length(metagens, names, 'Metagenome', 1000, 'Length')

    fig, ax = plt.subplots(3, sharey=True)

    boxplot_contig(ax[0], base_lengths, [''] * len(names), 'Contig length (bp)',
                   'Metagenome Contig Length Distribution')
    boxplot_contig(ax[1], base_lengths_g500, [''] * len(names), 'Contig length (bp)',
                   'Metagenome Contig Length Distribution (Contigs > 500 bp)')
    boxplot_contig(ax[2], base_lengths_g1000, names, 'Contig length (bp)',
                   'Metagenome Contig Length Distribution (Contigs > 1000bp)')

    fig.set_figheight(15)
    plt.yscale(value='log')
    plt.savefig('../../results/plots/Contig_Length_Distribution.png')
    plt.show()


def main():
    mypath = getcwd()
    onlydirs = [d for d in listdir(mypath) if isdir(join(mypath, d))]
    files = ['./{}/scaffolds.fa'.format(d) for d in onlydirs]
    sample_names = [n[7:] for n in onlydirs]
    parsed_names = ['./{}/{}_greater-1kb.fa'.format(d, d[7:]) for d in onlydirs]
    parse_files(files, parsed_names)
    metagens = make_df(parsed_names, sample_names)

    make_figure(metagens, names)
    print(len(parsed_names), len(directories), len(files))


if __name__ == '__main__':
    main()
