#! /bin/bash

# Create .tsv file of file name prefix and dataset (dataset = read file parent directory )
for i in $(ls Batch_0*/*.gz) ; do echo $i | awk -F'_R' '{print $1}' | awk -F'/' '{print $2 "\t" $1}' >> a.tsv ; done

# remove duplicates from R1/R2
sort -u a.tsv > b.tsv

# add header
echo -ne "sample\tdataset\n" | cat - b.tsv > samples_minimal.tsv

# rm temp files
rm a.tsv b.tsv

