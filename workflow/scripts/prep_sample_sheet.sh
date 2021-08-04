#! /bin/bash

rm samples.tsv
rm tmp.tsv

for i in $(ls Batch_0*/*.gz) ; do echo $i | awk -F'_R' '{print $1}' >> tmp.tsv ; done

while read p ; do echo $p ; awk -F'/' '{print $2 "\t" $1}' >> samples.tsv ; done < tmp.tsv

sed -i '1s/^/samples\tdataset\n/' samples.tsv

