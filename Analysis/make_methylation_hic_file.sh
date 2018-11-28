#!/bin/bash
## Creating hic file for use in methylation analysis

## In this example done for only one chromosome
## MethylDackel looks only at CpGs
MethylDackel perRead -o chr14_perRead -@ 4 /aidenlab/references/Homo_sapiens_assembly19.fasta chr14.bam 

## Must gather the reads since read names are printed more than once
awk '$5>0{a[$1]+=($4*$5/100); b[$1]+=$5}END{print "#perRead summary: readname, count methylated, total"; for (i in a){print i, a[i], b[i]}}' chr14_perRead  > chr14_perRead_compiled

## Grab chromosome 14 contacts (not intrafrag, min MAPQ > 0) from
## merged_nodups.txt file
awk '$2==14 && $6==14 && $4!=$8 && $9>0 &&$12>0' merged_nodups.txt > chr14_contacts.txt

## Create a merged_nodups containing only reads with methylation information.
## Attach that information to the readname and assign "m" or "u" based on it.
awk '$0!~/^#/ && (FNR==NR){percent=$2/$3;a[$1]=percent":"$3}(FNR!=NR)&&($(NF-1) in a)&&($NF in a){split(a[$(NF-1)],b,":"); if (b[1] >= .5){$2=$2"m"}else{$2=$2"u"} $(NF-1)=$(NF-1)":"a[$(NF-1)]; split(a[$NF],b,":");if (b[1] >= .5){$6=$6"m"}else{$6=$6"u"} $NF=$NF":"a[$NF]; print}' chr14_perRead_compiled chr14_contacts.txt > chr14_contacts_methylated.txt

## Sort by methylated/unmethylated chromosome, necessary for Hi-C file creation
sort -T. -k2,2d -k6,6d chr14_contacts_methylated.txt > chr14_contacts_methylated.txt.sort
mv chr14_contacts_methylated.txt.sort chr14_contacts_methylated.txt

## Special chrom.sizes file looks like this:
## 14m     107349540
## 14u     107349540
/aidenlab/juicebox pre -n chr14_contacts_methylated.txt chr14_contacts_methylated.hic meth.chrom.sizes
