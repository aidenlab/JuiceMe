#!/bin/bash
# Hi-Culfite script for making observed / expected matrices
# See README.md for more information

#Juicer tools, for running dump
juicer_tools=/aidenlab/juicebox

# Given a hic file that contains contacts divided by methylation status
# (via a simulated genome that repeats each chromosome; this allows us 
# to preserve asymmetric contacts between unmethylated and methylated)
# See make_methylation_hic_file.sh
hic_file=$1

chr=14

for res in 50 100 250 500 1000 2500
do
    # Create the M, U, and Y matrices
    val=$(($res*1000))
    $juicer_tools dump observed NONE $hic_file ${chr}m ${chr}m BP $val M_chr${chr}_${res}K.txt
    awk -v res=$val '{$1=int($1/res); $2=int($2/res); $3=int($3); print; if ($1!= $2){print $2,$1,$3}}' M_chr${chr}_${res}K.txt >  M_chr${chr}_${res}K.txt2
    mv  M_chr${chr}_${res}K.txt2  M_chr${chr}_${res}K.txt
    $juicer_tools dump observed NONE $hic_file ${chr}u ${chr}u BP $val U_chr${chr}_${res}K.txt
    awk -v res=$val '{$1=int($1/res); $2=int($2/res); $3=int($3); print; if ($1!= $2){print $2,$1,$3}}' U_chr${chr}_${res}K.txt >  U_chr${chr}_${res}K.txt2
    mv  U_chr${chr}_${res}K.txt2  U_chr${chr}_${res}K.txt
    $juicer_tools dump observed NONE $hic_file ${chr}m ${chr}u BP $val Y_chr${chr}_${res}K.txt
    awk -v res=$val '{$1=int($1/res); $2=int($2/res); $3=int($3); print}' Y_chr${chr}_${res}K.txt >  Y_chr${chr}_${res}K.txt2
    mv Y_chr${chr}_${res}K.txt2  Y_chr${chr}_${res}K.txt    

    # Make the methylation vector and O/E matrices
    python make_a.py $chr $res
    python make_oe.py $chr $res
done

