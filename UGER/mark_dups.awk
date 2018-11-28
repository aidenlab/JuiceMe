BEGIN{
    OFS="\t"; 
    fname=outputdir"/dups.txt"
    while ((getline < fname) > 0) {
        split($NF, a, "/"); 
        dups[a[1]]
    } 
    fname=outputdir"/opt_dups.txt"
    while ((getline < fname) > 0) {
        split($NF, a, "/"); 
        dups[a[1]]
    }
}
$1~/^@/{
    print
}
$1!~/^@/{
    split($1, a, "/"); 
    if (a[1] in dups) {
        $2 = or($2, 1024);
    }
    print
}
