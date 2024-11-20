models=(miseq hiseq novaseq nextseq)
pref="test_genomes"
genome=(Phix.fna Ecoli.fna Scer.fna Hsap.fna)
counts=(25000 25000 25000 25000)

# change files
for GN in ${genome[*]}; do
    sed -i "s/>/>$GN-/g" "$pref/$GN"
done


# generate

for MD in ${models[*]}; do
    for i in ${!genome[*]}; do
        GN=${genome[$i]}
        CN=${counts[$i]}
        iss generate --genomes $pref/$GN --model $MD --n_reads $CN -o $GN$MD --cpus 6
    done
    
    cat *${MD}_R1* > ${MD}_R1.fq
    cat *${MD}_R2* > ${MD}_R2.fq
    rm *.fastq
    rm *abundance*
done
