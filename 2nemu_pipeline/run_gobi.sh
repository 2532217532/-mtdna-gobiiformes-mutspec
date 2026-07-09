#!/bin/bash

PATH_TO_NEMU=/home/zengjl/gitclone/nemu-pipeline-nf/main.nf
PATH_TO_CONFIG=/home/zengjl/gitclone/nemu-pipeline-nf/nextflow.config
PATH_TO_INPUT=/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/nemu_input
PATH_TO_OUTPUT=/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/nemu_output

GENES=("Cytb" "CO1" "CO2" "CO3" "ND1" "ND2" "ND3" "ND4" "ND4L" "ND5" "ND6" "A6" "A8")

MAX_NJOBS=30
SLEEP_TIME=200  # secs
RUNLIMIT=20000
COUNTER=1
echo "Execute $RUNLIMIT genes"

for gene in ${GENES[@]}; do
    echo $gene
    INDIR=$PATH_TO_INPUT/$gene
    OUTDIR=$PATH_TO_OUTPUT/$gene

    for seq in $INDIR/*.fasta; do
        if [ $COUNTER -gt $RUNLIMIT ]; then
            exit 0
        fi
        
        gene_species=$(basename $seq .fasta)
        species=${gene_species#*__}

        workdir=$OUTDIR/$gene_species
        if [ -e $workdir ]; then 
            # echo "'$workdir' already exist: pass"
            continue
        fi
        let COUNTER++

        seq_abs=`realpath $seq`

        echo -e "workdir: $workdir"
        mkdir -p $workdir
        cd $workdir

        # -resume
        # -q -bg
        nextflow -bg -q -c $PATH_TO_CONFIG run $PATH_TO_NEMU -with-trace \
            --sequence $seq_abs --species_name $species \
            --DB /home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/ncbi_db/${gene}_BLAST --outdir .
        cd - >/dev/null
        sleep 1

        if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
            echo sleep $SLEEP_TIME secs
            sleep $SLEEP_TIME
        fi
    done
    echo -e '\n'
done

echo "Processed $COUNTER fastas"

# rm -rf $PATH_TO_OUTPUT/*/*/work
