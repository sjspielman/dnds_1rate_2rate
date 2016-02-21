#! /bin/bash

# SJS. Generate submission lines for dN/dS inference on Stampede. The stdout is sent to the file inference_commands


TYPE="asym"

SBATCH_RAW=raw_launcher_inference.slurm

#clean_paramfile $PARAMFILE

# create output directory as needed
REPODIR=$HOME/dnmu
OUTDIR=${REPODIR}/results
mkdir -p $OUTDIR
    
TYPE="gtr"
for REP in {1..50}
do
    # touch launcher file and paramfile to go along with it for this rep 
    LAUNCHFILE=launcher_rep${TYPE}.slurm
    PARAMFILE=inference_commands${TYPE}
    touch $LAUNCHFILE
    touch $PARAMFILE
    
    for N in n7 n8 n9 n10 n11
    do            
        for BL in bl0.0025 bl0.01 bl0.04 bl0.16 bl0.64
        do
            for METHOD in FUBAR1 FUBAR2 SLAC_GTR
            do
                    
                DATA=rep${REP}_${N}_${BL}_${TYPE}
                TREE=$REPODIR/data/trees/${N}_${BL}.tre 
                ALN=$REPODIR/data/alignments/${DATA}.fasta
                OUTFILE1=$OUTDIR/${DATA}_${METHOD}.txt        # dnds inference
                OUTFILE2=$OUTDIR/${DATA}_${METHOD}_nucfit.txt # nucleotide fit
                echo sh run_inference.sh $REPODIR $DATA $ALN $TREE $METHOD $OUTFILE1 $OUTFILE2 >> $PARAMFILE
 
            done                    
        done	
    done
    # Create launcher file for this paramfile and submit. 
    sed "s/PLACEHOLDER/$PARAMFILE/" $SBATCH_RAW > $LAUNCHFILE
    sbatch $LAUNCHFILE
done

