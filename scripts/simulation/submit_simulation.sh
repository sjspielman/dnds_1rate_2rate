#! /bin/bash

# SJS. Generate submission lines for simulation

TYPE="asym" # or, "nobias"


SBATCH_RAW=raw_launcher_simulation.slurm

# create output directory as needed
REPODIR=$HOME/dnmu
OUTDIR=${REPODIR}/data/alignments
mkdir -p $OUTDIR
    

for REP in 1
do
    # touch launcher file and paramfile to go along with it for this rep 
    LAUNCHFILE=launcher_rep${REP}.slurm
    PARAMFILE=inference_commands${REP}
    touch $LAUNCHFILE
    touch $PARAMFILE
    
    for N in n7 n8 n9 n10 n11
    do            
        for BL in bl0.0025 bl0.01 bl0.04 bl0.16 bl0.64
        do
                    
                DATA=rep${REP}_${N}_${BL}_${TYPE}
                TREEFILE=$REPODIR/data/trees/${N}_${BL}.tre 
                ALN1=$OUTDIR/${DATA}.fasta
                ALN2=$OUTDIR/${DATA}_withanc.fasta
                echo python simulate_alignments.py $TREEFILE $TYPE $ALN1 $ALN2 >> $PARAMFILE
 
            done                    
        done	
    done
    # Create launcher file for this paramfile and submit. 
    sed "s/PLACEHOLDER/$PARAMFILE/" $SBATCH_RAW > $LAUNCHFILE
    sbatch $LAUNCHFILE
done

