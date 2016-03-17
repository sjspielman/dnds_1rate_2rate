#! /bin/bash

# SJS. Generate submission lines for simulation

TYPE="bias_gtr" # OR:, "bias", "nobias"


SBATCH_RAW=raw_launcher_simulation.slurm

# create output directory as needed
REPODIR=$HOME/dnmu
OUTDIR=${REPODIR}/data/alignments
mkdir -p $OUTDIR
    

LAUNCHFILE=launcher_simulate_${TYPE}.slurm
PARAMFILE=simulation_commands_${TYPE}
touch $LAUNCHFILE
touch $PARAMFILE

for REP in {2..50}
do
    for N in n7 n8 n9 n10 n11
    do            
        for BL in bl0.0025 bl0.01 bl0.04 bl0.16 bl0.64
        do
                    
                DATA=rep${REP}_${N}_${BL}_${TYPE}
                TREEFILE=$REPODIR/data/trees/${N}_${BL}.tre 
                ALN1=${DATA}.fasta
                ALN2=${DATA}_withanc.fasta
                echo python simulate_alignments.py $TREEFILE $TYPE $ALN1 $ALN2 >> $PARAMFILE
 
        done	
    done
done

sed "s/PLACEHOLDER/$PARAMFILE/" $SBATCH_RAW > $LAUNCHFILE
#sbatch $LAUNCHFILE


