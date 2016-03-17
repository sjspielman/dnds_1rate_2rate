#! /bin/bash

# SJS. Generate submission lines for simulation using real trees


SBATCH_RAW=raw_launcher_simulation_realtrees.slurm

# create output directory as needed
REPODIR=$HOME/dnmu
OUTDIR=${REPODIR}/data/realtree_alignments
mkdir -p $OUTDIR


# vertrho is super fast so we do in idev session
# amine and h3 should be 3 hours each
# camelid and hivrt should be 1 hour each


for REP in {1..10}
do
    LAUNCHFILE=launcher_simulate_realtrees_${REP}.slurm
    PARAMFILE=simulation_commands_realtrees_${TYPE}
    touch $LAUNCHFILE
    touch $PARAMFILE
    
    for TYPE in gtr bias_gtr
    do
        for TREE in amine camelid hivrt h3
        do                                
            DATA=rep${REP}_${TREE}_${TYPE}
            TREEFILE=$REPODIR/data/trees/${TREE}.tre 
            ALN1=$OUTDIR/${DATA}.fasta
            ALN2=$OUTDIR/${DATA}_withanc.fasta
            echo python simulate_alignments.py $TREEFILE $TYPE $ALN1 $ALN2 >> $PARAMFILE
        done       
    done
    sed "s/PLACEHOLDER/$PARAMFILE/" $SBATCH_RAW > $LAUNCHFILE
    #sbatch $LAUNCHFILE	
done




