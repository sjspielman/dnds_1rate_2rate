#! /bin/bash

# SJS. Generate submission lines for simulation using real trees


SBATCH_RAW=raw_launcher_simulation_realtrees.slurm

# create output directory as needed
REPODIR=$HOME/dnmu
OUTDIR=${REPODIR}/data/alignments_realtrees
mkdir -p $OUTDIR


# vertrho is super fast so we do in idev session
# amine h3 should be 2 hours each
# camelid hivrt should be 30 min each


LAUNCHFILE=launcher_simulate_realtrees1.slurm
PARAMFILE=simulation_commands_realtrees1
touch $LAUNCHFILE
touch $PARAMFILE
    
for REP in {21..50}
do    
    for TYPE in gtr bias_gtr
    do
        for TREE in h3 amine
        do                                
            DATA=rep${REP}_${TREE}_${TYPE}
            TREEFILE=$REPODIR/data/trees/${TREE}.tre 
            ALN1=$OUTDIR/${DATA}.fasta
            ALN2=$OUTDIR/${DATA}_withanc.fasta
            echo python simulate_alignments.py $TREEFILE $TYPE $ALN1 $ALN2 >> $PARAMFILE
        done       
    done
done
sed "s/PLACEHOLDER/$PARAMFILE/" $SBATCH_RAW > $LAUNCHFILE
sbatch $LAUNCHFILE	




