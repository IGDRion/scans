#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=34G
#SBATCH --mail-user aurore.besson@univ-rennes.fr
#SBATCH --mail-type=ALL

##########################################################################
# Run lift-off from Shumate and Salzberg (2020)
# Developed by: Aurore Besson
# Created: 07/04/2023
##########################################################################

# input variable
PROGNAME=$(basename $0)
QUERY_GTF=$1
QUERY_FA=$2
TARGET_FA=$3
FLANK=$4
OUTFILE=$5
FEATURES=$6
OUTDIR=$7

# create usage function that will be called each time there is an ERROR
# >&2 	: print to STDERR
usage() {
    echo "#" >&2
    echo -e "# USAGE: $PROGNAME <QUERY_GTF> <QUERY_FA> <TARGET_FA> <FLANK> <OUTFILE> <FEATURES> <TEMP_DIR> " >&2
    exit 1;
}
##############################

# Test number of parameters
###########################
if [[ $# -ne 7 ]]; then
	echo "ERROR : wrong number of arguments">&2
	# call usage function
	usage;
fi

##########################
# CODE
###########################

# source environment
. /local/env/envconda3.sh
conda activate activate /home/genouest/cnrs_umr6290/abesson/conda_env/liftoff_env

# run liftoff

## run liftoff tool
echo "liftoff -g $QUERY_GTF -o $OUTFILE -p 8 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $OUTDIR -u $OUTDIR/unmapped_features.txt -copies"
liftoff -g $QUERY_GTF -o $OUTFILE -p 8 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $OUTDIR -u $OUTDIR/unmapped_features.txt -copies
