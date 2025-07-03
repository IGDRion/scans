#!/bin/bash
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

##########################################################################
# Run FEELnc_classifier
# Developed by: F. Degalez
# Update : A. Besson
##########################################################################

# input variable
PROGNAME=$(basename $0)
LNC=$1
MRNA=$2
LOG=$3
OUT=$4

# create usage function that will be called each time there is an ERROR
# >&2 	: print to STDERR
usage() {
    echo "#" >&2
    echo -e "# USAGE: $PROGNAME <lncRNA.gtf> <mRNA.gtf> <log_file> <FEELnc_outfile> " >&2
    exit 1;
}
##############################

# Test number of parameters
###########################
if [[ $# -ne 4 ]]; then
	echo "ERROR : wrong number of arguments">&2
	# call usage function
	usage;
fi

## Configuration of the FEELnc environment

. /local/env/envconda.sh
conda activate scans_env

FEELnc_classifier.pl -i $LNC -a $MRNA -l $LOG > $OUT 