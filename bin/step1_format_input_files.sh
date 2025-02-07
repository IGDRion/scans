#!/bin/bash
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user aurore.besson@univ-rennes.fr
#SBATCH --mail-type=ALL

##########################################################################
# Format input files to perform synteny analysis
# Developed by: F. Degalez
# Updated by : A. Besson
##########################################################################

###########################
# VARIABLES
###########################

# input variable
PROGNAME=$(basename $0)
CONFIG=""
WORKDIR=""
ORTHOLOGY=""

###########################
# CHECK / USAGE
###########################

# create usage function that will be called each time there is an ERROR
# >&2 	: print to STDERR
usage() {
    echo "#" >&2
    echo -e "# USAGE: $PROGNAME <CONFIG_FILE> <WORK_DIR> [-o|--orthology <ORTHOLOGY_DIR>]  " >&2
    echo -e "# Example: $PROGNAME config.txt results/ [-o|--orthology data/orthology/]  " >&2
    exit 1;
}


###########################
# GESTION ARGUMENTS
###########################

# mandatory arguments
CONFIG="$1"
WORKDIR="$2"
shift 2

# optional argument
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--orthology)
            if [[ -n $2 && -d $2 ]]; then
                ORTHOLOGY="$2"
                shift 2
            else
                echo "Error : Argument -o/--orthology needs a valid path." >&2
                exit 1
            fi
            ;;
        *)
            echo "Invalid argument: $1" >&2
            usage
            ;;
    esac
done

###########################
# CODE PRINCIPAL
###########################

mkdir ${WORKDIR}

echo "Config file analyzed: $CONFIG"
echo "Working directory: $WORKDIR"

# load conda env => en créer un spécifiquement pour ça ??
. /local/env/envconda.sh
conda activate /home/genouest/cnrs_umr6290/abesson/conda_env/jupyterR_env

#export PATH="/projects/dog/aurore/BrownLincs/module/bin:$PATH" # fonctionne pas -> à revoir pour appeler le script de partout

# launch R script to format gtf files
Rscript format_gtf.R $CONFIG $WORKDIR

# launch R script to format orthology if needed
if [[ -n $ORTHOLOGY ]]; then
    echo "Directory containing orthology files: $ORTHOLOGY"
    Rscript format_orthology.R $CONFIG $WORKDIR $ORTHOLOGY
else
    echo "Orthology option not used."
fi