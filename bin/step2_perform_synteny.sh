#!/bin/bash
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user aurore.besson@univ-rennes.fr
#SBATCH --mail-type=ALL

##########################################################################
# Perform synteny analysis 
# Developed by: F. Degalez
# Updated by : A. Besson
##########################################################################

###########################
# VARIABLES
###########################

# input variable
PROGNAME=$(basename $0)
CONFIG=$1
WORKDIR=$2
ORTHOLOGY=$3
shift 3 

# optional argument (by default value)
METHOD="method1"

###########################
# CHECK / USAGE
###########################

# create usage function that will be called each time there is an ERROR
# >&2 	: print to STDERR
usage() {
    echo "#" >&2
    echo -e "# USAGE: $PROGNAME <CONFIG_FILE> <WORK_DIR> <ORTHOLOGY_DIR> [-s|--synteny method1(default)] " >&2
    echo -e "# Example: $PROGNAME config.txt results/ orthology/ [-s|--synteny method1|method2|both] " >&2
    exit 1;
}

# Check number of mandatory arguments
if [ $# -lt 3 ]; then
    echo "ERROR : wrong number of mandatory arguments.">&2
    usage;
fi

if [ ! -d "$ORTHOLOGY" ]; then
    echo "$ORTHOLOGY is not a valid directory.">&2
    usage;
fi

###########################
# GESTION ARGUMENTS
###########################

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--synteny)
            if [[ $2 == "method1" || $2 == "method2" || $2 == "both" ]]; then
                METHOD=$2
                shift 2
            else
                echo "Error : wrong value for -s/--synteny option. Please choose: method1, method2 or both."
                exit 1
            fi
            ;;
        *)
            echo "Not recognized option : $1"
            usage;
            exit 1
            ;;
    esac
done


###########################
# CODE PRINCIPAL
###########################

echo "Config file analyzed: $CONFIG"
echo "Working directory: $WORKDIR"
echo "Directory containing homology files: $ORTHOLOGY"

# load conda env => en créer un spécifiquement pour ça ??
. /local/env/envconda.sh
conda activate /home/genouest/cnrs_umr6290/abesson/conda_env/jupyterR_env

# perform synteny according method argument
case "$METHOD" in
    method1)
        echo "Synteny method applied: method 1"
        Rscript perform_synteny_method1.R $CONFIG $WORKDIR $ORTHOLOGY
        ;;
    method2)
        echo "Synteny method applied: method 2"
        Rscript perform_synteny_method2.R $CONFIG $WORKDIR $ORTHOLOGY
        ;;
    both)
        echo "Synteny method applied: methods 1 and 2"
        Rscript perform_synteny_method1.R $CONFIG $WORKDIR $ORTHOLOGY
        Rscript perform_synteny_method2.R $CONFIG $WORKDIR $ORTHOLOGY
        ;;
esac
