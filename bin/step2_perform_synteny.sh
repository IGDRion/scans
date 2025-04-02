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

PROGNAME=$(basename $0)

###########################
# CHECK / USAGE
###########################

# create usage function that will be called each time there is an ERROR
# >&2 	: print to STDERR
usage() {
    echo "#" >&2
    echo -e "# USAGE: $PROGNAME --config <CONFIG_FILE> --workdir <WORK_DIR> --orthology <ORTHOLOGY_DIR> [options]"
    echo "Options:"
    echo "--synteny [method1|method2|both]   Synteny method to use (default: method1) "
    exit 1;
}

###########################
# VARIABLES
###########################

# input variable

CONFIG=""
WORKDIR=""
ORTHOLOGY=""
METHOD="method1"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG="$2"
            shift 2
            ;;
        --workdir)
            WORKDIR="$2"
            shift 2
            ;;
        --orthology)
            ORTHOLOGY="$2"
            shift 2
            ;;
        --synteny)
            METHOD="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done


# Check if mandatory arguments are provided
if [[ -z "$CONFIG" || -z "$WORKDIR" || -z "$ORTHOLOGY" ]]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

if [ ! -d "$ORTHOLOGY" ]; then
    echo "$ORTHOLOGY is not a valid directory.">&2
    usage;
fi

###########################
# CODE PRINCIPAL
###########################

echo "Config file analyzed: $CONFIG"
echo "Working directory: $WORKDIR"
echo "Directory containing homology files: $ORTHOLOGY"
echo "Synteny method used for analysis: $METHOD"

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
        echo "Synteny method 2: perform FEELnc classification."      
        bash getLNCclassFromFEELnc.sh "$CONFIG" "$WORKDIR"
        echo "bash getLNCclassFromFEELnc.sh $CONFIG $WORKDIR"

        echo "Synteny method 2: perform synteny."
        Rscript perform_synteny_method2.R "$CONFIG" "$WORKDIR" "$ORTHOLOGY"
        echo "bash getLNCclassFromFEELnc.sh $CONFIG $WORKDIR"

        ;;
    both)
        echo "Synteny method applied: method 1"
        Rscript perform_synteny_method1.R $CONFIG $WORKDIR $ORTHOLOGY

        echo "Synteny method 2: perform FEELnc classification."      
        bash getLNCclassFromFEELnc.sh "$CONFIG" "$WORKDIR"
        echo "bash getLNCclassFromFEELnc.sh $CONFIG $WORKDIR"

        echo "Synteny method 2: perform synteny."
        Rscript perform_synteny_method2.R "$CONFIG" "$WORKDIR" "$ORTHOLOGY"
        echo "bash getLNCclassFromFEELnc.sh $CONFIG $WORKDIR"
        ;;
esac
