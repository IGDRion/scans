#!/bin/bash

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
ORTHOLOGY="no_orthology"

###########################
# CHECK / USAGE
###########################

# create usage function that will be called each time there is an ERROR
# >&2 	: print to STDERR
usage() {
    echo "#" >&2
    echo -e "# USAGE: $PROGNAME --config <CONFIG_FILE> --output <WORK_DIR> [--orthology <ORTHOLOGY_DIR>]  " >&2
    echo -e "# Example: $PROGNAME --config config.txt -- output project/ [--orthology path/to/orthology/]  " >&2
    exit 1;
}


###########################
# GESTION ARGUMENTS
###########################

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="$2"; shift 2 ;;
        --output) WORKDIR="$2"; shift 2 ;;
        --orthology) ORTHOLOGY="$2"; shift 2 ;;
        *) usage ;;
    esac
done

# Check if mandatory arguments are provided
if [[ -z "$CONFIG" || -z "$WORKDIR" ]]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Arguments used
echo "Configuration file: $CONFIG"
echo "Path to results directory: $WORKDIR"

###########################
# CODE PRINCIPAL
###########################

mkdir ${WORKDIR}

echo "Config file analyzed: $CONFIG"
echo "Working directory: $WORKDIR"

# launch R script to format orthology if needed
echo "ANALYSIS IN PROGRESS: Format gtf files"
Rscript format_gtf.R $CONFIG $WORKDIR

if [[ "$ORTHOLOGY" == "no_orthology" ]] ; then
    echo "No orthology file reformatting"
else
    echo "ANALYSIS IN PROGRESS: Format orthology files"
    echo "Directory containing homology files: $ORTHOLOGY"
    Rscript format_orthology.R $CONFIG $WORKDIR $ORTHOLOGY
fi