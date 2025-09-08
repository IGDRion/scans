#!/bin/bash
#SBATCH --mem=34G
#SBATCH --cpus-per-task=16

##########################################################################
# Perform sequence alignment using Liftoff and Bedtools 
# Developed by: A. Besson
##########################################################################

PROGNAME=$(basename $0)

# Function to display usage
usage() {
    echo "Usage: $PROGNAME --config <CONFIG> --output <WORKDIR> [options]"
    echo "Options:"
    echo "  --biotype [mRNA|lncRNA]     Biotype to analyze (default: lncRNA)"
    echo "  --flank [0.0-1.0]           Flank fraction for liftoff (default: 0)"
    echo "  --coverage [0.0-1.0]        Coverage fraction for liftoff (default: 0.5)"
    echo "  --identity [0.0-1.0]        Sequence identity fraction for liftoff (default: 0)"
    exit 1
}

###########################
# VARIABLES
###########################

# Initialize variables
CONFIG=""
WORKDIR=""
BIOTYPE="lncRNA"
FLANK=0
COVERAGE=0.5
IDENTITY=0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="$2"; shift 2 ;;
        --output) WORKDIR="$2"; shift 2 ;;
        --biotype) BIOTYPE="$2"; shift 2 ;;
        --flank) FLANK="$2"; shift 2 ;;
        --coverage) COVERAGE="$2"; shift 2 ;;
        --identity) IDENTITY="$2"; shift 2 ;;
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
echo "Gene biotype to align: $BIOTYPE"
echo "Flank value for liftoff: $FLANK"
echo "Coverage cut-off for liftoff: $COVERAGE"
echo "Sequence identity cut-off for liftoff: $IDENTITY"

echo "Sequence alignment analysis - IN PROGRESS"

# read config file to retrieve each species pair
declare -A completeName shortName pathToGTF pathToFasta

while IFS=',' read -r name short gtf fasta; do
    completeName["$name"]="$name"
    shortName["$name"]="$short"
    pathToGTF["$name"]="$gtf"
    pathToFasta["$name"]="$fasta"
done < <(tail -n +2 "$CONFIG") 

# list of species to analyse
species=("${!completeName[@]}")

# generate all species pair to analyse
for ((i = 0; i < ${#species[@]}; i++)); do
    for ((j = 0; j < ${#species[@]}; j++)); do
        if [[ $i -ne $j ]]; then
            species1="${species[$i]}"
            species2="${species[$j]}"

            echo "Pair analysis : $species1 to $species2"
            echo "Sequence alignment starts : $(date +"%d.%m.%Y %H:%M")"

            # result files/directories
            ## directory results by species
            RESULTS="$WORKDIR"/scans_results/method3/"${shortName[$species1]}"
            echo "mkdir -p $RESULTS"

            ## files: subset annotation by biotype & feature database associated
            file_name=$(basename ${pathToGTF[$species1]} .gtf)_"$BIOTYPE".gtf
            QUERYDB="$RESULTS"/"$file_name"_db

            ## sub-directory for each pair analysis
            COMPDIR="$RESULTS"/"${shortName[$species1]}"_to_"${shortName[$species2]}"_"$BIOTYPE"
            
            # perform sequence aligment by pair
            ## check if output dir already exists
            if [ ! -d "$COMPDIR" ]; then
                ## check if feature database already exists
                if [ -f $QUERYDB ]; then
                    echo "Features database already exists: $QUERYDB."
                    echo "Liftoff process starts."
                    bash perform_liftoff_byPair.sh \
                        --query_sp ${shortName[$species1]} \
                        --query_fa ${pathToFasta[$species1]} \
                        --query_gtf ${pathToGTF[$species1]} \
                        --target_sp ${shortName[$species2]}\
                        --target_fa ${pathToFasta[$species2]} \
                        --target_gtf ${pathToGTF[$species2]} \
                        --workdir $RESULTS \
                        --biotype $BIOTYPE \
                        --flank $FLANK \
                        --coverage $COVERAGE \
                        --identity $IDENTITY \
                        --withDB TRUE
                    echo "Sequence alignment ends : $(date +"%d.%m.%Y %H:%M")"
                    echo " "
                else
                    echo "Features database doesn't exist."
                    echo "$QUERYDB"
                    echo "Liftoff process starts."
                    bash perform_liftoff_byPair.sh \
                        --query_sp ${shortName[$species1]} \
                        --query_fa ${pathToFasta[$species1]} \
                        --query_gtf ${pathToGTF[$species1]} \
                        --target_sp ${shortName[$species2]}\
                        --target_fa ${pathToFasta[$species2]} \
                        --target_gtf ${pathToGTF[$species2]} \
                        --workdir $RESULTS \
                        --biotype $BIOTYPE \
                        --flank $FLANK \
                        --coverage $COVERAGE \
                        --identity $IDENTITY \
                        --withDB FALSE
                    echo "Sequence alignment ends : $(date +"%d.%m.%Y %H:%M")"
                    echo " "
                fi
            else
                echo "Directory $COMPDIR already exists. Sequence analysis not performed."
                echo " "
            fi  
        fi
    done
done
echo "All sequence alignments done."
echo "$(date +"%d.%m.%Y %H:%M")"

# Generate global results by species
if [ $BIOTYPE = "lncRNA" ]; then
    echo "Generate global results for alignment analysis - IN PROGRESS"
    getGlobalResultMethod3.R $CONFIG $WORKDIR
else
    echo "mRNA biotype analyzed. No global results file generated."
fi
