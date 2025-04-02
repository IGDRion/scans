#!/bin/bash
#SBATCH --mem=34G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user aurore.besson@univ-rennes.fr
#SBATCH --mail-type=ALL

##########################################################################
# Perform sequence alignment using liftoff and bedtools 
# Developed by: A. Besson
# Updated by : 
##########################################################################

PROGNAME=$(basename $0)

# Function to display usage
usage() {
    echo "Usage: $PROGNAME --query_fa <QUERY_FA> --query_gtf <QUERY_GTF> --target_fa <TARGET_FA> --target_gtf <TARGET_GTF> --workdir <WORKDIR> [options]"
    echo "Options:"
    echo "  --biotype [mRNA|lncRNA|all]     Biotype to analyze (default: all)"
    echo "  --flank [0.0-1.0]               Flank fraction for liftoff (default: 0)"
    echo " --coverage [0.0-1.0]             Coverage fraction for liftoff (default: 0.5)"
    echo "  --identity [0.0-1.0]            Sequence identity fraction for liftoff (default: 0)"
    exit 1
}


###########################
# VARIABLES
###########################

# Initialize variables
QUERY_FA=""
QUERY_GTF=""
TARGET_FA=""
TARGET_GTF=""
WORKDIR=""
BIOTYPE="all"
FLANK=0
COVERAGE=0.5
IDENTITY=0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --query_fa)
            QUERY_FA="$2"
            shift 2
            ;;
        --query_gtf)
            QUERY_GTF="$2"
            shift 2
            ;;
        --target_fa)
            TARGET_FA="$2"
            shift 2
            ;;
        --target_gtf)
            TARGET_GTF="$2"
            shift 2
            ;;
        --workdir)
            WORKDIR="$2"
            shift 2
            ;;
        --biotype)
            BIOTYPE="$2"
            shift 2
            ;;
        --flank)
            FLANK="$2"
            shift 2
            ;;
        --coverage)
            COVERAGE="$2"
            shift 2
            ;;
        --identity)
            IDENTITY="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Check if mandatory arguments are provided
if [[ -z "$QUERY_FA" || -z "$QUERY_GTF" || -z "$TARGET_FA" || -z "$TARGET_GTF" || -z "$WORKDIR" ]]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Arguments used
echo "Query genome: $QUERY_FA"
echo "Query annotation: $QUERY_GTF"
echo "Target genome: $TARGET_FA"
echo "Target annotation: $TARGET_GTF"
echo "Path to results directory: $WORKDIR"
echo "Flank value for liftoff: $FLANK"
echo "Coverage cut-off for liftoff: $COVERAGE"
echo "Sequence identity cut-off for liftoff: $IDENTITY"

SPECIES1=$(basename $QUERY_GTF .gtf)
SPECIES2=$(basename $TARGET_GTF .gtf)

### STEP1 - liftoff
RESULTS="$WORKDIR"/work/method3/"$SPECIES1"_to_"$SPECIES2"
mkdir -p ${RESULTS}

echo "Step1: liftoff analysis - IN PROGRESS"
echo "Working directory: $WORKDIR"

# perform sequence alignment according BIOTYPE argument
echo "Gene biotype to analyze: $BIOTYPE"
case "$BIOTYPE" in
    all)
        ## select features to add to `-f`
        FEATURES=$RESULTS/all_features.txt
        grep -v '#' $QUERY_GTF | cut -f3 | sort | uniq > $FEATURES
        echo "Features file created: $FEATURES"

        ## create outfile and outdir
        LIFTOFF_DIR="$RESULTS"/liftoff_flank"$FLANK"
        LO_OUT=liftoff_"$SPECIES1"_to_"$SPECIES2"_flank"$FLANK".gtf

        ## launch liftoff script
        echo "Launch liftoff analysis: "
        ### source environment
        . /local/env/envconda3.sh
        conda activate activate /home/genouest/cnrs_umr6290/abesson/conda_env/liftoff_env

        ### run liftoff
        echo "liftoff -g $QUERY_GTF -o $LO_OUT -p 8 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $LIFTOFF_DIR -u $LIFTOFF_DIR/unmapped_features.txt -copies"
        liftoff -g $QUERY_GTF -o $LO_OUT -p 8 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $LIFTOFF_DIR -u $LIFTOFF_DIR/unmapped_features.txt -copies

        ;;
    lncRNA)
        echo "option not available yet"
        ;;
    mRNA)
        echo "option not available yet"
        ;;
esac

# move features file in liftoff dir
mv $FEATURES $LIFTOFF_DIR
mv $LO_OUT $RESULTS

LO_OUT2="$RESULTS"/"$LO_OUT"

### STEP2 - Filter liftoff results according cut-off arguments
echo "Step2: liftoff results filtering - IN PROGRESS"

# filtered file
LO_FILTER="$RESULTS"/$(basename $LO_OUT2 .gtf)_filtered.gtf

conda deactivate
conda activate /home/genouest/cnrs_umr6290/abesson/conda_env/jupyterR_env

Rscript filter_liftoffOutput.R $LO_OUT2 $LO_FILTER $COVERAGE $IDENTITY $RESULTS

### STEP3 - Bedtools intersect
echo "Step3: bedtools analysis - IN PROGRESS"
BEDTOOLS_DIR="$RESULTS"/bedtools_intersect
mkdir ${BEDTOOLS_DIR}

### create outfile for liftoff output, target annotation and bedtools output
LO_BED="$BEDTOOLS_DIR"/$(basename $LO_FILTER .gtf).bed
TARGET_BED="$BEDTOOLS_DIR"/"$SPECIES2".bed
FINAL_BED="$BEDTOOLS_DIR"/overlap_"$SPECIES1"_to_"$SPECIES2".bed

#### format gtf in bed 
. format_gtf2bed.sh $LO_FILTER > $LO_BED
echo "Bed file created: $LO_BED"

. format_gtf2bed.sh $TARGET_GTF > $TARGET_BED
echo "Bed file created: $TARGET_BED"

#### bedtools intersect - default fraction ~1bp
conda deactivate
. /local/env/envbedtools-2.27.1.sh

echo "Intersection analysis done by bedtools intersect."
bedtools intersect -a $LO_BED -b $TARGET_BED -wao -split > $FINAL_BED
echo "bedtools intersect -a $LO_BED -b $TARGET_BED -wao -split > $FINAL_BED"

### STEP4 - Sequence alignment analysis
echo "Step4: Sequence alignment analysis - IN PROGRESS"
ALIGN_DIR="$RESULTS"/alignment_analysis
mkdir ${ALIGN_DIR}

# load conda env => en créer un spécifiquement pour ça ??
conda deactivate
conda activate /home/genouest/cnrs_umr6290/abesson/conda_env/jupyterR_env

echo "Rscript seq_alignment_analysis.R $QUERY_GTF $TARGET_GTF $FINAL_BED $LIFTOFF_DIR/unmapped_features.txt $ALIGN_DIR"
Rscript seq_alignment_analysis.R $QUERY_GTF $TARGET_GTF $FINAL_BED $LIFTOFF_DIR/unmapped_features.txt $ALIGN_DIR
