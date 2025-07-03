#!/bin/bash

##########################################################################
# Perform sequence alignment using liftoff and bedtools 
# Developed by: A. Besson
# Updated by : 
##########################################################################

PROGNAME=$(basename $0)

# Function to display usage
usage() {
    echo "Usage: $PROGNAME --query_sp <QUERY_SPECIES> --query_fa <QUERY_FA> --query_gtf <QUERY_GTF> --target_sp <TARGET_SPECIES> --target_fa <TARGET_FA> --target_gtf <TARGET_GTF> --workdir <WORKDIR> [options]"
    echo "Options:"
    echo " --biotype [mRNA|lncRNA]     Biotype to analyze (default: lncRNA)"
    echo " --flank [0.0-1.0]               Flank fraction for liftoff (default: 0)"
    echo " --coverage [0.0-1.0]            Coverage fraction for liftoff (default: 0.5)"
    echo " --identity [0.0-1.0]            Sequence identity fraction for liftoff (default: 0)"
    echo " --withDB [logical]              to speed up the process, use feature database if already created"
    exit 1
}


###########################
# VARIABLES
###########################

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --query_sp) QUERY_SP="$2"; shift ;;
        --query_fa) QUERY_FA="$2"; shift ;;
        --query_gtf) QUERY_GTF="$2"; shift ;;
        --target_sp) TARGET_SP="$2"; shift ;;
        --target_fa) TARGET_FA="$2"; shift ;;
        --target_gtf) TARGET_GTF="$2"; shift ;;
        --workdir) WORKDIR="$2"; shift ;;
        --biotype) BIOTYPE="$2"; shift ;;
        --flank) FLANK="$2"; shift ;;
        --coverage) COVERAGE="$2"; shift ;;
        --identity) IDENTITY="$2"; shift ;;
        --withDB) DATABASE="$2"; shift ;;
        *) usage ;;
    esac
    shift
done

# sub-directory for each pair analysis
COMPDIR="$WORKDIR"/"$QUERY_SP"_to_"$TARGET_SP"_"$BIOTYPE"
mkdir -p ${COMPDIR}

# log file
LOG_FILE="$COMPDIR/alignment.log"
{

    # Arguments used
    echo "Query genome: $QUERY_FA"
    echo "Query annotation: $QUERY_GTF"
    echo "Target genome: $TARGET_FA"
    echo "Target annotation: $TARGET_GTF"
    echo "Path to results directory: $WORKDIR"
    echo "Gene biotype to analyse: $BIOTYPE"
    echo "Flank value for liftoff: $FLANK"
    echo "Coverage cut-off for liftoff: $COVERAGE"
    echo "Sequence identity cut-off for liftoff: $IDENTITY"
    echo "Features database available: $DATABASE"

    ### STEP1 - liftoff

    echo "Step1: liftoff analysis - $(date +"%d.%m.%Y %H:%M")"
    echo "Working directory: $COMPDIR"

    ## select features to add to `-f`
    FEATURES=$COMPDIR/all_features.txt
    grep -v '#' $QUERY_GTF | cut -f3 | sort | uniq > $FEATURES
    echo "Features file created: $FEATURES"

    ## create outfile and outdir
    LIFTOFF_DIR="$COMPDIR"/liftoff_flank"$FLANK"
    LO_OUT=liftoff_"$QUERY_SP"_to_"$TARGET_SP"_flank"$FLANK".gtf

    ## launch liftoff script

    ### source environment
    . /local/env/envconda.sh
    conda activate scans_env

    ### run liftoff
    file_name=$(basename $QUERY_GTF .gtf)_"$BIOTYPE".gtf
    file_path="$WORKDIR"/"$file_name"

    if [ "$DATABASE" = FALSE ]; then
        # subset annotation according biotype

        echo "Feature database not existing. Creation of annotation file according biotype."
        if [[ "$BIOTYPE" == "lncRNA" ]] ; then
            # list of biotype string corresponding to lncRNA
            lnc_list='"lncRNA" "lincRNA" "antisense" "sense_overlapping" "sense_intronic" "lnc_RNA"'
            grep -E "biotype ($(echo $lnc_list | sed 's/ /|/g'))" $QUERY_GTF > "$file_path"
        elif [[ "$BIOTYPE" == "mRNA" ]] ; then
            # list of biotype string corresponding to lncRNA
            mrna_list='"mRNA" "protein_coding"'
            grep -E "biotype ($(echo $mrna_list | sed 's/ /|/g'))" $QUERY_GTF > "$file_path"
        else
            echo "Error: wrong biotype to analyse"
        fi

        echo "Subset annotation by biotype created: "$file_path""

        echo "liftoff -g $file_path -o $LO_OUT -p 16 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $LIFTOFF_DIR -u $LIFTOFF_DIR/unmapped_features.txt -copies"
        liftoff -g $file_path -o $LO_OUT -p 16 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $LIFTOFF_DIR -u $LIFTOFF_DIR/unmapped_features.txt -copies

    elif [ "$DATABASE" = TRUE ]; then
        QUERY_DB="$file_path"_db
        echo "liftoff -db $QUERY_DB -o $LO_OUT -p 16 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $LIFTOFF_DIR -u $LIFTOFF_DIR/unmapped_features.txt -copies"
        liftoff -db $QUERY_DB -o $LO_OUT -p 16 $TARGET_FA $QUERY_FA -flank $FLANK -f $FEATURES -dir $LIFTOFF_DIR -u $LIFTOFF_DIR/unmapped_features.txt -copies
        
    else
        echo "Error for feature database (boolean value needed)"    
    fi


    # move features file in liftoff dir
    mv $FEATURES $LIFTOFF_DIR
    mv $LO_OUT $COMPDIR

    LO_OUT2="$COMPDIR"/"$LO_OUT"

    ### STEP2 - Filter liftoff results according cut-off arguments
    echo "Step2: liftoff results filtering - $(date +"%d.%m.%Y %H:%M")"

    # filtered file
    LO_FILTER="$COMPDIR"/$(basename $LO_OUT2 .gtf)_filtered.gtf

    Rscript filter_liftoffOutput.R $LO_OUT2 $LO_FILTER $COVERAGE $IDENTITY $COMPDIR

    ### STEP3 - Bedtools intersect
    echo "Step3: bedtools analysis - $(date +"%d.%m.%Y %H:%M")"
    BEDTOOLS_DIR="$COMPDIR"/bedtools_intersect
    mkdir ${BEDTOOLS_DIR}

    ### create outfile for liftoff output, target annotation and bedtools output
    LO_BED="$BEDTOOLS_DIR"/$(basename $LO_FILTER .gtf).bed
    TARGET_BED="$BEDTOOLS_DIR"/"$TARGET_SP".bed
    FINAL_BED="$BEDTOOLS_DIR"/overlap_"$QUERY_SP"_to_"$TARGET_SP".bed

    #### format gtf in bed 
    . format_gtf2bed.sh $LO_FILTER > $LO_BED
    echo "Bed file created: $LO_BED"

    . format_gtf2bed.sh $TARGET_GTF > $TARGET_BED
    echo "Bed file created: $TARGET_BED"

    #### bedtools intersect - default fraction ~1bp

    echo "Intersection analysis done by bedtools intersect."
    bedtools intersect -a $LO_BED -b $TARGET_BED -wao -split > $FINAL_BED
    echo "bedtools intersect -a $LO_BED -b $TARGET_BED -wao -split > $FINAL_BED"

    ### STEP4 - Sequence alignment analysis
    echo "Step4: Sequence alignment analysis - $(date +"%d.%m.%Y %H:%M")"
    ALIGN_DIR="$COMPDIR"/alignment_analysis
    mkdir ${ALIGN_DIR}

    echo "Rscript seq_alignment_analysis.R $QUERY_GTF $TARGET_GTF $FINAL_BED $LIFTOFF_DIR/unmapped_features.txt $ALIGN_DIR"
    Rscript seq_alignment_analysis.R $QUERY_GTF $TARGET_GTF $FINAL_BED $LIFTOFF_DIR/unmapped_features.txt $ALIGN_DIR

    echo "Sequence alignment done: $(date +"%d.%m.%Y %H:%M")"

} >> "$LOG_FILE" 2>&1