#!/bin/bash

##########################################################################
# Get lncRNA classification from FEELnc classifier 
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
    echo -e "# USAGE: $PROGNAME <CONFIG_FILE> <WORK_DIR> "
    exit 1;
}

###########################
# VARIABLES
###########################

# input variable
CONFIG=$1
WORKDIR=$2

# Test number of parameters
###########################
if [[ $# -ne 2 ]]; then
	echo "ERROR : wrong number of arguments">&2
	# call usage function
	usage;
fi

echo "Config file analyzed: $CONFIG"
echo "Working directory: $WORKDIR"

OUTDIR="$WORKDIR"/scans_results/method2/lncClassification
mkdir -p $OUTDIR


###########################
# CODE PRINCIPAL
###########################

## check that general biotype will not bug for ensembl files (keep tx and not gene)
## refseq data tx/exon have transcript_biotype and not gene_biotype
#regEx_lncRNA='gene_biotype\s\"(lncRNA|lincRNA|sense_intronic|sense_exonic|sense_overlapping|antisense|lnc_RNA)\"'
#regEx_mRNA='gene_biotype\s\"(protein_coding|mRNA)\"'
regEx_lncRNA='biotype\s\"(lncRNA|lincRNA|sense_intronic|sense_exonic|sense_overlapping|antisense|lnc_RNA)\"'
regEx_mRNA='biotype\s\"(protein_coding|mRNA)\"'

# Extraction PCG and lncRNA
sed 1d $CONFIG | while IFS=',' read -r completeName shortName pathToGTF pathToFasta;
do
    echo $pathToGTF ;
    grep -v "#" $pathToGTF | grep -P $regEx_lncRNA > ${OUTDIR}/${completeName}_LNCextracted.tmp.gtf;
    grep -v "#" $pathToGTF | grep -P $regEx_mRNA > ${OUTDIR}/${completeName}_mRNAextracted.tmp.gtf;
done

# Classification of lncRNA relative to mRNA
mkdir ${OUTDIR}/logs

while IFS=',' read -r completeName shortName pathToGTF pathToFasta; do
    output_name=$(basename "$pathToGTF" .gtf)
    LNC="${OUTDIR}/${completeName}_LNCextracted.tmp.gtf"
    MRNA="${OUTDIR}/${completeName}_mRNAextracted.tmp.gtf"
    LOG="${OUTDIR}/logs/${completeName}_feelncclassifier.log"
    OUT="${OUTDIR}/${completeName}_classes_feelncclassifier.txt"

    echo "Running FEELnc_classifier.sh for $completeName"
    bash FEELnc_classifier.sh "$LNC" "$MRNA" "$LOG" "$OUT" > "$LOG" 2>&1

    if [ $? -eq 0 ]; then
        echo "Job for $completeName completed successfully."
    else
        echo "Job for $completeName failed. Check log: $LOG"
    fi

done < <(sed 1d "$CONFIG")


# Extract gene level classification from transcript level classification
echo "All jobs completed. Launching classification at gene level."
echo "Rscript getGeneLevelClassification.R "${CONFIG}" "${OUTDIR}""
Rscript getGeneLevelClassification.R "${CONFIG}" "${OUTDIR}"