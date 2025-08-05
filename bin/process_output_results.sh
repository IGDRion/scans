#!/bin/bash

##########################################################################
# Process scans results: results summary , figures 
# Developed by: A. Besson
##########################################################################

PROGNAME=$(basename $0)

# Function to display usage
usage() {
    echo "Usage: $PROGNAME --summary/--upsetPlot [arguments specific to choosen mode]"
    echo " --summary -> to generate a global results file for one species on results obtained in the 3 methods:"
    echo "      --species         shortName species to analyze (relative to shortName in config file)"
    echo "      --resultsDir         path to `scans_results`"
    echo " --upsetPlot -> to generate upsetPlot (byPair: M1,M2,M3 with species1 vs species 2 OR byMethod: one species vs all in one method)"
    echo "      --species         shortName species to analyze (relative to shortName in config file)"
    echo "      --species2        target shortName species to perform comparison"
    echo "      --resultsDir         path to `scans_results`"
    echo "      --method [met1|met1|met3]   method choice to perform upsetPlot"
    exit 1
}

###########################
# VARIABLES
###########################

MODE=""
SPECIES=""
WORKDIR=""
SPECIES2=""

while [[ $# -gt 0 ]]; do
  key="$1"
  
  case $key in
    --summary)
      MODE="summary"
      shift 
      ;;
    --upsetPlot)
      MODE="upsetPlot"
      shift
      ;;
    --species)
      SPECIES="$2"
      shift; shift
      ;;
    --resultsDir)
      WORKDIR="$2"
      shift; shift
      ;;
    --species2)
      SPECIES2="$2"
      shift; shift
      ;;
    *)
      echo "Unknown argument : $1"
      exit 1
      ;;
  esac
done

# Check option mode and launch associated script
if [[ "$MODE" == "summary" ]]; then
  if [[ -z "$SPECIES" || -z "$WORKDIR" ]]; then
    echo "Missing arguments for --summary"
    usage
    exit 1
  fi
  echo "Global summary analysis of $SPECIES - IN PROGRESS."
  Rscript getGlobalSummary.R "$SPECIES" "$WORKDIR"

elif [[ "$MODE" == "upsetPlot" ]]; then
  if [[ -z "$SPECIES" || -z "$SPECIES2" || -z "$WORKDIR" ]]; then
    echo "Missing arguments for --upsetPlot"
    usage
    exit 1
  fi
  echo "Upset plot generation - IN PROGRESS."
  Rscript getUpsetPlot.R "$SPECIES" "$SPECIES2" "$WORKDIR"

else
  echo "You must use --summary or --upsetPlot"
  usage
  exit 1
fi
