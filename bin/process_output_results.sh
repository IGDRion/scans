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
    echo "      --resultsDir      path to \`scans_results\`"
    echo " --upsetPlot [byPair|byMethod] -> to generate upsetPlot (byPair: M1,M2,M3 with species1 vs species 2 OR byMethod: one species vs all in one method)"
    echo "      --species         shortName species to analyze (relative to shortName in config file) [requested for: byPair,byMethod]"
    echo "      --species2        target shortName species to perform comparison for byPair [requested for: byPair]"
    echo "      --method [met1|met1|met3]   method choice to perform upsetPlot [requested for: byMethod]"
    echo "      --resultsDir      path to \`scans_results\` [requested for: byPair,byMethod]"
    echo "      --config          path to config file used for analysis [requested for: byPair,byMethod]"

    exit 1
}

###########################
# VARIABLES
###########################

MODE=""
UPSETPLOT=""
SPECIES=""
WORKDIR=""
SPECIES2=""
METHOD=""
CONFIG=""

while [[ $# -gt 0 ]]; do
  key="$1"
  
  case $key in
    --summary)
      MODE="summary"
      shift 
      ;;
    --upsetPlot)
      MODE="upsetPlot"
      UPSETPLOT="$2"
      shift; shift
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
    --method)
      METHOD="$2"
      shift; shift
      ;;
    --config)
      CONFIG="$2"
      shift; shift
      ;;
    *)
      echo "Unknown argument : $1"
      usage
      ;;
  esac
done

# Check option mode and launch associated script
if [[ "$MODE" == "summary" ]]; then
  if [[ -z "$SPECIES" || -z "$WORKDIR" ]]; then
    echo "Missing arguments for --summary"
    usage
  fi
  echo "Global summary analysis of $SPECIES - IN PROGRESS."
  getGlobalSummary.R "$SPECIES" "$WORKDIR"

elif [[ "$MODE" == "upsetPlot" ]]; then
  if [[ -z "$UPSETPLOT" || -z "$SPECIES" || -z "$WORKDIR" || -z "$CONFIG" ]]; then
    echo "Missing arguments for --upsetPlot"
    usage
  fi
  echo "Upset plot generation - IN PROGRESS."
  if [[ "$UPSETPLOT" == "byPair" ]]; then
    if [[ -z "$SPECIES2" ]]; then
      echo "Missing requested argument for --upsetPlot byPair"
      usage
    fi
    getUpsetPlot.R "$UPSETPLOT" "$SPECIES" "$SPECIES2" "$WORKDIR" "$CONFIG"

  elif [[ "$UPSETPLOT" == "byMethod" ]]; then
    if [[ -z "$METHOD" ]]; then
      echo "Missing requested argument for --upsetPlot byMethod"
      usage
    fi

    case "$METHOD" in met1|met2|met3) ;; *)
      echo "Error: --method must be met1, met2 or met3."
      usage
      ;;
    esac
    getUpsetPlot.R "$UPSETPLOT" "$SPECIES" "$METHOD" "$WORKDIR" "$CONFIG"
  else
    echo "You must choose byPair or byMethod for --upsetPlot option"
    usage
  fi

else
  echo "You must use --summary or --upsetPlot"
  usage
fi
