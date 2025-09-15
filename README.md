# SCANS

## Overview
SCANS (Synteny Conservation Analysis for Non-coding sequences) is a tool allowing to study both lncRNA synteny and sequence conservation across species.
For each species to be compared, SCANS only requires a reference genome (.fa), a reference annotation (.gtf) containing gene biotype information (e.g., lncRNAs and protein-coding) and the protein-coding orthologs obtained for each pair of species under comparison. SCANS integrates two synteny-based methods that use one (permissive) or two (stringent) protein-coding anchors flanking lncRNAs in the query genome to search for their counterparts in the target genome ([Degalez, 2024](https://www.biorxiv.org/content/10.1101/2024.10.03.616473v1)). The tool also integrates a third complementary method, which directly aligns query lncRNAs to target genomes using a modified version of the Liftoff program ([Shumate and Salzberg, 2020](https://doi.org/10.1093/bioinformatics/btaa1016)).

![SCANS](https://github.com/user-attachments/assets/2a2b1625-7cb2-4bc3-8a68-3ee307c01abd)


## Installation

Clone scans repository
```
git clone https://github.com/IGDRion/scans.git
```
Go to scans repository and create conda environment
```
cd scans
bash create_SCANS_env.sh

```
Export SCANSPATH variable
```
export SCANSPATH=${PWD}
export PATH=$PATH:${SCANSPATH}/bin/
```

To test SCANS, a toy dataset is available in `toy_data` directory. Sequence files need to be decompressed before SCANS analysis. See the [tutorial](https://github.com/IGDRion/scans/wiki/Tutorial.md) for more details. 

## System requirements

SCANS requires 2.5GB of disk space for installation. By default, it requires 34 GB of RAM and 16 CPU threads to perform method 3 (sequence alignment) and only 4 CPU threads to perform method 2 (synteny). 

## Usage

SCANS is composed of 4 modules.

* `format_input_files.sh`: Format input files (annotation, orthology).
* `perform_synteny.sh`: Perform synteny analysis by method1, method2 or both.
* `perform_seq_alignment.sh`: Perform sequence alignment using liftoff.
* `process_output_results.sh`: Additional analysis (global results of all methods for one species + visualization) 

SCANS requires a configuration file as input containing the annotation (gtf format), the sequence (fasta format) and the name of each species to be compared:
```
completeName,shortName,pathToGTF,pathToFasta
Ectocarpus_sp7,ectosp7,/path/to/annotation/ectosp7.gtf,/path/to/sequence/Ectocarpus_sp7.fa
Ectocarpus_subulatus,esubulatus,/path/to/annotation/esubulatus.gtf,/path/to/sequence/Ectocarpus_subulatus.fa
Saccharina_latissima,slatissima,/path/to/annotation/slatissima.gtf,/path/to/sequence/Saccharina_latissima.fa
```

### Step1: format input files
```
usage: format_input_files.sh --config config --output WORKING_DIR [--orthology ORTHOLOGY_DIR]

Required input :
  --config                      configuration file with paths
  --output                      directory to save results

Options :
  --orthology ORTHOLOGY_DIR     path to directory containing orthology files to format for synteny analysis

```

**Note:**

Orthology files have to be concordant with annotation files (identical `gene_id`, same annotation version) and can be generated using tool as [`OrthoFinder`](https://github.com/davidemms/OrthoFinder) or directly downloaded from Ensembl database (Compara).

Orthology files in ORTHOLOGY_DIR must be formatted as  `<shortName_QUERY>_<shortName_TARGET>_homology.tsv` to be concordant with configuration file.

```
gene_id	<shortName_TARGET>_homolog	<shortName_TARGET>_homolog_orthology_type
Ec-05_002070	SL_01-10421	ortholog_one2one
Ec-05_001920	SL_01-106.1	ortholog_one2one
Ec-05_001910	SL_01-106.5	ortholog_one2one
Ec-05_001930	SL_01-107.5	ortholog_one2one
Ec-05_002080	SL_01-123.7	ortholog_one2one
Ec-05_002920	SL_01-12579	ortholog_one2many
```

### Step2: synteny analysis
```
usage: perform_synteny.sh --config config --output WORKING_DIR --orthology ORTHOLOGY_DIR [--synteny method1|method2|both]

Required input :
  --config config                     configuration file with paths
  --output WORKING_DIR                directory to save results
  --orthology ORTHOLOGY_DIR           directory containing formatted orthology files

Options :
  --synteny [method1|method2|both]    synteny method to perform (default: method1)

```

Feel free to see the [dedicated wiki page](https://github.com/IGDRion/scans/wiki/Understanding-how-methods-work#synteny-methods) for further explanations on methods used.

### Step3: sequence alignment analysis
```
usage: perform_seq_alignment.sh --config config --output WORKING_DIR

Required input :
  --config config              configuration file with paths
  --output WORKING_DIR         directory to save results

Options :
  --biotype [mRNA|lncRNA|all]  Biotype to analyze (default: lncRNA)
  --flank [0.0-1.0]            Flank fraction for liftoff (default: 0)
  --coverage [0.0-1.0]         Coverage fraction for liftoff (default: 0.5)
  --identity [0.0-1.0]         Sequence identity fraction for liftoff (default: 0)
```

Feel free to see the [dedicated wiki page](https://github.com/IGDRion/scans/wiki/Understanding-how-methods-work#sequence-alignment-method) for further explanations on methods used.

### Optional: additional output analysis
There are two optional output analysis available:
- perform summary analysis on the 3 method results
- visualize one-to-one orthologous lncRNAs results by species pair or by method for one species

##### Summary analysis option
```
usage: process_output_results.sh --summary --species shortName --resultsDir scans_results

Required input :
  --species shortName              Species shortName in config file
  --resultsDir scans_results       Directory path containing all output files from SCANS

```

##### Visualization option
```
usage: process_output_results.sh --upsetPlot [byPair|byMethod]

Required input for both `byPair|byMethod` options:
  --species shortName              Query shortName species in config file
  --resultsDir scans_results       Directory path containing all output files from SCANS
  --config config.txt              Path to config file used for analysis

Required input for `byPair` option (M1,M2,M3 with species1 vs species 2)
  --species2                       Target shortName species to perform comparison

Required input for `byMethod` option (one species vs all in one method)
  --method [met1|met1|met3]        Method choice to perform upsetPlot

```

Two visualization plots are available based on upset plot.

Using `byPair` option, you will obtain the number of one-to-one orthologous lncRNAs for one pair of species determined by each method separately or by several methods.

<img width="1000" height="800" alt="upsetPlot_ectosp7_vs_esubulatus_ortho1to1lnc_allMethods" src="https://github.com/user-attachments/assets/5d4d081d-e60d-4d21-849e-6c62fc444c00" />

Using `byMethod` option, you will obtain the number of one-to-one orthologous lncRNAs obtained by one choosen method for one species compared to all other species analyzed.

<img width="1500" height="1200" alt="upsetPlot_ectosp7_ortho1to1lnc_met1" src="https://github.com/user-attachments/assets/9146858e-7a91-4bd6-8571-5c83ce1bce5a" />

## Output

The output directory are presented as followed: 

```
WORKING_DIR
└── scans_results
    ├── input_data
    │   ├── allMerged_gnInfo.tsv
    │   ├── gnInfo
    │   │   └── <completeName>_gnInfo.tsv 
    │   └── orthology
    │       └── <shortNameQ>_<shortNameT>homology.tsv
    ├── method1
    │   ├── lncBetwPcg
    │   │   └── <completeName>_lncRNAbetweenPcg.tsv
    │   ├── mergedSyntenyBySpecies
    │   │   └── <completeName>_syntenyMerged.tsv
    │   └── syntenyByPair
    │       └── <shortNameQ>_<shortNameT>_synteny.tsv
    ├── method2
    │   ├── lncClassification
    │   │   └── <completeName>_lncConfiguration_feelncclassifier.tsv
    │   ├── mergedSyntenyBySpeciesFeelnc
    │   │   └── <completeName>_feelncMerged.tsv
    │   └── syntenyByPairFeelnc
    │       └── <shortNameQ>-<shortNameT>_lncConfigurationHomologyAggregated.tsv
    ├── method3
    │   ├── mergedseqAlignBySpecies
    │   │   └── <completeName>_seqAlignMerged.tsv
    │   └── <shortNameQ>
    │       └── <shortNameQ>_to_<shortNameT>
    │           ├── alignment_analysis
    │           │   ├── <shortNameQ>_to_<shortNameT>_mapped_knownGenes.txt
    │           │   ├── <shortNameQ>_to_<shortNameT>_mapped_unknownGenes.txt
    │           │   └── <shortNameQ>_to_<shortNameT>_unmapped_genes.txt
    │           ├── liftoff_<shortNameQ>_to_<shortNameT>_flank0.gtf
    │           ├── liftoff_<shortNameQ>_to_<shortNameT>_flank0_filtered.gtf
    │           └── liftoff_plot_coverage_seqID.png
    ├── plots
    └── results_summary
        ├── summaryByPair
        │   └── <shortNameQ>
        |       └── <shortNameQ>_<shortNameT>_globalResults.tsv
        └── <shortNameQ>_conservation.tsv

shortNameQ = query species shortName from config file (i.e. ectosp7)
shortNameT = target species shortName from config file (i.e. esubulatus)
completeName = species completeName from config file (i.e. Ectocarpus_sp7)
```

You will find additional information to understand output files on the [wiki page](https://github.com/IGDRion/scans/wiki/Understanding-output-files).


## License

This project is freely available under the GNU General Public License v3.0 (GPL-3.0). See the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.html) for details.

Developed by the IGDR (Institut de Génétique et Développement de Rennes) - UMR 6290 CNRS, University of Rennes
