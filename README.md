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
conda env create -f scans.yml
```
Export SCANSPATH variable
```
export SCANSPATH=${PWD}

export PATH=$PATH:${SCANSPATH}/bin/
```

## Usage

SCANS is composed of 3 modules.

* step1_format_input_files.sh	: Format input files (annotation, orthology).
* step2_perform_synteny.sh	: Perform synteny analysis by method1, method2 or both.
* step3_perform_seq_alignment.sh: Perform sequence alignment using liftoff.

SCANS requires a configuration file as input containing the annotation (gtf format), the sequence (fasta format) and the name of each species to be compared:
```
completeName,shortName,pathToGTF,pathToFasta
Ectocarpus_sp7,ectosp7,/path/to/annotation/ectosp7.gtf,/path/to/sequence/Ectocarpus_sp7.fa
Ectocarpus_subulatus,esubulatus,/path/to/annotation/esubulatus.gtf,/path/to/sequence/Ectocarpus_subulatus.fa
Saccharina_latissima,slatissima,/path/to/annotation/slatissima.gtf,/path/to/sequence/Saccharina_latissima.fa
```

### Step1: format input files
```
usage: step1_format_input_files.sh config WORKING_DIR [-o ORTHOLOGY_DIR]

Required input :
  config              configuration file with paths
  WORKING_DIR         directory to save results

Options :
  -o ORTHOLOGY_DIR    path to directory containing orthology files to format for synteny analysis

```
**Note:**
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
usage: step2_perform_synteny.sh --config config --workdir WORKING_DIR --orthology ORTHOLOGY_DIR [--synteny method1|method2|both]

Required input :
  --config config                     configuration file with paths
  --workdir WORKING_DIR               directory to save results
  --orthology ORTHOLOGY_DIR           directory containing formatted orthology files

Options :
  --synteny [method1|method2|both]    synteny method to perform (default: method1)

```
### Step3: sequence alignment analysis
```
usage: step3_perform_seq_alignment.sh --config config --workdir WORKING_DIR

Required input :
  --config config              configuration file with paths
  --workdir WORKING_DIR        directory to save results

Options :
  --biotype [mRNA|lncRNA|all]  Biotype to analyze (default: lncRNA)
  --flank [0.0-1.0]            Flank fraction for liftoff (default: 0)
  --coverage [0.0-1.0]         Coverage fraction for liftoff (default: 0.5)
  --identity [0.0-1.0]         Sequence identity fraction for liftoff (default: 0)
```

## Output

```
WORKING_DIR
└── work
    ├── input_data
    │   ├── allMerged_gnInfo.tsv
    │   └── gnInfo
    │       ├── <completeName>_gnInfo.tsv
    │       ├── Canis_lupus_familiaris_gnInfo.tsv
    │       ├── Homo_sapiens_gnInfo.tsv
    │       └── Mus_musculus_gnInfo.tsv
    ├── method1
    │   ├── lncBetwPcg
    │   │   ├── <completeName>_lncRNAbetweenPcg.tsv
    │   │   ├── Canis_lupus_familiaris_lncRNAbetweenPcg.tsv
    │   │   ├── Homo_sapiens_lncRNAbetweenPcg.tsv
    │   │   └── Mus_musculus_lncRNAbetweenPcg.tsv
    │   ├── mergedSyntenyBySpecies
    │   │   ├── <completeName>_syntenyMerged.tsv
    │   │   ├── Canis_lupus_familiaris_syntenyMerged.tsv
    │   │   ├── Homo_sapiens_syntenyMerged.tsv
    │   │   └── Mus_musculus_syntenyMerged.tsv
    │   └── syntenyByPair
    │       ├── <shortNameQ>_<shortNameT>_synteny.tsv
    │       ├── cfamiliaris_hsapiens_synteny.tsv
    │       ├── cfamiliaris_mmusculus_synteny.tsv
    │       ├── hsapiens_cfamiliaris_synteny.tsv
    │       ├── hsapiens_mmusculus_synteny.tsv
    │       ├── mmusculus_cfamiliaris_synteny.tsv
    │       └── mmusculus_hsapiens_synteny.tsv
    └── method3
        ├── hsap
        │    ├── hsap_to_cfam_lncRNA
        │    │    ├── alignment_analysis
        │    │    │   ├── hsap_to_cfam_mapped_knownGenes.txt
        │    │    │   ├── hsap_to_cfam_mapped_unknownGenes.txt
        │    │    │   └── hsap_to_cfam_unmapped_genes.txt
        │    │    ├── bedtools_intersect
        │    │    │   └── overlap_hsap_to_cfam.bed
        │    │    └── liftoff_flank0
        │    │        ├── liftoff_hsap_to_cfam_flank0.gtf
        │    │        └── unmapped_features.txt
        │    └── hsap_to_mmus_lncRNA
        └── mmus
            ├── mmus_to_cfam_lncRNA
            └── mmus_to_hsap_lncRNA
```

## License

This project is freely available under the GNU General Public License v3.0 (GPL-3.0). See the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.html) for details.

Developed by the IGDR (Institut de Génétique et Développement de Rennes) - UMR 6290 CNRS, University of Rennes
