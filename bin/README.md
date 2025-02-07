# How modules work

## VARIABLES

- CONFIG file format
```
completeName,shortName,pathToGTF
Homo_sapiens,hsapiens,/projects/dog/data/hg38_GRCh38/annotation/Ensembl112/Homo_sapiens.GRCh38.112.gtf
Mus_musculus,mmusculus,/projects/dog/data/mm10/annotation/Mus_musculus.GRCm39.112.gtf
Canis_lupus_familiaris,cfamiliaris,/projects/dog/data/ROS_Cfam_1.0/annotation/ensembl/Canis_lupus_familiaris.ROS_Cfam_1.0.112.gtf
```

- ORTHOLOGY directory with orthology files named as: `<shortNameQUERY>_<shortNameTARGET>_homology.tsv`
```
ensembl_gene_id	hsapiens_homolog_ensembl_gene	hsapiens_homolog_orthology_type
ENSCAFG00845000002	ENSG00000081913	ortholog_one2one
ENSCAFG00845000003		
ENSCAFG00845000004		
```

## LAUNCH SCRIPTS

```
step1_format_input_files.sh $CONFIG $WORKDIR [-o $ORTHOLOGY]

step2_perform_synteny.sh $CONFIG $WORKDIR $ORTHOLOGY [-s $METHOD]
# only method1 works for now (default)

step3_perform_seq_alignment.sh --query_fa $QUERY_FA --query_gtf $QUERY_GTF \
--target_fa $TARGET_FA --target_gtf $TARGET_GTF --workdir $WORKDIR \
[--biotype $BIOTYPE --flank $FLANK --identity $IDENTITY]
```


WORKDIR
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
        └── hsap_chr22_to_cfam
            ├── alignment_analysis
            │   ├── hsap_chr22_to_cfam_mapped_knownGenes.txt
            │   ├── hsap_chr22_to_cfam_mapped_knownTranscripts.txt
            │   ├── hsap_chr22_to_cfam_mapped_unknownGenes.txt
            │   ├── hsap_chr22_to_cfam_mapped_unknownTranscripts.txt
            │   └── hsap_chr22_to_cfam_unmapped_genes.txt
            ├── bedtools_intersect
            │   ├── cfam.bed
            │   ├── liftoff_hsap_chr22_to_cfam_flank0.bed
            │   └── overlap_hsap_chr22_to_cfam.bed
            └── liftoff_flank0
                ├── all_features.txt
                ├── liftoff_hsap_chr22_to_cfam_flank0.gtf
                ├── reference_all_genes.fa
                ├── reference_all_to_target_all.sam
                └── unmapped_features.txt

**Note:**
- For sequence alignment (method3), species names correspond to annotation file name
- could correspond to `shortName`

## TEST

```
# variables for steps 1 and 2
CONFIG=/projects/dog/aurore/BrownLincs/module/config.txt
WORKDIR=/projects/dog/aurore/BrownLincs/module/test1
ORTHOLOGY=/projects/dog/aurore/BrownLincs/module/test1/ortho_noformat

# scripts: format input + synteny method 1
sbatch step1_format_input_files.sh $CONFIG $WORKDIR -o $ORTHOLOGY
sbatch step2_perform_synteny.sh $CONFIG $WORKDIR $ORTHOLOGY

# variables for step 3 -> change input file names
QUERY_FA=$WORKDIR/data/hsap.fa
QUERY_GTF=$WORKDIR/data/hsap_chr22.gtf
TARGET_FA=$WORKDIR/data/cfam.fa
TARGET_GTF=$WORKDIR/data/cfam.gtf

# symbolic links to modify file name (!! this name will be used for output files !!)
ln -s /projects/dog/aurore/BrownLincs/data_GENCODE/GRCh38.primary_assembly.genome.fa $QUERY_FA
ln -s /projects/dog/aurore/LIFTOFF/MAMMALS/DATA/gencode.v47_chr22_PCG_lncRNA.gtf $QUERY_GTF
ln -s /projects/dog/data/canFam4/sequence/ensembl/Canis_lupus_familiaris-GCA_011100685.1-softmasked.fa $TARGET_FA
ln -s /projects/dog/data/canFam4/annotation/ensembl/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.108.gtf $TARGET_GTF

# script: sequence alignment
sbatch step3_perform_seq_alignment.sh --query_fa $QUERY_FA --query_gtf $QUERY_GTF --target_fa $TARGET_FA --target_gtf $TARGET_GTF --workdir $WORKDIR
```