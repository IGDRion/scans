# How modules work

## VARIABLES

- CONFIG =
```
completeName,shortName,pathToGTF
Homo_sapiens,hsapiens,/projects/dog/data/hg38_GRCh38/annotation/Ensembl112/Homo_sapiens.GRCh38.112.gtf
Mus_musculus,mmusculus,/projects/dog/data/mm10/annotation/Mus_musculus.GRCm39.112.gtf
Canis_lupus_familiaris,cfamiliaris,/projects/dog/data/ROS_Cfam_1.0/annotation/ensembl/Canis_lupus_familiaris.ROS_Cfam_1.0.112.gtf
```

- ORTHOLOGY = DIRECTORY with orthology files named as `<shortNameQUERY>_<shortNameTARGET>_homology.tsv
```
ensembl_gene_id	hsapiens_homolog_ensembl_gene	hsapiens_homolog_orthology_type
ENSCAFG00845000002	ENSG00000081913	ortholog_one2one
ENSCAFG00845000003		
ENSCAFG00845000004		
```

## LAUNCH SCRIPTS
step1_format_input_files.sh $CONFIG $WORKDIR [-o $ORTHOLOGY]

step2_perform_synteny.sh $CONFIG $WORKDIR $ORTHOLOGY [-s $METHOD]
### only method1 works for now (default)

step3_perform_seq_alignment.sh --query_fa $QUERY_FA --query_gtf $QUERY_GTF \
--target_fa $TARGET_FA --target_gtf $TARGET_GTF --workdir $WORKDIR \
[--biotype $BIOTYPE --flank $FLANK --identity $IDENTITY]


WORKDIR/
    |- work
        |- input_data
        |   |- allMerged_gnInfo.tsv
        |   |- gnInfo
        |       |- <completeName>_gnInfo.tsv
        |
        |- method1
        |   |- lncBetwPcg
        |   |   |- <completeName>_lncRNAbetweenPcg.tsv
        |   |
        |   |- mergedSyntenyBySpecies
        |   |   |- <completeName>_syntenyMerged.tsv
        |   |
        |   |- syntenyByPair
        |       |- <shortNameQ>_<shortNameT>_synteny.tsv
        |
        |- method3
        |



## tests

### variables for steps 1 and 2
CONFIG=/projects/dog/aurore/BrownLincs/module/config.txt
WORKDIR=/projects/dog/aurore/BrownLincs/module/test1
ORTHOLOGY=/projects/dog/aurore/BrownLincs/module/test1/ortho_noformat

### scripts: format input + synteny method 1
sbatch step1_format_input_files.sh $CONFIG $WORKDIR -o $ORTHOLOGY
sbatch step2_perform_synteny.sh $CONFIG $WORKDIR $ORTHOLOGY

### variables for step 3 -> change input file names
QUERY_FA=$WORKDIR/data/hsap.fa
QUERY_GTF=$WORKDIR/data/hsap_chr22.gtf
TARGET_FA=$WORKDIR/data/cfam.fa
TARGET_GTF=$WORKDIR/data/cfam.gtf

ln -s /projects/dog/aurore/BrownLincs/data_GENCODE/GRCh38.primary_assembly.genome.fa $QUERY_FA
ln -s /projects/dog/aurore/LIFTOFF/MAMMALS/DATA/gencode.v47_chr22_PCG_lncRNA.gtf $QUERY_GTF
ln -s /projects/dog/data/canFam4/sequence/ensembl/Canis_lupus_familiaris-GCA_011100685.1-softmasked.fa $TARGET_FA
ln -s /projects/dog/data/canFam4/annotation/ensembl/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.108.gtf $TARGET_GTF

### script: sequence alignment
sbatch step3_perform_seq_alignment.sh --query_fa $QUERY_FA --query_gtf $QUERY_GTF --target_fa $TARGET_FA --target_gtf $TARGET_GTF --workdir $WORKDIR
