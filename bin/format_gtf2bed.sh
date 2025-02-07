#!/bin/bash

#Format a gtf file to an .bed6

# MARCH 2016 :
# ADD option -ignoreGroupsWithoutExons
# for gencode file, wo this option, it resulted in:
# no exons defined for group ENSG00000223972.4, feature gene (perhaps try -ignoreGroupsWithoutExons)
#
# in : /home/genouest/umr6061/recomgen/tderrien/DATA/hg19/annotation/gencode19
# gtf2bed12.sh gencode.v19.annotation.wochr.gtf > gencode.v19.annotation.wochr.bed
# no exons defined for group ENSG00000223972.4, feature gene (perhaps try -ignoreGroupsWithoutExons)



# USAGE

usage() {
    echo "#" >&2
    echo -e "# USAGE: `basename $0` <file.gtf> \n OR\n    cat <my_file> | `basename $0`">&2
}

BINDIR=~tderrien/progs/UCSC/exec/


infile=$1

##################################

if [ -p /dev/stdin ];then

    #from STDIN
    cat /dev/stdin |  ${BINDIR}/gtfToGenePred -ignoreGroupsWithoutExons stdin stdout | ${BINDIR}/genePredToBed stdin stdout


elif [ $# -eq  1 ];then

     ${BINDIR}/gtfToGenePred -ignoreGroupsWithoutExons $infile stdout | ${BINDIR}/genePredToBed stdin stdout

else
    echo "#Error! no argument  file or file empty !"
    usage;

fi
