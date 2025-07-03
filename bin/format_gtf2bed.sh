#!/bin/bash

#Format a gtf file to an .bed6

# USAGE

usage() {
    echo "#" >&2
    echo -e "# USAGE: `basename $0` <file.gtf> \n OR\n    cat <my_file> | `basename $0`">&2
}

infile=$1

##################################

if [ -p /dev/stdin ];then

    #from STDIN
    cat /dev/stdin |  gtfToGenePred -ignoreGroupsWithoutExons stdin stdout | genePredToBed stdin stdout


elif [ $# -eq  1 ];then

    gtfToGenePred -ignoreGroupsWithoutExons $infile stdout | genePredToBed stdin stdout

else
    echo "#Error! no argument  file or file empty !"
    usage;

fi
