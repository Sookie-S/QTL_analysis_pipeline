#!/bin/sh

###################################################################
#Script Name    : pipeline.sh                                                                                          
#Description    : this script runs the pipeline to identify QTL regions and genes overlapping
#Args           : 5 arguments are required (LOD_threshold LOD_support_interval conditionName ChromosomePrefix GFF_file_Name)                                                                                                                                                            
#Author         : Soukaina Timouma                                             
#Email          : soukaina.timouma@manchester.ac.uk                                           
###################################################################

if [ "$1" = "" ]; then
    echo "Error: missing arguments:"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi
if [ "$2" = "" ]; then
    echo "Error: missing arguments"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi
if [ "$3" = "" ]; then
    echo "Error: missing arguments:"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi

if [ "$4" = "" ]; then
    echo "Error: missing arguments"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi
if [ "$5" = "" ]; then
    echo "Error: missing arguments"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi

if [ "$6" = "" ]; then
    echo "Error: missing arguments"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi
if [ "$7" = "" ]; then
    echo "Error: missing arguments"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi
if [ "$8" = "" ]; then
    echo "Error: missing arguments"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi
if [ "$9" = "" ]; then
    echo "Error: missing arguments"
    echo "./pipeline.sh LOD_threshold LOD_support_interval conditionName Name_Parent1 Name_parent2 ChromosomePrefix_Parent1 GFF_file_Name_Parent1 ChromosomePrefix_Parent2 GFF_file_Name_Parent2"
    exit 1
    break
fi


echo "\n-----Parameters selected:"
echo "- LOD threshold selected: $1"
echo "- LOD support interval selected: $2"
echo "- Condition selected: $3"
echo "- Name Parent 1: $4"
echo "- Name Parent 2: $5"
echo "- ChromosomePrefix Parent 1: $6"
echo "- GFF files selected Parent 1: $7"
echo "- ChromosomePrefix Parent 2: $8"
echo "- GFF files selected Parent 2: $9"
echo "-----\n"



mkdir -p ../Results/
mkdir -p ../Results/QTL_regions/
mkdir -p ../Results/QTL_regions/$3/
mkdir -p ../Results/Gene_overlap/
mkdir -p ../Results/Gene_overlap/$3/

echo "\n-----Extraction of QTL regions:"
./extract_QTL_regions.py --lod $1 --int $2 --dir ../Files/$3/ --out ../Results/QTL_regions/$3/ --nam1 $4 --nam2 $5 --pre1 $6 --pre2 $8
echo "-----\n"
echo "\n-----Research of genes overlapping QTL regions:"
./overlap_QTL_regions_genes.py --dir ../Results/QTL_regions/$3/ --lod $1 --out ../Results/Gene_overlap/$3/ --nam1 $4 --nam2 $5 --gff1 ../Files/Genomes/GFF_files/$7 --gff2 ../Files/Genomes/GFF_files/$9
echo "-----\n"


echo "####### END OF PIPELINE #######\n"