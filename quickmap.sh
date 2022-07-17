#!/bin/bash

# #########################################################################
# This script is likely to only be used for mapping mutations 
# in special cases
# It calls bwa from the bwa folder in the home directory and
# samtools is required and a permanent path must be set for 
# this script to work
# those programs are the known dependencies
#
# NOTE: scripts run in shell must be allowed to be read, write and 
# be executable
# use case: chmod +x scriptfilename
# Created by S. Dean Rider Jr., June 2021 for the Leffak Lab
# 
# if you are a novice, see:
# https://www.ebi.ac.uk/sites/ebi.ac.uk/files/content.ebi.ac.uk/materials/2014/140217_AgriOmics/dan_bolser_snp_calling.pdf
# as it will help you understand this script
#
# #########################################################################


# #########################################################################
#
# set path as a precaution and exit script if anything fails
#
# #########################################################################

export PATH=$PATH:$HOME/bwa
set -e

# #########################################################################
#
# assign filenames to arguments passed in from maplist
#
# #########################################################################

referencefile="$1"
readsfile="$2"

# #########################################################################
#
# continue script with input parameters of reference and reads file names
#
# #########################################################################

echo "Mapping Run Progress Log" > $referencefile.$readsfile.log
echo Reference file is $referencefile >> $referencefile.$readsfile.log
echo Reads file is $readsfile >> $referencefile.$readsfile.log
date >> $referencefile.$readsfile.log
echo " " >> $referencefile.$readsfile.log

# #########################################################################
#
# index and align reads to reference if not already done
#
# #########################################################################

echo -e "\033[1;37m ██████╗ ██╗    ██╗ █████╗ \033[0m";
echo -e "\033[1;37m ██╔══██╗██║    ██║██╔══██╗\033[0m";
echo -e "\033[1;37m ██████╔╝██║ █╗ ██║███████║\033[0m";
echo -e "\033[1;37m ██╔══██╗██║███╗██║██╔══██║\033[0m";
echo -e "\033[1;37m ██████╔╝╚███╔███╔╝██║  ██║\033[0m";
echo -e "\033[1;37m ╚═════╝  ╚══╝╚══╝ ╚═╝  ╚═╝\033[0m";
echo -e "\033[1;37m                           \033[0m";
echo
echo Checking for index files
echo Checking for index files >> $referencefile.$readsfile.log

if [ -e $referencefile.bwt ]
then
echo -e "\033[1;37m ███████╗██╗  ██╗██╗██████╗     ██╗███╗   ██╗██████╗ ███████╗██╗  ██╗\033[0m";
echo -e "\033[1;37m ██╔════╝██║ ██╔╝██║██╔══██╗    ██║████╗  ██║██╔══██╗██╔════╝╚██╗██╔╝\033[0m";
echo -e "\033[1;37m ███████╗█████╔╝ ██║██████╔╝    ██║██╔██╗ ██║██║  ██║█████╗   ╚███╔╝ \033[0m";
echo -e "\033[1;37m ╚════██║██╔═██╗ ██║██╔═══╝     ██║██║╚██╗██║██║  ██║██╔══╝   ██╔██╗ \033[0m";
echo -e "\033[1;37m ███████║██║  ██╗██║██║         ██║██║ ╚████║██████╔╝███████╗██╔╝ ██╗\033[0m";
echo -e "\033[1;37m ╚══════╝╚═╝  ╚═╝╚═╝╚═╝         ╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝\033[0m";
echo -e "\033[1;37m                                                                     \033[0m";
  echo $referencefile.bwt found, Skipping Indexing of Reference
  echo $referencefile.bwt found, Skipping Indexing of Reference >> $referencefile.$readsfile.log
else
echo -e "\033[1;37m ██╗███╗   ██╗██████╗ ███████╗██╗  ██╗\033[0m";
echo -e "\033[1;37m ██║████╗  ██║██╔══██╗██╔════╝╚██╗██╔╝\033[0m";
echo -e "\033[1;37m ██║██╔██╗ ██║██║  ██║█████╗   ╚███╔╝ \033[0m";
echo -e "\033[1;37m ██║██║╚██╗██║██║  ██║██╔══╝   ██╔██╗ \033[0m";
echo -e "\033[1;37m ██║██║ ╚████║██████╔╝███████╗██╔╝ ██╗\033[0m";
echo -e "\033[1;37m ╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝\033[0m";
echo -e "\033[1;37m                                      \033[0m";
  echo Indexing reference $referencefile with bwa
  echo Indexing reference $referencefile with bwa >> $referencefile.$readsfile.log
  echo
  bwa index $referencefile
  echo Indexing completed successfully >> $referencefile.$readsfile.log
  date >> $referencefile.$readsfile.log
echo " " >> $referencefile.$readsfile.log
  echo
fi

echo
echo Checking for existing .sam files

if [ -e $referencefile.sam ]
then
echo -e "\033[1;37m ███████╗██╗  ██╗██╗██████╗     ███╗   ███╗ █████╗ ██████╗ ██████╗ ██╗███╗   ██╗ ██████╗ \033[0m";
echo -e "\033[1;37m ██╔════╝██║ ██╔╝██║██╔══██╗    ████╗ ████║██╔══██╗██╔══██╗██╔══██╗██║████╗  ██║██╔════╝ \033[0m";
echo -e "\033[1;37m ███████╗█████╔╝ ██║██████╔╝    ██╔████╔██║███████║██████╔╝██████╔╝██║██╔██╗ ██║██║  ███╗\033[0m";
echo -e "\033[1;37m ╚════██║██╔═██╗ ██║██╔═══╝     ██║╚██╔╝██║██╔══██║██╔═══╝ ██╔═══╝ ██║██║╚██╗██║██║   ██║\033[0m";
echo -e "\033[1;37m ███████║██║  ██╗██║██║         ██║ ╚═╝ ██║██║  ██║██║     ██║     ██║██║ ╚████║╚██████╔╝\033[0m";
echo -e "\033[1;37m ╚══════╝╚═╝  ╚═╝╚═╝╚═╝         ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝     ╚═╝╚═╝  ╚═══╝ ╚═════╝ \033[0m";
echo -e "\033[1;37m                                                                                         \033[0m";
  echo $referencefile.sam found, Skipping Mapping of Reads
  echo $referencefile.sam found, Skipping Mapping of Reads >> $referencefile.$readsfile.log
else
echo -e "\033[1;37m ███╗   ███╗ █████╗ ██████╗ \033[0m";
echo -e "\033[1;37m ████╗ ████║██╔══██╗██╔══██╗\033[0m";
echo -e "\033[1;37m ██╔████╔██║███████║██████╔╝\033[0m";
echo -e "\033[1;37m ██║╚██╔╝██║██╔══██║██╔═══╝ \033[0m";
echo -e "\033[1;37m ██║ ╚═╝ ██║██║  ██║██║     \033[0m";
echo -e "\033[1;37m ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     \033[0m";
echo -e "\033[1;37m                            \033[0m";
  echo Mapping reads $readsfile onto reference $referencefile with bwa mem
  echo Mapping reads $readsfile onto reference $referencefile with bwa mem >> $referencefile.$readsfile.log
  echo
  bwa mem $referencefile $readsfile > $referencefile.sam  #./bwa changed to bwa
  echo Mapping completed successfully >> $referencefile.$readsfile.log
  date >> $referencefile.$readsfile.log
echo " " >> $referencefile.$readsfile.log
fi

# #########################################################################
#
# Generate .fai, .dict files 
#
# #########################################################################

echo -e "\033[1;37m ███████╗ █████╗ ██╗\033[0m";
echo -e "\033[1;37m ██╔════╝██╔══██╗██║\033[0m";
echo -e "\033[1;37m █████╗  ███████║██║\033[0m";
echo -e "\033[1;37m ██╔══╝  ██╔══██║██║\033[0m";
echo -e "\033[1;37m ██║     ██║  ██║██║\033[0m";
echo -e "\033[1;37m ╚═╝     ╚═╝  ╚═╝╚═╝\033[0m";
echo -e "\033[1;37m                    \033[0m";
echo
echo Generating $referencefile .fai file for future use
samtools fqidx $referencefile -o $referencefile.fai
echo

echo -e "\033[1;37m ██████╗ ██╗ ██████╗████████╗\033[0m";
echo -e "\033[1;37m ██╔══██╗██║██╔════╝╚══██╔══╝\033[0m";
echo -e "\033[1;37m ██║  ██║██║██║        ██║   \033[0m";
echo -e "\033[1;37m ██║  ██║██║██║        ██║   \033[0m";
echo -e "\033[1;37m ██████╔╝██║╚██████╗   ██║   \033[0m";
echo -e "\033[1;37m ╚═════╝ ╚═╝ ╚═════╝   ╚═╝   \033[0m";
echo -e "\033[1;37m                             \033[0m";
echo Making .dict file for future use
samtools dict $referencefile > $referencefile.dict

echo .fai, .dict files completed successfully >> $referencefile.$readsfile.log
date >> $referencefile.$readsfile.log
echo " " >> $referencefile.$readsfile.log

# #########################################################################
#
# while you are here, might as well get some stats collected
# the mpileup is using old samtools, will need piped through to
# bcftools for use in VCF conversions
#
# #########################################################################

echo
echo Generating Stats
echo -e "\033[1;37m ██████╗  █████╗ ███╗   ███╗\033[0m";
echo -e "\033[1;37m ██╔══██╗██╔══██╗████╗ ████║\033[0m";
echo -e "\033[1;37m ██████╔╝███████║██╔████╔██║\033[0m";
echo -e "\033[1;37m ██╔══██╗██╔══██║██║╚██╔╝██║\033[0m";
echo -e "\033[1;37m ██████╔╝██║  ██║██║ ╚═╝ ██║\033[0m";
echo -e "\033[1;37m ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝\033[0m";
echo -e "\033[1;37m                            \033[0m";
samtools view -S -b $referencefile.sam > $referencefile.bam
echo -e "\033[1;37m ███████╗ ██████╗ ██████╗ ████████╗\033[0m";
echo -e "\033[1;37m ██╔════╝██╔═══██╗██╔══██╗╚══██╔══╝\033[0m";
echo -e "\033[1;37m ███████╗██║   ██║██████╔╝   ██║   \033[0m";
echo -e "\033[1;37m ╚════██║██║   ██║██╔══██╗   ██║   \033[0m";
echo -e "\033[1;37m ███████║╚██████╔╝██║  ██║   ██║   \033[0m";
echo -e "\033[1;37m ╚══════╝ ╚═════╝ ╚═╝  ╚═╝   ╚═╝   \033[0m";
echo -e "\033[1;37m                                   \033[0m";
samtools sort $referencefile.bam -o $referencefile.sorted.bam
wait
echo -e "\033[1;37m ██████╗ ███████╗██████╗ ████████╗██╗  ██╗\033[0m";
echo -e "\033[1;37m ██╔══██╗██╔════╝██╔══██╗╚══██╔══╝██║  ██║\033[0m";
echo -e "\033[1;37m ██║  ██║█████╗  ██████╔╝   ██║   ███████║\033[0m";
echo -e "\033[1;37m ██║  ██║██╔══╝  ██╔═══╝    ██║   ██╔══██║\033[0m";
echo -e "\033[1;37m ██████╔╝███████╗██║        ██║   ██║  ██║\033[0m";
echo -e "\033[1;37m ╚═════╝ ╚══════╝╚═╝        ╚═╝   ╚═╝  ╚═╝\033[0m";
echo -e "\033[1;37m                                          \033[0m";
samtools depth -d 2000000 $referencefile.sorted.bam > $referencefile.DepthOfCoverage.out
wait
echo
echo Adding Headers to Data files and removing some unwanted files
echo reference position depth > DepthHeader
cat DepthHeader $referencefile.DepthOfCoverage.out > $referencefile.DepthOfCoverage.txt
rm DepthHeader
rm $referencefile.DepthOfCoverage.out
echo -e "\033[1;37m ██╗ ██████╗ ██╗   ██╗\033[0m";
echo -e "\033[1;37m ██║██╔════╝ ██║   ██║\033[0m";
echo -e "\033[1;37m ██║██║  ███╗██║   ██║\033[0m";
echo -e "\033[1;37m ██║██║   ██║╚██╗ ██╔╝\033[0m";
echo -e "\033[1;37m ██║╚██████╔╝ ╚████╔╝ \033[0m";
echo -e "\033[1;37m ╚═╝ ╚═════╝   ╚═══╝  \033[0m";
echo -e "\033[1;37m                      \033[0m";
samtools index $referencefile.sorted.bam
echo -e "\033[1;37m ███████╗████████╗ █████╗ ████████╗███████╗\033[0m";
echo -e "\033[1;37m ██╔════╝╚══██╔══╝██╔══██╗╚══██╔══╝██╔════╝\033[0m";
echo -e "\033[1;37m ███████╗   ██║   ███████║   ██║   ███████╗\033[0m";
echo -e "\033[1;37m ╚════██║   ██║   ██╔══██║   ██║   ╚════██║\033[0m";
echo -e "\033[1;37m ███████║   ██║   ██║  ██║   ██║   ███████║\033[0m";
echo -e "\033[1;37m ╚══════╝   ╚═╝   ╚═╝  ╚═╝   ╚═╝   ╚══════╝\033[0m";
echo -e "\033[1;37m                                           \033[0m";
samtools idxstats $referencefile.sorted.bam > $referencefile.ContigLengthNumreads.txt
echo -e "\033[1;37m ██████╗ ██╗██╗     ███████╗██╗   ██╗██████╗ \033[0m";
echo -e "\033[1;37m ██╔══██╗██║██║     ██╔════╝██║   ██║██╔══██╗\033[0m";
echo -e "\033[1;37m ██████╔╝██║██║     █████╗  ██║   ██║██████╔╝\033[0m";
echo -e "\033[1;37m ██╔═══╝ ██║██║     ██╔══╝  ██║   ██║██╔═══╝ \033[0m";
echo -e "\033[1;37m ██║     ██║███████╗███████╗╚██████╔╝██║     \033[0m";
echo -e "\033[1;37m ╚═╝     ╚═╝╚══════╝╚══════╝ ╚═════╝ ╚═╝     \033[0m";
echo -e "\033[1;37m                                             \033[0m";
samtools mpileup -f $referencefile $referencefile.sorted.bam > $referencefile.Mpileup.txt

echo Stats estimated successfully >> $referencefile.$readsfile.log
date >> $referencefile.$readsfile.log
echo " " >> $referencefile.$readsfile.log

# #########################################################################
#
# Moving files into package
#
# #########################################################################
echo -e "\033[1;37m ██████╗  █████╗  ██████╗██╗  ██╗ █████╗  ██████╗ ███████╗\033[0m";
echo -e "\033[1;37m ██╔══██╗██╔══██╗██╔════╝██║ ██╔╝██╔══██╗██╔════╝ ██╔════╝\033[0m";
echo -e "\033[1;37m ██████╔╝███████║██║     █████╔╝ ███████║██║  ███╗█████╗  \033[0m";
echo -e "\033[1;37m ██╔═══╝ ██╔══██║██║     ██╔═██╗ ██╔══██║██║   ██║██╔══╝  \033[0m";
echo -e "\033[1;37m ██║     ██║  ██║╚██████╗██║  ██╗██║  ██║╚██████╔╝███████╗\033[0m";
echo -e "\033[1;37m ╚═╝     ╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝\033[0m";
echo -e "\033[1;37m                                                          \033[0m";
mkdir package.$referencefile.$readsfile
mv $referencefile.* ./package.$referencefile.$readsfile/
echo -e "\033[1;37m ███╗   ██╗███████╗██╗  ██╗████████╗\033[0m";
echo -e "\033[1;37m ████╗  ██║██╔════╝╚██╗██╔╝╚══██╔══╝\033[0m";
echo -e "\033[1;37m ██╔██╗ ██║█████╗   ╚███╔╝    ██║   \033[0m";
echo -e "\033[1;37m ██║╚██╗██║██╔══╝   ██╔██╗    ██║   \033[0m";
echo -e "\033[1;37m ██║ ╚████║███████╗██╔╝ ██╗   ██║   \033[0m";
echo -e "\033[1;37m ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝   ╚═╝   \033[0m";
echo -e "\033[1;37m                                    \033[0m";
# #########################################################################
#
# The End
#
# #########################################################################

exit 0