#!/bin/bash -l

#Replace the next line by ech sample name #This should be replace for each of the samples for mouse adrenal gland: {SS2_17_343, SS2_17_344, SS2_18_111, SS2_18_114, SS2_18_116}
sample="SS2_17_343"
#Study case sufix
dirSfx="."
#Folder with bam files' folders generated with 1_map_reads_with_STAR and to output files with cell features
mapStarFldr="mapping.d/STAR.d"
#Folder with folders including fastaQC results
qcRdsFldr="qcReads.d"

python smmrzQC.py -s=$sample -p=$dirSfx -m=$mapStarFldr -q=$qcRdsFldr
