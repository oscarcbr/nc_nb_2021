#!/bin/bash -l

#Replace the next line by ech sample name #This should be replace for each of the samples for human adrenal gland: {SS2_19_182, SS2_19_184, SS2_19_186, SS2_18_197}
sample="SS2_19_182"
#Study case sufix
dirSfx=".noIntrns."
#Folder with bam files' folders generated with 1_map_reads_with_STAR and to output files with cell features
mapStarFldr="mapping.d/STAR.d"
#Folder with folders including fastaQC results
qcRdsFldr="qcReads.d"

python smmrzQC.py -s=$sample -p=$dirSfx -m=$mapStarFldr -q=$qcRdsFldr
