#!/bin/bash -l

#Replace the next line by ech sample name #This should be replace for each of the samples for neuroblastoma: {SS2_17_378, SS2_17_286, SS2_17_285, SS2_17_282, SS2_17_284, SS2_17_376, SS2_17_283, SS2_17_380, SS2_17_374, SS2_17_382, SS2_17_281} 
sample="SS2_17_378"
#Study case sufix
dirSfx=".noIntrns."
#Folder with bam files' folders generated with 1_map_reads_with_STAR and to output files with cell features
mapStarFldr="mapping.d/STAR.d"
#Folder with folders including fastaQC results
qcRdsFldr="qcReads.d"

python smmrzQC.py -s=$sample -p=$dirSfx -m=$mapStarFldr -q=$qcRdsFldr
