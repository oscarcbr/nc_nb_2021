#!/bin/bash -l

#Sample code
sample="SS2_17_378"#This should be replace for each of the samples for neuroblastoma: {SS2_17_378, SS2_17_286, SS2_17_285, SS2_17_282, SS2_17_284, SS2_17_376, SS2_17_283, SS2_17_380, SS2_17_374, SS2_17_382, SS2_17_281} 
#Study case sufix
dirSfx=".noIntrns."
#Folder with bam files' folders generated with 1_map_reads_with_STAR and files with cell features obtained with 3_summarize_quality_control_features
mapStarFldr="mapping.d/STAR.d"
#Output folder for expression files
exprssnFldr="quantification.d/HTseq.d"
#GTF file
gtfFile="hg38.genecodeV28Comp.ERCCeGFP.cfflinks.noIntrns.gnNms.biotyp.gtf"
#Species
species="Human"

python quantification.py -s=$sample -p=$dirSfx -m=$mapStarFldr -e=$exprssnFldr -g=$gtfFile -S=$species

