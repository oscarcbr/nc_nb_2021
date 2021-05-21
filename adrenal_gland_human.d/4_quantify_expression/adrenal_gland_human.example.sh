#!/bin/bash -l

#Sample code
sample="SS2_19_182"#This should be replace for each of the samples for human adrenal gland: {SS2_19_182, SS2_19_184, SS2_19_186, SS2_18_197}
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

