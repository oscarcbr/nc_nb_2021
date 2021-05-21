#!/bin/bash -l

#Sample code
sample="SS2_17_343"#This should be replace for each of the samples for mouse adrenal gland: {SS2_17_343, SS2_17_344, SS2_18_111, SS2_18_114, SS2_18_116}
#Study case sufix
dirSfx="."
#Folder with bam files' folders generated with 1_map_reads_with_STAR and files with cell features obtained with 3_summarize_quality_control_features
mapStarFldr="mapping.d/STAR.d"
#Output folder for expression files
exprssnFldr="quantification.d/HTseq.d"
#GTF file
gtfFile="mm10.genecodeV18Comp.ERCCYFPiCre.cfflinks.gnNms.biotyp.gtf"
#Species
species="Mouse"

python quantification.py -s=$sample -p=$dirSfx -m=$mapStarFldr -e=$exprssnFldr -g=$gtfFile -S=$species


