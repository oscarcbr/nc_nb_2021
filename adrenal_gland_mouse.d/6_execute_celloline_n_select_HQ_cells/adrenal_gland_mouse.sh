#!/bin/bash -l

#Input sample codes joined by '|'
smpls="SS2_17_343|SS2_17_344|SS2_18_111|SS2_18_114|SS2_18_116"
#GTF file
gtfFile="data/mm10.genecodeV18Comp.ERCCYFPiCre.cfflinks.gnNms.biotyp.gtf"
#Folder with files with cell features obtained with 3_summarize_quality_control_features
qcFldr="mapping.d/STAR.d"
#Species
species="Mouse"
#Study case sufix
dirSfx="."
#Study case name
stdyCsSfx="AG_mouse"
#File with cell cycle genes
cllCycl="data/melanomaCllCycl_mm10.tsv"
#Run celloline 
rCllOlnFtrs=TRUE
#Folder to write/read celloline features
ftrsFldr="qcCells.d/STAR.d"
#Folder with Fastq files generated with 0_execute_cutAdapt_and_FASTAQC
fastQMainFldr="qcReads.d"
#Output folder to write full output features (including those calculated by celloline/cellity) and list of cells included/excluded
outFldr="qcCells.d/HTseq.d"
#Expression output folder, will write the filtered expression matrices
exprssnFldr="quantification.d/HTseq.d"
#Shuffle input samples
shuffleSmpls=FALSE
#Minimum number of genes detected allowed
ftrMinTrshld=2000
#Maximum number of genes detected allowed
ftrMaxTrshld=8000


python smplFltrNSlctHQclls.py -s=$smpls -g=$gtfFile -e=$exprssnFldr -q=$qcFldr -S=$species -p=$dirSfx -c=$stdyCsSfx -C=$cllCycl -o=$outFldr -m=$ftrMinTrshld -M=$ftrMaxTrshld -T=$rCllOlnFtrs -F=$shuffleSmpls -t=$ftrsFldr -Q=$fastQMainFldr
