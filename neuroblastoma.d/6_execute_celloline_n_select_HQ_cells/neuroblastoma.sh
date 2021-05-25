#!/bin/bash -l

#Input sample codes joined by '|'
smpls="SS2_17_282|SS2_17_378|SS2_17_286|SS2_17_285|SS2_17_284|SS2_17_376|SS2_17_283|SS2_17_380|SS2_17_374|SS2_17_382|SS2_17_281"
#GTF file
gtfFile="data/hg38.genecodeV28Comp.ERCCeGFP.cfflinks.noIntrns.gnNms.biotyp.gtf"
#Folder with files with cell features obtained with 3_summarize_quality_control_features
qcFldr="mapping.d/STAR.d"
#Species
species="Human"
#Study case sufix
dirSfx=".noIntrns."
#Study case name
stdyCsSfx="NB"
#File with cell cycle genes
cllCycl="data/melanomaCllCycl.tsv"
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
