#!/bin/bash -l

#Sample codes joined by '|'
smpls="SS2_17_378|SS2_17_286|SS2_17_285|SS2_17_282|SS2_17_284|SS2_17_376|SS2_17_283|SS2_17_380|SS2_17_374|SS2_17_382|SS2_17_281"
#Study case sufix
dirSfx=".noIntrns."
#Folder with files with cell features obtained with 3_summarize_quality_control_features
qcFldr="mapping.d/STAR.d"
#Expression output folder, will write the formatted expression matrices
exprssnFldr="quantification.d/HTseq.d"
#File to convert gene symbols to ENSEMBL gene codes
pthCmmnGnNmToENSEMBLG="data/cmmnGnNmTOEnsmbl.hg38"
#Species
species="Human"
#Study case name
stdyCsSfx="NB"

python format_quantification.py -s=$smpls -p=$dirSfx -q=$qcFldr -e=$exprssnFldr -E=$pthCmmnGnNmToENSEMBLG -S=$species -c=$stdyCsSfx
