#!/bin/bash -l

#Sample codes joined by '|'
smpls="SS2_19_182|SS2_19_184|SS2_19_186|SS2_18_197"
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
stdyCsSfx="eAdrnl4v2"

python format_quantification.py -s=$smpls -p=$dirSfx -q=$qcFldr -e=$exprssnFldr -E=$pthCmmnGnNmToENSEMBLG -S=$species -c=$stdyCsSfx
