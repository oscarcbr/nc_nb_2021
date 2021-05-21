#!/bin/bash -l

#Sample codes joined by '|'
smpls="SS2_17_343|SS2_17_344|SS2_18_111|SS2_18_114|SS2_18_116"
#Study case sufix
dirSfx="."
#Folder with files with cell features obtained with 3_summarize_quality_control_features
qcFldr="mapping.d/STAR.d"
#Expression output folder, will write the formatted expression matrices
exprssnFldr="quantification.d/HTseq.d"
#File to convert gene symbols to ENSEMBL gene codes
pthCmmnGnNmToENSEMBLG="data/cmmnGnNmTOEnsmbl.mm10"
#Species
species="Mouse"
#Study case name
stdyCsSfx="AG_mouse"

python format_quantification.py -s=$smpls -p=$dirSfx -q=$qcFldr -e=$exprssnFldr -E=$pthCmmnGnNmToENSEMBLG -S=$species -c=$stdyCsSfx
