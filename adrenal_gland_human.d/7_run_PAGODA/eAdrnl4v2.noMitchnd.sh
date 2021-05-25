#!/bin/bash -l

#Sample codes joined by "|"
smpls="SS2_19_182|SS2_19_184|SS2_19_186|SS2_18_197"
#Gene expression matrix obtained with 6_execute_celloline_n_select_HQ_cells
exprssnFl='quantification.d/HTseq.d/hg38.noIntrns.eAdrnl4v2highQC.allGns.counts.htsq.clltyHQ.noMitchnd.csv'
#Output folder parental. Will write all PAGODA files and results and is amimed to further fold additional results.
outFldrPrnt='eAdrnl4v2.d'
#Seed for tSNE
cSeed=49897
#Species
species='Human'

python PAGODA_mstr.py -s=$smpls -e=$exprssnFl -S=$species -o=$outFldrPrnt -r=$cSeed
