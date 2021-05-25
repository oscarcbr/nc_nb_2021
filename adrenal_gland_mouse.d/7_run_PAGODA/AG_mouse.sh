#!/bin/bash -l

#Sample codes joined by "|"
smpls="SS2_17_343|SS2_17_344|SS2_18_111|SS2_18_114|SS2_18_116"
#Gene expression matrix obtained with 6_execute_celloline_n_select_HQ_cells
exprssnFl='quantification.d/HTseq.d/mm10.AG_mousehighQC.allGns.counts.htsq.clltyHQ.csv'
#Output folder parental. Will write all PAGODA files and results and is amimed to further fold additional results.
outFldrPrnt='AG_mouse.d'
#Seed for tSNE
cSeed=69804
#Species
species='Mouse'

python PAGODA_mstr.py -s=$smpls -e=$exprssnFl -S=$species -o=$outFldrPrnt -r=$cSeed
