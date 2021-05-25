#!/bin/bash -l

#Sample codes joined by "|"
smpls='SS2_17_281|SS2_17_282|SS2_17_283|SS2_17_284|SS2_17_285|SS2_17_286|SS2_17_374|SS2_17_376|SS2_17_378|SS2_17_380|SS2_17_382'
#Gene expression matrix obtained with 6_execute_celloline_n_select_HQ_cells
exprssnFl='quantification.d/HTseq.d/hg38.noIntrns.NBhighQC.allGns.counts.htsq.clltyHQ.noMitchnd.csv'
#Output folder parental. Will write all PAGODA files and results and is amimed to further fold additional results.
outFldrPrnt='NB.d'
#Seed for tSNE
cSeed=70777
#Species
species='Human'

python PAGODA_mstr.py -s=$smpls -e=$exprssnFl -S=$species -o=$outFldrPrnt -r=$cSeed
