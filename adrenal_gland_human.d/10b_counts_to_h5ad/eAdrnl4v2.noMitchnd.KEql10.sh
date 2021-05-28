#!/bin/bash -l

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

#Paths eAdrnl4v2.noMitchnd.KEql10
raw_cnts='data/eAdrnl4v2.noMitchnd.d/KEql10.d/eAdrnl4v2.noMitchnd.KEql10.cnts.csv'
clstrs_PAGODA='data/eAdrnl4v2.noMitchnd.d/KEql10.d/eAdrnl4v2.noMitchnd.KEql10.clstrsLbld.tsv'
depot_prefix='eAdrnl4v2.noMitchnd.KEql10'
#Files to make dictionaries
flSmplStg='data/index.smplsINSSstgs.txt'
flSmplOutcm='data/index.smplsOutcm.txt'
#Optional files to include tSNE from PAGODA
pagoda_tSNEFl='data/eAdrnl4v2.noMitchnd.d/KEql10.d/eAdrnl4v2.noMitchnd.KEql10.tsne.tsv'
flippingTsne=FALSE
flippingHrzntlTsne=TRUE


######################
#First methods to run#
######################
#If True will write a h5ad file with the basics counts and PAGODA embedding
w_bscs=TRUE

python3 counts_to_h5ad.py -b=$w_bscs -P=$pagoda_tSNEFl -F=$flippingTsne -r=$raw_cnts -c=$clstrs_PAGODA -d=$depot_prefix -s=$flSmplStg -O=$flSmplOutcm -H=$flippingHrzntlTsne
