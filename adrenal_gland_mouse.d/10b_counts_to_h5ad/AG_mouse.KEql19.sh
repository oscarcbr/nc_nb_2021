#!/bin/bash -l

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

#Paths AG_mouse.KEql19
raw_cnts='data/AG_mouse.d/KEql19.d/AG_mouse.KEql19.cnts.csv'
clstrs_PAGODA='data/AG_mouse.d/KEql19.d/AG_mouse.KEql19.clstrsLbld.tsv'
depot_prefix='AG_mouse.KEql19'
#Files to make dictionaries
flSmplStg='data/index.smplsINSSstgs.txt'
flSmplOutcm='data/index.smplsOutcm.txt'
#Optional files to include tSNE from PAGODA
pagoda_tSNEFl='data/AG_mouse.d/KEql19.d/AG_mouse.KEql19.tsne.tsv'
flippingTsne=TRUE
flippingHrzntlTsne=TRUE


######################
#First methods to run#
######################
#If True will write a h5ad file with the basics counts and PAGODA embedding
w_bscs=TRUE

python3 counts_to_h5ad.py -b=$w_bscs -P=$pagoda_tSNEFl -F=$flippingTsne -r=$raw_cnts -c=$clstrs_PAGODA -d=$depot_prefix -s=$flSmplStg -O=$flSmplOutcm -H=$flippingHrzntlTsne
