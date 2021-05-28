#!/bin/bash -l

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

##########################
#Paths AG_mouse
##########################
#Raw counts from HTSeq. Generate with 6_execute_celloline_n_select_HQ_cells
rawCnts_PAGODAnrmlzd='mm10.AG_mousehighQC.allGns.counts.htsq.clltyHQ.csv'
#Genes of interest to plot and annotation
flGnsToPlt='mm10_gnIntrstToPlt.txt'
#Genes of interest to plot and annotation
smplsGrpFl='smpl_AG_mouse.database.tsv'
#Sample column (0-based) in smplsGrpFl files
smplsClmn=0
#Group column (0-based) in smplsGrpFl files
grpsClmn=2
##########################
#Other parameters
##########################
#Study case sufix
stdyCsSfx="AG_mouse"
#Number of clusters
lKmins='19'
#Folder with cluster annotation file
clstrAnnttnsFldr="clstrAnntns.d"
#Parent folder with cluster information obtained with PAGODA
PAGODAFldrPrnt="AG_mouse.d"
#Output parental folder
outFldrPrnt='mm10_compile_h5ad.d'
#Flip vertically input PAGODA tSNE
flippingTsne=FALSE
#Flip horizontally input PAGODA tSNE
flippingHrzntlTsne=TRUE

######################
#First methods to run#
######################
#If True will return the normalized data expression in PAGODA
rtrnNrmlzExprssn=TRUE
# If True will return the distances between cells from wcords in PAGODA
rtrnDstncs=TRUE
# If True will return the tSNE for PAGODA
rtrnTsnes=TRUE
#If true will make adata data structures for each PAGODA result file
mkAdataObjcts=FALSE

python compile_h5ad.py -r=$rawCnts_PAGODAnrmlzd -f=$flippingTsne -L=$smplsGrpFl -s=$smplsClmn -g=$grpsClmn -N=$rtrnNrmlzExprssn -D=$rtrnDstncs -T=$rtrnTsnes -M=$mkAdataObjcts -S=$stdyCsSfx -A=$clstrAnnttnsFldr -F=$PAGODAFldrPrnt -k=$lKmins -O=$outFldrPrnt -p=$flGnsToPlt -H=$flippingHrzntlTsne

######################
#     Second run     #
######################
#If True will return the normalized data expression in PAGODA
rtrnNrmlzExprssn=FALSE
# If True will return the distances between cells from wcords in PAGODA
rtrnDstncs=FALSE
# If True will return the tSNE for PAGODA
rtrnTsnes=FALSE
#If true will make adata data structures for each PAGODA result file
mkAdataObjcts=TRUE

python3 compile_h5ad.py -r=$rawCnts_PAGODAnrmlzd -f=$flippingTsne -L=$smplsGrpFl -s=$smplsClmn -g=$grpsClmn -N=$rtrnNrmlzExprssn -D=$rtrnDstncs -T=$rtrnTsnes -M=$mkAdataObjcts -S=$stdyCsSfx -A=$clstrAnnttnsFldr -F=$PAGODAFldrPrnt -k=$lKmins -O=$outFldrPrnt -p=$flGnsToPlt -H=$flippingHrzntlTsne
