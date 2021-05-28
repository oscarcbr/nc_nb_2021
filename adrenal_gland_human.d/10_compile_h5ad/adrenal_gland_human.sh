#!/bin/bash -l

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

##########################
#Paths eAdrnl4v2
##########################
#Raw counts from HTSeq. Generate with 6_execute_celloline_n_select_HQ_cells
rawCnts_PAGODAnrmlzd='hg38.noIntrns.eAdrnl4v2highQC.allGns.counts.htsq.clltyHQ.noMitchnd.csv'
#Genes of interest to plot and annotation
flGnsToPlt='hg38_gnIntrstToPlt.txt'
#Genes of interest to plot and annotation
smplsGrpFl='hg38_eAdrnl4v2.tsv'
#Sample column (0-based) in smplsGrpFl files
smplsClmn=0
#Group column (0-based) in smplsGrpFl files
grpsClmn=2
##########################
#Other parameters
##########################
#Study case sufix
stdyCsSfx="eAdrnl4v2.noMitchnd"
#Number of clusters
lKmins='10'
#Folder with cluster annotation file
clstrAnnttnsFldr="clstrAnntns.d"
#Parent folder with cluster information obtained with PAGODA
PAGODAFldrPrnt="eAdrnl4v2.d"
#Output parental folder
outFldrPrnt='hg38_compile_h5ad.d'
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
