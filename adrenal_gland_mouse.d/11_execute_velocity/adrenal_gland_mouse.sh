#!/bin/bash -l

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

#Inputs for preprocessing
#Master DB file with columns: source, bam_root_fldr, and out_root_fldr
mstrDBfls='AG_mouse.KEql19.mstrDB.tsv'
#Repeat masker GTF file. Read scVelo tutorial for how to obtain it.
rmskFl='mm10_rmsk.gtf'
#Annotation GTF file
anntGTFFl='mm10.genecodeV18Comp.ERCCYFPiCre.cfflinks.gnNms.biotyp.gtf'
#Output loom merge
outLoomMrgd='AG_mouse.mrgd.loom'
#Prefix for the project
scPrfx='AG_mouse.KEql19.velocity'
#Sufix for the file name
sfxNm='.'
#h5ad file with adata including the embedding of interest
h5adEmbddngFl='AG_mouse.KEql19_cnts.h5ad'
#Minimum shared counts
min_shared_counts=3000
#Number of top genes to build velocity
n_top_genes=1000
#Number of genes to plot
n_genes_plot=100
#Species
species='Mouse'
#Reference genes file
file_genes_ref='gnsRef.mm10.txt'
#Interest genes file
file_gene_interest='gnsIntrst.mm10.txt'
#Precursor marker. To execute Palantir
bslGn='Erbb3'


######################
#First methods to run#
######################
# If True will run preprocessing
preprocs=TRUE
# If True will run scVelo
scvelo=TRUE

######################
#      Execute       #
######################
python3 velocyto_mstr.py -m=$mstrDBfls -r=$rmskFl -g=$anntGTFFl -l=$outLoomMrgd -p=$preprocs -v=$scvelo -i=$h5adEmbddngFl -s=$scPrfx -M=$n_top_genes -S=$min_shared_counts -P=$n_genes_plot -R=$file_genes_ref -I=$file_gene_interest -e=$species -a=$bslGn

