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
mstrDBfls='eAdrnl4v2.mstrDB.tsv'
#Repeat masker GTF file. Read scVelo tutorial for how to obtain it.
rmskFl='hg38_rmsk.gtf'
#Annotation GTF file
anntGTFFl='hg38.genecodeV28Comp.ERCCeGFP.cfflinks.noIntrns.gnNms.biotyp.gtf'
#Output loom merge file
outLoomMrgd='eAdrnl4v2.mrgd.loom'
#Prefix for the project
scPrfx='eAdrnl4v2.velocity'
#Sufix for the file name
sfxNm='noIntrns.'
#h5ad file with adata including the embedding of interest
h5adEmbddngFl='eAdrnl4v2.noMitchnd.KEql10_cnts.h5ad'
#Minimum shared counts
min_shared_counts=40
#Number of top genes to build velocity
n_top_genes=600
#Number of genes to plot
n_genes_plot=100
#Reference genes file
file_genes_ref='gnsRef.txt'
#Interest genes file
file_gene_interest='gnsIntrst.txt'
#Species
species='Human'
#Precursor marker. To execute Palantir
bslGn='ERBB3'


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
python3 velocyto_mstr.py -m=$mstrDBfls -r=$rmskFl -g=$anntGTFFl -l=$outLoomMrgd -p=$preprocs -v=$scvelo -i=$h5adEmbddngFl -s=$scPrfx -M=$n_top_genes -S=$min_shared_counts -P=$n_genes_plot -R=$file_genes_ref -I=$file_gene_interest -e=$species -a=$bslGn -G=TRUE
