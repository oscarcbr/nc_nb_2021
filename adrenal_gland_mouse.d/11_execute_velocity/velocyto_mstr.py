#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  velocyto_mstr.py
#  
#  Copyright 2019 Oscar C. Bedoya-Reina <oscarbed@ki.se>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

##################################
# Work on scanpy-based analysis  #                
##################################

import argparse,os
#Assert is running on python 3
import sys
assert sys.version_info.major>2#Assert is running on python 3


from numpy import array,inf,float32
from singlecell.velocyto import preprcss_mltplBAMFldrs,wrpr_scevlo

import scanpy as sc


#################################
#         Parse inputs          #                
#################################
def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def str2float(v):
	value = float32(v)
	return value
	
def str2int(v):
	value = int(v)
	return value
	


parser = argparse.ArgumentParser()
 
#input parameters for preprocessing 
parser.add_argument('-m','--mstrDBfls',help='Master DB file with columns: source, bam_root_fldr, and out_root_fldr',default=None)
parser.add_argument('-r','--rmskFl',help='Repeat masker GTF file',default=None)
parser.add_argument('-g','--anntGTFFl',help='Annotation GTF file',default=None)
parser.add_argument('-l','--outLoomMrgd',help='Output loom folder',default=None)
parser.add_argument('-i','--inFlh5ad',help='h5ad file with adata to plot',default=None)
parser.add_argument('-s','--scPrfx',help='prefix for the output plots',default=None)
parser.add_argument('-e','--species',help='Species',default=None)
#Programs to run
parser.add_argument('-p','--preprocs',help='If True will run preprocessing on BAM files',default=None,type=str2bool)
parser.add_argument('-v','--scvelo',help='If True will run scVelocity',default=None,type=str2bool)
#Optional parameters
parser.add_argument('-M', '--n_top_genes',help='If True will select the top number of genes to build the velocity',default=None,type=str2int)
parser.add_argument('-S', '--min_shared_counts',help='Minimumm number of share counts',default=2000,type=str2int)
parser.add_argument('-P', '--n_genes_plot',help='Number of genes to plot',default=100,type=str2int)
parser.add_argument('-R', '--file_genes_ref',help='Input file with genes of reference',default=None)
parser.add_argument('-I', '--file_gene_interest',help='Input file with genes of interest',default=None)
parser.add_argument('-a', '--mostAcntrlGn',help='Most ancestral gene names',default='RTTN')
#More additional parameters for bam files
parser.add_argument('-x', '--sfxNm',help='Sufix name',default='')
parser.add_argument('-A', '--sample',help='sample of interest',default=None)
parser.add_argument('-G', '--retain_genes',help='Switchg to retain genes in the list of genes of reference',default=False,type=str2bool)
#More additional extra parameters for embeddings
parser.add_argument('-E', '--embedding_basis',help='Embedding of interest',default='X_tsne_PAGODA_hc')
parser.add_argument('-B', '--color_basis',help='Clustering scheme of interest',default='PAGODA_hc')
parser.add_argument('-F', '--from_counts',help='Are data from counts?',default=True,type=str2bool)


args = parser.parse_args()

if args:
	#input parameters
	mstrDBfls = args.mstrDBfls
	rmskFl = args.rmskFl
	anntGTFFl = args.anntGTFFl
	outLoomMrgd = args.outLoomMrgd
	preprocs = args.preprocs
	scvelo = args.scvelo
	#input parameters for scvelo
	h5adEmbddngFl = args.inFlh5ad
	scPrfx = args.scPrfx
	#Filtering parameters
	n_top_genes = args.n_top_genes
	#Optional parameters
	min_shared_counts = args.min_shared_counts
	n_genes_plot = args.n_genes_plot
	file_genes_ref = args.file_genes_ref
	file_gene_interest = args.file_gene_interest
	species = args.species
	mostAcntrlGn = args.mostAcntrlGn
	#More additional parameters for bam files
	sfxNm = args.sfxNm
	sample = args.sample
	retain_genes = args.retain_genes
	#More additional extra parameters for embeddings
	clstr_basis = args.embedding_basis
	clr_basis = args.color_basis
	frmCnts = args.from_counts
	#Read files
	if file_genes_ref is not None:
		geneRefs=[l.splitlines()[0] for l in open(file_genes_ref,'r') if \
		l.strip()]
	else:
		geneRefs=['Mycn','Prrx1','Cldn11','Epas1','Mboat1','Ntrk2','Th', \
	'Muc3a', 'Pnmt', 'Sox6', 'Brca1', 'Sox5', 'Pias2','Ntrk1']
		if species=='Human':
			geneRefs=[g.upper() for g in geneRefs]
	if file_gene_interest is not None:
		lGnsIntrst=[l.splitlines()[0] for l in open(file_gene_interest, \
		'r') if l.strip()]
	else:
		lGnsIntrst=['Mycn','Prrx1','Cldn11','Epas1','Mboat1','Ntrk2','Th', \
	'Muc3a', 'Pnmt', 'Sox6', 'Brca1', 'Sox5', 'Pias2','Ntrk1']
		if species=='Human':
			lGnsIntrst=[g.upper() for g in lGnsIntrst]
	#retain_genes
	if retain_genes:
		retain_genes=sorted(set(lGnsIntrst).union(set(geneRefs)))
	else:
		retain_genes=None

print('Log info:')
print('\tmstrDBfls -> ',mstrDBfls)
print('\trmskFl -> ',rmskFl)
print('\tanntGTFFl -> ',anntGTFFl)
print('\toutLoomMrgd -> ',outLoomMrgd)
print('\tpreprocs -> ',preprocs)
print('\tscvelo -> ',scvelo)
print('\tinFlh5ad -> ',h5adEmbddngFl)
print('\tscPrfx -> ',scPrfx)
print('\tn_top_genes -> ',n_top_genes)
print('\tmin_shared_counts -> ',min_shared_counts)
print('\tn_genes_plot -> ',n_genes_plot)
print('\tfile_genes_ref -> ',file_genes_ref)
print('\tfile_gene_interest -> ',file_gene_interest)
print('\tgeneRefs -> ',geneRefs)
print('\tlGnsIntrst -> ',lGnsIntrst)
print('\tspecies -> ',species)
print('\tmostAcntrlGn -> ',mostAcntrlGn)
print('\tsfxNm -> ',sfxNm)
print('\tsample -> ',sample)
print('\tretain_genes -> ',retain_genes)
print('\tclstr_basis -> ',clstr_basis)
print('\tclr_basis -> ',clr_basis)

assert species in {'Human','Mouse'}

#################################################
#################################################
#################################################
#####   Make original links to BAM files   ######
#################################################
#################################################
#################################################


##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################
#Make a list of bam folders
lBAMFldrs,lPrfxNm,lOutFldr = [],[],[]
for l in open(mstrDBfls,'r'):
	if l.strip() and l[0]!='#':
		l=l.splitlines()[0].split('\t')
		source,bam_root_fldr,out_root_fldr,src_BAM_fdr,smpls = l
		if not os.path.exists(bam_root_fldr):
			os.mkdir(bam_root_fldr)
		#Make original links
		if not os.path.exists(bam_root_fldr) or not os.listdir(bam_root_fldr):
			smplsNms = sorted([s for s in smpls.split('|') if s.strip()])
			for smpl in smplsNms:
				mkLnks = False
				outBAMnFldrSmpl = os.path.join(bam_root_fldr,'%s.d'%smpl)
				#
				if not os.path.exists(outBAMnFldrSmpl):
					os.mkdir(outBAMnFldrSmpl)
					mkLnks = True
				elif not os.listdir(outBAMnFldrSmpl):
					mkLnks = True
				if mkLnks:
					srcBAMnFldrSmpl = os.path.join(src_BAM_fdr,'%s.%sd'% \
					(smpl,sfxNm))
					os.system('ln -sf %s/*.bam %s/.'%(srcBAMnFldrSmpl, \
					outBAMnFldrSmpl))
		#
		if sample is None:
			lBAMFldrNms = [f for f in os.listdir(bam_root_fldr)]
		else:
			lBAMFldrNms = [f for f in os.listdir(bam_root_fldr) if \
			f.find('%s.d'%sample)>-1]
		lPrfxNm.extend(['%s.%s'%(source,fldrNm.split('.d')[0]) for fldrNm in \
		lBAMFldrNms])
		lBAMFldrs.extend([os.path.join(bam_root_fldr,fldrNm) for fldrNm in \
		lBAMFldrNms])
		lOutFldr.extend([os.path.join(out_root_fldr,fldrNm) for fldrNm in \
		lBAMFldrNms])


#
#Make output folders
if not os.path.exists(out_root_fldr):
	os.mkdir(out_root_fldr)
#
for fldr in lOutFldr:
	if not os.path.exists(fldr):
		os.mkdir(fldr)
		
		

###########################
###########################
###########################
######   Switches   #######
###########################
###########################
###########################
#Scanpy's
rn_preprocs=0
if preprocs:
	rn_preprocs = 1 #If True will run veolcyto preprocessing on BAm files.

rn_scvelo = 0
if scvelo:
	rn_scvelo = 1#If True will scVelo 

##########################
##########################
##########################
######   Execute   #######
##########################
##########################
##########################


####################
####################
#Make adata set
####################
####################
if rn_preprocs:
	#Execute preprocessing
	preprcss_mltplBAMFldrs(lBAMFldrs,lPrfxNm,lOutFldr,rmskFl,anntGTFFl, \
	outLoomMrgd)


####################
####################
#Run scVelo
####################
####################
if rn_scvelo:
	#Execute preprocessing
	adata=sc.read(h5adEmbddngFl)
	if n_top_genes is None:
		wrpr_scevlo(adata,outLoomMrgd,scPrfx,mostAcntrlGn=mostAcntrlGn, \
		clstr_basis=clstr_basis,clr_basis=clr_basis,frmCnts=frmCnts)
	else:
		if geneRefs is not None:
			if lGnsIntrst is not None:
				wrpr_scevlo(adata,outLoomMrgd,scPrfx,min_shared_counts= \
				min_shared_counts,n_top_genes=n_top_genes,n_genes_plot=n_genes_plot, \
				geneRefs=geneRefs,lGnsIntrst=lGnsIntrst,mostAcntrlGn=mostAcntrlGn, \
				retain_genes=retain_genes,clstr_basis=clstr_basis, \
				clr_basis=clr_basis,frmCnts=frmCnts)
			else:
				wrpr_scevlo(adata,outLoomMrgd,scPrfx,min_shared_counts= \
				min_shared_counts,n_top_genes=n_top_genes,n_genes_plot=n_genes_plot, \
				geneRefs=geneRefs,mostAcntrlGn=mostAcntrlGn,retain_genes= \
				retain_genes,clstr_basis=clstr_basis,clr_basis=clr_basis, \
				frmCnts=frmCnts)
		else:
			if lGnsIntrst is not None:
				wrpr_scevlo(adata,outLoomMrgd,scPrfx,min_shared_counts= \
				min_shared_counts,n_top_genes=n_top_genes,n_genes_plot=n_genes_plot, \
				lGnsIntrst=lGnsIntrst,mostAcntrlGn=mostAcntrlGn,retain_genes= \
				retain_genes,clstr_basis=clstr_basis,clr_basis=clr_basis, \
				frmCnts=frmCnts)
			else:
				wrpr_scevlo(adata,outLoomMrgd,scPrfx,min_shared_counts= \
				min_shared_counts,n_top_genes=n_top_genes,n_genes_plot=n_genes_plot, \
				mostAcntrlGn=mostAcntrlGn,retain_genes=retain_genes, \
				clstr_basis=clstr_basis,clr_basis=clr_basis,frmCnts=frmCnts)
