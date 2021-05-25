#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  PAGODA_mstr.py
#  
#  Copyright 2018 Oscar C. Bedoya Reina <oscarbed@ki.se>
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

"""
This script is a wrapper to execute the PAGODA pipeline
"""


import argparse,os,sys

from singlecell.PAGODA_mod import runPAGODA,cntCllsTsneRDS

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

parser = argparse.ArgumentParser()
 
######################
#  Input parameters  #
######################
parser.add_argument('-s','--smpls',help='Input samples',default=None)
parser.add_argument('-e','--exprssnFl',help='Expression file',default=None)
parser.add_argument('-S','--species',help='Species',default=None)
parser.add_argument('-o','--outFldrPrnt',help='Output folder parental',default=None)
parser.add_argument('-k','--lKmins',help='Input number of K clusters',default='5|8|9|10|11|12|13')
parser.add_argument('-G','--smplsGrpFl',help='File with sample\tgroup (optional)',default=None)
parser.add_argument('-E','--ercc_file',help='ercc_file',default='ERCC_mix1.csv')
parser.add_argument('-c','--nCores',help='number of cores',default=1,type=int)
parser.add_argument('-F','--minNonfailed',help='minNonfailed PAGODA parameter',default=5,type=int)
parser.add_argument('-T','--minCountThreshold',help='minCountThreshold PAGODA parameter',default=2,type=int)
parser.add_argument('-P','--maxModelPlots',help='maxModelPlots PAGODA parameter',default=10,type=int)
parser.add_argument('-A','--maxAdjVar',help='maxAdjVar PAGODA parameter',default=5,type=int)
parser.add_argument('-C','--rnCllCycl_crrctn',help='Run hard cell cycle correction (regressing on GO:0007049 and GO:0051301 ontologies)',default=True,type=str2bool)
parser.add_argument('-r','--cSeed',help='Seed for tSNE',default=1,type=int)
parser.add_argument('-I','--lInptClrs',help='List of input colors',default=None,type=str)

args = parser.parse_args()

######################
######################
######################
#  Input parameters  #
######################
######################
######################
smpls = args.smpls
exprssnFl = args.exprssnFl
species = args.species
outFldrPrnt = args.outFldrPrnt
lKmins = args.lKmins
smplsGrpFl = args.smplsGrpFl
ercc_file = args.ercc_file
nCores = args.nCores
minNonfailed = args.minNonfailed
minCountThreshold = args.minCountThreshold
maxModelPlots = args.maxModelPlots
maxAdjVar = args.maxAdjVar
rnCllCycl_crrctn = args.rnCllCycl_crrctn
cSeed = args.cSeed
lInptClrs = args.lInptClrs


"""
Set up initial parameters
"""

#number of clusters
kMins=[2]
kMins.extend([int(v) for v in lKmins.split('|') if v.strip()])
smplsNms = sorted([s for s in smpls.split('|') if s.strip()])

#list of colors
if lInptClrs is not None:
	lInptClrs = lInptClrs.split('|')

#Test for samples consistency and define species names
assert len(set(smplsNms))==len(smplsNms)
assert os.path.exists(exprssnFl)
assert os.path.exists(outFldrPrnt) 
assert os.path.exists(ercc_file)
assert species in {'Human','Mouse'}
if species=='Mouse':
	sppPrfx='mm10'
else:
	sppPrfx='hg38'


#Define a dictionary of names and samples to make figures and applications
if smplsGrpFl is not None:
	dGrpSmpl = {}
	for l in open(smplsGrpFl,'r'):
		if l.strip():
			smpl,grp = l.splitlines()[0].split('\t')
			if dGrpSmpl.has_key(grp):
				dGrpSmpl[grp].append(smpl)
			else:
				dGrpSmpl[grp] = [smpl]
				


#########################################################
################### Execute analysis ####################
#########################################################

#Run PAGODA and cluster cells in kMins clusters
for k_minSlhtts in kMins:
	if k_minSlhtts==2:
		outFldr = os.path.join(outFldrPrnt,'KMin%s.d'%k_minSlhtts)
	else:
		outFldr = os.path.join(outFldrPrnt,'KEql%s.d'%k_minSlhtts)
		#For cell cycle corrected
		if rnCllCycl_crrctn:
			outFldr_cc_prnt = os.path.join(outFldrPrnt,'cllCyclCrrtd.d')
			outFldr_cc = os.path.join(outFldr_cc_prnt,'KEql%s.d'% \
			k_minSlhtts)
			#
			for oFldr in [outFldr_cc_prnt,outFldr_cc]:
				if not os.path.exists(oFldr):
					os.mkdir(oFldr)
	#Create folder
	if not os.path.exists(outFldr):
		os.mkdir(outFldr)
	#Set standard output files
	errMdlPdf = os.path.join(outFldr,'errMdl.pdf')#Error model pdf
	errMdlRds = os.path.join(outFldr,'errMdl.rds')#Error model rsd file
	varNormPdf = os.path.join(outFldr,'varNorm.pdf')#Normalized variance plots
	varNormDepthRds = os.path.join(outFldr,'varNormDepth.rds')#Normalized variance and depth rds file
	go_envRds = os.path.join(outFldr,'go_envRds.rds')#GO gene set environment rds file
	pwpcaRds  = os.path.join(outFldr,'pwpca.rds')#Weighted first PC magnitudes for each GO gene set
	ovrDsprsnPdf = os.path.join(outFldr,'ovrDsprsn.pdf')#Plot of observed overdispersion for each GO gene set
	ovrDsprsnCsv = os.path.join(outFldr,'ovrDsprsn.rds')#Observed overdispersion for each GO gene set rds file
	ovrDsprsnDeNovoPdf = os.path.join(outFldr,'ovrDsprsnDeNovo.pdf')#Plot of observed overdispersion for each de novo gene set
	ovrDsprsnDeNovoRds = os.path.join(outFldr,'ovrDsprsnDeNovo.rds')#Observed overdispersion for each de novo gene set rds file
	clpcaRds = os.path.join(outFldr,'clpca.rds')#Background variance magnitude file for de novo gene set rds file
	clpcaPdf = os.path.join(outFldr,'clpca.pdf')#Plot of background variance magnitude file for de novo gene sets
	vslzRsltsPdf = os.path.join(outFldr,'vslzRslts.pdf')#Plot file to visualize results
	hcRds = os.path.join(outFldr,'hc.rds')#Hierarchical clustering rds file
	hcPdf = os.path.join(outFldr,'hc.pdf'))#Plot of the hierarchical clustering
	silPdf = os.path.join(outFldr,'silhouettes.pdf')#Average silhouettes plot for different K from 1 to k_minSlhtts clusters
	silWidthPdf = os.path.join(outFldr,'silhouettes_width.pdf')#Silhouettes width for the selected k_minSlhtts
	tamrRds = os.path.join(outFldr,'tamr.rds')#Reduce redundant aspects first step rds file
	tamr2Rds = os.path.join(outFldr,'tamr2.rds')#Reduce redundant aspects second step rds file
	wcordRds = os.path.join(outFldr,'wcord.rds')#Distance matrix between cells rds file
	cellClusteringRds = os.path.join(outFldr,'cellClustering.rds')#Overall cell clustering rds file
	tSNEPagodaRDS = os.path.join(outFldr,'tSNEPagoda.rds')#tSNE rds file
	clstrAllSmplsPdf = os.path.join(outFldr,'clstrAllSmpls.pdf')#Plots to visualize the hierarchical clustering and tSNE
	allSmplsPdfClstrs = os.path.join(outFldr,'allSmpls_sample.pdf')#Plots to visualize the tSNE
	allSmplPdf = os.path.join(outFldr,'allSmpls_hc.pdf')#Plot to visualize the samples in the tSNE embedding
	appRds = os.path.join(outFldr,'all_app.rds')#PAGODA application file with all the results
	appRdsAllSmpl = os.path.join(outFldr,'allSmpls_app.rds')#PAGODA application file with all the results (copy)
	#Create symbolic links from PAGODA processed files for K==2
	if k_minSlhtts!=2:
		srcFldr = os.path.join(outFldrPrnt,'KMin2.d')
		for nmFl in ['clpca.rds','clpca.pdf','errMdl.pdf','errMdl.rds', \
		'go_envRds.rds','hc.rds','ovrDsprsnDeNovo.pdf','ovrDsprsnDeNovo.rds', \
		'ovrDsprsn.pdf','ovrDsprsn.rds','pwpca.rds','tamr2.rds', \
		'tamr.rds','varNormDepth.rds','varNorm.pdf','varNorm.rds', \
		'wcord.rds']:
			if not os.path.exists(os.path.join(outFldr,nmFl)):
				os.symlink(os.path.join(srcFldr,nmFl),os.path.join(outFldr, \
				nmFl))
	#Make lists of samples for input
	l_lNames,l_clstrsPdf,l_clstrSmplPdf,l_appRds = [smplsNms], \
	[os.path.join(outFldr,'all_hc.pdf')],[os.path.join(outFldr, \
	'all_sample.pdf')],[os.path.join(outFldr,'all_app.rds')]
	for smpl in smplsNms:
		l_lNames.append(smpl)
		l_clstrsPdf.append(os.path.join(outFldr,'%s_hc.pdf'%smpl))
		l_clstrSmplPdf.append(os.path.join(outFldr,'%s_sample.pdf'%smpl))
		l_appRds.append(os.path.join(outFldr,'%s_app.rds'%smpl))
	if smplsGrpFl is not None:
		for name in dGrpSmpl.keys():
			smpls = dGrpSmpl[name]
			l_lNames.append(smpls)
			l_clstrsPdf.append(os.path.join(outFldr,'%s_hc.pdf'% \
			name))
			l_clstrSmplPdf.append(os.path.join(outFldr,'%s_sample.pdf'% \
			name))
			l_appRds.append(os.path.join(outFldr,'%s_app.rds'%name))
	#Run PAGODA
	if k_minSlhtts==2:#For K==2
		rnPGDinKMin2=False#switch to runs PAGODA in KMin2
		for nmFl in ['clpca.rds','clpca.pdf','errMdl.pdf','errMdl.rds', \
		'go_envRds.rds','hc.rds','ovrDsprsnDeNovo.pdf','ovrDsprsnDeNovo.rds', \
		'ovrDsprsn.pdf','ovrDsprsn.rds','pwpca.rds','tamr2.rds', \
		'tamr.rds','varNormDepth.rds','varNorm.pdf','varNorm.rds', \
		'wcord.rds']:
			if not os.path.exists(os.path.join(outFldrPrnt,'KMin2.d',nmFl)):
				rnPGDinKMin2=True
		#
		if rnPGDinKMin2:
			runPAGODA(exprssnFl,ercc_file,nCores,minCountThreshold, \
			minNonfailed,maxModelPlots,errMdlRds,errMdlPdf,maxAdjVar,varNormPdf, \
			varNormRds,varNormDepthRds,go_envRds,pwpcaRds,ovrDsprsnPdf, \
			ovrDsprsnCsv,ovrDsprsnDeNovoPdf,clpcaRds,ovrDsprsnDeNovoRds, \
			vslzRsltsPdf,hcRds,tamrRds,tamr2Rds,wcordRds,hcPdf,silPdf, \
			silWidthPdf,cellClusteringRds,tSNEPagodaRDS,clstrAllSmplsPdf, \
			allSmplsPdfClstrs,allSmplPdf,appRdsAllSmpl,l_lNames,l_clstrsPdf, \
			l_clstrSmplPdf,l_appRds,slhttsHght=False,clpcaPdf=clpcaPdf, \
			nCoresNoise=1,species=species.lower(),k_minSlhtts=k_minSlhtts, \
			cSeed=cSeed,lInptClrs=lInptClrs)
		#run with cell cycle correction. This will correct with genes in 
		#'GO:0007049' and 'GO:0051301' gene ontologies
		if rnCllCycl_crrctn:
			rnPGDinKMin2_cc=False#switch to runs PAGODA in KMin2
			for nmFl in ['clpca.cc.rds','clpca.cc.pdf','tSNEPagoda.cc.rds', \
			'hc.cc.rds','ovrDsprsnDeNovo.cc.pdf','ovrDsprsnDeNovo.cc.rds', \
			'ovrDsprsn.cc.pdf','ovrDsprsn.cc.rds','pwpca.cc.rds','tamr2.cc.rds', \
			'tamr.cc.rds','varNormDepth.cc.rds','varNorm.cc.pdf','varNorm.rds', \
			'wcord.cc.rds','cellClustering.cc.rds','silhouettes.cc.pdf']:
				if not os.path.exists(os.path.join(outFldrPrnt,'KMin2.d',nmFl)):
					rnPGDinKMin2_cc=True
			if rnPGDinKMin2_cc:
				runPAGODA(exprssnFl,ercc_file,nCores,minCountThreshold, \
				minNonfailed,maxModelPlots,errMdlRds,errMdlPdf,maxAdjVar, \
				varNormPdf,varNormRds,varNormDepthRds,go_envRds,pwpcaRds, \
				ovrDsprsnPdf,ovrDsprsnCsv,ovrDsprsnDeNovoPdf,clpcaRds, \
				ovrDsprsnDeNovoRds,vslzRsltsPdf,hcRds,tamrRds,tamr2Rds, \
				wcordRds,hcPdf,silPdf,silWidthPdf,cellClusteringRds, \
				tSNEPagodaRDS,clstrAllSmplsPdf,allSmplsPdfClstrs, \
				allSmplPdf,appRdsAllSmpl,l_lNames,l_clstrsPdf, \
				l_clstrSmplPdf,l_appRds,slhttsHght=False,clpcaPdf=clpcaPdf, \
				nCoresNoise=1,species=species.lower(),k_minSlhtts=k_minSlhtts, \
				crrct_by_fctr=True,fctrsToCrrct=['GO:0007049','GO:0051301'], \
				cSeed=cSeed,lInptClrs=lInptClrs)
	else:#For K>2
		runPAGODA(exprssnFl,ercc_file,nCores,minCountThreshold, \
		minNonfailed,maxModelPlots,errMdlRds,errMdlPdf,maxAdjVar,varNormPdf, \
		varNormRds,varNormDepthRds,go_envRds,pwpcaRds,ovrDsprsnPdf, \
		ovrDsprsnCsv,ovrDsprsnDeNovoPdf,clpcaRds,ovrDsprsnDeNovoRds, \
		vslzRsltsPdf,hcRds,tamrRds,tamr2Rds,wcordRds,hcPdf,silPdf, \
		silWidthPdf,cellClusteringRds,tSNEPagodaRDS,clstrAllSmplsPdf, \
		allSmplsPdfClstrs,allSmplPdf,appRdsAllSmpl,l_lNames,l_clstrsPdf, \
		l_clstrSmplPdf,l_appRds,slhttsHght=False,clpcaPdf=clpcaPdf, \
		nCoresNoise=1,species=species.lower(),k_minSlhtts=k_minSlhtts, \
		k_maxSlhtts=k_minSlhtts,cSeed=cSeed,lInptClrs=lInptClrs)
		#run with cell cycle correction
		if rnCllCycl_crrctn:
			for trgtNmFl,refNmFl in [('wcord.rds','wcord.cc.rds'), \
			('vslzRslts.pdf','vslzRslts.cc.pdf'),('varNorm.rds', \
			'varNormDepth.cc.rds'),('varNormDepth.rds', \
			'varNormDepth.cc.rds'),('tamr.rds','tamr.cc.rds'), \
			('tamr2.rds','tamr2.cc.rds'),('pwpca.rds','pwpca.cc.rds'), \
			('ovrDsprsn.rds','ovrDsprsn.cc.rds'),('ovrDsprsn.pdf', \
			'ovrDsprsn.cc.pdf'),('ovrDsprsnDeNovo.rds', \
			'ovrDsprsnDeNovo.cc.rds'),('ovrDsprsnDeNovo.pdf', \
			'ovrDsprsnDeNovo.cc.pdf'),('hc.rds','hc.cc.rds'),('go_envRds.rds', \
			'go_envRds.rds'),('errMdl.rds','errMdl.rds'),('errMdl.pdf', \
			'errMdl.pdf'),('clpca.rds','clpca.cc.rds'),('clpca.pdf', \
			'clpca.cc.pdf')]:
				trgtFl = os.path.join(outFldr_cc,trgtNmFl)
				if not os.path.exists(trgtFl):
					refFl = os.path.join(outFldrPrnt,'KMin2.d',refNmFl)
					assert os.path.exists(refFl)
					os.symlink(refFl,trgtFl)
			#Make lists of samples for input
			l_lNames,l_clstrsPdf,l_clstrSmplPdf,l_appRds = [smplsNms], \
			[os.path.join(outFldr_cc,'all_hc.pdf')],[os.path.join(outFldr_cc, \
			'all_sample.pdf')],[os.path.join(outFldr_cc,'all_app.rds')]
			for smpl in smplsNms:
				l_lNames.append(smpl)
				l_clstrsPdf.append(os.path.join(outFldr_cc, \
				'%s_hc.pdf'%smpl))
				l_clstrSmplPdf.append(os.path.join(outFldr_cc, \
				'%s_sample.pdf'%smpl))
				l_appRds.append(os.path.join(outFldr_cc,'%s_app.rds'%smpl))
			if smplsGrpFl is not None:
				for name in dGrpSmpl.keys():
					smpls = dGrpSmpl[name]
					l_lNames.append(smpls)
					l_clstrsPdf.append(os.path.join(outFldr_cc, \
					'%s_hc.pdf'%name))
					l_clstrSmplPdf.append(os.path.join(outFldr_cc,'%s_sample.pdf'% \
					name))
					l_appRds.append(os.path.join(outFldr_cc,'%s_app.rds'% \
					name))
			#####################
			#Redefine file names for cell-cycle corrected results
			#####################
			errMdlPdf = os.path.join(outFldr_cc,'errMdl.pdf')#Error model pdf
			errMdlRds = os.path.join(outFldr_cc,'errMdl.rds')#Error model rsd file
			varNormPdf = os.path.join(outFldr_cc,'varNorm.pdf')#Normalized variance plots
			varNormRds = os.path.join(outFldr_cc,'varNorm.rds')#Normalized variance rds file
			varNormDepthRds = os.path.join(outFldr_cc,'varNormDepth.rds')#Normalized variance and depth rds file
			ovrDsprsnPdf = os.path.join(outFldr_cc,'ovrDsprsn.pdf')#Plot of observed overdispersion for each GO gene set
			ovrDsprsnCsv = os.path.join(outFldr_cc,'ovrDsprsn.rds')#Observed overdispersion for each GO gene set rds file
			pwpcaRds  = os.path.join(outFldr_cc,'pwpca.rds')#Weighted first PC magnitudes for each GO gene set
			go_envRds = os.path.join(outFldr_cc,'go_envRds.rds')#GO gene set environment rds file
			ovrDsprsnDeNovoPdf = os.path.join(outFldr_cc,'ovrDsprsnDeNovo.pdf')#Plot of observed overdispersion for each de novo gene set
			ovrDsprsnDeNovoRds = os.path.join(outFldr_cc,'ovrDsprsnDeNovo.rds')#Observed overdispersion for each de novo gene set rds file
			clpcaRds = os.path.join(outFldr_cc,'clpca.rds')#Background variance magnitude file for de novo gene set rds file
			clpcaPdf = os.path.join(outFldr_cc,'clpca.pdf')#Plot of background variance magnitude file for de novo gene sets
			vslzRsltsPdf = os.path.join(outFldr_cc,'vslzRslts.pdf')#Plot file to visualize results
			hcRds = os.path.join(outFldr_cc,'hc.rds')#Hierarchical clustering rds file
			hcPdf = os.path.join(outFldr_cc,'hc.pdf')#Plot of the hierarchical clustering
			silPdf = os.path.join(outFldr_cc,'silhouettes.pdf')#Average silhouettes plot for different K from 1 to k_minSlhtts clusters
			silWidthPdf = os.path.join(outFldr_cc,'silhouettes_width.pdf')#Silhouettes width for the selected k_minSlhtts
			tamrRds = os.path.join(outFldr_cc,'tamr.rds')#Reduce redundant aspects first step rds file
			tamr2Rds = os.path.join(outFldr_cc,'tamr2.rds')#Reduce redundant aspects second step rds file
			wcordRds = os.path.join(outFldr_cc,'wcord.rds')#Distance matrix between cells rds file
			cellClusteringRds = os.path.join(outFldr_cc,'cellClustering.rds')#Overall cell clustering rds file
			tSNEPagodaRDS = os.path.join(outFldr_cc,'tSNEPagoda.rds')#tSNE rds file
			clstrAllSmplsPdf = os.path.join(outFldr_cc,'clstrAllSmpls.pdf')#Plots to visualize the hierarchical clustering and tSNE
			allSmplsPdfClstrs = os.path.join(outFldr_cc,'allSmpls_sample.pdf')#Plots to visualize the tSNE
			allSmplPdf = os.path.join(outFldr_cc,'allSmpls_hc.pdf')#Plot to visualize the samples in the tSNE embedding
			appRds = os.path.join(outFldr_cc,'all_app.rds')#PAGODA application file with all the results
			appRdsAllSmpl = os.path.join(outFldr_cc,'allSmpls_app.rds')#PAGODA application file with all the results (copy)
			#Run PAGODA for cell cycle corrected
			runPAGODA(exprssnFl,ercc_file,nCores,minCountThreshold, \
			minNonfailed,maxModelPlots,errMdlRds,errMdlPdf,maxAdjVar, \
			varNormPdf,varNormRds,varNormDepthRds,go_envRds,pwpcaRds, \
			ovrDsprsnPdf,ovrDsprsnCsv,ovrDsprsnDeNovoPdf,clpcaRds, \
			ovrDsprsnDeNovoRds,vslzRsltsPdf,hcRds,tamrRds,tamr2Rds, \
			wcordRds,hcPdf,silPdf,silWidthPdf,cellClusteringRds, \
			tSNEPagodaRDS,clstrAllSmplsPdf,allSmplsPdfClstrs,allSmplPdf, \
			appRdsAllSmpl,l_lNames,l_clstrsPdf,l_clstrSmplPdf,l_appRds, \
			slhttsHght=False,clpcaPdf=clpcaPdf,nCoresNoise=1, \
			species=species.lower(),k_minSlhtts=k_minSlhtts, \
			k_maxSlhtts=k_minSlhtts,cSeed=cSeed,lInptClrs=lInptClrs)
	#
print 'Done!'


#report the numbers of cells in each cluster in each tSNE
lTsneFls = []
for k_minSlhtts in kMins:
	if rnCllCycl_crrctn:
		lFldrsIntrst = [outFldrPrnt,outFldr_cc_prnt]
	else:
		lFldrsIntrst = [outFldrPrnt]
	for outFldrDad in lFldrsIntrst:
		outFldr = os.path.join(outFldrDad,'KEql%s.d'%k_minSlhtts)
		tSNEPagodaRDS = os.path.join(outFldr,'tSNEPagoda.rds')
		if os.path.exists(tSNEPagodaRDS):
			lTsneFls.append(tSNEPagodaRDS)
#report the numbers of cells in each cluster, execute
for f in lTsneFls:
	try:
		nClls = cntCllsTsneRDS(f)
		print 'For file %s, tSNE has %s cells'%(f,nClls)
	except:
		pass
