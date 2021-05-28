#!/usr/bin/env R
# -*- coding: utf-8 -*-
#
#  compile_h5ad.py
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


import matplotlib
import sys
import argparse,os

from numpy import array,inf,isinf,where

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
 
#input parameters  
parser.add_argument('-r','--rawCnts_PAGODAnrmlzd', default=None)
parser.add_argument('-f','--flippingTsne', default=False, type=str2bool)
#Files to make dictionaries
parser.add_argument('-E','--lExcldSmpls',help='List of samples to exclude: excldSmpl1|excldSmpl2|...',default=None)
parser.add_argument('-L','--smplsGrpFl',help='File with sample, group and outcome in smplsGrpFl',default=None)
parser.add_argument('-s','--smplsClmn',help='Column in each smplsGrpFl with the sample name',default=0,type=int)
parser.add_argument('-g','--grpsClmn',help='Column in each smplsGrpFl with the group name',default=None,type=int)
parser.add_argument('-u','--outcmClmn',help='Column in each smplsGrpFl with the outcome name',default=None,type=int)
parser.add_argument('-p','--flTypGnsToPlt', default=None)
#Methods to run
parser.add_argument('-N','--rtrnNrmlzExprssn', default=False, type=str2bool)
parser.add_argument('-D','--rtrnDstncs', default=False, type=str2bool)
parser.add_argument('-T','--rtrnTsnes', default=False, type=str2bool)
parser.add_argument('-M','--mkAdataObjcts', default=False, type=str2bool)
#Additional
parser.add_argument('-S','--stdyCsSfx',help='Study case sufix',default=None)
parser.add_argument('-A','--clstrAnnttnsFldr',help='Input folder with all files of clusters annotations (e.g. stdyCsSfx.KEqlX.tsv, with #Cluster_number\tCluster_annotation\n...)',default=None)
parser.add_argument('-n','--clstrNmbrClmn',help='Cluster number column in cluster file',default=0,type=int)
parser.add_argument('-a','--clstrAnntnClmn',help='Cluster annotation column in cluster file',default=1,type=int)
parser.add_argument('-F','--PAGODAFldrPrnt',help='Input folder with all PAGODA results',default=None)
parser.add_argument('-k','--lKmins',help='Input number of K clusters',default=None)
parser.add_argument('-O','--outFldrPrnt',help='output folder',default=None)
parser.add_argument('-H','--flippingHrzntlTsne',help='Flipping tSNE horizontally',default=False, type=str2bool)

###################################
##########    Switches  ###########
###################################


args = parser.parse_args()

if args:
	#input paramenters
	rawCnts_PAGODAnrmlzd = args.rawCnts_PAGODAnrmlzd
	flippingTsne = args.flippingTsne
	#Files to make dictionaries
	lExcldSmpls = args.lExcldSmpls
	smplsClmn = args.smplsClmn
	grpsClmn = args.grpsClmn
	outcmClmn = args.outcmClmn
	flTypGnsToPlt = args.flTypGnsToPlt
	#Methods to run
	rtrnNrmlzExprssn = args.rtrnNrmlzExprssn# If True will return the normalized data expression in PAGODA
	rtrnDstncs = args.rtrnDstncs# If True will return the distances between cells from wcords in PAGODA
	rtrnTsnes = args.rtrnTsnes# If True will return the tSNE for PAGODA
	mkAdataObjcts = args.mkAdataObjcts#If true will make adata data structures for each PAGODA result file
	#Additional
	smplsGrpFl = args.smplsGrpFl
	stdyCsSfx = args.stdyCsSfx
	clstrAnnttnsFldr = args.clstrAnnttnsFldr
	clstrNmbrClmn = args.clstrNmbrClmn
	clstrAnntnClmn = args.clstrAnntnClmn
	PAGODAFldrPrnt = args.PAGODAFldrPrnt
	lKmins = args.lKmins
	outFldrPrnt = args.outFldrPrnt
	flippingHrzntlTsne = args.flippingHrzntlTsne


#Print variables
print('Log info:')
print('rawCnts_PAGODAnrmlzd -> ',rawCnts_PAGODAnrmlzd)
print('flippingTsne -> ',flippingTsne)	
print('flippingHrzntlTsne -> ',flippingHrzntlTsne)	
print('lExcldSmpls -> ',lExcldSmpls)	
print('smplsGrpFl -> ',smplsGrpFl)
print('smplsClmn -> ',smplsClmn)
print('grpsClmn -> ',grpsClmn)
print('outcmClmn -> ',outcmClmn)	
print('flTypGnsToPlt -> ',flTypGnsToPlt)	
print('rtrnNrmlzExprssn -> ',rtrnNrmlzExprssn)	
print('rtrnDstncs -> ',rtrnDstncs)	
print('rtrnTsnes -> ',rtrnTsnes)	
print('mkAdataObjcts -> ',mkAdataObjcts)
print('stdyCsSfx -> ',stdyCsSfx)
print('clstrAnnttnsFldr -> ',clstrAnnttnsFldr)
print('clstrNmbrClmn -> ',clstrNmbrClmn)
print('clstrAnntnClmn -> ',clstrAnntnClmn)
print('PAGODAFldrPrnt -> ',PAGODAFldrPrnt)	
print('lKmins -> ',lKmins)	
print('outFldrPrnt -> ',outFldrPrnt)	

#Set output folder and files
l_KNum = [int(v) for v in lKmins.split('|') if v.strip()]
assert len(l_KNum)==1
kNum = l_KNum[0]
assert outFldrPrnt is not None 
outFldrSmpl = os.path.join(outFldrPrnt,'%s.d'%stdyCsSfx)
if kNum==2:
	outFldr = os.path.join(outFldrSmpl,'KMin%s.d'%kNum)
else:
	outFldr = os.path.join(outFldrSmpl,'KEql%s.d'%kNum)

if not os.path.exists(outFldrPrnt):
	os.mkdir(outFldrPrnt)
if not os.path.exists(outFldrSmpl):
	os.mkdir(outFldrSmpl)
if not os.path.exists(outFldr):
	os.mkdir(outFldr)

#Define output file names
if kNum==2:
	dstncMtrxRdsInput = os.path.join(PAGODAFldrPrnt,'KMin%s.d'%kNum, \
	'cellClustering.rds')
	rdsPAGODAAppInput = os.path.join(PAGODAFldrPrnt,'KMin%s.d'%kNum, \
	'all_app.rds')
	dstncMtrxRds = os.path.join(outFldr,'%s.KMin%s.pagodaDm.rds'% \
	(stdyCsSfx,kNum))
	rdsPAGODAApp = os.path.join(outFldr,'%s.KMin%s.pagodaApp.rds'% \
	(stdyCsSfx,kNum))
	clstrs_PAGODAnrmlzd = os.path.join(outFldr,'%s.KMin%s.clstrsLbld.tsv'% \
	(stdyCsSfx,kNum))
	outCnts_PAGODAnrmlzd = os.path.join(outFldr,'%s.KMin%s.pagodaApp.nmrlzd.tsv'% \
	(stdyCsSfx,kNum))
	outDm_PAGODAnrmlzd = os.path.join(outFldr,'%s.KMin%s.pagodaDm.tsv'% \
	(stdyCsSfx,kNum))
	pagoda_tSNEFl = os.path.join(outFldr,'%s.KMin%s.tsne.tsv'% \
	(stdyCsSfx,kNum))
	depot = os.path.join(outFldr,'%s.KMin%s.h5ad'%(stdyCsSfx,kNum))
	outRawCnts_PAGODAnrmlzd = os.path.join(outFldr,'%s.KMin%s.cnts.csv'% \
	(stdyCsSfx,kNum))
	flClstrNmOrd = os.path.join(outFldr,'%s.KMin%s.pagodaApp.clstrNmOrd.txt'% \
	(stdyCsSfx,kNum))
else:
	dstncMtrxRdsInput = os.path.join(PAGODAFldrPrnt,'KEql%s.d'%kNum, \
	'cellClustering.rds')
	rdsPAGODAAppInput = os.path.join(PAGODAFldrPrnt,'KEql%s.d'%kNum, \
	'all_app.rds')
	dstncMtrxRds = os.path.join(outFldr,'%s.KEql%s.pagodaDm.rds'% \
	(stdyCsSfx,kNum))
	rdsPAGODAApp = os.path.join(outFldr,'%s.KEql%s.pagodaApp.rds'% \
	(stdyCsSfx,kNum))
	clstrs_PAGODAnrmlzd = os.path.join(outFldr,'%s.KEql%s.clstrsLbld.tsv'% \
	(stdyCsSfx,kNum))
	outCnts_PAGODAnrmlzd = os.path.join(outFldr,'%s.KEql%s.pagodaApp.nmrlzd.tsv'% \
	(stdyCsSfx,kNum))
	outDm_PAGODAnrmlzd = os.path.join(outFldr,'%s.KEql%s.pagodaDm.tsv'% \
	(stdyCsSfx,kNum))
	pagoda_tSNEFl = os.path.join(outFldr,'%s.KEql%s.tsne.tsv'% \
	(stdyCsSfx,kNum))
	depot = os.path.join(outFldr,'%s.KEql%s.h5ad'%(stdyCsSfx,kNum))
	outRawCnts_PAGODAnrmlzd = os.path.join(outFldr,'%s.KEql%s.cnts.csv'% \
	(stdyCsSfx,kNum))
	flClstrNmOrd = os.path.join(outFldr,'%s.KEql%s.pagodaApp.clstrNmOrd.txt'% \
	(stdyCsSfx,kNum))
		


#Assert is running on python with the correct version
if mkAdataObjcts:
	assert sys.version_info.major>2#Assert is running on python 3
	from singlecell.scanpy_mod import rtrnClstrsPAGODARdsUnq, \
	rtrnClstrsPAGODARdsRuslanUnq,rtrnClstrsPAGODARdsRuslanUnqV2,addRawDt, \
	rtrnExprssnNms,rnBscPrcss
	import scanpy as sc
elif rtrnNrmlzExprssn or rtrnDstncs or rtrnTsnes:
	assert sys.version_info.major==2#Assert is running on python 2
	from singlecell.formats import mkTabFromArray,rtrnExprssnNms,getTsne, \
	readRDdata,wrprRtrnNrmlzdPAGODA,rtrnClstrsPAGODARds

assert os.path.exists(rawCnts_PAGODAnrmlzd)

##################################
##################################
##################################
#####    Get input files    ######
##################################
##################################
##################################

#Build links to files
if not os.path.exists(rdsPAGODAApp):
	assert os.path.exists(rdsPAGODAAppInput)
	os.symlink(rdsPAGODAAppInput,rdsPAGODAApp)

if not os.path.exists(dstncMtrxRds):
	if os.path.exists(dstncMtrxRdsInput):
		os.symlink(dstncMtrxRdsInput,dstncMtrxRds)
	else:
		dstncMtrxRds = None
	
if not os.path.exists(outRawCnts_PAGODAnrmlzd):
	assert os.path.exists(rawCnts_PAGODAnrmlzd)
	os.symlink(rawCnts_PAGODAnrmlzd,outRawCnts_PAGODAnrmlzd)


#Build lists of output folders, samples and number of clusters
if not os.path.exists(clstrs_PAGODAnrmlzd):
	assert clstrAnnttnsFldr is not None and os.path.exists(clstrAnnttnsFldr)
	assert stdyCsSfx is not None
	assert lKmins is not None
	#
	dClstrNms = {}
	if kNum==2:
		clstrAnnttnsFl = os.path.join(clstrAnnttnsFldr,'%s.KMin%s.tsv'% \
		(stdyCsSfx,kNum))
	else:
		clstrAnnttnsFl = os.path.join(clstrAnnttnsFldr,'%s.KEql%s.tsv'% \
		(stdyCsSfx,kNum))
	#
	assert os.path.exists(clstrAnnttnsFl)
	sClstrAnntn = set()
	for l in open(clstrAnnttnsFl,'r'):
		if l.strip() and l[0]!='#':
			l = l.splitlines()[0].split('\t')
			clstrNmbr = l[clstrNmbrClmn]
			clstrAnntn = l[clstrAnntnClmn]
			assert not dClstrNms.has_key(clstrNmbr)
			dClstrNms[clstrNmbr] = clstrAnntn
			assert clstrAnntn not in sClstrAnntn
			sClstrAnntn.add(clstrAnntn)
	#write cluster names
	if not os.path.exists(flClstrNmOrd):
		opndFlClstrNmOrd = open(flClstrNmOrd,'w')
		opndFlClstrNmOrd.write('\n'.join(sorted(sClstrAnntn)))
		opndFlClstrNmOrd.close()
		print('Please re-order cluster names in file %s after first run...'% \
		flClstrNmOrd)
	#Select genes and cells from PAGODA server and retrieve the counts
	nrmlzdCountsMmus,smplNms,cllNmsArray,ar_gnNms = \
	wrprRtrnNrmlzdPAGODA(rdsPAGODAApp)
	ar_clstrs = rtrnClstrsPAGODARds(rdsPAGODAApp)
	#clusters
	ar_nmdClstrs = array([dClstrNms[str(c)] for c in array(ar_clstrs)])
	mkTabFromArray(cllNmsArray,array(['cluster']),array([ar_nmdClstrs]).T, \
	clstrs_PAGODAnrmlzd,frstClmnNm='cell')

		
##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

dCll_clstrs_PAGODAnrmlzd = dict([(l.split()[0],l.split()[1]) for l in \
open(clstrs_PAGODAnrmlzd).read().splitlines()[1:] if l.strip()])

#Stages, expected outcome, cell cycle genes
lPAGODAAppClstrNmOrd = [l.splitlines()[0] for l in open(flClstrNmOrd,'r') \
if l.strip()]

#Make dictionary of genes to plot
if flTypGnsToPlt is not None:
	dTypGnsToPlt = {}
	for l in open(flTypGnsToPlt,'r'):
		if l.strip():
			gnNm,cllTyp = l.splitlines()[0].split('\t')
			if cllTyp in dTypGnsToPlt:
				dTypGnsToPlt[cllTyp].append(gnNm)
			else:
				dTypGnsToPlt[cllTyp]=[gnNm]


#Test for samples consistency and define species names
if lExcldSmpls is not None:
	sSmplsTExcld = set([s for s in lExcldSmpls.split('|') if s.strip()])
else:
	sSmplsTExcld = set()

#Define a dictionary of names and samples to make figures and applications
dSmplStg,dSmplOutcm = {},{}
for l in open(smplsGrpFl,'r'):
	if l.strip() and l[0]!='#':
		l = l.splitlines()[0].split('\t')
		smpl = l[smplsClmn]
		if smpl not in sSmplsTExcld:
			if grpsClmn is not None:
				grp = l[grpsClmn]
				assert smpl not in set(dSmplStg.keys())
				dSmplStg[smpl] = grp
			else:
				dSmplStg[smpl] = 'Unk'
			#
			if outcmClmn is not None:
				outcm = l[outcmClmn]
				assert smpl not in set(dSmplOutcm.keys())
				dSmplOutcm[smpl] = outcm
			else:
				dSmplOutcm[smpl] = 'Unk'


###################################
####### Load SC references ########
###################################
#Run in python 2.7x
if rtrnNrmlzExprssn:
	#Load counts and reference
	cnts_nrmlzd_PAGODA,gnNms_nrmlzd_PAGODA = rtrnExprssnNms(rdsPAGODAApp)
	gnNms_nrmlzd_PAGODA = array(cnts_nrmlzd_PAGODA.rownames)
	cllNms_nrmlzd_PAGODA = array(cnts_nrmlzd_PAGODA.colnames)
	cnts_nrmlzd_PAGODA = array(cnts_nrmlzd_PAGODA)
	mkTabFromArray(gnNms_nrmlzd_PAGODA,cllNms_nrmlzd_PAGODA, \
	cnts_nrmlzd_PAGODA,outCnts_PAGODAnrmlzd,sep=',')

#Run in python 2.7x
if rtrnDstncs:
	#Save the PAGODA-normalized distances to tsv
	#Load distances
	if dstncMtrxRds is not None:
		dm_nrmlzd_PAGODA = readRDdata(dstncMtrxRds)
		cllNms_nrmlzd_PAGODA_v1 = array(dm_nrmlzd_PAGODA.rownames)
		cllNms_nrmlzd_PAGODA_v2 = array(dm_nrmlzd_PAGODA.colnames)
		dm_nrmlzd_PAGODA = array(dm_nrmlzd_PAGODA)
		mkTabFromArray(cllNms_nrmlzd_PAGODA_v1,cllNms_nrmlzd_PAGODA_v2, \
		dm_nrmlzd_PAGODA,outDm_PAGODAnrmlzd,sep=',')


#Run in python 2.7x
if rtrnTsnes:
	#Save the PAGODA-normalized counts to tsv
	#Load distances
	tSNEFl_nrmlzd_PAGODA = getTsne(rdsPAGODAApp)
	gnNms_nrmlzd_PAGODA = array(tSNEFl_nrmlzd_PAGODA.rownames)
	cllNms_nrmlzd_PAGODA = array(['tSNE1','tSNE2'])
	tSNEFl_nrmlzd_PAGODA = array(tSNEFl_nrmlzd_PAGODA)
	mkTabFromArray(gnNms_nrmlzd_PAGODA,cllNms_nrmlzd_PAGODA, \
	tSNEFl_nrmlzd_PAGODA,pagoda_tSNEFl,sep=',')




###################################
####### Make scanpy objects #######
###################################
#Run in Python 3.6X
if mkAdataObjcts:
	#Write adata
	dCll_clstrs_PAGODAnrmlzd = dict([(l.split()[0],l.split()[1]) \
	for l in open(clstrs_PAGODAnrmlzd).read().splitlines()[1:] if \
	l.strip()])
	#Correct for missing data
	sSmplNames = set(dCll_clstrs_PAGODAnrmlzd.keys())
	while sSmplNames:
		smpl = sSmplNames.pop()
		smpl = '.'.join(smpl.split('.')[:-1])
		if smpl not in dSmplStg:
			dSmplStg[smpl] = 'Unk'
		if smpl not in dSmplOutcm:
			dSmplOutcm[smpl] = 'Unk'
	#
	adata_cllCltr= addRawDt(outCnts_PAGODAnrmlzd,rawCnts_PAGODAnrmlzd, \
	dCll_clstrs_PAGODAnrmlzd,pagoda_tSNEFl,flippingTsne=flippingTsne, \
	dSmplStg=dSmplStg,dSmplOutcm=dSmplOutcm,flippingHrzntlTsne= \
	flippingHrzntlTsne)
	#Add colors
	if rdsPAGODAApp.find('Wnt1_tail_E105.')>-1:
		rtrnClstrsPAGODARdsRuslanUnq(rdsPAGODAApp, \
		adata_cllCltr,lPAGODAAppClstrNmOrd= \
		lPAGODAAppClstrNmOrd)
		vmin,vmax=0,0.8
	elif rdsPAGODAApp.find('Wnt1_tail_E105_nc.')>-1:
		rtrnClstrsPAGODARdsRuslanUnqV2(rdsPAGODAApp, \
		adata_cllCltr,lPAGODAAppClstrNmOrd=  \
		lPAGODAAppClstrNmOrd)
		vmin,vmax=0,0.8
	else:
		rtrnClstrsPAGODARdsUnq(rdsPAGODAApp,adata_cllCltr, \
		lPAGODAAppClstrNmOrd=lPAGODAAppClstrNmOrd)
		vmin,vmax=-7,7
	#
	if flTypGnsToPlt is not None:
		sRef=set(adata_cllCltr.var_names)
		dTypGnsToPlt = dict([(k,list(set(v).intersection(sRef))) \
		for k,v in dTypGnsToPlt.items()])
		print('Drawing non-raw data...')
		#
		if 'absent_N' in set(adata_cllCltr.obs['PAGODA_hc']):
			adata_cllCltrTmp = adata_cllCltr[adata_cllCltr.obs['PAGODA_hc']!= \
			'absent_N']
			###Re-do
			#get order of cluster names
			srtdCat = adata_cllCltrTmp.obs['PAGODA_hc'].cat.categories
			clrCtgSrtd = array([adata_cllCltrTmp.obs[adata_cllCltrTmp. \
			obs['PAGODA_hc']==cat]['cluster_color'][0] for cat in srtdCat])
			adata_cllCltrTmp.uns['PAGODA_hc_colors']=clrCtgSrtd
			###
			sc.pl.heatmap(adata_cllCltrTmp,dTypGnsToPlt,groupby='PAGODA_hc', \
			dendrogram=True,swap_axes=True,use_raw=False,save=True, \
			figsize=(10,10),**{'cmap':'seismic','vmin':vmin,'vmax':vmax})
			sc.pl.tsne(adata_cllCltrTmp,color=['PAGODA_hc'], \
			legend_loc='on data',use_raw=False,save=True)
		else:
			sc.pl.heatmap(adata_cllCltr,dTypGnsToPlt,groupby='PAGODA_hc', \
			dendrogram=True,swap_axes=True,use_raw=False,save=True, \
			figsize=(10,10),**{'cmap':'seismic','vmin':vmin,'vmax':vmax})
			sc.pl.tsne(adata_cllCltr,color=['PAGODA_hc'], \
			legend_loc='on data',use_raw=False,save=True)
	adata_cllCltr.write(depot, compression='gzip')




