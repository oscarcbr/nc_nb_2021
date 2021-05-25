#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  db_pop_marker_enrichment.py
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
Determine marker enrichments for independent cell populations
"""

from singlecell.between_samples import clcEnrchmntGrps
from singlecell.plots import pltCllTypEnrchmnt
from singlecell.formats import rtrnClrsOrdrFrmH5ad
from numpy import array,mean,float32
from numpy.random import shuffle
from itertools import combinations
from string import upper

import argparse,os
import sys


"""
Analysis conducted to comparisons for counts and expression between 
clusters and samples.
"""

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
 
#input parameters for preprocessing 
parser.add_argument('-s','--src',help='File	name',default=None)
parser.add_argument('-g','--clmn_gn',help='Gene name column',default=None)
parser.add_argument('-r','--clmn_grp',help='Annotation group column',default=None)
parser.add_argument('-f','--src_fl',help='Enrichment file location',default=None)
parser.add_argument('-t','--trgtSpp',help='Target species',default=None)
#Optional parameters
parser.add_argument('-d','--dataFldr',help='Input folder with all annotations',default='data')
parser.add_argument('-D','--rn_inAllDB',help='If True will run all databases',default=True,type=str2bool)


args = parser.parse_args()

if args:
	#input parameters
	src = args.src
	clmn_gn = args.clmn_gn
	clmn_grp = args.clmn_grp
	src_fl = args.src_fl
	trgtSpp = args.trgtSpp
	#optional
	rn_inAllDB = args.rn_inAllDB
	dataFldr = args.dataFldr


##############
#Switches for all data
##############
tstClstrEncrhmnt=True#if True will test cell type enrichment for each db
pltClstrEncrhmnt=True#if True will plot cell type enrichment for each db
sExcldTrgtDB=set()#avoid these as target databases (stil used them as references)



#Print log information
print('Log info:')
print('src -> ',src)
print('clmn_gn -> ',clmn_gn)
print('clmn_grp -> ',clmn_grp)
print('src_fl -> ',src_fl)
print('trgtSpp -> ',trgtSpp)
print('rn_inAllDB -> ',rn_inAllDB)
print('dataFldr -> ',dataFldr)
print('tstClstrEncrhmnt -> ',tstClstrEncrhmnt)
print('pltClstrEncrhmnt -> ',pltClstrEncrhmnt)
print('sExcldTrgtDB -> ',sExcldTrgtDB)


##############
##############
##############
#Input paths
##############
##############
fl1to1Ort = os.path.join(dataFldr,'hg38_mm10_orthlgs_1to1.txt')
dClmnTrgtClmnRef = {'Human':(0,1),'Mouse':(1,0)}
##############
#Enrichment DB paths
##############
#mm10
anntEnsmblFl_mm10 = os.path.join(dataFldr, \
'GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv')
ar_gnNms_mm10 = array([l.splitlines()[0].split('\t')[1] for l in \
open(anntEnsmblFl_mm10,'r') if l.strip() and l[0]!='#'])
CllAssFl_GOBP_mm10 = os.path.join(dataFldr,'GO_biological_process.mm10.tsv')
CllAssFl_GOMF_mm10 = os.path.join(dataFldr,'GO_molecular_function.mm10.tsv')
CllAssFl_GOCC_mm10 = os.path.join(dataFldr,'GO_cellular_component.mm10.tsv')
CllAssFl_Krtype_mm10 = os.path.join(dataFldr,'karyotype_band.mm10.tsv')
lAll_CllAssFl_mm10 = [CllAssFl_GOBP_mm10,CllAssFl_GOMF_mm10, \
CllAssFl_GOCC_mm10,CllAssFl_Krtype_mm10]
lAll_DbType_mm10 = ['GO_BP','GO_MF','GO_CC','Karyotype']
lAll_Dbspp_mm10 = ['Mouse','Mouse','Mouse','Mouse']
lAll_clmnGn_mm10 = [0,0,0,0]
lAll_clmnGrp_mm10 = [1,1,1,3]
lAll_outFldrs_mm10 = ['GO_BP','GO_MF','GO_CC','karyotype']
lAll_h5ad_mm10 = [None,None,None,None]
#hg38
anntEnsmblFl_hg38 = os.path.join(dataFldr, \
'GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv')
ar_gnNms_hg38 = array([l.splitlines()[0].split('\t')[1] for l in \
open(anntEnsmblFl_hg38,'r') if l.strip() and l[0]!='#'])
CllAssFl_GOBP_hg38 = os.path.join(dataFldr,'GO_biological_process.tsv')
CllAssFl_GOMF_hg38 = os.path.join(dataFldr,'GO_molecular_function.tsv')
CllAssFl_GOCC_hg38 = os.path.join(dataFldr,'GO_cellular_component.tsv')
CllAssFl_Krtype_hg38 = os.path.join(dataFldr,
'karyotype_band.hg38.tsv')
lAll_CllAssFl_hg38 = [CllAssFl_GOBP_hg38,CllAssFl_GOMF_hg38, \
CllAssFl_GOCC_hg38,CllAssFl_Krtype_hg38]
lAll_DbType_hg38 = ['GO_BP','GO_MF','GO_CC','Karyotype']
lAll_Dbspp_hg38 = ['Human','Human','Human','Human']
lAll_clmnGn_hg38 = [0,0,0,0]
lAll_clmnGrp_hg38 = [1,1,1,3]
lAll_outFldrs_hg38 = ['GO_BP','GO_MF','GO_CC','karyotype']
lAll_h5ad_hg38 = [None,None,None,None]
#Make dictionary of ar_gnNms and target species
dAr_gnNms = {'Human':ar_gnNms_hg38,'Mouse':ar_gnNms_mm10}

##############
#Obtain database for enrichment based on file
##############
indxDBFl = os.path.join(dataFldr,'index.enrchmnts.txt')
for l in open(indxDBFl,'r'):
	if l.strip() and l[0]!='#':
		l=l.splitlines()[0].split('\t')
		n=lAll_outFldrs_mm10.append(l[1]),lAll_outFldrs_hg38.append(l[1])
		n=lAll_clmnGn_mm10.append(int(l[2])), \
		lAll_clmnGn_hg38.append(int(l[2]))
		n=lAll_clmnGrp_mm10.append(int(l[3])), \
		lAll_clmnGrp_hg38.append(int(l[3]))
		n=lAll_CllAssFl_mm10.append(l[4]),lAll_CllAssFl_hg38.append(l[4])
		n=lAll_Dbspp_mm10.append(l[5]),lAll_Dbspp_hg38.append(l[5])
		n=lAll_DbType_mm10.append('ref'),lAll_DbType_hg38.append('ref')
		h5adDB = l[6].strip()
		if h5adDB=='NA':
			h5adDB = None
		n=lAll_h5ad_mm10.append(h5adDB),lAll_h5ad_hg38.append(h5adDB)


#list files and datatypes
dAll_CllAssctnTbls = {'Human':lAll_CllAssFl_hg38,'Mouse': \
lAll_CllAssFl_mm10}
dAll_DbType = {'Human':lAll_DbType_hg38,'Mouse':lAll_DbType_mm10}
dAll_Dbspp = {'Human':lAll_Dbspp_hg38,'Mouse':lAll_Dbspp_mm10}
dAll_clmnGn = {'Human':lAll_clmnGn_hg38,'Mouse':lAll_clmnGn_mm10}
dAll_clmnGrp = {'Human':lAll_clmnGrp_hg38,'Mouse':lAll_clmnGrp_mm10}
dAll_outFldrs = {'Human':lAll_outFldrs_hg38,'Mouse':lAll_outFldrs_mm10}
dAll_h5ad = {'Human':lAll_h5ad_hg38,'Mouse':lAll_h5ad_mm10}
##############
#Set other required input variables
##############
varNormDepthRds = ''#dummy
gnsEncrhdClstrFldr = ''#dummy
clmnrefDBSpInFl = 2


########################################################
########################################################
########################################################
########################################################
########################################################
################    Enrichment tests    ################
########################################################
########################################################
########################################################
########################################################
########################################################
"""
Test for enrichment of different cell types in each cluster for DBs
"""
if tstClstrEncrhmnt:
	#Run cells enrichment in each cluster
	print "Test for enrichment of different cell types in each cluster for DBs (upregulated)..."
	if rn_inAllDB:
		lindxDBFl = [l for l in open(indxDBFl,'r') if l.strip()]
		shuffle(lindxDBFl)
	else:
		lindxDBFl = ['\t'.join(['NA',src,str(clmn_gn),str(clmn_grp), \
		src_fl,trgtSpp,'NA'])]
	#
	while lindxDBFl:
		lInfo = lindxDBFl.pop()
		if lInfo.strip() and lInfo[0]!='#':
			#Parse require information
			cttn,src,clmn_gn,clmn_grp,src_fl,trgtSpp,h5ad_db = \
			lInfo.splitlines()[0].split('\t')
			if src in sExcldTrgtDB:
				continue
			clmn_gn,clmn_grp = int(clmn_gn),int(clmn_grp)
			#Obtain input info
			dClstrsGns_inpt = {}#{group lists:gene list}
			l_gnsEncrhdClstrFls,ar_gnNms_input = set(),dAr_gnNms[trgtSpp]
			for lInf in open(src_fl).read().splitlines():
				if lInf.strip() and lInf[0]!='#':
					lInf=lInf.split('\t')
					gnNm,grp = lInf[clmn_gn],lInf[clmn_grp]
					l_gnsEncrhdClstrFls.add(grp)
					if dClstrsGns_inpt.has_key(grp):
						dClstrsGns_inpt[grp].add(gnNm)
					else:
						dClstrsGns_inpt[grp]=set([gnNm])
			#Obtain complementary input info
			lAll_CllAssFl,lAll_DbType,lAll_Dbspp,lAll_clmnGn,lAll_clmnGrp, \
			lAll_outFldrs = dAll_CllAssctnTbls[trgtSpp],dAll_DbType[trgtSpp], \
			dAll_Dbspp[trgtSpp],dAll_clmnGn[trgtSpp],dAll_clmnGrp[trgtSpp], \
			dAll_outFldrs[trgtSpp]
			clmnTrgtIn1Ort1Fl,clmnRefIn1Ort1Fl = dClmnTrgtClmnRef[trgtSpp]
			print clmnRefIn1Ort1Fl,clmnTrgtIn1Ort1Fl,fl1to1Ort
			#Output folder
			rsltFldr = os.path.split(src_fl)[0]
			enrchmntAnlysFldr = os.path.join(rsltFldr,'cllTypeEnrchmnts.d')
			if not os.path.exists(enrchmntAnlysFldr):
				os.mkdir(enrchmntAnlysFldr)
			enrchmntAnlysFldr = os.path.join(rsltFldr, \
			'cllTypeEnrchmnts.d','%s.d'%src)
			if not os.path.exists(enrchmntAnlysFldr):
				os.mkdir(enrchmntAnlysFldr)
			#Calculate the enrichments for each annotation database		
			nCllAssctnTbls = len(lAll_CllAssFl)
			for q in xrange(nCllAssctnTbls):
				refSpp,cllAssctnTblFl,dbType,fldrNm,clmnGnInFl, \
				clmnCllTypInFl = lAll_Dbspp[q],lAll_CllAssFl[q], \
				lAll_DbType[q],lAll_outFldrs[q],lAll_clmnGn[q], \
				lAll_clmnGrp[q]
				#Test enrichment for each cluster independently
				if dbType is None:
					outFl = os.path.join(enrchmntAnlysFldr, \
					'SINCERA.tsv')
				else:
					outFl = os.path.join(enrchmntAnlysFldr,'%s.tsv'% \
					fldrNm)
				if not os.path.exists(outFl):
					clcEnrchmntGrps(varNormDepthRds,gnsEncrhdClstrFldr, \
					cllAssctnTblFl,outFl,l_gnsEncrhdClstrFls,dbType=dbType, \
					dClstrsGns_inpt=dClstrsGns_inpt,clmnGnInFl=clmnGnInFl, \
					clmnCllTypInFl=clmnCllTypInFl,clmnTrgtIn1Ort1Fl= \
					clmnTrgtIn1Ort1Fl,clmnRefIn1Ort1Fl=clmnRefIn1Ort1Fl, \
					ar_gnNms_input=ar_gnNms_input,trgtDBSp=trgtSpp, \
					refDBSp=refSpp,fl1to1Ort=fl1to1Ort,clmnrefDBSpInFl= \
					clmnrefDBSpInFl)
	print 'Done!'				
	



"""
Plot cell type enrichment for each cell cluster
"""
if pltClstrEncrhmnt:
	#Plot cells types enrichment in each cluster
	print "Plotting enrichment of different cell types in each cluster (upregulated)..."
	if rn_inAllDB:
		lindxDBFl = [l for l in open(indxDBFl,'r') if l.strip()]
		shuffle(lindxDBFl)
	else:
		lindxDBFl = ['\t'.join(['NA',src,str(clmn_gn),str(clmn_grp), \
		src_fl,trgtSpp,'NA'])]
	#		
	while lindxDBFl:
		lInfo = lindxDBFl.pop()
		if lInfo.strip() and lInfo[0]!='#':		
			cttn,src,clmn_gn,clmn_grp,src_fl,trgtSpp,h5ad_fl = \
			lInfo.splitlines()[0].split('\t')
			#Return colors in databases
			if h5ad_fl=='NA':
				ar_setClstrs,ar_setClstrsClrs=None,None
			else:
				ar_setClstrs,ar_setClstrsClrs=rtrnClrsOrdrFrmH5ad(h5ad_fl)
			#Obtain complementary input info
			lAll_CllAssFl,lAll_DbType,lAll_Dbspp,lAll_clmnGn, \
			lAll_clmnGrp,lAll_outFldrs = dAll_CllAssctnTbls[trgtSpp], \
			dAll_DbType[trgtSpp],dAll_Dbspp[trgtSpp], \
			dAll_clmnGn[trgtSpp],dAll_clmnGrp[trgtSpp], \
			dAll_outFldrs[trgtSpp]
			lAll_h5ad = dAll_h5ad[trgtSpp]
			#Output folder
			rsltFldr = os.path.split(src_fl)[0]
			enrchmntAnlysFldr = os.path.join(rsltFldr, \
			'cllTypeEnrchmnts.d','%s.d'%src)
			print rsltFldr
			#Plot the enrichments for each annotation database		
			nCllAssctnTbls = len(lAll_CllAssFl)
			for qq in xrange(nCllAssctnTbls):
				refSpp,cllAssctnTblFl,dbType,fldrNm,clmnGnInFl, \
				clmnCllTypInFl = lAll_Dbspp[qq],lAll_CllAssFl[qq], \
				lAll_DbType[qq],lAll_outFldrs[qq],lAll_clmnGn[qq], \
				lAll_clmnGrp[qq]
				#Test enrichment for each cluster independently
				for trshld in [0.05,0.01,0.001]:
					if dbType is None:
						encrhmntFl = os.path.join(enrchmntAnlysFldr, \
						'SINCERA.tsv')
					else:
						encrhmntFl = os.path.join(enrchmntAnlysFldr, \
						'%s.tsv'%fldrNm)
						if trshld==0.05:
							encrhmntPlt = os.path.join(enrchmntAnlysFldr, \
							'%s.svg'%fldrNm)
						elif trshld==0.01:
							encrhmntPlt = os.path.join(enrchmntAnlysFldr, \
							'%s_01.svg'%fldrNm)
						elif trshld==0.001:
							encrhmntPlt = os.path.join(enrchmntAnlysFldr, \
							'%s_001.svg'%fldrNm)
					if os.path.exists(encrhmntFl) and \
					((not os.path.exists(encrhmntPlt) or \
					os.path.getsize(encrhmntPlt)<300)  or h5ad_fl is not None):
						pltCllTypEnrchmnt(encrhmntFl,encrhmntPlt, \
						trshld=trshld,grtrThn=1)
				encrhmntFl = os.path.join(enrchmntAnlysFldr, \
				'%s.tsv'%fldrNm)
				encrhmntPlt = os.path.join(enrchmntAnlysFldr, \
				'%s_top10.svg'%fldrNm)
				if os.path.exists(encrhmntFl) and \
				((not os.path.exists(encrhmntPlt) or \
				os.path.getsize(encrhmntPlt)<300) or h5ad_fl is not None):
					pltCllTypEnrchmnt(encrhmntFl,encrhmntPlt, \
					trshld=trshld,grtrThn=1,topK=10,ar_setClstrs_in=ar_setClstrs, \
					ar_setClstrsClrs_in=ar_setClstrsClrs)
	print 'Done!'

