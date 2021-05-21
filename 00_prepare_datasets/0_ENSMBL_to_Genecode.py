#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ENSMBL_to_Genecode.py
#  
#  Copyright 2020 Oscar C. Bedoya Reina <oscarbed@ki.se>
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
Script to compare the GRCh38_p12/GRCm38_p6 with genecodeV28Comp/genecodeV18Comp respectively
"""

import os

######################
######################
#####  Switches  #####
######################
######################
cmprAnntns = True # If True will compare the GRCh38_p12/GRCm38_p6 with genecodeV28Comp/genecodeV18Comp respectively

inFl_hg38p12 = 'GRCh38_p12_ENSMBLgnNm.tsv'
inFl_mm10p6 = 'GRCm38_p6_ENSEMBLgnNmEdtd.txt'
inFl_hg38anntn = 'cmmnGnNmTOEnsmbl.hg38'
inFl_mm10anntn = 'cmmnGnNmTOEnsmbl.mm10'

sGns_hg38p12 = set([l.splitlines()[0].split('\t')[1] for l in \
open(inFl_hg38p12,'r') if l.strip() and l.split('\t')[0]!='Gene name' \
and l[0]!='#'])
sGns_mm10p6 = set([l.splitlines()[0].split('\t')[1] for l in \
open(inFl_mm10p6,'r') if l.strip() and l.split('\t')[0]!='Gene name' \
and l[0]!='#'])
dSppGns = {'Human':sGns_hg38p12,'Mouse':sGns_mm10p6}



"""
The following is to make  a comparison between the ENSEMBL's GRCh38_p12 and
GRCm38_p6 annotations, and genecodeV28Comp and genecodeV18Comp annotations,
respectively. 
"""
######################
######################
if cmprAnntns:
	#Look for codes with redundant names and correct
	def cntGns(inf,clGnNm=1,clnEnsmbl=0,maxNumbr=1,sGnsIntrst=None):
		"""
		Method to return a dictionary of ENSEMBL codes and gene names.
		Output: dOutFnl is a dictionary of gene names as keys and ENSEMBL
		codes as values. sEnsmblMltpl is a set of ensembl codes with
		redundancies.
		"""
		dOut={}
		for l in open(inf,'r'):
			if l.strip() and l[0]!='#':
				l=l.splitlines()[0].split('\t')
				ensmbl,gnNm=l[clnEnsmbl],l[clGnNm]
				if sGnsIntrst is None: 
					if dOut.has_key(gnNm):
						dOut[gnNm].add(ensmbl)
					else:
						dOut[gnNm]=set([ensmbl])
				elif ensmbl.split('.')[0] in sGnsIntrst:
					if dOut.has_key(gnNm):
						dOut[gnNm].add(ensmbl)
					else:
						dOut[gnNm]=set([ensmbl])
		#
		dOutFnl={}
		sEnsmblMltpl = set()
		for k,v in dOut.items():
			if len(v)>maxNumbr:
				dOutFnl[k]=v
				sEnsmblMltpl.update(v)
		return dOutFnl,sEnsmblMltpl
	#
	#p12/p6 versions to genecodeV28Comp/genecodeV18Comp files
	outFlCnvrt_mm10 = 'GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv'
	outFlCnvrt_hg38 = 'GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv'
	flNonAltrntv_mm10 = 'GRCm38_p6_ENSMBLGnNmChr_gnNm.noAltrntvChrms.tsv'
	flNonAltrntv_hg38 = 'GRCh38_p12_ENSMBLGnNmChr_gnNm.noAltrntvChrms.tsv'
	sNonAltrntvGnNms_hg38 = set([l.split('\t')[0] for l in \
	open(flNonAltrntv_hg38,'r') if l.strip() and l[0]!='#'])
	sNonAltrntvGnNms_mm10 = set([l.split('\t')[0] for l in \
	open(flNonAltrntv_mm10,'r') if l.strip() and l[0]!='#'])
	#Look for duplicated and ambiguous annotations to exclude
	dhg38p12NmsCnts,sToDel = cntGns(inFl_hg38p12)
	dmm10p6NmsCnts,sToDel = cntGns(inFl_mm10p6)
	dhg38anntnNmsCnts,sToDel = cntGns(inFl_hg38anntn,clnEnsmbl=1,clGnNm=0)
	dmm10anntnNmsCnts,sToDel = cntGns(inFl_mm10anntn,clnEnsmbl=1,clGnNm=0)
	#In comparison with ENSEMBL genes only present in non-alternative chromosome
	dhg38p12NmsCnts,sEnsmblMltplhg38p12NmsCnts = cntGns(inFl_hg38p12, \
	sGnsIntrst=sNonAltrntvGnNms_hg38)
	dmm10p6NmsCnts,sEnsmblMltplmm10p6NmsCnts = cntGns(inFl_mm10p6, \
	sGnsIntrst=sNonAltrntvGnNms_mm10)
	dhg38anntnNmsCnts,sEnsmblMltplhg38anntnNmsCnts = cntGns(inFl_hg38anntn, \
	clnEnsmbl=1,clGnNm=0,sGnsIntrst=sNonAltrntvGnNms_hg38)
	dmm10anntnNmsCnts,sEnsmblMltplmm10anntnNmsCnts = cntGns(inFl_mm10anntn, \
	clnEnsmbl=1,clGnNm=0,sGnsIntrst=sNonAltrntvGnNms_mm10)
	##########
	####NOTE: SOD2 and DIABLO for being an important MES marker were 
	####included (even if they had redundant ENSEMBL code) after manual 
	####curation under the ensembl codes ENSG00000112096 and
	####ENSG00000184047 respectively.
	sEnsmblMltplhg38p12NmsCnts.remove('ENSG00000112096')
	sEnsmblMltplhg38p12NmsCnts.remove('ENSG00000184047')
	##########
	#Compare the ENSEMBL codes
	#Get ENSEMBL code for the p12/p6 annotations
	dENSM_hg38p12 = dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(inFl_hg38p12).read().splitlines() if l.strip() and l.split('\t')[0]!='Gene name' \
	and l[0]!='#' and l.split('\t')[0] not in sEnsmblMltplhg38p12NmsCnts])
	dENSM_mm10p6 =  dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(inFl_mm10p6).read().splitlines() if l.strip() and l.split('\t')[0]!='Gene name'\
	and l[0]!='#' and l.split('\t')[0] not in sEnsmblMltplmm10p6NmsCnts])
	#Get ENSEMBL code for the genecodeV28Comp/genecodeV18Comp annotations
	dENSM_hg38anntn = dict([(l.split('\t')[1],l.split('\t')[0]) for l in \
	open(inFl_hg38anntn).read().splitlines() if l.strip() and l.split('\t')[0]!='Gene name' \
	and l[0]!='#' and l.split('\t')[1] not in sEnsmblMltplhg38anntnNmsCnts])
	dENSM_mm10anntn = dict([(l.split('\t')[1],l.split('\t')[0]) for l in \
	open(inFl_mm10anntn).read().splitlines() if l.strip() and l.split('\t')[0]!='Gene name'\
	and l[0]!='#' and l.split('\t')[1] not in sEnsmblMltplmm10anntnNmsCnts])
	#
	"""
	All these previous dictionaries have duplicated keys in the databases
	so, they need to be excluded in further comparisons
	"""
	sENSM_hg38p12=set(dENSM_hg38p12.keys())
	sENSM_mm10p6=set(dENSM_mm10p6.keys())
	sENSM_hg38anntn=set(dENSM_hg38anntn.keys())
	sENSM_mm10anntn=set(dENSM_mm10anntn.keys())
	sENSM_hg38Intrsctn=sENSM_hg38anntn.intersection(sENSM_hg38p12)
	sENSM_mm10Intrsctn=sENSM_mm10anntn.intersection(sENSM_mm10p6)
	#sorted(sENSM_hg38anntn.difference(sENSM_hg38p12)) and 
	#sorted(sENSM_mm10anntn.difference(sENSM_mm10p6)) were inspected 
	#manually. All of the codes in these lists are deprecated markers in
	#p12/p6 versions or duplicated/ambigously annotated.
	###########
	#Now datasets to convert p12/p6 versions to genecodeV28Comp/genecodeV18Comp
	#are going to be created
	ooutFlCnvrt_mm10 = open(outFlCnvrt_mm10,'w')
	ooutFlCnvrt_mm10.write('#GRCm38_p6 gene name\tgenecodeV18Comp gene name\tENSEMBL code\n')
	for ENSM in sorted(sENSM_mm10Intrsctn):
		ooutFlCnvrt_mm10.write('%s\n'%'\t'.join([dENSM_mm10p6[ENSM], \
		dENSM_mm10anntn[ENSM],ENSM]))
	ooutFlCnvrt_mm10.close()
	#
	ooutFlCnvrt_hg38 = open(outFlCnvrt_hg38,'w')
	ooutFlCnvrt_hg38.write('#GRCh38_p12 gene name\tgenecodeV18Comp gene name\tENSEMBL code\n')
	for ENSM in sorted(sENSM_hg38Intrsctn):
		ooutFlCnvrt_hg38.write('%s\n'%'\t'.join([dENSM_hg38p12[ENSM], \
		dENSM_hg38anntn[ENSM],ENSM]))
	ooutFlCnvrt_hg38.close()
	#
