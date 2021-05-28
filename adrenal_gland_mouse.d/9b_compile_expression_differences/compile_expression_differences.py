#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  compile_expression_differences.py
#  
#  Copyright 2021 Oscar C. Bedoya Reina <oscarbed@ki.se>
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
These scripts aim to merge the difference expression of genes in single
files.
"""

import os

###########################
# Input folders and files #
###########################

#File with databases of interest to compile
inFldIndx='dpdtClstAnntnDff.tsv'
lIndx=[l.splitlines()[0].split('\t') for l in open(inFldIndx).read(). \
splitlines() if l.strip() and l[0]!='#']


#Folder and files to transform GRCm38_p6 gene symbols to genecodeV18Comp
srcFldrAnnttn='$PATH_to_"00_prepare_datasets"'
lAnnttnFls = ['GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv', \
'GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv']


###########
# Methods #
###########
#Method to make symbolic links
def mkSLnkDffExprsn(srcFldr,dptFldr,lSrcF = ['dffExprssn.d', \
'dffExprssn.negative.d','dffExprssnAvrg.d','dffExprssnAvrg.negative.d']):
	"""
	This method makes symbolic links from a source folder to a depot 
	folder.
	Input: srcFldr is the source folder. dptFldr is the output folder. 
	lSrcF is optionally a list of files to make the symbolic link.
	"""
	for srcF in lSrcF:
		if not os.path.exists('%s/%s'%(dptFldr,srcF)):
			os.system('ln -sf %s/%s %s/.'%(srcFldr,srcF,dptFldr))
	return 0
	
#To create a dictionary from the annotation
def mkDctrFrmAnntn(inAnnttn):
	"""
	Method to return a dictionary from a 2-columns inAnnttn file.
	Input: inAnnttn is the input annotation file.
	Output: dClstrNms is the output dictionary with keys in the first
	column and values in the second.
	"""
	dClstrNms = dict([l.split('\t') for l in open(inAnnttn).read(). \
	splitlines() if l.strip() and l[0]!='#'])
	return dClstrNms

#To merge the average expression in the databases
def rtrnd(infl,clmnK=0,clmnV=1):
	"""
	Method to merge the average gene expression of each cluster a single
	dictionary.
	Input: infl is the input file with average expression, clmnK is the 
	key (cluster) column position (0-based), and clmnV is the value 
	(cluster) column position (0-based).
	Output: dKV is the output dictionary.
	"""
	dKV={}
	for l in open(infl,'r'):
		if l.strip() and l.split('\t')[0]!='gene':
			l=l.splitlines()[0]
			k=l.split('\t')[clmnK]
			v=l.split('\t')[clmnV]
			if dKV.has_key(k):
				raise Exception('Error %s,%s'%(k,v))
				dKV[k].add(v)
			else:
				dKV[k]=set([v])
	dKV=dict([(k,'|'.join(v)) for k,v in dKV.items()])
	return dKV

#Write specific gene signature files
def wrtEnrchmnt(dClstrNms,outFlUpODwnrgltd,dptFldr):
	"""
	Method to write specific gene signature files.
	Input: dClstrNms is a dictionary of clusters and names. 
	outFlUpODwnrgltd is the output file. dptFldr is the path to the 
	folder with expression differences folders.
	"""
	#
	outFlUpODwnrgltd=open(outFlUpODwnrgltd,'w')
	outFlUpODwnrgltd.write('#Gene\tCell_type\n')
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssn.d','%s_allPairs.tsv'%clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_uprgltd\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssn.negative.d','%s_allPairs.tsv'% \
		clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_dwnrgltd\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssn.d','%s_allPairs.tsv'%clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_upOdwn\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssn.negative.d','%s_allPairs.tsv'% \
		clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_upOdwn\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	outFlUpODwnrgltd.close()
	#
	return 0

#Write genes significantly expressed in each cluster
def wrtEnrchmntAvrg(dClstrNms,outFlUpODwnrgltd,dptFldr):
	"""
	Method to write genes significantly expressed in each cluster.
	Input: dClstrNms is a dictionary of clusters and names. 
	outFlUpODwnrgltd is the output file. dptFldr is the path to the 
	folder with expression differences folders.
	"""
	outFlUpODwnrgltd=open(outFlUpODwnrgltd,'w')
	outFlUpODwnrgltd.write('#Gene\tCell_type\n')
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssnAvrg.d', \
		'clstr_%s_avrgAll.ovrlClstr.tsv'%clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#' and l.split('\t')[0]!='gene':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_uprgltd\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssnAvrg.negative.d', \
		'clstr_%s_avrgAll.ovrlClstr.tsv'%clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#' and l.split('\t')[0]!='gene':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_dwnrgltd\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssnAvrg.d', \
		'clstr_%s_avrgAll.ovrlClstr.tsv'%clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#' and l.split('\t')[0]!='gene':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_upOdwn\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	for clstr in dClstrNms.keys():
		inFl=os.path.join(dptFldr,'dffExprssnAvrg.negative.d', \
		'clstr_%s_avrgAll.ovrlClstr.tsv'%clstr)
		for l in open(inFl,'r'):
			if l.strip() and l[0]!='#' and l.split('\t')[0]!='gene':
				gnNm=l.split('\t')[0]
				outFlUpODwnrgltd.write('%s\t%s_upOdwn\n'%(gnNm, \
				dClstrNms[clstr]))
	#
	outFlUpODwnrgltd.close()
	return 0

#Write average gene expression in each cluster
def wrtCrrltn(dClstrNms,outFlUpODwnrgltd,dptFldr):
	"""
	Method to write the average gene expression in each cluster.
	Input: dClstrNms is a dictionary of clusters and names. 
	outFlUpODwnrgltd is the output file. dptFldr is the path to the 
	folder with expression differences folders.
	"""
	dGnNmsClstrExprssn={}
	sGnNms=set()
	sClstrs=set()
	for clstr in dClstrNms.keys():
			sClstrs.add(clstr)
			cllnNm=dClstrNms[clstr]
			dGnNmsExprssn=rtrnd(os.path.join(dptFldr,'dffExprssnAvrg.d', \
			'clstr_%s_avrgAll.tsv'%clstr))
			for gnNm in dGnNmsExprssn.keys():
					sGnNms.add(gnNm)
					if dGnNmsClstrExprssn.has_key(gnNm):
						dGnNmsClstrExprssn[gnNm][clstr]=dGnNmsExprssn[gnNm]
					else:
						dGnNmsClstrExprssn[gnNm]={clstr:dGnNmsExprssn[gnNm]}
	#
	sGnNms=sorted(sGnNms)
	sClstrs=sorted(sClstrs)
	outFl=open(outFlUpODwnrgltd,'w')
	outFl.write('#Gene\t%s\n'%'\t'.join([dClstrNms[c] for c in sClstrs]))
	for gnNm in sGnNms:
			outFl.write('%s\t%s\n'%(gnNm,'\t'.join([dGnNmsClstrExprssn[gnNm][c] \
			for c in sClstrs])))
	outFl.close()
	return 0



###########
# Execute #
###########
for dptFldr,inAnnttn,srcFldr,prfx in lIndx:
	dptFldrPrnt = os.path.split(dptFldr)[0]
	outFlUpODwnrgltd = os.path.join(dptFldrPrnt,'%s.upOdwn.tsv'%prfx)
	outFlUpODwnrgltdAvrg = os.path.join(dptFldrPrnt,'%s.upOdwn.avrg.tsv'% \
	prfx)
	outFlUpODwnrgltdCrrltn = os.path.join(dptFldrPrnt,'%s.crrltn.tsv'%prfx)
	#
	if not os.path.exists(dptFldr):
		os.system('mkdir -p %s'%dptFldr)
	#Make symbolic links to annotation files
	mkSLnkDffExprsn(srcFldr,dptFldr)
	#Make symbolic links to expression folders
	mkSLnkDffExprsn(srcFldrAnnttn,dptFldr,lAnnttnFls)
	#Make dictionary from annotation
	dClstrNms = mkDctrFrmAnntn(inAnnttn)
	#Run differential expression compilation -specific gene signature-
	wrtEnrchmnt(dClstrNms,outFlUpODwnrgltd,dptFldr)
	#Run differential expression compilation
	wrtEnrchmntAvrg(dClstrNms,outFlUpODwnrgltdAvrg,dptFldr)
	#Run average expression compilation
	wrtCrrltn(dClstrNms,outFlUpODwnrgltdCrrltn,dptFldr)


