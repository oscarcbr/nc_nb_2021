#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  gnExprssnDffncs_mstr.py
#  
#  Copyright 2018 Oscar C. Bedoya Reina <oscar@oscar-J53kOiSr>
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
Determine marker enrichments for independent populations
"""
import argparse,os,sys

from singlecell.between_samples import rtrnCntsSmplINClstrs, \
clcPrWsFDRforEchClstr,clcPrWsFDRforClstrAvrg

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
 
######################
#  Input parameters  #
######################
parser.add_argument('-o','--outFldrPrnt',help='Output folder parental',default=None)
parser.add_argument('-N','--clcCllNmbrsInClstrs',help='If True will calculate the cell number differences between clusters and samples',default=True,type=str2bool)
parser.add_argument('-e','--prWsFDRForEchClstr',help='If True will calculate the genes that characterize each cluster. This is the intersection of all upregulated in each cluster vs. all others.',default=True,type=str2bool)
parser.add_argument('-E','--prWsFDRForEchClstrAvrg',help='If True will calculate the genes that characterize each cluster. This is the intersection of all upregulated in each cluster vs. all others by average.',default=True,type=str2bool)
parser.add_argument('-d','--prWsFDRForDepltdClstr',help='If True will calculate the genes that characterize each cluster. This is the intersection of all DOWNregulated in each cluster vs. all others.',default=True,type=str2bool)
parser.add_argument('-D','--prWsFDRForDepltdClstrAvrg',help='If True will calculate the genes that characterize each cluster. This is the intersection of all downregulated in each cluster vs. all others by average.',default=True,type=str2bool)
parser.add_argument('-c','--clcSgnfcntCmbtnOvrlpng',help='If True will calculate the intersection for every combination of target and background population.',default=False,type=str2bool)
parser.add_argument('-t','--clcUprgltnDffNBStgClstr',help='If True will calculate the different upregulation of genes for each NB stage and cluster.',default=False,type=str2bool)
parser.add_argument('-u','--clcUprgltnDffNBOutcmClstr',help='If True will Calculate the different upregulation of genes for each NB stage outcome and cluster.',default=False,type=str2bool)
parser.add_argument('-k','--lKmins',help='Input number of K clusters',default='5|8|9|10|11|12|13')
parser.add_argument('-g','--lSmplsGrpFl',help='File 1 with sample\tgroup\tstage\t...|File 2 with sample\tgroup\tstage\t...|...',default=None)
parser.add_argument('-s','--smplsClmn',help='Column in each smplsGrpFl with the sample name',default=0,type=int)
parser.add_argument('-G','--grpsClmn',help='Column in each smplsGrpFl with the group name',default=None,type=int)
parser.add_argument('-O','--outcmClmn',help='Column in each smplsGrpFl with the outcome name',default=None,type=int)
parser.add_argument('-x','--lExcldSmpls',help='List of samples to exclude: excldSmpl1|excldSmpl2|...',default=None)
parser.add_argument('-C','--complDBs',help='#If true will compile databases for files',default=True,type=str2bool)


args = parser.parse_args()

######################
######################
######################
#  Input parameters  #
######################
######################
######################
#Input requirements
outFldrPrnt = args.outFldrPrnt
clcCllNmbrsInClstrs  = args.clcCllNmbrsInClstrs
prWsFDRForEchClstr  = args.prWsFDRForEchClstr
prWsFDRForEchClstrAvrg = args.prWsFDRForEchClstrAvrg
prWsFDRForDepltdClstr = args.prWsFDRForDepltdClstr
prWsFDRForDepltdClstrAvrg = args.prWsFDRForDepltdClstrAvrg
clcSgnfcntCmbtnOvrlpng = args.clcSgnfcntCmbtnOvrlpng
clcUprgltnDffNBStgClstr = args.clcUprgltnDffNBStgClstr
clcUprgltnDffNBOutcmClstr = args.clcUprgltnDffNBOutcmClstr
lKmins = args.lKmins
lSmplsGrpFl = args.lSmplsGrpFl
smplsClmn = args.smplsClmn
grpsClmn = args.grpsClmn
outcmClmn = args.outcmClmn
lExcldSmpls = args.lExcldSmpls
complDBs = args.complDBs

assert os.path.exists(outFldrPrnt)
assert lKmins is not None
assert lSmplsGrpFl is not None

#Build lists of output folders, samples and number of clusters
l_rsltFldrs=[]
l_KNum=[]
l_KNum.extend([int(v) for v in lKmins.split('|') if v.strip()])
for kNum in l_KNum:
	if kNum==2:
		outFldr = os.path.join(outFldrPrnt,'KMin%s.d'%kNum)
	else:
		outFldr = os.path.join(outFldrPrnt,'KEql%s.d'%kNum)
	assert os.path.exists(outFldr)
	l_rsltFldrs.append(outFldr)


#Test for samples consistency and define species names
assert lSmplsGrpFl is not None
lSmplsGrpFl = sorted([s for s in lSmplsGrpFl.split('|') if s.strip()])
if lExcldSmpls is not None:
	sSmplsTExcld = set([s for s in lExcldSmpls.split('|') if s.strip()])
else:
	sSmplsTExcld = set()
	
	

#Define a dictionary of names and samples to make figures and applications
strdSmplsRef,dSmplGrp,dSmplOutcome = [],{},{}
#
for smplsGrpFl in lSmplsGrpFl:
	for l in open(smplsGrpFl,'r'):
		if l.strip() and l[0]!='#':
			l = l.splitlines()[0].split('\t')
			smpl = l[smplsClmn]
			if smpl not in sSmplsTExcld:
				strdSmplsRef.append(smpl)
				if grpsClmn is not None:
					grp = l[grpsClmn]
					assert not dSmplGrp.has_key(smpl)
					dSmplGrp[smpl] = grp
				else:
					dSmplGrp[smpl] = 'U'
				#
				if outcmClmn is not None:
					outcm = l[outcmClmn]
					assert not dSmplOutcome.has_key(smpl)
					dSmplOutcome[smpl] = outcm
				else:
					dSmplOutcome[smpl] = 'U'



strdSmplsRef.sort()

#Print log information
print('Log info:')
print('outFldrPrnt -> ',outFldrPrnt)
print('clcCllNmbrsInClstrs -> ',clcCllNmbrsInClstrs)
print('prWsFDRForEchClstr -> ',prWsFDRForEchClstr)
print('prWsFDRForEchClstrAvrg -> ',prWsFDRForEchClstrAvrg)
print('prWsFDRForDepltdClstr -> ',prWsFDRForDepltdClstr)
print('prWsFDRForDepltdClstrAvrg -> ',prWsFDRForDepltdClstrAvrg)
print('clcSgnfcntCmbtnOvrlpng -> ',clcSgnfcntCmbtnOvrlpng)
print('clcUprgltnDffNBStgClstr -> ',clcUprgltnDffNBStgClstr)
print('clcUprgltnDffNBOutcmClstr -> ',clcUprgltnDffNBOutcmClstr)
print('lKmins -> ',lKmins)
print('lSmplsGrpFl -> ',lSmplsGrpFl)
print('smplsClmn -> ',smplsClmn)
print('grpsClmn -> ',grpsClmn)
print('outcmClmn -> ',outcmClmn)
print('strdSmplsRef -> ',strdSmplsRef)
print('dSmplGrp -> ',dSmplGrp)
print('dSmplOutcome -> ',dSmplOutcome)
print('complDBs -> ',complDBs)


########################################################################
########################################################################
########################################################################
#########################                      #########################
######################## Calculate cell numbers ########################
#########################                      #########################
########################################################################
########################################################################
########################################################################

"""
First let's see the number of cells called for each cluster by the 
hierarchical clustering method used.
"""
if clcCllNmbrsInClstrs:
	print "Calculating cluster's cell numbers enrichment..."
	#Run analysis for the various folders
	for p in xrange(len(l_rsltFldrs)):
		fldrAnlys = os.path.join(l_rsltFldrs[p],'cllNmbrs.d')
		if not os.path.exists(fldrAnlys):
			os.mkdir(fldrAnlys)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		outFlIndvls = os.path.join(fldrAnlys,'indvls.tsv')
		outFlSmpls = os.path.join(fldrAnlys,'smpls.tsv')
		outPltIndvls = os.path.join(fldrAnlys,'indvls.svg')
		outPltSmpls = os.path.join(fldrAnlys,'smpls.svg')
		varNormDepthRds = os.path.join(l_rsltFldrs[p],'varNormDepth.rds')
		hcRds = os.path.join(l_rsltFldrs[p],'hc.rds')
		#make dictionary of clusters, samples and cell numbers
		ar_ClstrSmpl,ar_clstrs,ar_smpls = rtrnCntsSmplINClstrs \
		(varNormDepthRds,hcRds,dSmplGrp,outFlSmpls,outFlIndvls, \
		outPltSmpls,outPltIndvls,strdSmplsRef,k_num=k_num,skpErrs=True)
	print 'Done!'


########################################################################
########################################################################
########################################################################
#########################                      #########################
####################### Specific gene signature ########################
#########################                      #########################
########################################################################
########################################################################
########################################################################
"""
Let's calculate the genes that characterize each cluster. This is the
intersection of all upregulated in each cluster vs. all others.
"""
if prWsFDRForEchClstr:
	#Run analysis for the various folders
	print "Calculating cluster's cell numbers enrichment in each sample..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p],'dffExprssn.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		varNormDepthRds = os.path.join(l_rsltFldrs[p],'varNormDepth.rds')
		hcRds = os.path.join(l_rsltFldrs[p],'hc.rds')
		#list of clusters of interest
		lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR enrichment
		clcPrWsFDRforEchClstr(varNormDepthRds,hcRds,lPoptns, \
		enrchmntAnlysFldr,k_num=k_num,rplcFls=True)
	print 'Done!'

"""
Let's calculate the genes that characterize each cluster. This is the
intersection of all downregulated in each cluster vs. all others.
"""
if prWsFDRForDepltdClstr:
	#Run analysis for the various folders
	print "Calculating cluster's cell numbers enrichment in each sample..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p], \
		'dffExprssn.negative.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		varNormDepthRds = os.path.join(l_rsltFldrs[p],'varNormDepth.rds')
		hcRds = os.path.join(l_rsltFldrs[p],'hc.rds')
		#list of clusters of interest
		lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR depletion
		clcPrWsFDRforEchClstr(varNormDepthRds,hcRds,lPoptns, \
		enrchmntAnlysFldr,k_num=k_num,rplcFls=True,grtr=False)
	print 'Done!'


########################################################################
########################################################################
########################################################################
#########################                      #########################
######################### Average calculations #########################
#########################                      #########################
########################################################################
########################################################################
########################################################################
"""
Let's calculate the genes that characterize each cluster. This is the
intersection of all upregulated in each cluster vs. all others by average
"""
if prWsFDRForEchClstrAvrg:
	#Run analysis for the various folders
	print "Calculating cluster's significance difference by average in each sample..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p],'dffExprssnAvrg.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		varNormDepthRds = os.path.join(l_rsltFldrs[p],'varNormDepth.rds')#l_appFls[p])
		hcRds = os.path.join(l_rsltFldrs[p],'hc.rds')
		#list of clusters of interest
		lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR enrichment
		clcPrWsFDRforClstrAvrg(varNormDepthRds,hcRds,lPoptns, \
		enrchmntAnlysFldr,k_num=k_num,rplcFls=False)
	print 'Done!'
		
"""
Let's calculate the genes that characterize each cluster. This is the
intersection of all downregulated in each cluster vs. all others by average
"""
if prWsFDRForDepltdClstrAvrg:
	#Run analysis for the various folders
	print "Calculating cluster's significance difference by average in each sample..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p], \
		'dffExprssnAvrg.negative.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		varNormDepthRds = os.path.join(l_rsltFldrs[p],'varNormDepth.rds')
		hcRds = os.path.join(l_rsltFldrs[p],'hc.rds')
		#list of clusters of interest
		lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR enrichment
		clcPrWsFDRforClstrAvrg(varNormDepthRds,hcRds,lPoptns, \
		enrchmntAnlysFldr,k_num=k_num,rplcFls=False,grtr=False)
	print 'Done!'


