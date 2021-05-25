#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  smplFltrNSlctHQclls.py
#  
#  Copyright 2018 Oscar C Bedoya-Reina <oscarbed@ki.se>
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
The idea is to run a simple filtering on number of genes and collect
other features
"""

import argparse,os,sys

from singlecell.formats import htseqLFlsToGRM
from singlecell.statistics import clcClltyFeatrsDbgd
from singlecell.cellity_mod import cmptCllnFtrs
from numpy import array,random,where
from datetime import datetime

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
parser.add_argument('-e','--exprssnFldr',help='Expression folder',default=None)
parser.add_argument('-q','--qcFldr',help='Folder with map STAR results',default=None)
parser.add_argument('-g','--gtfFile',help='GTF file',default=None)
parser.add_argument('-S','--species',help='Species',default=None)
parser.add_argument('-p','--dirSfx',help='sufix',default=None)
parser.add_argument('-c','--stdyCsSfx',help='study case',default=None)
parser.add_argument('-Q','--fastQMainFldr',help='fastQMainFldr',default=None)
parser.add_argument('-C','--cllCycl',help='cllCycl',default=None)
parser.add_argument('-o','--outFldr',help='outFldr',default=None)
parser.add_argument('-l','--slctBySmpl',help='slctBySmpl',default=False,type=str2bool)
parser.add_argument('-Y','--rCllFltr',help='run_simple filtering',default=True,type=str2bool)
parser.add_argument('-T','--rCllOlnFtrs',help='celloline on features',default=True,type=str2bool)
parser.add_argument('-t','--ftrsFldr',help='features folder',default=None)
parser.add_argument('-F','--shuffleSmpls',help='shuffle samples on celloline features',default=True,type=str2bool)
parser.add_argument('-E','--mkExprssnExcldLQclls',help='mk_exprssn_excludelowHQcells',default=True,type=str2bool)
parser.add_argument('-m','--ftrMinTrshld',help='Minimum number of genes for HQ cells',default=2000,type=int)
parser.add_argument('-M','--ftrMaxTrshld',help='Maximum number of genes for HQ cells',default=8000,type=int)

args = parser.parse_args()

######################
######################
######################
#  Input parameters  #
######################
######################
######################
smpls = args.smpls
exprssnFldr = args.exprssnFldr
qcFldr = args.qcFldr
ftrsFldr = args.ftrsFldr
gtfFile = args.gtfFile
species = args.species
dirSfx = args.dirSfx
stdyCsSfx = args.stdyCsSfx
fastQMainFldr = args.fastQMainFldr
cllCycl = args.cllCycl
outFldr = args.outFldr
slctBySmpl = args.slctBySmpl
rCllFltr = args.rCllFltr#run_celloline_htSeq_hg38
rCllOlnFtrs = args.rCllOlnFtrs#celloline on features
shuffleSmpls = args.shuffleSmpls#celloline on features
mkExprssnExcldLQclls = args.mkExprssnExcldLQclls#mk_exprssn_excludelowHQcells
ftrMinTrshld = args.ftrMinTrshld#mk_exprssn_excludelowHQcells
ftrMaxTrshld = args.ftrMaxTrshld#mk_exprssn_excludelowHQcells

"""
Let's run celloline.
"""

smplsNms = sorted([s for s in smpls.split('|') if s.strip()])
if shuffleSmpls:
	#
	crrntTme = datetime.now()
	trgtTme = datetime(crrntTme.year,crrntTme.month,crrntTme.day, \
	crrntTme.hour-1,crrntTme.minute,crrntTme.second,crrntTme.microsecond)
	#
	sSmplsDn = set(['.'.join(f.split('.')[:-2]) for f in os.listdir(ftrsFldr) \
	if os.path.getsize(os.path.join(ftrsFldr,f)) and (f.find('.P')>-1 or \
	datetime.fromtimestamp(os.path.getmtime(os.path.join(ftrsFldr,f)))> \
	trgtTme) and f.find('.allGns.')==-1])
	smplsNms = [s for s in smplsNms if s not in sSmplsDn]
	random.shuffle(smplsNms)
	print 'Shuffle set in order %s'%', '.join(smplsNms)


assert len(set(smplsNms))==len(smplsNms)
assert os.path.exists(exprssnFldr) and os.path.exists(qcFldr) and \
os.path.exists(gtfFile)
assert species in {'Human','Mouse'}
if species=='Mouse':
	sppPrfx='mm10'
else:
	sppPrfx='hg38'


cllNms=sorted(['A1','A19','A6','B15','B24','C11','C20','C8','D17','D4','E13', \
'E22','F1','F19','F6','G15','G24','H11','H20','H8','I17','I4','J13', \
'J22','K1','K19','K6','L15','L24','M11','M20','M8','N17','N4','O13', \
'O22','P1','P19','P6','A10','A2','A7','B16','B3','C12','C21','C9','D18', \
'D5','E14','E23','F10','F2','F7','G16','G3','H12','H21','H9','I18','I5', \
'J14','J23','K10','K2','K7','L16','L3','M12','M21','M9','N18','N5', \
'O14','O23','P10','P2','P7','A11','A20','A8','B17','B4','C13','C22','D1', \
'D19','D6','E15','E24','F11','F20','F8','G17','G4','H13','H22','I1','I19', \
'I6','J15','J24','K11','K20','K8','L17','L4','M13','M22','N1','N19','N6', \
'O15','O24','P11','P20','P8','A12','A21','A9','B18','B5','C14','C23', \
'D10','D2','D7','E16','E3','F12','F21','F9','G18','G5','H14','H23','I10', \
'I2','I7','J16','J3','K12','K21','K9','L18','L5','M14','M23','N10','N2', \
'N7','O16','O3','P12','P21','P9','A13','A22','B1','B19','B6','C15','C24', \
'D11','D20','D8','E17','E4','F13','F22','G1','G19','G6','H15','H24', \
'I11','I20','I8','J17','J4','K13','K22','L1','L19','L6','M15','M24', \
'N11','N20','N8','O17','O4','P13','P22','A14','A23','B10','B2','B7', \
'C16','C3','D12','D21','D9','E18','E5','F14','F23','G10','G2','G7','H16', \
'H3','I12','I21','I9','J18','J5','K14','K23','L10','L2','L7','M16','M3', \
'N12','N21','N9','O18','O5','P14','P23','A15','A24','B11','B20','B8', \
'C17','C4','D13','D22','E1','E19','E6','F15','F24','G11','G20','G8','H17', \
'H4','I13','I22','J1','J19','J6','K15','K24','L11','L20','L8','M17','M4', \
'N13','N22','O1','O19','O6','P15','P24','A16','A3','B12','B21','B9','C18', \
'C5','D14','D23','E10','E2','E7','F16','F3','G12','G21','G9','H18','H5', \
'I14','I23','J10','J2','J7','K16','K3','L12','L21','L9','M18','M5','N14', \
'N23','O10','O2','O7','P16','P3','A17','A4','B13','B22','C1','C19','C6', \
'D15','D24','E11','E20','E8','F17','F4','G13','G22','H1','H19','H6','I15', \
'I24','J11','J20','J8','K17','K4','L13','L22','M1','M19','M6','N15','N24', \
'O11','O20','O8','P17','P4','A18','A5','B14','B23','C10','C2','C7','D16', \
'D3','E12','E21','E9','F18','F5','G14','G23','H10','H2','H7','I16','I3', \
'J12','J21','J9','K18','K5','L14','L23','M10','M2','M7','N16','N3','O12', \
'O21','O9','P18','P5'])

"""
Process to obtain subsets of genes: cell cycling, mitochondrial, and 
only protein-coding
"""
#make a set of all mitochondrial genes in annotation and PCGs and no cell cycle
#Only protein coding
sExcldGnNmPCG = set([l.split('gene "')[1].split('"')[0] for l in open(gtfFile,'r') \
if l.strip() and l.split('gene_biotype "')[1].split('"')[0]!='protein_coding' and \
l.find('gene "ERCC-')==-1])
sGnNmIncld = set([l.split('gene "')[1].split('"')[0] for l in open(gtfFile,'r') \
if l.strip() and l.split('gene_biotype "')[1].split('"')[0]=='protein_coding' and \
l.find('gene "ERCC-')==-1])
for gnNm in sGnNmIncld.intersection(sExcldGnNmPCG):
	sExcldGnNmPCG.remove(gnNm)
#No mitochondrial
sExcldGnNmNOMito=set([l.split('gene "')[1].split('"')[0] for l in \
open(gtfFile,'r') if l.strip() and l.split()[0]=='chrM'])
#No cell cycle
sExcldGnNmNOcc=set([l.split()[0] for l in open(cllCycl,'r') if l.strip() \
and l[0]!='#'  and l.split()[1]=='G2/M'])
#No PCG and no mitochondrial
snMnPCG = sExcldGnNmPCG.union(sExcldGnNmNOMito)
#No cell cycle and no mitochondrial
snMnCC = sExcldGnNmNOcc.union(sExcldGnNmNOMito)
#No cell cycle and no mitochondrial and no PCG
snMnCCnPCG = sExcldGnNmNOcc.union(sExcldGnNmNOMito).union(sExcldGnNmPCG)

print 'slctBySmpl ==>',slctBySmpl


def smplCllFltrng(counts_fl,outFlFtrs,outFlClssfy,ftrMinTrshld=2000, \
	ftrMaxTrshld=8000,ftrNm='#Detected genes'):
	"""
	Method to run simple filtering based on the number of gene counts
	"""
	print 'Simple filtering by gene counts...'
	ar_ftrs = array([l.split('\t') for l in open(outFlFtrs).read().splitlines() if \
	l.strip()])
	posFtrNm = where(ar_ftrs[0,:]==ftrNm)[0]
	outL = ['cell\tquality']
	for cllFtrs in ar_ftrs[1:,]:
		if ftrMinTrshld<=int(cllFtrs[posFtrNm][0])<=ftrMaxTrshld:
			outL.append('%s\t%s'%(cllFtrs[0],1))
		else:
			outL.append('%s\t%s'%(cllFtrs[0],0))
	oOutFlClssfy = open(outFlClssfy,'w')
	oOutFlClssfy.write('\n'.join(outL))
	oOutFlClssfy.close()
	return 0
	


"""
The first step is to select only cells that FAILED at most three FastaQC
tests, and extract the features required to run celloline 
"""
#run_celloline_features_htSeq
if rCllOlnFtrs:
	ar_cllNm_all,ar_bamFlPths_all,ar_fastQflPths_all,ar_smplNm_all, \
	ar_outFls_all = [],[],[],[],[]	
	assert len(set(smplsNms))==len(smplsNms)
	#
	if shuffleSmpls:
		smplsReady = set()
		sSmplsDn = set(['.'.join(f.split('.')[:-2]) for f in os.listdir(ftrsFldr) \
		if os.path.getsize(os.path.join(ftrsFldr,f)) and (f.find('.P')>-1 or \
		datetime.fromtimestamp(os.path.getmtime(os.path.join(ftrsFldr,f)))> \
		trgtTme) and f.find('.allGns.')==-1])
		smplsNms = [s for s in smplsNms if s not in sSmplsDn. \
		union(smplsReady)]
		random.shuffle(smplsNms)
		print 'Shuffle set in order %s'%', '.join(smplsNms)
	#
	while smplsNms:
		smplsIn=smplsNms.pop()
		print 'Processing sample %s...'%smplsIn
		if shuffleSmpls:
			smplsReady.add(smplsIn)
		##
		ar_cllNm,ar_bamFlPths,ar_fastQflPths,ar_smplNm,ar_outFls = [], \
		[],[],[],[]
		smplFldr = os.path.join(qcFldr,'%s%sd'%(smplsIn,dirSfx))
		qcFl = os.path.join(qcFldr,'%s%sstats.txt'%(smplsIn,dirSfx))
		dCllNmsNmbrFls = dict([(l.split('\t')[0],l.count('FAIL')) for l \
		in open(qcFl).read().splitlines()[1:] if l.strip()])
		for cllNm in cllNms:
			nmbrQCFldTst = dCllNmsNmbrFls[cllNm]
			if nmbrQCFldTst<3:
				inpuSAM = os.path.join(smplFldr,'%s.bam'%cllNm)
				fstqFldr = os.path.join(fastQMainFldr,smplsIn,'rawdata', \
				cllNm)
				fastQfl = os.path.join(fstqFldr,'%s_R1_trimmed.fq.gz'% \
				cllNm)
				ftrsFl = os.path.join(ftrsFldr,'%s.%s.ftrs'%(smplsIn, \
				cllNm))
				ar_cllNm.append(cllNm)
				ar_smplNm.append(smplsIn)
				ar_bamFlPths.append(inpuSAM)
				ar_fastQflPths.append(fastQfl)
				ar_outFls.append(ftrsFl)
		#
		ar_cllNm_all.extend(ar_cllNm)
		ar_bamFlPths_all.extend(ar_bamFlPths)
		ar_fastQflPths_all.extend(ar_fastQflPths)
		ar_smplNm_all.extend(ar_smplNm)
		ar_outFls_all.extend(ar_outFls)
		#Build one feature file for each sample
		if slctBySmpl:
			clltyFtrsFl = os.path.join(ftrsFldr, \
			'%s.highQC.allGns.noKIFIsfrms.ftrs.cllty.csv'%smplsIn)
			clcClltyFeatrsDbgd(array(ar_bamFlPths),array(ar_fastQflPths), \
			gtfFile,clltyFtrsFl,array(ar_smplNm),array(ar_cllNm), \
			array(ar_outFls),incldERCC=True)
		#
		if shuffleSmpls:
			sSmplsDn = set(['.'.join(f.split('.')[:-2]) for f in \
			os.listdir(ftrsFldr) if os.path.getsize(os.path.join(ftrsFldr, \
			f)) and (f.find('.P')>-1 or datetime.fromtimestamp(os.path. \
			getmtime(os.path.join(ftrsFldr,f)))>trgtTme) and f.find \
			('.allGns.')==-1])
			smplsNms = [s for s in smplsNms if s not in sSmplsDn. \
			union(smplsReady)]
			random.shuffle(smplsNms)
			print 'Shuffle set in order %s'%', '.join(smplsNms)
	#
	if not slctBySmpl:
		ar_bamFlPths,ar_fastQflPths = array(ar_bamFlPths_all), \
		array(ar_fastQflPths_all)
		ar_smplNm,ar_cllNm,ar_outFls = array(ar_smplNm_all), \
		array(ar_cllNm_all),array(ar_outFls_all)
		#
		clltyFtrsFl = os.path.join(ftrsFldr, \
		'%s%s%shighQC.allGns.noKIFIsfrms.ftrs.cllty.csv'%(sppPrfx,dirSfx, \
		stdyCsSfx))
		clcClltyFeatrsDbgd(ar_bamFlPths,ar_fastQflPths,gtfFile,clltyFtrsFl, \
		ar_smplNm,ar_cllNm,ar_outFls,incldERCC=True)
	#
	if shuffleSmpls:
		raise Exception('Shuffle needs to be turned off in order to continue...')
	else:
		smplsNms = sorted([s for s in smpls.split('|') if s.strip()])


"""
The first step is to select only cells that FAILED at most three FastaQC
tests, and extract the features required to run simple cell filtering 
"""
#run_smplCllFltrng
if rCllFltr:
	if slctBySmpl:
		for smplsIn in smplsNms:
			counts_fl = os.path.join(exprssnFldr, \
			'%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.cllty.csv'% \
			(smplsIn,dirSfx))
			stats_fl = os.path.join(ftrsFldr, \
			'%s%shighQC.allGns.noKIFIsfrms.ftrs.cllty.csv'%(smplsIn, \
			dirSfx))
			outFlFtrs = os.path.join(outFldr, \
			'%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.allFtrs.cllty.csv'% \
			(smplsIn,dirSfx))
			outFlClssfy = os.path.join(outFldr, \
			'%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.clsfy.cllty.csv'% \
			(smplsIn,dirSfx))
			#Compute features and filter
			if not os.path.exists(outFlFtrs):
				cmptCllnFtrs(counts_fl,stats_fl,outFlFtrs,organism= \
				species.lower())
			smplCllFltrng(counts_fl,outFlFtrs,outFlClssfy,ftrMinTrshld, \
			ftrMaxTrshld)
	#Merged
	#
	else:
		counts_fl = os.path.join(exprssnFldr, \
		'%s%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.cllty.csv'% \
		(sppPrfx,dirSfx,stdyCsSfx))
		stats_fl = os.path.join(ftrsFldr, \
		'%s%s%shighQC.allGns.noKIFIsfrms.ftrs.cllty.csv'%(sppPrfx, \
		dirSfx,stdyCsSfx))
		outFlFtrs = os.path.join(outFldr, \
		'%s%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.allFtrs.cllty.csv'% \
		(sppPrfx,dirSfx,stdyCsSfx))
		outFlClssfy = os.path.join(outFldr, \
		'%s%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.clsfy.cllty.csv'% \
		(sppPrfx,dirSfx,stdyCsSfx))
		#Compute features and filter
		if not os.path.exists(outFlFtrs):
			cmptCllnFtrs(counts_fl,stats_fl,outFlFtrs,organism= \
			species.lower())
		smplCllFltrng(counts_fl,outFlFtrs,outFlClssfy,ftrMinTrshld, \
		ftrMaxTrshld)


#mk_exprssn_excludelowHQcells
if mkExprssnExcldLQclls:
	lFlsHtseqAll = []
	for smplsIn in smplsNms:
		print smplsIn
		exprssnFldrCnts = os.path.join(exprssnFldr,smplsIn)
		qcFl = os.path.join(qcFldr,'%s%sstats.txt'%(smplsIn,dirSfx))
		dCllNmsNmbrFls = dict([(l.split('\t')[0],l.count('FAIL')) for l \
		in open(qcFl).read().splitlines()[1:] if l.strip()])
		lFlsHtseq = []
		for cll in cllNms:
			#Check if the cell has not passed more than 2 quality control checks
			nmbrQCFldTst = dCllNmsNmbrFls[cll]
			if nmbrQCFldTst<3:
				smplCllFl = os.path.join(exprssnFldrCnts, \
				'%s%scounts.htsq'%(cll,dirSfx))
				lFlsHtseqAll.append(smplCllFl)
	####
	#Execute
	####
	for sGnsToExcld,flSfx in ((None,'.'),(sExcldGnNmNOcc,'.noCllCycl.'), \
	(sExcldGnNmNOMito,'.noMitchnd.'),(snMnPCG,'.noMitchnd.onlyPCG.'), \
	(snMnCC,'.noMitchnd.noCllCycl.'),(snMnCCnPCG, \
	'.noMitchnd.noCllCycl.onlyPCG.')):
		#Return low quality cell names	
		#Test if the test should be by sample or by joined samples
		if slctBySmpl:
			sCllsExcldLib = set()
			for smplsIn in smplsNms:
				clltyFl = os.path.join(outFldr, \
				'%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.clsfy.cllty.csv'% \
				(smplsIn,dirSfx))
				sCllsExcldLib.update(set([l.split()[0] for l in open(clltyFl,'r') if \
				l.strip() and l.splitlines()[0].split()[1]=='0']))
			outGRMFlGnsCnts = os.path.join(exprssnFldr, \
			'%s%s%shighQC.allGns.counts.htsq.clltyHQ.smplBsd%scsv'%(sppPrfx,dirSfx, \
			stdyCsSfx,flSfx))
			outGRMFlSpkInsCnts = os.path.join(exprssnFldr, \
			'%s%s%shighQC.spkIns.counts.htsq.clltyHQ.smplBsd%scsv'%(sppPrfx,dirSfx, \
			stdyCsSfx,flSfx))
			#
		else:
			clltyFl_allTiss = os.path.join(outFldr, \
			'%s%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.clsfy.cllty.csv'% \
			(sppPrfx,dirSfx,stdyCsSfx))
			sCllsExcldLib = set([l.split()[0] for l in open(clltyFl_allTiss,'r') 
			if l.strip() and l.splitlines()[0].split()[1]=='0'])
			outGRMFlGnsCnts = os.path.join(exprssnFldr, \
			'%s%s%shighQC.allGns.counts.htsq.clltyHQ%scsv'%(sppPrfx,dirSfx, \
			stdyCsSfx,flSfx))
			outGRMFlSpkInsCnts = os.path.join(exprssnFldr, \
			'%s%s%shighQC.spkIns.counts.htsq.clltyHQ%scsv'%(sppPrfx,dirSfx, \
			stdyCsSfx,flSfx))
		#
		#Format expressions
		if not os.path.exists(outGRMFlGnsCnts) or not \
		os.path.exists(outGRMFlSpkInsCnts):
			htseqLFlsToGRM(lFlsHtseqAll,outGRMFlGnsCnts,outGRMFlSpkInsCnts, \
			nmbrSpkIns=0,excldNmPrfx='%scounts.htsq'%dirSfx, \
			sExcldSmples=sCllsExcldLib,sExcldGnNm=sGnsToExcld)
