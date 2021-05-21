#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  format_quantification.py
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
The idea is to format the expression of genes of interest.
"""

import argparse,os,sys

from singlecell.formats import htseqLFlsToGRM

#################################
#         Parse inputs          #                
#################################
parser = argparse.ArgumentParser()
 
######################
#  Input parameters  #
######################
parser.add_argument('-s','--smpls',help='Input sample',default=None)
parser.add_argument('-e','--exprssnFldr',help='Expression folder',default=None)
parser.add_argument('-q','--qcFldr',help='Folder with map STAR results',default=None)
parser.add_argument('-E','--pthCmmnGnNmToENSEMBLG',help='File to convert gene symbols to ENSEMBL gene codes',default=None)
parser.add_argument('-S','--species',help='Species',default=None)
parser.add_argument('-p','--dirSfx',help='Additional sufix',default=None)
parser.add_argument('-c','--stdyCsSfx',help='Study case sufix',default=None)

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
pthCmmnGnNmToENSEMBLG = args.pthCmmnGnNmToENSEMBLG
species = args.species
dirSfx = args.dirSfx
stdyCsSfx = args.stdyCsSfx

smplsNms = sorted([s for s in smpls.split('|') if s.strip()])

print('Log info:')
print('\tsmpls -> ',smpls)
print('\texprssnFldr -> ',exprssnFldr)
print('\tpthCmmnGnNmToENSEMBLG -> ',pthCmmnGnNmToENSEMBLG)
print('\tqcFldr -> ',qcFldr)
print('\tspecies -> ',species)
print('\tdirSfx -> ',dirSfx)
print('\tstdyCsSfx -> ',stdyCsSfx)

assert os.path.exists(exprssnFldr) and os.path.exists(qcFldr) and \
os.path.exists(pthCmmnGnNmToENSEMBLG)
assert species in {'Human','Mouse'}
if species=='Mouse':
	sppPrfx='mm10'
else:
	sppPrfx='hg38'

#List of sorted cells
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
sExcldGnNm = set(['KIF1bAlpha','KIF1bBeta'])#Exclude KIF1B isoforms from output ENSEMBL counts

#################
#################
#################
###  Execute  ###
#################
#################
#################
#
lFlsHtseqAll = []
for smplssIn in smplsNms:
	print smplssIn
	exprssnFldrCnts = os.path.join(exprssnFldr,smplssIn)
	qcFl = os.path.join(qcFldr,'%s%sstats.txt'%(smplssIn,dirSfx))
	dCllNmsNmbrFls = dict([(l.split('\t')[0],l.count('FAIL')) for l \
	in open(qcFl).read().splitlines()[1:] if l.strip()])
	lFlsHtseq = []
	for cll in cllNms:
		#Check if the cell has not passed more than 2 quality control checks
		nmbrQCFldTst = dCllNmsNmbrFls[cll]
		if nmbrQCFldTst<3:
			smplsCllFl = os.path.join(exprssnFldrCnts, \
			'%s%scounts.htsq'%(cll,dirSfx))
			lFlsHtseqAll.append(smplsCllFl)
			lFlsHtseq.append(smplsCllFl)
	#Write quantification for each sample
	outClltyFlGnsCnts = os.path.join(exprssnFldr, \
	'%s%shighQC.allGns.counts.htsq.csv'%(smplssIn,dirSfx))
	outClltyFlGnsCntsNoKIFIsfrms = os.path.join(exprssnFldr, \
	'%s%shighQC.allGns.noKIFIsfrms.counts.htsq.csv'%(smplssIn,dirSfx))
	outClltyFlSpkInsCnts = os.path.join(exprssnFldr, \
	'%s%shighQC.spkIns.counts.htsq.csv'%(smplssIn,dirSfx))
	outClltyFlGnsCntsEnsmbl = os.path.join(exprssnFldr, \
	'%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.cllty.csv'% \
	(smplssIn,dirSfx))
	outClltyFlSpkInsCntsEnsmbl = os.path.join(exprssnFldr, \
	'%s%shighQC.spkIns.noKIFIsfrms.ensmbl.counts.htsq.cllty.csv'% \
	(smplssIn,dirSfx))
	#
	if (not os.path.exists(outClltyFlGnsCnts) or not \
	os.path.exists(outClltyFlGnsCntsNoKIFIsfrms) or not \
	os.path.exists(outClltyFlGnsCntsEnsmbl)): 
		htseqLFlsToGRM(lFlsHtseq,outClltyFlGnsCntsNoKIFIsfrms, \
		outClltyFlSpkInsCnts,nmbrSpkIns=0,rowsNm='',spltChrc='\t', \
		sExcldGnNm=sExcldGnNm,excldNmPrfx='%scounts.htsq'%dirSfx, \
		excldERCCfrmCnts=True)
		htseqLFlsToGRM(lFlsHtseq,outClltyFlGnsCnts,outClltyFlSpkInsCnts, \
		nmbrSpkIns=0,rowsNm='',spltChrc='\t',excldNmPrfx='%scounts.htsq'%dirSfx, \
		excldERCCfrmCnts=True)
		htseqLFlsToGRM(lFlsHtseq,outClltyFlGnsCntsEnsmbl, \
		outClltyFlSpkInsCntsEnsmbl,nmbrSpkIns=0,rowsNm='',spltChrc='\t', \
		sExcldGnNm=sExcldGnNm,excldNmPrfx='%scounts.htsq'%dirSfx, \
		pthCmmnGnNmToENSEMBLG=pthCmmnGnNmToENSEMBLG,nrmlzDpth=True, \
		excldERCCfrmCnts=True)

#Write quantification for all samples
outGRMFlGnsCntsNoKIFIsfrms = os.path.join(exprssnFldr, \
'%s%s%shighQC.allGns.noKIFIsfrms.counts.htsq.csv'%(sppPrfx,dirSfx, \
stdyCsSfx))
outGRMFlGnsCnts = os.path.join(exprssnFldr, \
'%s%s%shighQC.allGns.counts.htsq.csv'%(sppPrfx,dirSfx,stdyCsSfx))
outGRMFlSpkInsCnts = os.path.join(exprssnFldr, \
'%s%s%shighQC.spkIns.counts.htsq.csv'%(sppPrfx,dirSfx,stdyCsSfx))
htseqLFlsToGRM(lFlsHtseqAll,outGRMFlGnsCntsNoKIFIsfrms, \
outGRMFlSpkInsCnts,nmbrSpkIns=0,rowsNm='',spltChrc='\t', \
sExcldGnNm=sExcldGnNm,excldNmPrfx='%scounts.htsq'%dirSfx, \
excldERCCfrmCnts=True)
htseqLFlsToGRM(lFlsHtseqAll,outGRMFlGnsCnts,outGRMFlSpkInsCnts, \
nmbrSpkIns=0,excldNmPrfx='%scounts.htsq'%dirSfx,excldERCCfrmCnts=True)
outGRMFlGnsCntsEnsmbl = os.path.join(exprssnFldr, \
'%s%s%shighQC.allGns.noKIFIsfrms.ensmbl.counts.htsq.cllty.csv'% \
(sppPrfx,dirSfx,stdyCsSfx))
outGRMFlSpkInsCntsEnsmbl = os.path.join(exprssnFldr, \
'%s%s%shighQC.spkIns.noKIFIsfrms.ensmbl.counts.htsq.cllty.csv'% \
(sppPrfx,dirSfx,stdyCsSfx))
htseqLFlsToGRM(lFlsHtseqAll,outGRMFlGnsCntsEnsmbl, \
outGRMFlSpkInsCntsEnsmbl,excldNmPrfx='%scounts.htsq'%dirSfx,rowsNm='', \
nmbrSpkIns=0,spltChrc='\t',sExcldGnNm=sExcldGnNm, \
pthCmmnGnNmToENSEMBLG=pthCmmnGnNmToENSEMBLG,nrmlzDpth=True, \
excldERCCfrmCnts=True)

