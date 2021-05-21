#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  quantification.py
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
The idea is to quantify the expression of genes of interest
"""

import argparse,os,sys

from singlecell.quantification import wrpHtSeq
from singlecell.formats import htseqLFlsToGRM



#################################
#         Parse inputs          #                
#################################
parser = argparse.ArgumentParser()
 
######################
#  Input parameters  #
######################
parser.add_argument('-s','--smpl',help='Input sample',default=None)
parser.add_argument('-e','--exprssnFldr',help='Expression folder',default=None)
parser.add_argument('-m','--mapStarFldr',help='Folder with map STAR results',default=None)
parser.add_argument('-g','--gtfFile',help='GTF file',default=None)
parser.add_argument('-S','--species',help='Species',default=None)
parser.add_argument('-p','--dirSfx',help='sufix',default=None)

args = parser.parse_args()

######################
######################
######################
#  Input parameters  #
######################
######################
######################
smpl = args.smpl
exprssnFldr = args.exprssnFldr
mapStarFldr = args.mapStarFldr
gtfFile = args.gtfFile
species = args.species
dirSfx = args.dirSfx

smplNms = [smpl]


print('Log info:')
print('\tsmpl -> ',smpl)
print('\texprssnFldr -> ',exprssnFldr)
print('\tgtfFile -> ',gtfFile)
print('\tmapStarFldr -> ',mapStarFldr)
print('\tspecies -> ',species)
print('\tdirSfx -> ',dirSfx)

assert os.path.exists(exprssnFldr) and os.path.exists(mapStarFldr) and os.path.exists(gtfFile)
assert species in {'Human','Mouse'}

#The following dictionary will allow to quantify KIF1B isoforms
dSppdTplChrIntrvlGnNm = {'Mouse':{('chr4',((149233734,149238483),)):'KIF1bAlpha',('chr4',((149176319,149180810),(149181484,149181602),(149181809,149182001),(149182268,149182417),(149184309,149184430),(149186103,149186174),(149187591,149187830),(149187999,149188144),(149189723,149189784),(149191149,149191282),(149192573,149192687),(149198432,149198537),(149199258,149199342),(149202508,149202574),(149204186,149204294),(149206718,149206773),(149207849,149207967),(149210042,149210132),(149213296,149213458),(149213623,149213752),(149214079,149214164),(149214907,149215025),(149220542,149220790),(149222227,149222364),(149223230,149223408),(149225090,149225238),(149226263,149226356))):'KIF1bBeta'},'Human':{('chr1',((10303163,10304647),)):'KIF1bAlpha',('chr1',((10320043,10320136),(10321709,10321857),(10323884,10324062),(10324758,10324895),(10326111,10326359),(10334520,10334638),(10336657,10336742),(10337074,10337203),(10337371,10337533),(10339769,10339859),(10342050,10342168),(10343232,10343287),(10345845,10345953),(10347761,10347827),(10348649,10348733),(10352664,10352736),(10355370,10355396),(10360929,10361043),(10361692,10361825),(10363283,10363344),(10365100,10365245),(10365409,10365648),(10368467,10368538),(10371141,10371262),(10374316,10374465),(10374854,10375046),(10375255,10375373),(10376545,10376587),(10378339,10378402))):'KIF1bBeta'}}
dTplChrIntrvlGnNm = dSppdTplChrIntrvlGnNm[species]

#List of cell names
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

#################
#################
#################
###  Execute  ###
#################
#################
#################

for smpl in smplNms:
	algnFldr = os.path.join(mapStarFldr,'%s%sd'%(smpl,dirSfx))
	for cllNm in cllNms:
		BAMfl = os.path.join(algnFldr,'%s.bam'%cllNm)
		ouFldr = os.path.join(exprssnFldr,smpl)
		if not os.path.exists(ouFldr):
			os.system('mkdir %s'%ouFldr)
		outFl = os.path.join(ouFldr,'%s%scounts.htsq'%(cllNm,dirSfx))
		#quantify with HTSeq
		wrpHtSeq(BAMfl,gtfFile,outFl,splcAwr=False,splcAwrWExcptns=True, \
		geneID='gene',legend='##Results for sample %s and cell %s'%(smpl, \
		cllNm),dTplChrIntrvlGnNm=dTplChrIntrvlGnNm,uniquelyMap=False, \
		rnDbgd=True)
