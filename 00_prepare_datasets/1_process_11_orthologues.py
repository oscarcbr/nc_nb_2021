#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  process_11_orthologues.py
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
This code aims to return a list of 1:1 orthologues between human and mouse
"""

import os

######################
######################
#####  Switches  #####
######################
######################
crrctOrthlgsFl = True# If True will correct the orthologous file by converting GRCh38_p12/GRCm38_p6 with genecodeV28Comp/genecodeV18Comp respectively

inFl_hg38p12 = 'GRCh38_p12_ENSMBLgnNm.tsv'
inFl_mm10p6 = 'GRCm38_p6_ENSEMBLgnNmEdtd.txt'

sGns_hg38p12 = set([l.splitlines()[0].split('\t')[1] for l in \
open(inFl_hg38p12,'r') if l.strip() and l.split('\t')[0]!='Gene name' \
and l[0]!='#'])
sGns_mm10p6 = set([l.splitlines()[0].split('\t')[1] for l in \
open(inFl_mm10p6,'r') if l.strip() and l.split('\t')[0]!='Gene name' \
and l[0]!='#'])
dSppGns = {'Human':sGns_hg38p12,'Mouse':sGns_mm10p6}

######################
######################
if crrctOrthlgsFl:
	#Ortholgous files
	inFlOrthlg = 'hg38_mm10_orthlgs_1to1.txt.ori'
	outFlOrthlg = 'hg38_mm10_orthlgs_1to1.txt'
	#p12/p6 versions to genecodeV28Comp/genecodeV18Comp files
	outFlCnvrt_mm10 = 'GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv'
	outFlCnvrt_hg38 = 'GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv'
	#Get ENSEMBL code for the p12/p6 annotations
	dhg38p12_gnCdV28 = dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(outFlCnvrt_hg38).read().splitlines() if l.strip() and l.split('\t')[0]!='Gene name' \
	and l[0]!='#'])
	dmm10p6_gnCdV18 =  dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(outFlCnvrt_mm10).read().splitlines() if l.strip() and l.split('\t')[0]!='Gene name'\
	and l[0]!='#'])
	sSelf_mm10 = set(dmm10p6_gnCdV18.keys())
	sSelf_hg38 = set(dhg38p12_gnCdV28.keys())
	#Get orthologues
	#Next line was added to bipass the Mycn+Mycs paralogue status. 
	#Mycn is a major player in neuroblastoma.
	lOoutFlOrthlg = ['%s\n'%'\t'.join(['MYCN\tMycn\t100\t100.00\t1\tortholog_one2many'])]
	#
	for l in open(inFlOrthlg,'r'):
		if l.strip():
			if l[0]=='#':
				hdr = l
			else:
				l=l.splitlines()[0].split('\t')
				hgGnNm,mmGnNm=l[0],l[1]
				if hgGnNm in sSelf_hg38 and mmGnNm in sSelf_mm10:
					l[0]=dhg38p12_gnCdV28[hgGnNm]
					l[1]=dmm10p6_gnCdV18[mmGnNm]
					lOoutFlOrthlg.append('%s\n'%'\t'.join(l))
	#Write file
	ooutFlOrthlg = open(outFlOrthlg,'w')
	ooutFlOrthlg.write(hdr)
	ooutFlOrthlg.write(''.join(sorted(lOoutFlOrthlg)))
	ooutFlOrthlg.close()
