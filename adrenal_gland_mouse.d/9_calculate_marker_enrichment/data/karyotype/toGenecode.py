#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  toGenecode.py
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
Script to convert GRCm38 gene symbols to genecode
"""

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
sSelf_mm10 = set(dmm10p6_gnCdV18.values())
sSelf_hg38 = set(dhg38p12_gnCdV28.values())
#
dSppGnsCnvrt={'Human':dhg38p12_gnCdV28,'Mouse':dmm10p6_gnCdV18}
dSppSself = {'Human':sSelf_hg38,'Mouse':sSelf_mm10}
#Convert
for inFl,spp in (('karyotype_band.hg38.tsv.ori','Human'),('karyotype_band.mm10.tsv.ori','Mouse')):
	outFl = inFl.replace('.ori','')
	oOutFl = open(outFl,'w')
	dGns=dSppGnsCnvrt[spp]
	for l in open(inFl,'r'):
		if l.strip() and l[0]!='#' and l.split()[0]!='Gene':
			l=l.splitlines()[0].split('\t')
			gnNm = l[0]
			if dGns.has_key(gnNm):
				l[0]=dGns[gnNm]
				oOutFl.write('%s\n'%'\t'.join(l))
			else:
				oOutFl.write(l)
	oOutFl.close()
