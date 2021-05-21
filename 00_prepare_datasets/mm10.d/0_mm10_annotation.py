#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  mm10_annotation.py
#  
#  Copyright 2018 Oscar C. Bedoya Reina <oscar@ki.se>
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
The idea is to annotate the mm10 genome using the latest version of 
GENECODE database
"""

from numpy import array
import gzip
import os

from numpy import array
from singlecell.annotation import addBioType,addGnNm,appndEGFP,appndERCC, \
appndICRE,mkDctnryRefSeqGnNm,mkGTFfromUCSCAnntFrmt

"""
The goal is to manage specific annotation formats for mm10.

The idea is to generate expanded annotations with the upgraded version 
of Gencode (v18) and included the annotations for iCre.
"""


#########
#Analysis
#########
make_gtf_from_genecodev18_anntn = True #If true will make a gtf annotation from a UCSC genecodev18 table: 1) suitable for cufflinks, 2) excluding alternative chromosomes, 3) calling gene names in X and Y chromosomes differentially, 4) excluding intron annotation
mrgChrms = True#If True will merge the sequence from each chromosomes into a single reference
addBioType_genecodev18_anntn = True#If True will add biotype annotations for each transcript

dataFldr = 'mm10.d'
outFldr = 'mm10.d'


"""
iCre sequence was included for interest of the research into the 
annotations.
I downloaded the iCre sequence from the GenBank reference AY056050.1.
I confirmed that the sequence was the same as the one reported by
Shimshek et al. 2002. and saved in the file chrICRE.fa
"""
"""
>chrICRE
AAGCTTGTCCACCATGGTGCCCAAGAAGAAGAGGAAAGTCTCCAACCTGC
TGACTGTGCACCAAAACCTGCCTGCCCTCCCTGTGGATGCCACCTCTGAT
GAAGTCAGGAAGAACCTGATGGACATGTTCAGGGACAGGCAGGCCTTCTC
TGAACACACCTGGAAGATGCTCCTGTCTGTGTGCAGATCCTGGGCTGCCT
GGTGCAAGCTGAACAACAGGAAATGGTTCCCTGCTGAACCTGAGGATGTG
AGGGACTACCTCCTGTACCTGCAAGCCAGAGGCCTGGCTGTGAAGACCAT
CCAACAGCACCTGGGCCAGCTCAACATGCTGCACAGGAGATCTGGCCTGC
CTCGCCCTTCTGACTCCAATGCTGTGTCCCTGGTGATGAGGAGAATCAGA
AAGGAGAATGTGGATGCTGGGGAGAGAGCCAAGCAGGCCCTGGCCTTTGA
ACGCACTGACTTTGACCAAGTCAGATCCCTGATGGAGAACTCTGACAGAT
GCCAGGACATCAGGAACCTGGCCTTCCTGGGCATTGCCTACAACACCCTG
CTGCGCATTGCCGAAATTGCCAGAATCAGAGTGAAGGACATCTCCCGCAC
CGATGGTGGGAGAATGCTGATCCACATTGGCAGGACCAAGACCCTGGTGT
CCACAGCTGGTGTGGAGAAGGCCCTGTCCCTGGGGGTTACCAAGCTGGTG
GAGAGATGGATCTCTGTGTCTGGTGTGGCTGATGACCCCAACAACTACCT
GTTCTGCCGGGTCAGAAAGAATGGTGTGGCTGCCCCTTCTGCCACCTCCC
AACTGTCCACCCGGGCCCTGGAAGGGATCTTTGAGGCCACCCACCGCCTG
ATCTATGGTGCCAAGGATGACTCTGGGCAGAGATACCTGGCCTGGTCTGG
CCACTCTGCCAGAGTGGGTGCTGCCAGGGACATGGCCAGGGCTGGTGTGT
CCATCCCTGAAATCATGCAGGCTGGTGGCTGGACCAATGTGAACATAGTG
ATGAACTACATCAGAAACCTGGACTCTGAGACTGGGGCCATGGTGAGGCT
GCTCGAGGATGGGGACTGAAACTGAGTCGA
"""

"""
The comprenhensive set of genecode v18 annotated genes was obtained from 
the UCSC broswer tables on 14.01.2019 and saved to the path:
genecodeV18Comp.txt.
Further gtf files  were created for make explicitly aware of introns 
and exons in CDS and UTRs. Is important that the gtf file to consider 
transcripts with different chromosome location and XY, with the next 
script.
"""
if make_gtf_from_genecodev18_anntn:
	anntnID='mm10_ERCC_eGFP_iCre'
	method = 'prseBedToStndrdGTF'	
	UCSCAnntFrmt = os.path.join(outFldr,'genecodeV18Comp.txt')
	gflOut = os.path.join(outFldr, \
	'mm10.genecodeV18Comp.ERCCeGFPiCre.cfflinks.gtf')
	mkGTFfromUCSCAnntFrmt(UCSCAnntFrmt,gflOut,anntnID,method, \
	stblCfflinks=True,excldIntrns=True,excldAltChrms=True,cnsdrYX=True)
	dRefSeqGnNm = mkDctnryRefSeqGnNm(UCSCAnntFrmt)
	gflOutGnNm = os.path.join(outFldr,\
	'mm10.genecodeV18Comp.ERCCeGFPiCre.cfflinks.gnNms.gtf')
	addGnNm(gflOut,dRefSeqGnNm,gflOutGnNm,anntHasDot=True, \
	stblCfflinks=True)
	ERCCgtf = os.path.join(dataFldr,'ERCC92.gtf')
	appndERCC(gflOutGnNm,ERCCgtf,anntnID)
	appndEGFP(gflOutGnNm,anntnID)
	appndICRE(gflOutGnNm,anntnID)



"""
I will merge the chromosomes in the following order for later analysis.
The fasta sequences for the chromosomes were obtained from 
http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr*.fa on
14.01.2019.
Note that only the chromosomes in the following "lChrms" list were 
obtained. 
The ERCC sequences and annotation files were obtained from the vendor 
website.
"""
if mrgChrms:
	lChrms=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9', \
	'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18', \
	'chr19','chrM','chrX','chrY','chrEGFP','chrICRE','ERCC-00002', \
	'ERCC-00003','ERCC-00004','ERCC-00009','ERCC-00012','ERCC-00013', \
	'ERCC-00014','ERCC-00016','ERCC-00017','ERCC-00019','ERCC-00022', \
	'ERCC-00024','ERCC-00025','ERCC-00028','ERCC-00031','ERCC-00033', \
	'ERCC-00034','ERCC-00035','ERCC-00039','ERCC-00040','ERCC-00041', \
	'ERCC-00042','ERCC-00043','ERCC-00044','ERCC-00046','ERCC-00048', \
	'ERCC-00051','ERCC-00053','ERCC-00054','ERCC-00057','ERCC-00058', \
	'ERCC-00059','ERCC-00060','ERCC-00061','ERCC-00062','ERCC-00067', \
	'ERCC-00069','ERCC-00071','ERCC-00073','ERCC-00074','ERCC-00075', \
	'ERCC-00076','ERCC-00077','ERCC-00078','ERCC-00079','ERCC-00081', \
	'ERCC-00083','ERCC-00084','ERCC-00085','ERCC-00086','ERCC-00092', \
	'ERCC-00095','ERCC-00096','ERCC-00097','ERCC-00098','ERCC-00099', \
	'ERCC-00104','ERCC-00108','ERCC-00109','ERCC-00111','ERCC-00112', \
	'ERCC-00113','ERCC-00116','ERCC-00117','ERCC-00120','ERCC-00123', \
	'ERCC-00126','ERCC-00130','ERCC-00131','ERCC-00134','ERCC-00136', \
	'ERCC-00137','ERCC-00138','ERCC-00142','ERCC-00143','ERCC-00144', \
	'ERCC-00145','ERCC-00147','ERCC-00148','ERCC-00150','ERCC-00154', \
	'ERCC-00156','ERCC-00157','ERCC-00158','ERCC-00160','ERCC-00162', \
	'ERCC-00163','ERCC-00164','ERCC-00165','ERCC-00168','ERCC-00170', \
	'ERCC-00171']
	inFldr='mm10.d/chroms'
	outFl = 'mm10.d/mrgdChroms.d/mm10.fa'
	ooutFl = open(outFl,'w')
	for inFl in lChrms:
		inFst = '%s.fa'%inFl
		opndInFst = open(os.path.join(inFldr,inFst)).read()
		if opndInFst[-1:]=='\n':
			ooutFl.write(opndInFst)
		else:
			ooutFl.write('%s\n'%opndInFst)

	ooutFl.close()

"""
In order to include a broader quality analysis by QoRTs I will use the
annotation of "gene type" downloaded from Ensembl for mm10 with the 
following fields: "Transcript stable ID" and "Transcript type". The 
values are going to be included in the GTF file using the field 
"gene_biotype". The file downloaded from Ensembl was saved in 
mm10.genecodeV18Comp.geneType.tsv.
The reference was mm10/GRCm38.p6 obtained in 14.01.2019.
"""
if addBioType_genecodev18_anntn:
	anntFl = os.path.join(outFldr,"mm10.genecodeV18Comp.geneType.tsv")
	gflInGnNm = os.path.join(outFldr,\
	'mm10.genecodeV18Comp.ERCCeGFPiCre.cfflinks.gnNms.gtf')
	gflOutGnNmAnntd = os.path.join(outFldr,\
	'mm10.genecodeV18Comp.ERCCeGFPiCre.cfflinks.gnNms.biotyp.gtf')
	dRefSeqGnNm = mkDctnryRefSeqGnNm(anntFl,twoClmn=True)
	addBioType(gflInGnNm,dRefSeqGnNm,gflOutGnNmAnntd)

