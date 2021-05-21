#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hg38_annotation.py
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
The idea is to annotate the hg38 genome using the latest version of 
GENECODE database
"""

from numpy import array
from singlecell.annotation import addBioType,addGnNm,appndEGFP,appndERCC, \
mkDctnryRefSeqGnNm,mkGTFfromUCSCAnntFrmt

import gzip
import os

"""
The goal is to manage specific annotation formats for hg38.

The idea is to generate expanded annotations with the upgraded version 
of Gencode (v28) and included the annotations for iCre.
"""


#########
#Analysis
#########
make_gtf_from_genecodeV28_anntn = True #If true will make a gtf annotation from a UCSC genecodeV28 table: 1) suitable for cufflinks, 2) excluding alternative chromosomes, 3) calling gene names in X and Y chromosomes differentially, 4) excluding intron annotation
mrgChrms = True#If True will merge the sequence from each chromosomes into a single reference
addBioType_genecodeV28_anntn = True#If True will add biotype annotations for each transcript

dataFldr = os.getcwd()
outFldr = os.getcwd()

"""
The comprenhensive set of genecode V28 annotated genes was obtained from 
the UCSC broswer tables on 25.02.2019 and saved to the path:
genecodeV28Comp.txt.
Further gtf files were created for make explicitly aware of introns and 
exons in CDS and UTRs. Is important that the gtf file to consider 
transcripts with different chromosome location and XY, with the next 
script.
"""
if make_gtf_from_genecodeV28_anntn:
	anntnID='hg38_ERCC_eGFP'
	method = 'prseBedToStndrdGTF'	
	UCSCAnntFrmt = os.path.join(dataFldr,'genecodeV28Comp.txt')
	gflOut = os.path.join(outFldr, \
	'hg38.genecodeV28Comp.ERCCeGFP.cfflinks.noIntrns.gtf')
	mkGTFfromUCSCAnntFrmt(UCSCAnntFrmt,gflOut,anntnID,method, \
	stblCfflinks=True,excldIntrns=True,excldAltChrms=True,cnsdrYX=True)
	dRefSeqGnNm = mkDctnryRefSeqGnNm(UCSCAnntFrmt)
	gflOutGnNm = os.path.join(outFldr,\
	'hg38.genecodeV28Comp.ERCCeGFP.cfflinks.noIntrns.gnNms.gtf')
	addGnNm(gflOut,dRefSeqGnNm,gflOutGnNm,anntHasDot=True, \
	stblCfflinks=True)
	ERCCgtf = os.path.join(dataFldr,'ERCC92.gtf')
	appndERCC(gflOutGnNm,ERCCgtf,anntnID)
	appndEGFP(gflOutGnNm,anntnID)


"""
I will merge the chromosomes in the following order for later analysis.
The fasta sequences for the chromosomes were obtained from 
http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr*.fa. 
Note that only the chromosomes in the following "lChrms" list were 
obtained. 
The ERCC sequences and annotation files were obtained from the vendor 
website.
"""
if mrgChrms:
	lChrms=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9', \
	'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18', \
	'chr19','chr20','chr21','chr22','chrM','chrX','chrY','chrEGFP','ERCC-00002', \
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
	inFldr=os.path.join(dataFldr,'chroms')
	outFl = os.path.join(dataFldr,'mrgdChroms.d','hg38.fa')
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
annotation of "gene type" downloaded from Ensembl for hg38 with the 
following fields: "Transcript stable ID" and "Transcript type". The 
values are going to be included in the GTF file using the field 
"gene_biotype". The file downloaded from Ensembl was saved in 
hg38.genecodeV28Comp.geneType.tsv. 
The reference was hg38/GRCh38.p12 obtained in 25.02.2019.
"""
if addBioType_genecodeV28_anntn:
	anntFl = os.path.join(outFldr,"hg38.genecodeV28Comp.geneType.tsv")
	gflInGnNm = os.path.join(outFldr,\
	'hg38.genecodeV28Comp.ERCCeGFP.cfflinks.noIntrns.gnNms.gtf')
	gflOutGnNmAnntd = os.path.join(outFldr,\
	'hg38.genecodeV28Comp.ERCCeGFP.cfflinks.noIntrns.gnNms.biotyp.gtf')
	dRefSeqGnNm = mkDctnryRefSeqGnNm(anntFl,twoClmn=True)
	addBioType(gflInGnNm,dRefSeqGnNm,gflOutGnNmAnntd)

