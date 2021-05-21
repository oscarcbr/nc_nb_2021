#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  smmrzQC.py
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
This script summarizes various quality features calculated in the previous steps
"""

import zipfile
import argparse,os,sys

from numpy import array,object,zeros,random
from singlecell.formats import rtrn_ArCllNms_ArVlsIntrst,rtrn_ArVlsIntrst_FastQC, \
rtrn_ArPrcntl_ArCovrg,rtrn_Stst,rtrn_ArBtype_ArTtlCnts
from singlecell.plots import pltHzntlBarPlotCmltv,pltPointPlot
from scipy.stats import pearsonr



#################################
#         Parse inputs          #                
#################################
parser = argparse.ArgumentParser()
 
######################
#  Input parameters  #
######################
parser.add_argument('-s','--smpl',help='Input sample',default=None)
parser.add_argument('-p','--dirSfx',help='Files sufix',default=None)
parser.add_argument('-m','--mapStarFldr',help='Folder with map STAR results',default=None)
parser.add_argument('-q','--qcRdsFldr',help='Folder with map QC results',default=None)

args = parser.parse_args()

######################
######################
######################
#  Input parameters  #
######################
######################
######################
smpl = args.smpl
dirSfx = args.dirSfx
mapStarFldr = args.mapStarFldr
qcRdsFldr = args.qcRdsFldr

lSmpls = [smpl]


print('Log info:')
print('\tsmpl -> ',smpl)
print('\tdirSfx -> ',dirSfx)
print('\tqcRdsFldr -> ',qcRdsFldr)
print('\tmapStarFldr -> ',mapStarFldr)

assert os.path.exists(qcRdsFldr) and os.path.exists(mapStarFldr)

lCllNms=['A1','A19','A6','B15','B24','C11','C20','C8','D17','D4','E13', \
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
'O21','O9','P18','P5']

#################
#################
#################
###  Execute  ###
#################
#################
#################

def smmrzQCStats(lSmpls,lCllNms,dirSfx,mapStarFldr,qcRdsFldr,pltHzntlBarPlot=True):
	"""
	Method to summarize and plot the QC results.
	"""
	#
	lCllNms.sort()
	dNmsInFlStstcs={'READ_PAIR_OK':'Number_reads_uniquely_mapped', \
	'TOTAL_READ_PAIRS':'Number_total_reads_HQ', \
	'DROPPED_NOT_UNIQUE_ALIGNMENT':'Number_reads_multimapped', \
	'minObservedReadLength':'minimum_read_length', \
	'maxObservedReadLength':'maximum_read_length','PREALIGNMENT_READ_CT': \
	'Number_total_reads','ReadPairs_AmbigGene': \
	'Number_reads_mapped_ambigously','ReadPairs_UniqueGene': \
	'Number_reads_mapped_unambigously','ReadPairs_UniqueGene_CDS': \
	'Number_reads_mapped_unambigously_CDSs','ReadPairs_UniqueGene_UTR': \
	'Number_reads_mapped_unambigously_UTRs','ReadPairs_NoGene': \
	'Number_reads_not_mapped_to_exons','ReadPairs_NoGene_Intron': \
	'Number_reads_mapped_to_introns','ReadPairs_NoGene_OneKbFromGene': \
	'Number_reads_mapped_to_1KB_from_gene','ReadPairs_NoGene_TenKbFromGene': \
	'Number_reads_mapped_to_10KB_from_gene', \
	'ReadPairs_NoGene_MiddleOfNowhere': \
	'Number_reads_mapped_to_more_10KB_from_gene','Genes_WithNonzeroCounts': \
	'Number_genes_detected','AVG_GC':'Average_GC_content','deletionLoci': \
	'Number_deletion_loci','insertionLoci':'Number_insertion_loci', \
	'deletionEventCt':'Number_times_deletions_in_a_read','insertionEventCt': \
	'Number_times_insertions_in_a_read','highCoverageDeletionLoci': \
	'Number_high-coverage_deletion_loci','highCoverageInsertionLoci': \
	'Number_high-coverage_insertion_loci','PAIR_CONTAINS_INDEL': \
	'Number_reads_mapped_unambigously_with_INDELS', \
	'HAS_REF_BASE_SWAP_R1':'Number_reads_mapped_unambigously_with_substitutions'}
	sBtypNms=set(['3prime_overlapping_ncRNA','IG_C_gene','IG_D_gene', \
	'IG_J_gene','IG_V_gene','Mt_rRNA','Mt_tRNA','TEC','TR_C_gene', \
	'TR_D_gene','TR_J_gene','TR_V_gene','TR_V_pseudogene','UNK','antisense', \
	'bidirectional_promoter_lncRNA','lincRNA','macro_lncRNA','miRNA', \
	'misc_RNA','non_coding','non_stop_decay','nonsense_mediated_decay', \
	'polymorphic_pseudogene','processed_pseudogene','processed_transcript', \
	'protein_coding','rRNA','retained_intron','ribozyme','sRNA','scRNA', \
	'scaRNA','sense_intronic','sense_overlapping','snRNA','snoRNA', \
	'transcribed_unitary_pseudogene','unk','unprocessed_pseudogene', \
	'vaultRNA'])#replace unk to Spike-In or other
	addStst=set(['eveness_of_coverage_(R2)'])
	dNmsInFlReadsQC={'Basic Statistics':'FastQC_Basic_statistics', \
	'Per base sequence quality':'FastQC_Per_base_sequence_quality', \
	'Per sequence quality scores':'FastQC_Per_sequence_quality_scores', \
	'Per base sequence content':'FastQC_Per_base_sequence_content', \
	'Per sequence GC content':'FastQC_Per_sequence_GC_content', \
	'Per base N content':'FastQC_Per_base_N_content', \
	'Sequence Length Distribution':'FastQC_Sequence_length_distribution', \
	'Sequence Duplication Levels':'FastQC_Sequence_duplication_levels', \
	'Overrepresented sequences':'FastQC_Overrepresented_sequences', \
	'Adapter Content':'FastQC_Adapter_content','Kmer Content': \
	'FastQC_Kmer_content'}
	#Sorted Keys
	dMrgDctnrs = dict([(k,dNmsInFlStstcs[k]) for k in dNmsInFlStstcs.keys()])
	dMrgDctnrs.update(dict([(k,'RNA_type_%s'%k) for k in sBtypNms]))
	dMrgDctnrs['unk']='RNA_type_Spike-in'
	dMrgDctnrs.update(dict([(k,k) for k in addStst]))
	dMrgDctnrs.update(dNmsInFlReadsQC)
	srtdKeys = ['PREALIGNMENT_READ_CT','TOTAL_READ_PAIRS','READ_PAIR_OK', \
	'DROPPED_NOT_UNIQUE_ALIGNMENT','ReadPairs_AmbigGene','ReadPairs_NoGene', \
	'ReadPairs_NoGene_Intron','ReadPairs_NoGene_OneKbFromGene', \
	'ReadPairs_NoGene_TenKbFromGene','ReadPairs_NoGene_MiddleOfNowhere', \
	'ReadPairs_UniqueGene','ReadPairs_UniqueGene_CDS', \
	'ReadPairs_UniqueGene_UTR','AVG_GC','minObservedReadLength', \
	'maxObservedReadLength','HAS_REF_BASE_SWAP_R1','PAIR_CONTAINS_INDEL', \
	'highCoverageDeletionLoci','highCoverageInsertionLoci','insertionLoci', \
	'insertionEventCt','deletionLoci','deletionEventCt', \
	'Genes_WithNonzeroCounts','eveness_of_coverage_(R2)','Basic Statistics',  \
	'Per base sequence quality','Per sequence quality scores', \
	'Per base sequence content','Per sequence GC content', \
	'Per base N content','Sequence Length Distribution', \
	'Sequence Duplication Levels','Overrepresented sequences', \
	'Adapter Content','Kmer Content','3prime_overlapping_ncRNA','IG_C_gene', \
	'IG_D_gene','IG_J_gene','IG_V_gene','Mt_rRNA','Mt_tRNA','TEC','TR_C_gene', \
	'TR_D_gene','TR_J_gene','TR_V_gene','TR_V_pseudogene','UNK','antisense', \
	'bidirectional_promoter_lncRNA','lincRNA','macro_lncRNA','miRNA', \
	'misc_RNA','non_coding','non_stop_decay','nonsense_mediated_decay', \
	'polymorphic_pseudogene','processed_pseudogene','processed_transcript', \
	'protein_coding','rRNA','retained_intron','ribozyme','sRNA','scRNA', \
	'scaRNA','sense_intronic','sense_overlapping','snRNA','snoRNA', \
	'transcribed_unitary_pseudogene','unk','unprocessed_pseudogene', \
	'vaultRNA']
	#Make results		
	for smpl in lSmpls:
		outFl = os.path.join(mapStarFldr,'%s%sstats.txt'%(smpl,dirSfx))
		if not os.path.exists(outFl):
		#~ if True:
			ooutFl = open(outFl,'w')
			ooutFl.write('Cell\t%s\n'%'\t'.join([dMrgDctnrs[k] for k in \
			srtdKeys]))
			genDir= \
			os.path.join(mapStarFldr,'%s%sd'%(smpl,dirSfx))
			for cllNm in lCllNms:
				fastaQCDir= os.path.join(qcRdsFldr,'%s/rawdata/%s'% \
				(smpl,cllNm))
				dCllAttrbtsStst=dict([(k,'nan') for k in srtdKeys])
				qortsFldr = os.path.join(genDir,'%s.qorts.d'%cllNm)
				flQCGnExprsnnBodyCvrg=os.path.join(qortsFldr, \
				'QC.geneBodyCoverage.by.expression.level.txt.gz')
				flQCSmmry=os.path.join(qortsFldr,'QC.summary.txt')
				#Return attributes for QoRTs
				dAttrbtsVal = rtrn_Stst(flQCSmmry,set(srtdKeys))
				dCllAttrbtsStst.update(dAttrbtsVal)	
				ar_Prcntl,ar_Covrg=rtrn_ArPrcntl_ArCovrg \
				(flQCGnExprsnnBodyCvrg)
				R,twoTlP = pearsonr(ar_Prcntl,ar_Covrg)
				dCllAttrbtsStst['eveness_of_coverage_(R2)']=str(R)
				#Return biotype results
				flQCBioType=os.path.join(qortsFldr, \
				'QC.biotypeCounts.txt.gz')
				ar_Btype,ar_Rslts=rtrn_ArBtype_ArTtlCnts(flQCBioType)
				dBtypTtl = dict([(ar_Btype[p],ar_Rslts[p]) for p in \
				xrange(len(ar_Btype))])
				dCllAttrbtsStst.update(dBtypTtl)
				#Return FastQC stats
				fastqcZip=os.path.join(fastaQCDir, \
				'%s_R1_trimmed_fastqc.zip'%cllNm)
				cmprssdFl=zipfile.ZipFile(fastqcZip)
				dFstQCtstRslts=dict([(v.split('\t')[1],v.split('\t')[0]) \
				for v in cmprssdFl.read('%s_R1_trimmed_fastqc/summary.txt'%cllNm). \
				splitlines() if v.strip()])
				dCllAttrbtsStst.update(dFstQCtstRslts)
				#Write results
				ooutFl.write('%s\t%s\n'%(cllNm,'\t'.join([dCllAttrbtsStst[k] \
				for k in srtdKeys])))
			#
			ooutFl.close()	
	"""
	Next I will make plots to summarize the quality values for each of the
	cells.
	"""
	if pltHzntlBarPlot:
		lReadsStatsIIn = ['Cell','Number_total_reads','Number_total_reads_HQ', \
		'Number_reads_uniquely_mapped','Number_reads_multimapped']
		lReadsStatsIOut = ['Cell','Number_reads_uniquely_mapped', \
		'Number_reads_multimapped','Number_reads_low_quality']
		lReadsStatsII = ['Cell','Number_reads_mapped_ambigously', \
		'Number_reads_mapped_to_more_10KB_from_gene', \
		'Number_reads_mapped_to_10KB_from_gene', \
		'Number_reads_mapped_to_1KB_from_gene', \
		'Number_reads_mapped_to_introns', \
		'Number_reads_mapped_unambigously_UTRs', \
		'Number_reads_mapped_unambigously_CDSs']
		lRNAStats=['Cell','RNA_type_3prime_overlapping_ncRNA','RNA_type_IG_C_gene', \
		'RNA_type_IG_D_gene','RNA_type_IG_J_gene','RNA_type_IG_V_gene', \
		'RNA_type_Mt_rRNA','RNA_type_Mt_tRNA','RNA_type_TEC','RNA_type_TR_C_gene', \
		'RNA_type_TR_D_gene','RNA_type_TR_J_gene','RNA_type_TR_V_gene', \
		'RNA_type_TR_V_pseudogene','RNA_type_UNK','RNA_type_antisense', \
		'RNA_type_bidirectional_promoter_lncRNA','RNA_type_lincRNA', \
		'RNA_type_macro_lncRNA','RNA_type_miRNA','RNA_type_misc_RNA', \
		'RNA_type_non_coding','RNA_type_non_stop_decay', \
		'RNA_type_nonsense_mediated_decay','RNA_type_polymorphic_pseudogene', \
		'RNA_type_processed_pseudogene','RNA_type_processed_transcript', \
		'RNA_type_protein_coding','RNA_type_rRNA','RNA_type_retained_intron', \
		'RNA_type_ribozyme','RNA_type_sRNA','RNA_type_scRNA','RNA_type_scaRNA', \
		'RNA_type_sense_intronic','RNA_type_sense_overlapping','RNA_type_snRNA', \
		'RNA_type_snoRNA','RNA_type_transcribed_unitary_pseudogene', \
		'RNA_type_Spike-in','RNA_type_unprocessed_pseudogene','RNA_type_vaultRNA']
		lFastaQCStats=['Cell','FastQC_Basic_statistics','FastQC_Per_base_sequence_quality', \
		'FastQC_Per_sequence_quality_scores','FastQC_Per_base_sequence_content', \
		'FastQC_Per_sequence_GC_content','FastQC_Per_base_N_content', \
		'FastQC_Sequence_length_distribution','FastQC_Sequence_duplication_levels', \
		'FastQC_Overrepresented_sequences','FastQC_Adapter_content', \
		'FastQC_Kmer_content']
		lFastaQCStatsOut=['Cell','Number_failed_tests','Number_warning_tests', \
		'Number_passed_tests']
		lNmbrGnsDtctd=['Cell','Number_genes_detected']
		lXYvrblsGrps=[['Average_GC_content','eveness_of_coverage_(R2)'], \
		['Number_times_deletions_in_a_read','Number_times_insertions_in_a_read'], \
		['Number_high-coverage_deletion_loci','Number_high-coverage_insertion_loci'], \
		['Number_reads_mapped_unambigously_with_INDELS', \
		'Number_reads_mapped_unambigously_with_substitutions']]
		lXYlbls=['Value','Count','Count','Count']
		#Make arrays to plot bar and XY plots
		for smpl in lSmpls:
			smmryFl = \
			os.path.join(mapStarFldr,'%s%sstats.txt'%(smpl,dirSfx))
			#Barplots Reads Stats I 
			outPltReadsStatsI = \
			os.path.join(mapStarFldr,'%s%sstats.readStatsI.svg'%(smpl,dirSfx))
			ar_vlsStatsIIn=rtrn_ArCllNms_ArVlsIntrst(smmryFl,lReadsStatsIIn)
			ar_vlsStatsIOut=zeros((ar_vlsStatsIIn.shape[0],len(lReadsStatsIOut)), \
			dtype=object)
			ar_vlsStatsIOut[:,0]=ar_vlsStatsIIn[:,0]
			ar_vlsStatsIOut[:,3]=ar_vlsStatsIIn[:,1]-ar_vlsStatsIIn[:,2]
			ar_vlsStatsIOut[:,1]=ar_vlsStatsIIn[:,3]
			ar_vlsStatsIOut[:,2]=ar_vlsStatsIIn[:,4]
			#Sort by Number_reads_uniquely_mapped
			lNmbrUnqlyMpdNCllNm = [(ar_vlsStatsIOut[pos,1], \
			ar_vlsStatsIOut[pos,0]) for pos in xrange(ar_vlsStatsIIn.shape[0])]
			ordrdClls=[v[1] for v in lNmbrUnqlyMpdNCllNm]
			lNmbrUnqlyMpdNCllNm.sort()
			lNmbrUnqlyMpdNCllNm.reverse()
			print lNmbrUnqlyMpdNCllNm,ar_vlsStatsIOut.shape
			ar_cllNmsOrdrd = array([ordrdClls.index(v[1]) for v in \
			lNmbrUnqlyMpdNCllNm])
			ordrdClls=array(ordrdClls)
			#Finish to plot
			if not os.path.exists(outPltReadsStatsI):
				ar_vlsStatsIOut=ar_vlsStatsIOut[ar_cllNmsOrdrd,:]
				pltHzntlBarPlotCmltv(array(lReadsStatsIOut),ar_vlsStatsIOut, \
				outPltReadsStatsI)
				del(ar_vlsStatsIOut)
			del(ar_vlsStatsIIn)
			del(lNmbrUnqlyMpdNCllNm)
			#Barplots Reads Stats II
			outPltReadsStatsII = \
			os.path.join(mapStarFldr,'%s%sstats.readStatsII.svg'%(smpl,dirSfx))
			if not os.path.exists(outPltReadsStatsII):
				ar_vlsStatsII=rtrn_ArCllNms_ArVlsIntrst(smmryFl,lReadsStatsII)
				assert (ordrdClls==ar_vlsStatsII[:,0]).all()
				ar_vlsStatsII=ar_vlsStatsII[ar_cllNmsOrdrd,:]
				pltHzntlBarPlotCmltv(array(lReadsStatsII),ar_vlsStatsII, \
				outPltReadsStatsII)
				del(ar_vlsStatsII)
			#Barplots RNA Stats
			outPltRNAStats = \
			os.path.join(mapStarFldr,'%s%sstats.RNAStats.svg'%(smpl,dirSfx))
			if not os.path.exists(outPltRNAStats):
				ar_vlsRNAStats=rtrn_ArCllNms_ArVlsIntrst(smmryFl,lRNAStats)
				assert (ordrdClls==ar_vlsRNAStats[:,0]).all()
				ar_vlsRNAStats=ar_vlsRNAStats[ar_cllNmsOrdrd,:]
				pltHzntlBarPlotCmltv(array(lRNAStats),ar_vlsRNAStats, \
				outPltRNAStats)
				del(ar_vlsRNAStats)
			#Barplots FastQC Stats
			outPltFastQCtats = \
			os.path.join(mapStarFldr,'%s%sstats.fstQCStats.svg'%(smpl,dirSfx))
			if not os.path.exists(outPltFastQCtats):
				ar_vlsFastaQCStats=rtrn_ArVlsIntrst_FastQC(smmryFl,lFastaQCStats)
				assert (ordrdClls==ar_vlsFastaQCStats[:,0]).all()
				ar_vlsFastaQCStats=ar_vlsFastaQCStats[ar_cllNmsOrdrd,:]
				pltHzntlBarPlotCmltv(array(lFastaQCStatsOut),ar_vlsFastaQCStats, \
				outPltFastQCtats)
				del(ar_vlsFastaQCStats)
			#Number of genes detected
			outPltGnsDtctd = \
			os.path.join(mapStarFldr,'%s%sstats.nmbrGnsDtctd.svg'%(smpl,dirSfx))
			if not os.path.exists(outPltGnsDtctd):
				ar_vlsGnsDtctd=rtrn_ArCllNms_ArVlsIntrst(smmryFl,lNmbrGnsDtctd)
				assert (ordrdClls==ar_vlsGnsDtctd[:,0]).all()
				ar_vlsGnsDtctd=ar_vlsGnsDtctd[ar_cllNmsOrdrd,:]
				pltHzntlBarPlotCmltv(array(lNmbrGnsDtctd),ar_vlsGnsDtctd, \
				outPltGnsDtctd)
				del(ar_vlsGnsDtctd)
			#XY plot
			for xyVrblGrpPos in xrange(len(lXYvrblsGrps)):
				outPltXYVrblStat = \
				os.path.join(mapStarFldr,'%s.%s%sstats.XYStats.svg'%(smpl,xyVrblGrpPos,dirSfx))
				if not os.path.exists(outPltXYVrblStat):
					xyVrblGrp=[v for v in lXYvrblsGrps[xyVrblGrpPos]]
					xyVrblGrp.insert(0,'Cell')
					yLabelXYVar=lXYlbls[xyVrblGrpPos]
					ar_vlsXYVarStats=rtrn_ArCllNms_ArVlsIntrst(smmryFl,xyVrblGrp)
					assert (ordrdClls==ar_vlsXYVarStats[:,0]).all()
					ar_vlsXYVarStats=ar_vlsXYVarStats[ar_cllNmsOrdrd,:]
					pltPointPlot(array(xyVrblGrp),ar_vlsXYVarStats,outPltXYVrblStat, \
					yLabel=yLabelXYVar)


smmrzQCStats(lSmpls,lCllNms,dirSfx,mapStarFldr,qcRdsFldr)
