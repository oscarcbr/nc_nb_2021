#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  process_furlan.py
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
Determine marker enrichments for independent populations
"""

from singlecell.between_samples import clcPrWsFDRforEchClstr, \
clcPrWsFDRforClstrAvrg
from numpy import array,mean,float32
from itertools import combinations
from string import upper

import os

"""
Analysis conducted to comparisons for counts and expression between 
clusters and samples, specifically for Furlan et al. 2018.
"""

##############
#Switches for all data
##############
prWsFDRForEchClstrFurlans=True#if True will calculate the genes that characterize each cluster. This is the intersection of all upregulated in each cluster vs. all others (specific signature) in Furlan's data
prWsFDRForEchClstrFurlansDpltn=True#if True will calculate the genes that characterize each cluster. This is the intersection of all downregulated in each cluster vs. all others in Furlan's data
prWsFDRForEchClstrFurlansAvrg=True#if True will calculate the genes that characterize each cluster. This is the intersection of all upregulated in each cluster vs. all others in Furlan's data (average)
prWsFDRForEchClstrFurlansDpltnAvrg=True#if True will calculate the genes that characterize each cluster. This is the intersection of all downregulated in each cluster vs. all others in Furlan's data (average)
complDBFurlans=True#If True will compile the characteristic genes of each cluster into a database of gene names and cluster/cell type
mrgDBs=True#If True will merge the results into single output files


########################################################
########################################################
########################################################
########################################################
########################################################
##############   Switches for all data    ##############
########################################################
########################################################
########################################################
########################################################
########################################################


#######################################
#Obtain results for Furlan's neural crest results
#######################################

#Input requirements
mnFldrNoIntrs = os.getcwd()
outFldr = os.getcwd()
#List of folders with results
l_rsltFldrs = [os.path.join(mnFldrNoIntrs,'E13.d'), \
os.path.join(mnFldrNoIntrs,'E13.KEql6.d'), \
os.path.join(mnFldrNoIntrs,'E12.d')]
l_KNum=[6,6,5] #list of number of clusters.


"""
Let's calculate the genes that characterize each cluster. This is the
intersection of all upregulated in each cluster vs. all others (e.g. 
specific gene signature) in Furlans' data .
"""
if prWsFDRForEchClstrFurlans:
	#Run analysis for the various folders
	print "Calculating cluster's specific signature in Furlans' data (upregulated)..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p],'dffExprssn.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		PAGODA_Rds = os.path.join(l_rsltFldrs[p],'app.final.rds')
		#list of clusters of interest
		if p==0:#added to exclude clusters 4 and 6 (gray and yellow)
			lPoptns = [x+1 for x in xrange(k_num) if x+1 not in {4,6}]
		else:
			lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR enrichment
		hcRds = ''
		clcPrWsFDRforEchClstr(PAGODA_Rds,hcRds,lPoptns, \
		enrchmntAnlysFldr,k_num=k_num,rplcFls=False,frmPAGODAappFl=True, \
		FDRtrhld=0.01)
	print 'Done!'
		

"""
Now the gene significantly upregulated for each cluster.
"""
if prWsFDRForEchClstrFurlansAvrg:
	#Run analysis for the various folders
	print "Calculating cluster's expression significance in Furlans' by average in each sample..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p],'dffExprssnAvrg.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		PAGODA_Rds = os.path.join(l_rsltFldrs[p],'app.final.rds')
		#list of clusters of interest
		if p==0:#added to exclude clusters 4 and 6 (gray and yellow)
			lPoptns = [x+1 for x in xrange(k_num) if x+1 not in {4,6}]
		else:
			lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR enrichment
		hcRds = ''
		clcPrWsFDRforClstrAvrg(PAGODA_Rds,hcRds,lPoptns, \
		enrchmntAnlysFldr,k_num=k_num,rplcFls=False,frmPAGODAappFl=True, \
		FDRtrhld=0.01)
	print 'Done!'


"""
Let's calculate the genes that characterize each cluster. This is the
intersection of all downregulated in each cluster vs. all others in Furlans'
data.
"""
if prWsFDRForEchClstrFurlansDpltn:
	#Run analysis for the various folders
	print "Calculating cluster's clusters uniquely downregulated in Furlans' data..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p],'dffExprssn.negative.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		PAGODA_Rds = os.path.join(l_rsltFldrs[p],'app.final.rds')
		#list of clusters of interest
		if p==0:#added to exclude clusters 4 and 6 (gray and yellow)
			lPoptns = [x+1 for x in xrange(k_num) if x+1 not in {4,6}]
		else:
			lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR enrichment
		hcRds = ''
		clcPrWsFDRforEchClstr(PAGODA_Rds,hcRds,lPoptns,enrchmntAnlysFldr, \
		k_num=k_num,rplcFls=False,frmPAGODAappFl=True,FDRtrhld=0.01, \
		grtr=False)
	print 'Done!'


"""
Now the gene significantly downregulated for each cluster.
"""
if prWsFDRForEchClstrFurlansDpltnAvrg:
	#Run analysis for the various folders
	print "Calculating cluster's expression significance in Furlans' (downregulated) by average in each sample..."
	for p in xrange(len(l_rsltFldrs)):
		enrchmntAnlysFldr = os.path.join(l_rsltFldrs[p],'dffExprssnAvrg.negative.d')
		if not os.path.exists(enrchmntAnlysFldr):
			os.mkdir(enrchmntAnlysFldr)
		k_num = l_KNum[p]
		#calculate cell number differences between clusters and samples
		PAGODA_Rds = os.path.join(l_rsltFldrs[p],'app.final.rds')
		#list of clusters of interest
		if p==0:#added to exclude clusters 4 and 6 (gray and yellow)
			lPoptns = [x+1 for x in xrange(k_num) if x+1 not in {4,6}]
		else:
			lPoptns = [x+1 for x in xrange(k_num)]
		#Calculate pairwise FDR enrichment
		hcRds = ''
		clcPrWsFDRforClstrAvrg(PAGODA_Rds,hcRds,lPoptns,enrchmntAnlysFldr, \
		k_num=k_num,rplcFls=False,frmPAGODAappFl=True,FDRtrhld=0.01, \
		grtr=False)


"""
Compile the characteristic genes of each cluster into a database of
gene names and cluster/cell type
"""
if complDBFurlans:
	print 'Compiling list into databases...'
	hdr='#The current database was created from PAGODA E12 and E13 mouse neural crest development\n#Gene\tCell_type\n'
	lInFldr = [os.path.join(mnFldrNoIntrs,'E13.d','dffExprssn.d'), \
	os.path.join(mnFldrNoIntrs,'E13.d','dffExprssnAvrg.d'), \
	os.path.join(mnFldrNoIntrs,'E13.d','dffExprssn.negative.d'), \
	os.path.join(mnFldrNoIntrs,'E13.d','dffExprssnAvrg.negative.d')]
	loutFlDB = [os.path.join(mnFldrNoIntrs,'E13.d','E13_NC.tsv'), \
	os.path.join(mnFldrNoIntrs,'E13.d','E13_NC.avrg.tsv'), \
	os.path.join(mnFldrNoIntrs,'E13.d','E13_NC.negative.tsv'), \
	os.path.join(mnFldrNoIntrs,'E13.d','E13_NC.negative.avrg.tsv')]
	for p in xrange(len(lInFldr)):
		inFldr = lInFldr[p]
		outFlDB = loutFlDB[p]
		dClustrNumbrName = {'1':'SCP','2':'Chromaffin_cells','3': \
		'Sympathoblast','5':'Bridge_cells'}
		opndFl = open(outFlDB,'w')
		opndFl.write(hdr)
		for flPrfx,cllTyp in dClustrNumbrName.items():
			if os.path.exists(os.path.join(inFldr,'%s_allPairs.tsv'%flPrfx)):
				inFlGns = os.path.join(inFldr,'%s_allPairs.tsv'%flPrfx)
			else:
				inFlGns = os.path.join(inFldr,'clstr_%s_avrgAll.ovrlClstr.tsv'%flPrfx)
			if not os.path.exists(inFlGns):
				print 'File %s do not exists!'%inFlGns
				opndFl.write('\t\n')
			else:
				for l in open(inFlGns,'r'):
					if l.strip() and l.split()[0]!='gene':
						gnNm = l.split('\t')[0]
						opndFl.write('%s\t%s\n'%(gnNm,cllTyp))
		opndFl.close()
	print 'Done!'

"""
Compile the characteristic genes of each cluster into a database of
gene names and cluster/cell type
"""
if complDBFurlans:
	print 'Compiling list into databases...'
	hdr='#The current database was created from PAGODA E12 and E13 mouse neural crest development\n#Gene\tCell_type\n'
	lInFldr = [os.path.join(mnFldrNoIntrs,'E13.KEql6.d','dffExprssn.d'), \
	os.path.join(mnFldrNoIntrs,'E13.KEql6.d','dffExprssnAvrg.d'), \
	os.path.join(mnFldrNoIntrs,'E13.KEql6.d','dffExprssn.negative.d'), \
	os.path.join(mnFldrNoIntrs,'E13.KEql6.d','dffExprssnAvrg.negative.d')]
	loutFlDB = [os.path.join(mnFldrNoIntrs,'E13.KEql6.d','E13_NC.KEql6.tsv'), \
	os.path.join(mnFldrNoIntrs,'E13.KEql6.d','E13_NC.KEql6.avrg.tsv'), \
	os.path.join(mnFldrNoIntrs,'E13.KEql6.d','E13_NC.KEql6.negative.tsv'), \
	os.path.join(mnFldrNoIntrs,'E13.KEql6.d','E13_NC.KEql6.negative.avrg.tsv')]	
	for p in xrange(len(lInFldr)):
		inFldr = lInFldr[p]
		outFlDB = loutFlDB[p]
		dClustrNumbrName = {'1':'SCP_blue_1','2':'Chromaffin_cells','3': \
		'Sympathoblast','4':'SCP_yellow_4','5':'Bridge_cells','6':'gray_6'}
		opndFl = open(outFlDB,'w')
		opndFl.write(hdr)
		for flPrfx,cllTyp in dClustrNumbrName.items():
			if os.path.exists(os.path.join(inFldr,'%s_allPairs.tsv'%flPrfx)):
				inFlGns = os.path.join(inFldr,'%s_allPairs.tsv'%flPrfx)
			else:
				inFlGns = os.path.join(inFldr,'clstr_%s_avrgAll.ovrlClstr.tsv'% \
				flPrfx)
			if not os.path.exists(inFlGns):
				print 'File %s do not exists!'%inFlGns
				opndFl.write('\t\n')
			else:
				for l in open(inFlGns,'r'):
					if l.strip() and l.split()[0]!='gene':
						gnNm = l.split('\t')[0]
						opndFl.write('%s\t%s\n'%(gnNm,cllTyp))
		opndFl.close()
	print 'Done!'

"""
Compile the characteristic genes of each cluster into a database of
gene names and cluster/cell type
"""
if complDBFurlans:
	print 'Compiling list into databases...'
	hdr='#The current database was created from PAGODA E12 and E13 mouse neural crest development\n#Gene\tCell_type\n'
	lInFldr = [os.path.join(mnFldrNoIntrs,'E12.d','dffExprssn.d'), \
	os.path.join(mnFldrNoIntrs,'E12.d','dffExprssnAvrg.d'), \
	os.path.join(mnFldrNoIntrs,'E12.d','dffExprssn.negative.d'), \
	os.path.join(mnFldrNoIntrs,'E12.d','dffExprssnAvrg.negative.d')]
	loutFlDB = [os.path.join(mnFldrNoIntrs,'E12.d','E12_NC.tsv'), \
	os.path.join(mnFldrNoIntrs,'E12.d','E12_NC.avrg.tsv'), \
	os.path.join(mnFldrNoIntrs,'E12.d','E12_NC.negative.tsv'), \
	os.path.join(mnFldrNoIntrs,'E12.d','E12_NC.negative.avrg.tsv')]
	for p in xrange(len(lInFldr)):
		inFldr = lInFldr[p]
		outFlDB = loutFlDB[p]
		dClustrNumbrName = {'1':'Bridge_cells_red_1','2':'Yellow_2','3': \
		'Chromaffin_cells_green_3','4':'SCP_blue_4','5':'Sympathoblast_purple_5'}
		opndFl = open(outFlDB,'w')
		opndFl.write(hdr)
		for flPrfx,cllTyp in dClustrNumbrName.items():
			if os.path.exists(os.path.join(inFldr,'%s_allPairs.tsv'%flPrfx)):
				inFlGns = os.path.join(inFldr,'%s_allPairs.tsv'%flPrfx)
			else:
				inFlGns = os.path.join(inFldr,'clstr_%s_avrgAll.ovrlClstr.tsv'% \
				flPrfx)
			if not os.path.exists(inFlGns):
				print 'File %s do not exists!'%inFlGns
				opndFl.write('\t\n')
			else:
				for l in open(inFlGns,'r'):
					if l.strip() and l.split()[0]!='gene':
						gnNm = l.split('\t')[0]
						opndFl.write('%s\t%s\n'%(gnNm,cllTyp))
		opndFl.close()
	print 'Done!'



#The following dictionary was cured manually to covert gene 
#names/symbols in Furlan's 2018 dataset to GRCm38.p6 gene symbols.

dOriNwName={'0610007N19Rik':'Snhg18','0610007P14Rik':'Erg28','0610011F06Rik':'Mettl26', \
'0610031J06Rik':'Glmp','0610038L08Rik':'Tbc1d22bos','1110001A16Rik':'Cebpzos', \
'1110001J03Rik':'Fmc1','1110007C09Rik':'Card19','1110008F13Rik':'Rab5if', \
'1110008J03Rik':'Rita1','1110034G24Rik':'Shld1','1110035M17Rik':'Zfp652os', \
'1110037F02Rik':'Virma','1110054M08Rik':'Lppos','1110058L19Rik':'Sdhaf4', \
'1190002F15Rik':'Lockd','1190002N15Rik':'Dipk2a','1200011I18Rik':'Gpalpp1', \
'1200014J11Rik':'Ncbp3','1300018J18Rik':'Selenoo','1500011K16Rik':'Mtln', \
'1500012F01Rik':'Zfas1','1500016L03Rik':'Lhx1os','1600002H07Rik':'Tedc2', \
'1600016N20Rik':'Lmntd2','1600029D21Rik':'Plet1','1700009P17Rik':'Cfap126', \
'1700011H14Rik':'ccdc198','1700011J10Rik':'Sp3os','1700021K19Rik':'Rubcn', \
'1700024F13Rik':'2810429I04Rik','1700024P03Rik':'Ccdc34os','1700026L06Rik': \
'Spaca9','1700034F02Rik':'Clhc1','1700040L02Rik':'Cabcoco1','1700045I19Rik': \
'Rnf138rt1','1700052N19Rik':'Armt1','1700066J24Rik':'Fam219aos','1700085B03Rik': \
'Tmem132cos','1700094D03Rik':'4933434E20Rik','1700106J16Rik':'Ccdc182', \
'1700109F18Rik':'Hoxd3os1','1700112E06Rik':'Lrmda','1810011O10Rik':'Tcim', \
'1810019J16Rik':'Kdf1','1810041L15Rik':'Shisal1','1810043H04Rik':'Ndufaf8', \
'2010001M06Rik':'Abhd11os','2010015L04Rik':'Cfap74','2010107G23Rik':'Fam241b', \
'2010111I01Rik':'Aopep','2210013O21Rik':'Nbdy','2210015D19Rik':'Ube2d-ps', \
'2210404O07Rik':'Smim24','2210416O15Rik':'Cuedc1','2210417K05Rik':'Hk1os', \
'2310008H04Rik':'Spidr','2310008H09Rik':'Knop1','2310014L17Rik':'Rnf225', \
'2310035C23Rik':'Relch','2310045N01Rik':'Borcs8','2310047M10Rik':'Borcs6', \
'2310067B10Rik':'Tmem94','2310069G16Rik':'A930017M01Rik','2410004N09Rik': \
'Epb41l4aos','2410015M20Rik':'Micos13','2410016O06Rik':'Riox1','2410057H14Rik': \
'Tspan2os','2410066E13Rik':'Mturn','2410089E03Rik':'Cplane1','2410127L17Rik': \
'Carnmt1','2410137F16Rik':'Srrm4os','2510003E04Rik':'Kif1bp','2510049J12Rik': \
'Mkrn2os','2610015P09Rik':'Ccdc191','2610017I09Rik':'Pantr1','2610018G03Rik': \
'Stk26','2610019F03Rik':'Tdrp','2610034B18Rik':'Arpin','2610034M16Rik':'Pdzph1', \
'2610203C20Rik':'3110039I08Rik','2610204G22Rik':'Atad3aos','2610305D13Rik': \
'Zfp979','2610524H06Rik':'1500011B03Rik','2700029M09Rik':'Hpf1','2700060E02Rik': \
'Rtraf','2700070H01Rik':'None','2700081L22Rik':'Smc2os','2700086A05Rik':'Hoxaas3', \
'2700094K13Rik':'Selenoh','2810008D09Rik':'Snhg20','2810011L19Rik':'Tunar', \
'2810055G20Rik':'Mir99ahg','2810408M09Rik':'Trp53rka','2810417H13Rik':'Pclaf', \
'2810442I21Rik':'Eldr','2900008C10Rik':'Bcor','2900056M20Rik':'Kantr','2900093L17Rik': \
'Frmpd1os','3010026O09Rik':'Mrnip','3110002H16Rik':'Rmc1','3110007F17Rik': \
'Cldn34c1','3110021A11Rik':'None','3110043O21Rik':'C9orf72','3110047P20Rik': \
'Nwd2','3110052M02Rik':'Zfp983','3110057O12Rik':'Abhd18','3110062M04Rik':'Cyren', \
'3632451O06Rik':'Armh4','4632415K11Rik':'Meak7','4632428N05Rik':'Vsir', \
'4632434I11Rik':'Ddias','4732415M23Rik':'Spata33','4732456N10Rik':'Krt90', \
'4833424O15Rik':'Plppr5','4921507L20Rik':'Alkbh3os1','4921530D09Rik':'Ccdc183', \
'4921530L18Rik':'Stamos','4922501C03Rik':'Cep162','4922501L14Rik':'Erich3', \
'4930422G04Rik':'Zgrf1','4930427A07Rik':'Tedc1','4930429B21Rik':'None', \
'4930444A02Rik':'Pomk','4930455F23Rik':'Ccdc181','4930465A12Rik':'Csmd2os', \
'4930469G21Rik':'Tex50','4930487H11Rik':'None','4930506M07Rik':'Shtn1', \
'4930524L23Rik':'Flicr','4930529F24Rik':'Nrg3os','4930529M08Rik':'Cfap61', \
'4930538K18Rik':'Tmem269','4930564G21Rik':'Grip1os2','4931406H21Rik':'None', \
'4931408D14Rik':'None','4931417G12Rik':'Fancd2os','4931429I11Rik':'Jhy', \
'4931430N09Rik':'None','4932418E24Rik':'Ccdc187','4933401P06Rik':'Ppp2r2cos', \
'4933403G14Rik':'Map10','4933411K20Rik':'Cfap97','4933416C03Rik':'Taf7l2', \
'4933426D04Rik':'Prkag2os1','4933426M11Rik':'Susd6','4933436C20Rik':'Crnde', \
'5031426D15Rik':'None','5033411D12Rik':'Sugct','5330426P16Rik':'Dubr', \
'5430417L22Rik':'Inafm2','5730408K05Rik':'None','5730422E09Rik':'Mm2pr', \
'5730508B09Rik':'Fam241a','5730559C18Rik':'Inava','5730577I03Rik':'Fbxl12os', \
'5830415F09Rik':'Trmo','5830418K08Rik':'Cep295','5930403L14Rik':'None', \
'6030419C18Rik':'Insyn1','6330403M23Rik':'Igip','6330408A02Rik':'Zswim9', \
'6330416G13Rik':'Tmem268','6430573F11Rik':'Trmt9b','6720401G13Rik':'Firre', \
'6720456H20Rik':'Tmem260','8430403D17Rik':'Bcas3os1','8430408G22Rik':'Depp1', \
'8430410A17Rik':'Hmces','8430419L09Rik':'Fam234b','8430427H17Rik':'Nol4l', \
'9030624J02Rik':'Vps35l','9130011E15Rik':'Armh3','9130017N09Rik':'Misp', \
'9230110C19Rik':'Cfap300','9230115E21Rik':'Arhgap20os','9330133O14Rik':'None', \
'9430008C03Rik':'Snhg17','9430020K01Rik':'Jcad','9430076C15Rik':'Creb5', \
'9430083A17Rik':'None','9530048J24Rik':'Bcas3os2','9530091C08Rik':'None', \
'9930013L23Rik':'Cemip','A030009H04Rik':'Rnf227','A130049A11Rik':'Phtf1os', \
'A230046K03Rik':'Washc4','A230070E04Rik':'None','A230073K19Rik':'None', \
'A330021E22Rik':'Cfap69','A330050B17Rik':'Mccc1os','A430071A18Rik':'Rptoros', \
'A430105I19Rik':'Ccdc9b','A730017C20Rik':'Minar2','A730085A09Rik':'Kank4os', \
'A830010M20Rik':'Btbd8','A830080D01Rik':'Bclaf3','A930011O12Rik':'Mir124a-1hg', \
'AA415398':'Frg2f1','Aaed1':'Prxl2c','AB099516':'Methig1','Acpl2':'Pxylp1', \
'Adc':'Azin2','Adck3':'Coq8a','Adck4':'Coq8b','Adrbk1':'Grk2','Adrbk2':'Grk3', \
'AF357355':'Gm26922','AF357359':'None','AF357426':'None','Agpat6':'Gpat4', \
'Agphd1':'Hykk','AI118078':'Tmem266','AI314180':'Ecpas','AI316807':'Smim19', \
'AI414108':'Igsf9b','AI462493':'Uqcc3','AI464131':'Myorg','AI480653':'Xxylt1', \
'AI597468':'Tmem263','AI836003':'Ccdc184','AI854517':'Mir9-3hg','Aim1l':'Crybg2', \
'AK010878':'Gon7','AK129341':'Cep126','Ankrd32':'Slf1','Apitd1':'Cenps', \
'Apoa1bp':'Naxe','Athl1':'Pgghg','Atp5s':'Dmac2l','Atp5sl':'Dmac2','Atpbd4': \
'Dph6','AU023871':'Mpig6b','AV051173':'None','AW549542':'Tbx3os1','Azi1':'Cep131', \
'B130006D01Rik':'None','B230120H23Rik':'Map3k20','B230378P21Rik':'Plxna4os1', \
'B3galtl':'B3glct','B430319G15Rik':'Gm38391','B930003M22Rik':'None', \
'B930041F14Rik':'Fndc10','Baat1':'Brat1','Bai1':'Adgrb1','Bai2':'Adgrb2', \
'Bai3':'Adgrb3','BC005764':'Plppr3','BC016423':'Tasor2','BC017612':'Smco4', \
'BC017643':'Cybc1','BC018242':'Plppr2','BC022687':'Clba1','BC023829':'Tmem185a', \
'BC024582':'Bach2os','BC026585':'Cryzl2','BC026590':'Abitram','BC029214':'Paxx', \
'BC030307':'Ttc41','BC030336':'Mosmo','BC037034':'Map11','BC039771':'Rin2', \
'BC051628':'Fndc11','BC056474':'Wdr83os','BC068157':'Prr36','BC089491':'Selenov', \
'BC100451':'Cep295nl','Beta-s':'None','Bloc1s2a':'Bloc1s2','Bre':'Babam2', \
'Bzrap1':'Tspoap1','C030039L03Rik':'Zfp607b','C130030K03Rik':'Grik2', \
'C230052I12Rik':'Faap24','C230081A13Rik':'Peak1','C230091D08Rik':'Snhg14', \
'C330006A16Rik':'Tmem250-ps','C330027C09Rik':'Cip2a','C630020P19Rik':'Pcsk2os1', \
'C77370':'Nexmif','Casc5':'Knl1','Ccbl1':'Kyat1','Ccbl2':'Kyat3','Ccdc101':'Sgf29', \
'Ccdc104':'Cfap36','Ccdc108':'Cfap65','Ccdc109b':'Mcub','Ccdc11':'Cfap53', \
'Ccdc111':'Primpol','Ccdc164':'Drc1','Ccdc176':'Bbof1','Ccdc19':'Cfap45','Ccdc23': \
'Svbp','Ccdc37':'Cfap100','Ccdc41':'Cep83','Ccdc64':'Bicdl1','Ccrn4l':'Noct', \
'Cd97':'Adgre5','Cdk3-ps':'Cdk3','Cecr5':'Hdhd5','Cecr6':'Tmem121b','Cep110': \
'Cntrl','Chi3l1':'Chil1','Clcn4-2':'Clcn4','Cldn25':'Cldnd1','Cml1':'Nat8f1', \
'Crxos1':'Crxos','Csda':'Ybx3','Csrp2bp':'Kat14','Ctage5':'Mia2','Ctgf':'Ccn2', \
'Cxcr7':'Ackr3','Cxx1a':'Rtl8a','Cxx1b':'Rtl8b','Cxx1c':'Rtl8c','Cyb5':'Cyb5a', \
'Cybasc3':'Cyb561a3','Cyr61':'Ccn1','D10Bwg1379e':'Arfgef3','D10Jhu81e':'Gatd3a', \
'D10Wsu52e':'Rtcb','D14Abb1e':'Tasor','D15Ertd621e':'Fam91a1','D17Ertd648e':'None', \
'D17Wsu104e':'Mydgf','D17Wsu92e':'Ilrun','D19Bwg1357e':'Pum3','D19Ertd737e': \
'Fam204a','D330045A20Rik':'Radx','D330050I16Rik':'None','D3Bwg0562e':'Plppr4', \
'D4Ertd617e':'None','D4Wsu53e':'Rsrp1','D630037F22Rik':'Tbc1d32','D630041G03Rik': \
'None','D7Ertd715e':'Snhg14','D8Ertd82e':'Prag1','D930015E06Rik':'Tmem131l', \
'Darc':'Ackr1','Dbc1':'Ccar2','Deb1':'Ss18l2','Dfna5':'Gsdme','Dfnb59':'Pjvk', \
'Diap2':'Diaph2','Diap3':'Diaph3','Dirc2':'Slc49a4','Dlx6as1':'Dlx6os1','Dnahc1': \
'Dnah1','Dnahc10':'Dnah10','Dnahc11':'Dnah11','Dnahc2':'Dnah2','Dnahc7a':'Dnah7a', \
'Dnahc7b':'Dnah7b','Dnahc8':'Dnah8','Dnahc9':'Dnah9','Dnalc1':'Dnal1','Dopey1': \
'Dop1a','Dopey2':'Dop1b','Dos':'Cbarp','Dpcd':'Gm17018','Dscr3':'Vps26c','Dyx1c1': \
'Dnaaf4','E030011O05Rik':'None','E130012A19Rik':'Epop','E130112N10Rik':'Vamp1', \
'E130309D14Rik':'Ccdc92b','E130309F12Rik':'Plppr1','E330033B04Rik':'None', \
'Efcab4a':'Cracr2b','Efcab4b':'Cracr2a','Efha2':'Micu3','Eif2c1':'Ago1','Eif2c2': \
'Ago2','Eif2c3':'Ago3','Eltd1':'Adgrl4','Emr1':'Adgre1','Epb4.1':'Epb41','Epb4.1l1': \
'Epb41l1','Epb4.1l2':'Epb41l2','Epb4.1l3':'Epb41l3','Epb4.1l4a':'Epb41l4a', \
'Epb4.1l4b':'Epb41l4b','Epb4.1l5':'Epb41l5','Epb4.9':'Dmtn','Ept1':'Selenoi', \
'Erbb2ip':'Erbin','Etohi1':'Zfp971','F930015N05Rik':'None','Fam103a1':'Ramac', \
'Fam105a':'Otulinl','Fam108a':'Abhd17a','Fam108c':'Abhd17c','Fam109a':'Pheta1', \
'Fam109b':'Pheta2','Fam115a':'Tcaf1','Fam115c':'Tcaf2','Fam132a':'C1qtnf12', \
'Fam132b':'Erfe','Fam134a':'Retreg2','Fam134b':'Retreg1','Fam154b':'Saxo2', \
'Fam159a':'Shisal2a','Fam173b':'Atpsckmt','Fam175a':'Abraxas1','Fam178a':'Slf2', \
'Fam188a':'Mindy3','Fam188b':'Mindy4','Fam194a':'Erich6','Fam195a':'Mcrip2', \
'Fam196a':'Insyn2a','Fam198a':'Gask1a','Fam198b':'Gask1b','Fam19a1':'Tafa1', \
'Fam19a2':'Tafa2','Fam19a5':'Tafa5','Fam203a':'Hgh1','Fam21':'Washc2','Fam211a': \
'Lrrc75a','Fam211b':'Lrrc75b','Fam212a':'Inka1','Fam213a':'Prxl2a','Fam213b': \
'Prxl2b','Fam35a':'Shld2','Fam46a':'Tent5a','Fam58b':'Ccnq','Fam5b':'Brinp2', \
'Fam5c':'Brinp3','Fam60a':'Sinhcaf','Fam63a':'Mindy1','Fam63b':'Mindy2','Fam64a': \
'Pimreg','Fam65a':'Ripor1','Fam65b':'Ripor2','Fam69a':'Dipk1a','Fam69b':'Dipk1b', \
'Fam86':'Eef2kmt','Fam96a':'Ciao2a','Fam96b':'Ciao2b','Fbxo18':'Fbh1','Figf': \
'Vegfd','Ftsjd1':'Cmtr2','Garem':'Garem1','Gatsl2':'Castor2','Gatsl3':'Castor1', \
'Gbas':'Nipsnap2','Glt25d2':'Colgalt2','Gltscr1':'Bicra','Gltscr1l':'Bicral', \
'Gm10536':'None','Gm10677':'None','Gm10768':'None','Gm10789':'None','Gm11202': \
'Rab11fip4os1','Gm11435':'Heatr9','Gm11974':'Snhg15','Gm129':'Ciart','Gm12942': \
'None','Gm13152':'Zfp982','Gm13157':'Zfp984','Gm13399':'Dbhos','Gm13476':'Zeb2os', \
'Gm13704':'Chn1os3','Gm14378':'Tgfbr3l','Gm14492':'Prr33','Gm14872':'Tomm6os', \
'Gm1564':'Meioc','Gm15800':'Hectd4','Gm16039':'Umad1','Gm16065':'Tbx3os2', \
'Gm16119':'Uckl1os','Gm16197':'Trp53cor1','Gm16515':'Natd1','Gm16516':'Map2k3os', \
'Gm166':'Ccdc189','Gm1661':'Armh1','Gm17296':'Tarbp1','Gm19784':'None','Gm19897': \
'None','Gm20748':'Bvht','Gm2382':'Mthfsl','Gm2518':'Taf6l','Gm3230':'None', \
'Gm3833':'Mnd1-ps','Gm4944':'Zfp994','Gm4980':'Tpbgl','Gm514':'Ubap1l','Gm5506': \
'Eno1b','Gm5607':'Sox1ot','Gm561':'Smim26','Gm5803':'Hnrnpa1l2-ps2','Gm5918': \
'Ccpg1os','Gm6484':'Angptl8','Gm6642':'None','Gm7120':'Tmem267','Gm7173':'Cfap47', \
'Gm7325':'Mymx','Gm872':'Cfap54','Gm9054':'None','Gm9079':'None','Gm996':'Ajm1', \
'Gpr123':'Adgra1','Gpr124':'Adgra2','Gpr125':'Adgra3','Gpr126':'Adgrg6','Gpr30': \
'Gper1','Gpr56':'Adgrg1','Gpr64':'Adgrg2','Gpr98':'Adgrv1','Gsg2':'Haspin', \
'Gtdc2':'Pomgnt2','Gtl3':'Cfap20','Gucy1a3':'Gucy1a1','Gucy1b3':'Gucy1b1', \
'Gyk':'Gk','H2-Ke2':'Pfdn6','Hdgfrp2':'Hdgfl2','Hdgfrp3':'Hdgfl3','Heatr2': \
'Dnaaf5','Hfe2':'Hjv','Hiatl1':'Mfsd14b','Hmga2-ps1':'None','Hmha1':'Arhgap45', \
'Hn1':'Jpt1','Hn1l':'Jpt2','Hnrpdl':'Hnrnpdl','Hnrpll':'Hnrnpll','Hras1':'Hras', \
'Hrsp12':'Rida','Igf2as':'Igf2os','Ikbkap':'Elp1','Inadl':'Patj','Ipw':'Snhg14', \
'Ispd':'Crppa','Itfg3':'Fam234a','Jhdm1d':'Kdm7a','Klhl17':'Plekhn1','l7Rn6': \
'Hikeshi','Lace1':'Afg1l','Large':'Large1','Ldoc1l':'Rtl6','Lect1':'Cnmd', \
'Lepre1':'P3h1','Leprel1':'P3h2','Leprel2':'P3h3','Leprel4':'P3h4','Lnp':'Lnpk', \
'LOC100503496':'None','LOC100504703':'None','LOC100861615':'Gm3411','LOC106740': \
'None','LOC171588':'None','Loh12cr1':'Borcs5','Lphn1':'Adgrl1','Lphn2':'Adgrl2', \
'Lphn3':'Adgrl3','Lrdd':'Pidd1','Lrrc16a':'Carmil1','Lrrc16b':'Carmil3','Lrrc33': \
'Nrros','Lrrc48':'Drc3','Lsmd1':'Naa38','Lyrm5':'Etfrf1','Mesdc1':'Tlnrd1', \
'Mesdc2':'Mesd','Mettl13':'Eef1aknmt','Mettl20':'Etfbkmt','Mfi2':'Meltf','Mfsd7b': \
'Flvcr1','Mfsd7c':'Flvcr2','Micalcl':'Mical2','Mina':'Riox2','Mir1196':'None', \
'Mir143hg':'Carmn','Mir496':'Mir496a','Mir5115':'None','Mira':'None','Mkl2': \
'Mrtfb','Mlf1ip':'Cenpu','Mll1':'Kmt2a','Mll2':'Kmt2b','Mll3':'Kmt2c','Mll5': \
'Kmt2e','Mllt4':'Afdn','Mnf1':'Uqcc2','Mrp63':'Mrpl57','Mtap7d3':'Map7d3', \
'Mterfd2':'Mterf4','Mterfd3':'Mterf2','Mtl5':'Tesmin','Mum1l1':'Pwwp3b','Murc': \
'Cavin4','Mut':'Msh5','Mycl1':'Mycl','Myeov2':'Cops9','N28178':'Phf24','Nadkd1': \
'Nadk2','Nat6':'Naa80','Ncrna00085':'Spaca6','Ncrna00086':'Smim10l2a','Ngfrap1': \
'Bex3','Nhp2l1':'Snu13','Nim1':'Nim1k','Nrp':'Alkbh1','Oraov1':'LTO1','Papd4': \
'Tent2','Park2':'Prkn','Pcdha10':'Pcdha11','Pcnxl2':'Pcnx2','Pcnxl4':'Pcnx4', \
'Pddc1':'Gatd1','Peo1':'Twnk','Phf15':'Jade2','Phf16':'Jade3','Phf17':'Jade1', \
'Pion':'Gsap','Plk1s1':'Kiz','Ppap2a':'Plpp1','Ppap2b':'Plpp3','Ppap2c':'Plpp2', \
'Ppapdc1b':'Plpp5','Ppapdc2':'Plpp6','Ppapdc3':'Plpp7','Ppp2r4':'Ptpa','Prkcdbp': \
'Cavin3','Prkrir':'Thap12','Prr24':'Inafm1','Ptchd2':'Disp3','Ptpla':'Hacd1', \
'Ptplad1':'Hacd3','Ptplad2':'Hacd4','Ptplb':'Hacd2','Ptrf':'Cavin1','Pvrl1': \
'Nectin1','Pvrl2':'Nectin2','Pvrl3':'Nectin3','Rab1':'Rab1a','Rab7l1':'Rab29', \
'Rabl5':'Ift22','Rgag4':'Rtl5','Ric8':'Ric8a','Rltpr':'Carmil2','Rp2h':'Rp2', \
'Rqcd1':'Cnot9','Rsg1':'Cplane2','Rtdr1':'Rsph14','Sc4mol':'Msmo1','Scn2a1': \
'Scn2a','Sdccag3':'Entr1','Sdpr':'Cavin2','Selk':'Selenok','Selm':'Selenom', \
'Selrc1':'Coa7','Sep15':'Selenof','Sepn1':'Selenon','Sepp1':'Selenop','Setd8': \
'Kmt5a','Sfrs18':'Pnisr','Sgol1':'Sgo1','Sgol2':'Sgo2a','Skiv2l2':'Mtrex', \
'Slc22a13b-ps':'Slc22a13b','Slc24a6':'Slc8b1','Slmo1':'Prelid3a','Slmo2': \
'Prelid3b','Smcr7':'Mief2','Smek1':'Ppp4r3a','Smek2':'Ppp4r3b','Snhg7':'None', \
'Soga2':'Mtcl1','Spin2':'Spin2c','Sqrdl':'Sqor','Srcrb4d':'Ssc4d','Ssfa2': \
'Itprid2','Stk30':'Mok','Stra13':'Bhlhe40','Stra13':'Cenpx','Stra13':'Bhlhe40', \
'Stxbp3a':'Stxbp3','Stxbp3b':'Stxbp3-ps','Suv420h1':'Kmt5b','Suv420h2':'Kmt5c', \
'Tbrg3':'None','Tceb1':'Eloc','Tceb3':'Eloa','Tenc1':'Tns2','Tmem110':'Stimate', \
'Tmem180':'Mfsd13a','Tmem194':'Nemp1','Tmem2':'Cemip2','Tmem35':'Tmem35a', \
'Tmem48':'Ndc1','Tmem5':'Rxylt1','Tmem55a':'Pip4p2','Tmem55b':'Pip4p1','Tmem57': \
'Maco1','Tmem66':'Saraf','Tssc1':'Eipr1','Ttc18':'Cfap70','Ubl4':'Ubl4a', \
'Ufd1l':'Ufd1','Uqcc':'Uqcc1','Usmg5':'Atp5md','Utp11l':'Utp11','Vimp':'Selenos', \
'Vwa9':'Ints14','Wbp5':'Tceal9','Wbscr16':'Rcc1l','Wbscr17':'Galnt17','Wbscr27': \
'Mettl27','Wdr52':'Cfap44','Wdr67':'Tbc1d31','Wdr96':'Cfap43','Whsc1':'Nsd2', \
'Zbtbd6':'Kbtbd6','Zcchc11':'Tut4','Zcchc16':'Rtl4','Zcchc6':'Tut7','Zfhx2as': \
'Zfhx2os','Zfp607':'Zfp607a','Znf512b':'Zfp512b','Zufsp':'Zup1'}

########################################
########################################
########################################
###  Merge databases for expression  ###
########################################
########################################
########################################
#If True next will merge up and downregulated genes
if mrgDBs:
	ordrdSfxPrfx = [('','_uprgltd'),('.negative','_dwnrgltd'),('','_upOdwn'), \
	('.negative','_upOdwn')]	
	#List of folders with results
	l_rsltFldrs = [os.path.join(mnFldrNoIntrs,'E13.d'), \
	os.path.join(mnFldrNoIntrs,'E13.KEql6.d'), \
	os.path.join(mnFldrNoIntrs,'E12.d')]
	for fldr in l_rsltFldrs:
		sfx=os.path.split(fldr)[1].replace('KEql6','K6').replace('E13','E13_NC'). \
		replace('E12','E12_NC').split('.d')[0]
		sGnNms=set()
		outFl = open(os.path.join(outFldr,'%s.upOdwn.tsv'%sfx),'w')
		outFl.write('#Gene\tCell_type\n')
		for prfxFl,sfxFl in ordrdSfxPrfx:
			inFl=os.path.join(fldr,'%s%s.tsv'%(sfx.replace('.K6','.KEql6'),prfxFl))
			for l in open(inFl,'r'):
				if l.strip() and l[0]!='#':
					l=l.splitlines()[0].split('\t')
					gnNm=dOriNwName.get(l[0],l[0])
					l[1]='%s%s'%(l[1],sfxFl)
					if gnNm is not None and (gnNm,l[1]) not in sGnNms:
						l[0]=gnNm
						sGnNms.add((gnNm,l[1]))
						outFl.write('%s\n'%'\t'.join(l))
		#
		sGnNms=set()
		outFl = open(os.path.join(outFldr,'%s.upOdwn.avrg.tsv'%sfx),'w')
		outFl.write('#Gene\tCell_type\n')
		for prfxFl,sfxFl in ordrdSfxPrfx:
			inFl=os.path.join(fldr,'%s%s.avrg.tsv'%(sfx.replace('.K6', \
			'.KEql6'),prfxFl))
			for l in open(inFl,'r'):
				if l.strip() and l[0]!='#':
					l=l.splitlines()[0].split('\t')
					gnNm=dOriNwName.get(l[0],l[0])
					l[1]='%s%s'%(l[1],sfxFl)
					if gnNm is not None and (gnNm,l[1]) not in sGnNms:
						l[0]=gnNm
						sGnNms.add((gnNm,l[1]))
						outFl.write('%s\n'%'\t'.join(l))
		outFl.close()
		



##########################################################
##########################################################
##########################################################
###  Compile average gene expression for each cluster  ###
##########################################################
##########################################################
##########################################################

#To merge the average expression in the databases

E12_NC = True#If True will merge the average expression for E12
E13_NC = True#If True will merge the average expression for E13

def rtrnd(infl,clmnK=0,clmnV=1):
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



###########
#For E12_NC
###########
if E12_NC:
	exprssAvrgFldr=os.path.join(mnFldrNoIntrs,'E12.d','dffExprssnAvrg.d')
	lFls=[f for f in os.listdir(exprssAvrgFldr) if f.find('clstr_')>-1 and \
	f.find('_avrgAll.tsv')>-1]
	dClstrNms={'1':'Bridge_cells_red_1','2':'Yellow_2','3':'Chromaffin_cells_green_3', \
	'4':'SCP_blue_4','5':'Sympathoblast_purple_5'}
	dGnNmsClstrExprssn={}
	sGnNms=set()
	sClstrs=set()
	for clstr in dClstrNms.keys():
			sClstrs.add(clstr)
			cllnNm=dClstrNms[clstr]
			dGnNmsExprssn=rtrnd(os.path.join(exprssAvrgFldr, \
			'clstr_%s_avrgAll.tsv'%clstr))
			for gnNm in dGnNmsExprssn.keys():
					sGnNms.add(gnNm)
					if dGnNmsClstrExprssn.has_key(gnNm):
							dGnNmsClstrExprssn[gnNm][clstr]= \
							dGnNmsExprssn[gnNm]
					else:
							dGnNmsClstrExprssn[gnNm]={clstr: \
							dGnNmsExprssn[gnNm]}
	#
	#
	sGnNms=sorted(sGnNms)
	sClstrs=sorted(sClstrs)
	outFl=open(os.path.join(outFldr,'E12_NC.crrltn.tsv'),'w')
	outFl.write('#Gene\t%s\n'%'\t'.join([dClstrNms[c] for c in sClstrs]))
	for gnNm in sGnNms:
		newGnNm = dOriNwName.get(gnNm,gnNm)
		if newGnNm is not None:
			outFl.write('%s\t%s\n'%(newGnNm,'\t'.join([dGnNmsClstrExprssn[gnNm][c] \
			for c in sClstrs])))
	#
	outFl.close()



###########
#For E13_NC
###########
if E13_NC:
	exprssAvrgFldr=os.path.join(mnFldrNoIntrs,'E13.KEql6.d', \
	'dffExprssnAvrg.d')
	lFls=[f for f in os.listdir(exprssAvrgFldr) if f.find('clstr_')>-1 and \
	f.find('_avrgAll.tsv')>-1]
	dClstrNms={'1':'SCP','2':'Chromaffin_cells','3':'Sympathoblast','4': \
	'SCP_yellow_4','5':'Bridge_cells','6':'gray_6'}
	dGnNmsClstrExprssn={}
	sGnNms=set()
	sClstrs=set()
	for clstr in dClstrNms.keys():
		sClstrs.add(clstr)
		cllnNm=dClstrNms[clstr]
		dGnNmsExprssn=rtrnd(os.path.join(exprssAvrgFldr,'clstr_%s_avrgAll.tsv'% \
		clstr))
		for gnNm in dGnNmsExprssn.keys():
			sGnNms.add(gnNm)
			if dGnNmsClstrExprssn.has_key(gnNm):
				dGnNmsClstrExprssn[gnNm][clstr]=dGnNmsExprssn[gnNm]
			else:
				dGnNmsClstrExprssn[gnNm]={clstr:dGnNmsExprssn[gnNm]}
	#
	#
	sGnNms=sorted(sGnNms)
	sClstrs=sorted(sClstrs)
	outFl=open(os.path.join(outFldr,'E13_NC.crrltn.tsv'),'w')
	outFl.write('#Gene\t%s\n'%'\t'.join([dClstrNms[c] for c in sClstrs]))
	for gnNm in sGnNms:
		newGnNm = dOriNwName.get(gnNm,gnNm)
		if newGnNm is not None:
			outFl.write('%s\t%s\n'%(newGnNm,'\t'.join([dGnNmsClstrExprssn[gnNm][c] \
			for c in sClstrs])))
	#
	outFl.close()

