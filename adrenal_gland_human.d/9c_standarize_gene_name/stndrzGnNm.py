#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  stndzGnNm.py
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
These script aim to standarized gene names/symbols of databases and 
data sources other than those generated in this study, to genecode gene
names/symbols.
"""

import os

######################
######################
#####  Switches  #####
######################
######################
tstIndxEnrchmnts = True # If True will test if the files for the index enrichments exists (for other databases different than ours)
tstIndxCrrltns = True # If True will test if the files for the index correlations exists (for other databases different than ours)
autmzOldToNwNm = True # If True will make a dictionary to patch gene names conversion to GRCm38_p6 names
crrctIndxEnrchmntsFls = True # If True will correct the indexed enrichments to convert GRCh38_p12/GRCm38_p6 with genecodeV28Comp/genecodeV18Comp, respectively
crrctIndxCrrltnsFls = True # If True will correct the indexed correlations to convert GRCh38_p12/GRCm38_p6 with genecodeV28Comp/genecodeV18Comp, respectively
tstIndxEnrchmntsFls = True# If True will test if the differences in between the depot and the final file are corresponding to genes not present in GRCh38_p12/GRCm38_p6 but present in genecodeV28Comp/genecodeV18Comp (enrichments)
tstIndxCrrltnsFls = True# If True will test if the differences in between the depot and the final file are corresponding to genes not present in GRCh38_p12/GRCm38_p6 but present in genecodeV28Comp/genecodeV18Comp (correlations)

inFl_hg38p12 = 'data/GRCh38_p12_ENSMBLgnNm.tsv'
inFl_mm10p6 = 'data/GRCm38_p6_ENSEMBLgnNmEdtd.txt'
inFl_hg38anntn = 'data/cmmnGnNmTOEnsmbl.hg38'
inFl_mm10anntn = 'data/cmmnGnNmTOEnsmbl.mm10'

sGns_hg38p12 = set([l.splitlines()[0].split('\t')[1] for l in \
open(inFl_hg38p12,'r') if l.strip() and l.split('\t')[0]!='Gene name' \
and l[0]!='#'])
sGns_mm10p6 = set([l.splitlines()[0].split('\t')[1] for l in \
open(inFl_mm10p6,'r') if l.strip() and l.split('\t')[0]!='Gene name' \
and l[0]!='#'])
dSppGns = {'Human':sGns_hg38p12,'Mouse':sGns_mm10p6}


######################
######################
if tstIndxEnrchmnts:
	#The sets of genes in sGnsHg and sGnsMm are known to not be present
	#in GRCh38_p12 / GRCm38_p6 annotations, respectively. Manually,
	#the equivalent name for these genes in GRCh38_p12 / GRCm38_p6 
	#annotations (if available) must be included in the input databases.
	#The details of these genes are commented next.
	sGnsHg=set(['KIRREL', 'MRC2', 'ZD76G03', 'NOTCH2NL', 'HLA+AC0-C', \
	'NKX2+AC0-1', 'CMAH', 'C1ORF54', 'BAGE3', 'PTRF', 'LHFP', 'C4ORF32', \
	'SERTAD4+AC0-AS1', 'EPR1', 'CYR61', 'FAM46A', 'KIAA1462'])
	sGnsMm=set(['Cirh1a', 'C920025E04Rik', 'Gm16532', 'Mb21d1', \
	'6330403A02Rik', 'Bcmo1', '4930415F15Rik', 'Dlx6as2', 'Fam46c', \
	'Gyltl1b', 'Gm684', 'Fam101a', '4932443I19Rik', 'Gm13845', \
	'1500002O10Rik'])
	"""
	For human GRCh38_p12_ENSMBLgnNm the dictionary is the following:
	dOldGnNmNwGnNmHg = {'BAGE3':None,'C1ORF54':'C1orf54','C4ORF32': \
	'FAM241A','CMAH':'CMAHP','CYR61':'CCN1','EPR1':None,'FAM46A':'TENT5A', \
	'HLA+AC0-C':'HLA-C','KIAA1462':'JCAD','KIRREL':'KIRREL1','LHFP':None, \
	'MRC2':'AC080038.1','NKX2+AC0-1':'NKX2-1','NOTCH2NL':'NOTCH2NLA', \
	'PTRF':'CAVIN1','SERTAD4+AC0-AS1':'SERTAD4-AS1','ZD76G03':None}
	#
	For mouse GRCm38_p6_ENSMBLgnNm the dictionary is the following:
	dOldGnNmNwGnNmMm = {'1500002O10Rik':'Gad1os','4930415F15Rik': \
	'Spata48','4932443I19Rik':'Cfap97d2','6330403A02Rik':'Stum','Bcmo1': \
	'Bco1','C920025E04Rik':'H2-T-ps','Cirh1a':'Utp4','Dlx6as2':'Dlx6os2', \
	'Fam101a':'Rflna','Fam46c':'Tent5c','Gm13845':'Mkln1os','Gm16532': \
	'Smim17','Gm684':'Colca2','Gyltl1b':'Large2','Mb21d1':'Cgas'}
	"""
	dSppNtPrsntGns={'Human':set(),'Mouse':set()}
	inFlEnrchmnt = 'index.enrchmnts.txt'
	for l in open(inFlEnrchmnt,'r'):
		if l.strip() and l[0]!='#':
			l=l.splitlines()[0]
			infl = l.split('\t')[4]
			gnClmn = int(l.split('\t')[2])
			spp = l.split('\t')[5]
			sGns = dSppGns[spp]
			if not os.path.exists(infl):
				fldrNm,flNm = os.path.split(infl)
				infl=os.path.join(fldrNm,'depot.d',flNm)
				try:
					assert os.path.exists(infl)
					for lne in open(infl,'r'):
						if lne.strip() and lne[0]!='#':
							lne = lne.splitlines()[0].split('\t')
							if lne[gnClmn] not in sGns:
								if lne[gnClmn] in sGnsHg or lne[gnClmn] in sGnsMm:
									print '%s %s'%(lne[gnClmn],infl)
								dSppNtPrsntGns[spp].add(lne[gnClmn])
				except:
					print 'Exception: %s'%infl
	#dSppNtPrsntGns will hold the gene symbols/names that are not present
	#in the ENSEMBL gene annotations GRCh38_p12/GRCm38_p6
	"""
	Gene names/symbols in dSppNtPrsntGns need to be corrected in the 
	original files. So the original files should be moved to files with 
	the ".ori" sufix, the gene symbols corrected and saved in files 
	without the ".ori" sufix, and finally located in a "depot.d" within
	the "location" folder included in the "index.enrchmnts.txt" file. 

	After correcting the names, a second run of the script must have as
	result: dSppNtPrsntGns=={'Mouse': set([]),'Human': set([])}
	"""

######################
######################
if tstIndxCrrltns:
	dSppNtPrsntGns={'Human':set(),'Mouse':set()}
	inFlCrrltn = 'index.crrltns.txt'
	#The sets of genes in sGnsHg and sGnsMm are known to not be present
	#in GRCh38_p12 / GRCm38_p6 annotations, respectively. Manually,
	#the equivalent name for these genes in GRCh38_p12 / GRCm38_p6 
	#annotations (if available) must be included in the input databases.
	#The details of these genes are commented next.
	sGnsHg=set()
	sGnsMm=set(['Selm','1700049G17Rik','2700050L05Rik','Gm101','Selk','Ssfa2','Adck4','BC033916', \
	'2310042D19Rik','H2-Ke2','BC061194','3010026O09Rik','F930015N05Rik','C77370','Gm13476','Ngfrap1', \
	'Fam132b','4632415K11Rik','Fam179b','1500016L03Rik','Fam115c','BC056474','2210417K05Rik','Mfsd7b', \
	'2010107G23Rik','A530054K11Rik','Dak','Fam73a','Fam73b','1700106J16Rik','Lphn1','Hmha1','Nov', \
	'Fam213a','Wbp7','Wdr96','Wbp5','4930429B21Rik','2510003E04Rik','1700019G17Rik','Lyrm5','Rab7l1', \
	'2900002K06Rik','2810442I21Rik','Gm6531','Adck3','Gnb2l1','4930415F15Rik','Sqrdl','Ankrd32','Gm561', \
	'AK129341','Fam101b','Fam101a','Gm10474','Garem','Acpl2','LOC106740','1700012B15Rik','Gm3833', \
	'5730559C18Rik','Stk30','BC020402','6720401G13Rik','2210404J11Rik','Ptchd2','4930471I20Rik', \
	'AI464131','Gm13845','Rab1','Vimp','Gm4532','Wbscr17','Dus2l','Mllt4','Gm16119','2700089E24Rik', \
	'1600016N20Rik','Fam35a','Ptrf','Fam213b','2410127L17Rik','9130011E15Rik','Gm6402','4933426M11Rik', \
	'Fam86','Igj','Fam203a','AI848285','AI414108','Lins','Ept1','2210018M11Rik','Raet1a','1200014J11Rik', \
	'4931408D14Rik','Gpr98','Gm3230','Gbas','Spink3','1700025E21Rik','Uqcc','Ppap2c','Fam65b','Fam65a', \
	'Gm13251','Gm13154','Gm13157','Cxx1b','Gm13152','AI118078','Fam115a','1700009P17Rik','5033411D12Rik', \
	'Xrcc6bp1','C130030K03Rik','Gm11974','Phf16','Skiv2l2','C230081A13Rik','4732415M23Rik','Whsc1l1', \
	'5031426D15Rik','Tex40','Tmem180','0610009D07Rik','Papd7','Papd4','Papd5','Diap3','Diap2','Diap1', \
	'1600002H07Rik','Fam188b','Wdr20a','BC021891','Zfp85-rs1','Naprt1','Cybasc3','D4Wsu53e','Fbxo18', \
	'Ccdc37','Ctgf','Ccrn4l','B630005N14Rik','AI314180','C330006A16Rik','l7Rn6','B3gnt1','Sepp1', \
	'4932415G12Rik','Solh','1700085B03Rik','Gm14378','Mkl1','Hnrpll','Rtdr1','Agpat6','Aaed1','Ifi27l1', \
	'Gm4975','Agpat9','Wdr16','Dom3z','D2Wsu81e','4921530D09Rik','Gm16039','2010001M06Rik','Dscr3', \
	'Fam46a','3110035E14Rik','3110002H16Rik','Fam46b','C030016D13Rik','D830031N03Rik','Efcab4b', \
	'6330408A02Rik','9330133O14Rik','Fam134c','Fam134b','Fam134a','Ccdc64','Ccdc67','Pcnxl4','Pcnxl3', \
	'A330050B17Rik','5730408K05Rik','2610524H06Rik','2610034B18Rik','4930455F23Rik','Fam212b','BC039771', \
	'1700022P22Rik','Gm5506','Tmem194b','1700024P03Rik','4930451C15Rik','Oraov1','BB283400','Gm11985', \
	'Agphd1','1700003M02Rik','4921525B02Rik','Fam194a','2010107E04Rik','Ccdc129','Fam46c','Gm6251', \
	'2310045N01Rik','Dlx6as2','Gm14872','Gm14873','Mir5115','Gyltl1b','Figf','C030039L03Rik', \
	'1700007G11Rik','Aim1l','4632428N05Rik','Ctage5','C920025E04Rik','Dlx6as1','Gpr125','AV051173', \
	'Fam188a','Wapal','4921506M07Rik','Gpr115','Prmt10','D19Ertd737e','Smcr7','Lsmd1','Gm5803','Tmem35', \
	'BC021785','Gtpbp5','D630032N06Rik','Pvrl3','Plk1s1','Soga2','Pvrl4','Minos1','8430427H17Rik', \
	'Klhdc5','Mir3473','4933436C20Rik','6330416G13Rik','4930469G21Rik','BC003331','Wdr85','Slc24a6', \
	'E030011O05Rik','Snhg7','9130017N09Rik','Gm4980','3110043O21Rik','4921507L20Rik','D330050I16Rik', \
	'Dirc2','2700081L22Rik','Gm9047','Tmem66','Dnahc2','Dnahc1','AI836003','Ptplb','Dnahc5','A230073K19Rik', \
	'Fam105a','4930487H11Rik','Tmem57','Mki67ip','D330045A20Rik','2900092D14Rik','4833403I15Rik', \
	'1810043H04Rik','BC017643','1600029D21Rik','Gm7325','Fam211a','Gm4371','Gm10007','D15Ertd621e', \
	'Fam109b','Fam109a','Gm12942','4930427A07Rik','Gm20199','5830415F09Rik','6720456H20Rik','D10Jhu81e', \
	'BC089491','4732456N10Rik','Gm17801','Gm14420','BC023829','9430076C15Rik','AA415398','Gm16516', \
	'Mtl5','2700029M09Rik','B430319G15Rik','AI846148','D630045M09Rik','5330426P16Rik','F630111L10Rik', \
	'Gm6642','Wdr20b','Mll2','Mrp63','Ccdc23','Gpr64','Gm16197','Obfc1','Gpr111','Chi3l1','Ccdc90a', \
	'Dopey2','5830418K08Rik','Slc22a13b-ps','Mesdc1','Mesdc2','Usmg5','2310061J03Rik','BC018507', \
	'A830010M20Rik','Wisp2','AU019823','Wisp1','Fam173b','Zcchc11','Gm4636','Hiat1','AF357355', \
	'Tsga14','Dopey1','AF357359','Gm17365','Gltpd1','0610038L08Rik','AF357426','D3Bwg0562e', \
	'Fert2','2210013O21Rik','Zbtbd6','Sep15','Eif2c3','4930506M07Rik','B130006D01Rik','Gm996', \
	'4930526I15Rik','Gm5126','4930441O14Rik','2410057H14Rik','B230216G23Rik','Rp2h','2610305D13Rik', \
	'3110057O12Rik','Loh12cr1','1190002F15Rik','Gatsl3','Gatsl2','1700052N19Rik','AI597468','Gtl3', \
	'Enthd2','Trove2','9430020K01Rik','3632451O06Rik','Gm20187','Ccdc42b','3110062M04Rik','1300002K09Rik', \
	'Fam196a','2310003H01Rik','Gpr126','2210404O09Rik','Gpr124','2610019F03Rik','LOC100504703','Gpr123', \
	'2310036O22Rik','1110057K04Rik','Wdr67','Wdr65','BC026585','Agxt2l2','2310014L17Rik','Hmga2-ps1', \
	'Ncrna00085','1110058L19Rik','A930011O12Rik','Hn1l','D14Abb1e','Athl1','Wbscr27','Dos','AA388235', \
	'Stxbp3a','2310067B10Rik','Spin2','Mycl1','Rfwd2','Fam64a','1700013G23Rik','Ncrna00086','Wbscr16', \
	'Tceb1','1500011K16Rik','4933426D04Rik','1700109F18Rik','2810008D09Rik','5730508B09Rik', \
	'2310047M10Rik','Mettl10','Ccbl2','C230052I12Rik','Shfm1','Gm608','Prune','Srcrb4d','2810055G20Rik', \
	'Glt25d1','1810033B17Rik','Fam26f','Dgcr14','Dnahc17','Fam5b','Tssc1','Nrp','Rqcd1','D7Ertd715e', \
	'1700101E01Rik','2010003O02Rik','Lnp','Gm10536','A630066F11Rik','Mina','Cecr5','Rtfdc1','D17Wsu104e', \
	'2310069G16Rik','Prosc','Gm10789','Qtrtd1','Gm514','Ccdc109b','Zfp828','Gsg2','1700015E13Rik', \
	'Mettl13','Inadl','Mum1l1','Gtdc2','1110034G24Rik','Gm129','Gm4944','Gm2382','1700034F02Rik', \
	'Pin1-ps1','Rltpr','Ccdc53','Ccdc55','Cecr6','2010107G12Rik','4930568K20Rik','BC037034', \
	'4931406H21Rik','1110054M08Rik','D330022K07Rik','2610034M16Rik','Leprel2','Leprel1','Mlf1ip', \
	'1110004E09Rik','C030046E11Rik','Nt5c3l','Hdgfrp2','Hdgfrp3','Ccdc132','Abp1','BC017612','Mterf', \
	'4930414L22Rik','Siglec5','Gltscr1l','5430435G22Rik','1700094D03Rik','2900055J20Rik','Nhp2l1', \
	'Dnahc11','Fam159b','1110007C09Rik','Ftsjd1','Fam159a','Nup62-il4i1','2700086A05Rik','Igf2as', \
	'1300018J18Rik','Dfnb59','Cml1','Casc5','Gm13826','5430428K19Rik','Mir219-2','2900008C10Rik', \
	'Mfsd7c','9830001H06Rik','Dfna5','Itfg3','E130309D14Rik','Gm20748','2610100L16Rik','Baat1', \
	'AA474331','Fam26e','Bzrap1','Dnahc8','Ubl4','Dnahc9','E130309F12Rik','Gm10677','Epb4.1l4b', \
	'A930013F10Rik','2010015L04Rik','Nadkd1','Adrbk1','Adrbk2','9230115E21Rik','Tceb3','Tceb2', \
	'Dnalc4','Ccdc164','Gyk','9430083A17Rik','Dnalc1','Hras1','Eftud1','Myeov2','2310008H09Rik', \
	'LOC100504608','Msx1as','2610203C20Rik','Dyx1c1','Rabl5','4632434I11Rik','A730085A09Rik','Smek2', \
	'Gm9079','Pion','Sgol2','Fam195b','Pacrgl','Gareml','Fam195a','Clcn4-2','Large','2900056M20Rik', \
	'2810474O19Rik','Gm15800','Gm9897','Stxbp3b','E130012A19Rik','Ldoc1l','Leprel4','Tenc1','AI450353', \
	'Apoa1bp','Gucy1b3','BC027231','Gltscr1','Gltscr2','Zfp191','Mettl20','6330403M23Rik','Zfp607', \
	'BC125332','2810408M09Rik','Lrrc48','Zufsp','C78339','2810011L19Rik','Lrrc16a','Gm2518','Ric8', \
	'Rgag4','A630020A06','Ddx26b','Gm5480','Phf15','Cd97','4931429I11Rik','AI316807','D10Bwg1379e', \
	'Fam103a1','Prkrir','1700021K19Rik','Beta-s','Ccdc19','Ppap2b','Zfyve20','Mum1','C330027C09Rik', \
	'0610011F06Rik','Hrsp12','4930405P13Rik','1200011I18Rik','Sdccag3','4932443I19Rik','AU023871', \
	'Grlf1','2410076I21Rik','Znf512b','1190003J15Rik','E430025E21Rik','6430573F11Rik','Gm684','Phf17', \
	'1100001G20Rik','Dnahc7b','A730017C20Rik','Emr4','Ccdc176','N6amt2','Dbc1','4922501C03Rik','Mgea5', \
	'Tmem2','Vprbp','A330021E22Rik','C230091D08Rik','Gm1045','2410012M07Rik','Narg2','BC068157','Fam108a', \
	'Fam108b','Mettl21d','Zfp862','Gm14005','D630037F22Rik','Prr24','A730098P11Rik','1700011J10Rik', \
	'2700094K13Rik','A230070E04Rik','Whsc1','4921530L18Rik','Bcmo1','Tmem194','1190002N15Rik','Gm11602', \
	'Csrp2bp','Ptpla','Wdr52','Cyb5','Peo1','E330033B04Rik','Gm3258','Mfsd4','1700011H14Rik', \
	'4933427G17Rik','2410015M20Rik','B230120H23Rik','2310008H04Rik','1110008J03Rik','Ispd', \
	'5730457N03Rik','Fam63a','Fam63b','BC026590','2610301G19Rik','Rsg1','Wibg','Gm1564','Gm11744','Ndnl2', \
	'2210404O07Rik','BC018242','Mtap7d3','5730577I03Rik','2210416O15Rik','Hfe2','Ppapdc1b','Gm10336', \
	'Selrc1','Wash','Azi1','Cxx1a','Cxx1c','4930528A17Rik','LOC100503496','5430417L22Rik','2610015P09Rik', \
	'Utp11l','Mterfd3','9230110C19Rik','Murc','2610018G03Rik','Jhdm1d','Gm16515','Aim1','2410004N09Rik', \
	'Hn1','4932416H05Rik','A630007B06Rik','Cldn25','6330407A03Rik','B3galtl','3110047P20Rik', \
	'2610002J02Rik','2810403A07Rik','2700060E02Rik','AK010878','Klhl17','Pex11c','Gm5595','Apitd1', \
	'Zfp71-rs1','6030419C18Rik','Fam108c','Pddc1','Gm13704','Cnih','1500032L24Rik','2610019E17Rik', \
	'Rnmtl1','1500012F01Rik','3110052M02Rik','8430403D17Rik','1110001A16Rik','9830147E19Rik', \
	'E030019B06Rik','2310044G17Rik','Cep110','A230046K03Rik','Zfml','C1s','2610204G22Rik','Ict1', \
	'4931417G12Rik','Gm13718','Erbb2ip','Cirh1a','Epb4.1l2','Epb4.1l3','Gm16532','Epb4.1l1','BC029214', \
	'Ccdc41','Epb4.1l5','9430008C03Rik','Zcchc5','Gm7120','Suv420h1','Suv420h2','Ccdc101','Epb4.9','Mnf1', \
	'Ccdc104','0610007N19Rik','3230401D17Rik','2510049J12Rik','D930015E06Rik','8430419L09Rik', \
	'2610020H08Rik','2310010M20Rik','Mb21d1','Smek1','Has2as','Gm10560','Bai1','Bai2','Bai3', \
	'9030624J02Rik','Tbrg3','Gm3925','Gm17296','Epb4.1','Gucy1a3','2310035C23Rik','1810011O10Rik', \
	'4933413G19Rik','D17Wsu92e','Lrrc16b','Fam154b','Cyr61','Asun','Wbscr22','9630033F20Rik','Ccdc111', \
	'Gm5820','AI462493','A130049A11Rik','Cxcr7','D10Wsu52e','Prkcdbp','Gm10664','Gm19424','Fam178a', \
	'D19Bwg1357e','1700012D01Rik','B930041F14Rik','D4Ertd617e','Atp5s','BC030307','Lace1','Trp53rk', \
	'A430071A18Rik','Sepn1','Scn2a1','BC049635','4930422G04Rik','Csda','Efha1','Mterfd2','Mterfd1', \
	'Efha2','Rps4y2','BC068281','Gm12060','9130206I24Rik','3110021A11Rik','0610007P14Rik','1110001J03Rik', \
	'2410089E03Rik','Ccdc94','Dnahc7a','D8Ertd82e','Gm5918','1810026J23Rik','Fam212a','Tmem48','AW549542', \
	'Gm5088','Ufd1l','Gm13247','Fam198a','Fam19a4','Fam198b','Gm6484','4930529M08Rik','BC016423', \
	'Fam132a','4930529F24Rik','9030617O03Rik','2210015D19Rik','AI480653','1110037F02Rik','Pet112','Mut', \
	'Slmo2','Ttc18','Carkd','Slmo1','Mkl2','2410016O06Rik','Gm17769','Gm17762','Mir143hg','B930003M22Rik', \
	'Tmem27','2810428I15Rik','Taf4a','Gm5607','Ppp2r4','BC048502','Acn9','Lepre1','2010012O05Rik','Mira', \
	'D630013N20Rik','BC005764','Sepw1','Dpcr1','Hiatl1','2410066E13Rik','Fam96a','Fam96b','Gm711','Lect1', \
	'Ppapdc1a','Pvrl1','Deb1','Gm19897','5930403L14Rik','1700026L06Rik','Fam21','Lphn3','Gm166','Park2', \
	'Lphn2','Fam175a','Fam175b','Ftsj2','Fam60a','AA987161','Sgol1','Mll3','8430410A17Rik','Mll1', \
	'Tmem55b','9430016H08Rik','Adc','Mll5','Nim1','8430408G22Rik','3110007F17Rik','Tmem110', \
	'5730422E09Rik','0610037L13Rik','1700016L04Rik','Gm3414','4833424O15Rik','4930524L23Rik','Ccrl1', \
	'0610009O20Rik','Ccdc79','4933403G14Rik','BC053749','1700019L03Rik','Wisp3','Nat6','1700045I19Rik', \
	'Ptplad2','1500002O10Rik','2410018M08Rik','Sdpr','Heatr2','Hnrpdl','2810417H13Rik','1810041L15Rik', \
	'Ppapdc2','Pvrl2','1700024F13Rik','Fam150a','Fam105b','4930444A02Rik','Fam150b','Mir5109','Ccbl1', \
	'Epb4.1l4a','Wdr8','Zfp259','Sfrs18','1700084C01Rik','1810013D10Rik','A030009H04Rik','A830080D01Rik', \
	'AI854517','C030030A07Rik','Ppap2a','AI314831','4921511H03Rik','Rn45s','Ftsjd2','2810407C02Rik', \
	'C630020P19Rik','Efcab4a','Tmem5','1110035M17Rik','LOC100861615','0610031J06Rik','Gm4827', \
	'6330403A02Rik','Bloc1s2a','Lrdd','2410137F16Rik','Ptplad1','Gpr44','3110001D03Rik','Glt25d2', \
	'Bet3l','4933411K20Rik','E130112N10Rik','Setd8','Gm7173','Fam58b','Spt1','Prosapip1','BC030336', \
	'4930538K18Rik','Atp5sl','BC048609','Darc','Ikbkap','Zcchc6','4930500J02Rik','Cpsf3l','Gm9199', \
	'Lrrc33','Smcr7l','Cml2','Gm9054','Vwa9','BC031361','Etohi1','1110018G07Rik','Cdk3-ps', \
	'9530091C08Rik','Gm216','2010111I01Rik','Fam19a5','BC022687','Sc4mol','Gm6194','Bre','Dpcd', \
	'4933422H20Rik','0610012H03Rik','2610017I09Rik','Ppapdc3','Eif2c4','A230056J06Rik','Atpbd4','Eif2c1', \
	'1110008F13Rik','Eif2c2','Fam69c','Fam69b','Fam69a','Tmem55a','Gm1661','1700112E06Rik', \
	'9430023L20Rik','Stra13','Gpr56'])
	for l in open(inFlCrrltn,'r'):
		if l.strip() and l[0]!='#':
			l=l.splitlines()[0]
			infl = l.split('\t')[3]
			gnClmn = 0
			spp = l.split('\t')[4]
			sGns = dSppGns[spp]
			if not os.path.exists(infl):
				fldrNm,flNm = os.path.split(infl)
				infl=os.path.join(fldrNm,'depot.d',flNm)
				try:
					assert os.path.exists(infl)
					for lne in open(infl,'r'):
						if lne.strip() and lne[0]!='#':
							lne = lne.splitlines()[0].split('\t')
							if lne[gnClmn] not in sGns:
								if lne[gnClmn] in sGnsHg or lne[gnClmn] \
								in sGnsMm:
									print '%s %s'%(lne[gnClmn],infl)
								dSppNtPrsntGns[spp].add(lne[gnClmn])
				except:
					print 'Exception: %s'%infl
	"""
	Due to the large number of genes not present in the GRCm38_p6 
	annotation, the "autmzOldToNwNm" script was used to obtain the 
	GRCm38_p6_ENSMBLgnNm  versions of the gene symbols/names. The names 
	were then replace in the original database using the file 
	"ptchToCnvrtToGRCm38_p6.tsv" generated by the script.
	
	Gene names/symbols in "ptchToCnvrtToGRCm38_p6.tsv" need to be 
	corrected in the original files. So the original files should be 
	moved to files with the ".ori" sufix, the gene symbols corrected and 
	saved in files without the ".ori" sufix, and finally located in a 
	"depot.d" within the "location" folder included in the 
	"index.enrchmnts.txt" file. 

	After correcting the names, a second run of the script must have as
	result: dSppNtPrsntGns=={'Mouse': set([]),'Human': set([])}
	"""

######################
######################
if autmzOldToNwNm:
	#
	import requests
	server = "https://jul2019.rest.ensembl.org"#jul2019 database of GRCm38_p6
	#	
	def rtrnd(infl,clmnK,clmnV):
		dKV={}
		for l in open(infl,'r'):
			if l.strip():
				l=l.splitlines()[0].split('\t')
				k=l[clmnK]
				v=l[clmnV]
				if dKV.has_key(k):
					dKV[k].add(v)
				else:
					dKV[k]=set([v])
		dKV=dict([(k,'|'.join(v)) for k,v in dKV.items()])
		return dKV
	#
	outFlOldNew='data/ptchToCnvrtToGRCm38_p6.tsv'
	dENSMBLtoGnNm = dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(inFl_mm10p6).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene stable ID'])
	#
	sGnsMm=set(['Selm','1700049G17Rik','2700050L05Rik','Gm101','Selk','Ssfa2','Adck4','BC033916', \
	'2310042D19Rik','H2-Ke2','BC061194','3010026O09Rik','F930015N05Rik','C77370','Gm13476','Ngfrap1', \
	'Fam132b','4632415K11Rik','Fam179b','1500016L03Rik','Fam115c','BC056474','2210417K05Rik','Mfsd7b', \
	'2010107G23Rik','A530054K11Rik','Dak','Fam73a','Fam73b','1700106J16Rik','Lphn1','Hmha1','Nov', \
	'Fam213a','Wbp7','Wdr96','Wbp5','4930429B21Rik','2510003E04Rik','1700019G17Rik','Lyrm5','Rab7l1', \
	'2900002K06Rik','2810442I21Rik','Gm6531','Adck3','Gnb2l1','4930415F15Rik','Sqrdl','Ankrd32','Gm561', \
	'AK129341','Fam101b','Fam101a','Gm10474','Garem','Acpl2','LOC106740','1700012B15Rik','Gm3833', \
	'5730559C18Rik','Stk30','BC020402','6720401G13Rik','2210404J11Rik','Ptchd2','4930471I20Rik', \
	'AI464131','Gm13845','Rab1','Vimp','Gm4532','Wbscr17','Dus2l','Mllt4','Gm16119','2700089E24Rik', \
	'1600016N20Rik','Fam35a','Ptrf','Fam213b','2410127L17Rik','9130011E15Rik','Gm6402','4933426M11Rik', \
	'Fam86','Igj','Fam203a','AI848285','AI414108','Lins','Ept1','2210018M11Rik','Raet1a','1200014J11Rik', \
	'4931408D14Rik','Gpr98','Gm3230','Gbas','Spink3','1700025E21Rik','Uqcc','Ppap2c','Fam65b','Fam65a', \
	'Gm13251','Gm13154','Gm13157','Cxx1b','Gm13152','AI118078','Fam115a','1700009P17Rik','5033411D12Rik', \
	'Xrcc6bp1','C130030K03Rik','Gm11974','Phf16','Skiv2l2','C230081A13Rik','4732415M23Rik','Whsc1l1', \
	'5031426D15Rik','Tex40','Tmem180','0610009D07Rik','Papd7','Papd4','Papd5','Diap3','Diap2','Diap1', \
	'1600002H07Rik','Fam188b','Wdr20a','BC021891','Zfp85-rs1','Naprt1','Cybasc3','D4Wsu53e','Fbxo18', \
	'Ccdc37','Ctgf','Ccrn4l','B630005N14Rik','AI314180','C330006A16Rik','l7Rn6','B3gnt1','Sepp1', \
	'4932415G12Rik','Solh','1700085B03Rik','Gm14378','Mkl1','Hnrpll','Rtdr1','Agpat6','Aaed1','Ifi27l1', \
	'Gm4975','Agpat9','Wdr16','Dom3z','D2Wsu81e','4921530D09Rik','Gm16039','2010001M06Rik','Dscr3', \
	'Fam46a','3110035E14Rik','3110002H16Rik','Fam46b','C030016D13Rik','D830031N03Rik','Efcab4b', \
	'6330408A02Rik','9330133O14Rik','Fam134c','Fam134b','Fam134a','Ccdc64','Ccdc67','Pcnxl4','Pcnxl3', \
	'A330050B17Rik','5730408K05Rik','2610524H06Rik','2610034B18Rik','4930455F23Rik','Fam212b','BC039771', \
	'1700022P22Rik','Gm5506','Tmem194b','1700024P03Rik','4930451C15Rik','Oraov1','BB283400','Gm11985', \
	'Agphd1','1700003M02Rik','4921525B02Rik','Fam194a','2010107E04Rik','Ccdc129','Fam46c','Gm6251', \
	'2310045N01Rik','Dlx6as2','Gm14872','Gm14873','Mir5115','Gyltl1b','Figf','C030039L03Rik', \
	'1700007G11Rik','Aim1l','4632428N05Rik','Ctage5','C920025E04Rik','Dlx6as1','Gpr125','AV051173', \
	'Fam188a','Wapal','4921506M07Rik','Gpr115','Prmt10','D19Ertd737e','Smcr7','Lsmd1','Gm5803','Tmem35', \
	'BC021785','Gtpbp5','D630032N06Rik','Pvrl3','Plk1s1','Soga2','Pvrl4','Minos1','8430427H17Rik', \
	'Klhdc5','Mir3473','4933436C20Rik','6330416G13Rik','4930469G21Rik','BC003331','Wdr85','Slc24a6', \
	'E030011O05Rik','Snhg7','9130017N09Rik','Gm4980','3110043O21Rik','4921507L20Rik','D330050I16Rik', \
	'Dirc2','2700081L22Rik','Gm9047','Tmem66','Dnahc2','Dnahc1','AI836003','Ptplb','Dnahc5','A230073K19Rik', \
	'Fam105a','4930487H11Rik','Tmem57','Mki67ip','D330045A20Rik','2900092D14Rik','4833403I15Rik', \
	'1810043H04Rik','BC017643','1600029D21Rik','Gm7325','Fam211a','Gm4371','Gm10007','D15Ertd621e', \
	'Fam109b','Fam109a','Gm12942','4930427A07Rik','Gm20199','5830415F09Rik','6720456H20Rik','D10Jhu81e', \
	'BC089491','4732456N10Rik','Gm17801','Gm14420','BC023829','9430076C15Rik','AA415398','Gm16516', \
	'Mtl5','2700029M09Rik','B430319G15Rik','AI846148','D630045M09Rik','5330426P16Rik','F630111L10Rik', \
	'Gm6642','Wdr20b','Mll2','Mrp63','Ccdc23','Gpr64','Gm16197','Obfc1','Gpr111','Chi3l1','Ccdc90a', \
	'Dopey2','5830418K08Rik','Slc22a13b-ps','Mesdc1','Mesdc2','Usmg5','2310061J03Rik','BC018507', \
	'A830010M20Rik','Wisp2','AU019823','Wisp1','Fam173b','Zcchc11','Gm4636','Hiat1','AF357355', \
	'Tsga14','Dopey1','AF357359','Gm17365','Gltpd1','0610038L08Rik','AF357426','D3Bwg0562e', \
	'Fert2','2210013O21Rik','Zbtbd6','Sep15','Eif2c3','4930506M07Rik','B130006D01Rik','Gm996', \
	'4930526I15Rik','Gm5126','4930441O14Rik','2410057H14Rik','B230216G23Rik','Rp2h','2610305D13Rik', \
	'3110057O12Rik','Loh12cr1','1190002F15Rik','Gatsl3','Gatsl2','1700052N19Rik','AI597468','Gtl3', \
	'Enthd2','Trove2','9430020K01Rik','3632451O06Rik','Gm20187','Ccdc42b','3110062M04Rik','1300002K09Rik', \
	'Fam196a','2310003H01Rik','Gpr126','2210404O09Rik','Gpr124','2610019F03Rik','LOC100504703','Gpr123', \
	'2310036O22Rik','1110057K04Rik','Wdr67','Wdr65','BC026585','Agxt2l2','2310014L17Rik','Hmga2-ps1', \
	'Ncrna00085','1110058L19Rik','A930011O12Rik','Hn1l','D14Abb1e','Athl1','Wbscr27','Dos','AA388235', \
	'Stxbp3a','2310067B10Rik','Spin2','Mycl1','Rfwd2','Fam64a','1700013G23Rik','Ncrna00086','Wbscr16', \
	'Tceb1','1500011K16Rik','4933426D04Rik','1700109F18Rik','2810008D09Rik','5730508B09Rik', \
	'2310047M10Rik','Mettl10','Ccbl2','C230052I12Rik','Shfm1','Gm608','Prune','Srcrb4d','2810055G20Rik', \
	'Glt25d1','1810033B17Rik','Fam26f','Dgcr14','Dnahc17','Fam5b','Tssc1','Nrp','Rqcd1','D7Ertd715e', \
	'1700101E01Rik','2010003O02Rik','Lnp','Gm10536','A630066F11Rik','Mina','Cecr5','Rtfdc1','D17Wsu104e', \
	'2310069G16Rik','Prosc','Gm10789','Qtrtd1','Gm514','Ccdc109b','Zfp828','Gsg2','1700015E13Rik', \
	'Mettl13','Inadl','Mum1l1','Gtdc2','1110034G24Rik','Gm129','Gm4944','Gm2382','1700034F02Rik', \
	'Pin1-ps1','Rltpr','Ccdc53','Ccdc55','Cecr6','2010107G12Rik','4930568K20Rik','BC037034', \
	'4931406H21Rik','1110054M08Rik','D330022K07Rik','2610034M16Rik','Leprel2','Leprel1','Mlf1ip', \
	'1110004E09Rik','C030046E11Rik','Nt5c3l','Hdgfrp2','Hdgfrp3','Ccdc132','Abp1','BC017612','Mterf', \
	'4930414L22Rik','Siglec5','Gltscr1l','5430435G22Rik','1700094D03Rik','2900055J20Rik','Nhp2l1', \
	'Dnahc11','Fam159b','1110007C09Rik','Ftsjd1','Fam159a','Nup62-il4i1','2700086A05Rik','Igf2as', \
	'1300018J18Rik','Dfnb59','Cml1','Casc5','Gm13826','5430428K19Rik','Mir219-2','2900008C10Rik', \
	'Mfsd7c','9830001H06Rik','Dfna5','Itfg3','E130309D14Rik','Gm20748','2610100L16Rik','Baat1', \
	'AA474331','Fam26e','Bzrap1','Dnahc8','Ubl4','Dnahc9','E130309F12Rik','Gm10677','Epb4.1l4b', \
	'A930013F10Rik','2010015L04Rik','Nadkd1','Adrbk1','Adrbk2','9230115E21Rik','Tceb3','Tceb2', \
	'Dnalc4','Ccdc164','Gyk','9430083A17Rik','Dnalc1','Hras1','Eftud1','Myeov2','2310008H09Rik', \
	'LOC100504608','Msx1as','2610203C20Rik','Dyx1c1','Rabl5','4632434I11Rik','A730085A09Rik','Smek2', \
	'Gm9079','Pion','Sgol2','Fam195b','Pacrgl','Gareml','Fam195a','Clcn4-2','Large','2900056M20Rik', \
	'2810474O19Rik','Gm15800','Gm9897','Stxbp3b','E130012A19Rik','Ldoc1l','Leprel4','Tenc1','AI450353', \
	'Apoa1bp','Gucy1b3','BC027231','Gltscr1','Gltscr2','Zfp191','Mettl20','6330403M23Rik','Zfp607', \
	'BC125332','2810408M09Rik','Lrrc48','Zufsp','C78339','2810011L19Rik','Lrrc16a','Gm2518','Ric8', \
	'Rgag4','A630020A06','Ddx26b','Gm5480','Phf15','Cd97','4931429I11Rik','AI316807','D10Bwg1379e', \
	'Fam103a1','Prkrir','1700021K19Rik','Beta-s','Ccdc19','Ppap2b','Zfyve20','Mum1','C330027C09Rik', \
	'0610011F06Rik','Hrsp12','4930405P13Rik','1200011I18Rik','Sdccag3','4932443I19Rik','AU023871', \
	'Grlf1','2410076I21Rik','Znf512b','1190003J15Rik','E430025E21Rik','6430573F11Rik','Gm684','Phf17', \
	'1100001G20Rik','Dnahc7b','A730017C20Rik','Emr4','Ccdc176','N6amt2','Dbc1','4922501C03Rik','Mgea5', \
	'Tmem2','Vprbp','A330021E22Rik','C230091D08Rik','Gm1045','2410012M07Rik','Narg2','BC068157','Fam108a', \
	'Fam108b','Mettl21d','Zfp862','Gm14005','D630037F22Rik','Prr24','A730098P11Rik','1700011J10Rik', \
	'2700094K13Rik','A230070E04Rik','Whsc1','4921530L18Rik','Bcmo1','Tmem194','1190002N15Rik','Gm11602', \
	'Csrp2bp','Ptpla','Wdr52','Cyb5','Peo1','E330033B04Rik','Gm3258','Mfsd4','1700011H14Rik', \
	'4933427G17Rik','2410015M20Rik','B230120H23Rik','2310008H04Rik','1110008J03Rik','Ispd', \
	'5730457N03Rik','Fam63a','Fam63b','BC026590','2610301G19Rik','Rsg1','Wibg','Gm1564','Gm11744','Ndnl2', \
	'2210404O07Rik','BC018242','Mtap7d3','5730577I03Rik','2210416O15Rik','Hfe2','Ppapdc1b','Gm10336', \
	'Selrc1','Wash','Azi1','Cxx1a','Cxx1c','4930528A17Rik','LOC100503496','5430417L22Rik','2610015P09Rik', \
	'Utp11l','Mterfd3','9230110C19Rik','Murc','2610018G03Rik','Jhdm1d','Gm16515','Aim1','2410004N09Rik', \
	'Hn1','4932416H05Rik','A630007B06Rik','Cldn25','6330407A03Rik','B3galtl','3110047P20Rik', \
	'2610002J02Rik','2810403A07Rik','2700060E02Rik','AK010878','Klhl17','Pex11c','Gm5595','Apitd1', \
	'Zfp71-rs1','6030419C18Rik','Fam108c','Pddc1','Gm13704','Cnih','1500032L24Rik','2610019E17Rik', \
	'Rnmtl1','1500012F01Rik','3110052M02Rik','8430403D17Rik','1110001A16Rik','9830147E19Rik', \
	'E030019B06Rik','2310044G17Rik','Cep110','A230046K03Rik','Zfml','C1s','2610204G22Rik','Ict1', \
	'4931417G12Rik','Gm13718','Erbb2ip','Cirh1a','Epb4.1l2','Epb4.1l3','Gm16532','Epb4.1l1','BC029214', \
	'Ccdc41','Epb4.1l5','9430008C03Rik','Zcchc5','Gm7120','Suv420h1','Suv420h2','Ccdc101','Epb4.9','Mnf1', \
	'Ccdc104','0610007N19Rik','3230401D17Rik','2510049J12Rik','D930015E06Rik','8430419L09Rik', \
	'2610020H08Rik','2310010M20Rik','Mb21d1','Smek1','Has2as','Gm10560','Bai1','Bai2','Bai3', \
	'9030624J02Rik','Tbrg3','Gm3925','Gm17296','Epb4.1','Gucy1a3','2310035C23Rik','1810011O10Rik', \
	'4933413G19Rik','D17Wsu92e','Lrrc16b','Fam154b','Cyr61','Asun','Wbscr22','9630033F20Rik','Ccdc111', \
	'Gm5820','AI462493','A130049A11Rik','Cxcr7','D10Wsu52e','Prkcdbp','Gm10664','Gm19424','Fam178a', \
	'D19Bwg1357e','1700012D01Rik','B930041F14Rik','D4Ertd617e','Atp5s','BC030307','Lace1','Trp53rk', \
	'A430071A18Rik','Sepn1','Scn2a1','BC049635','4930422G04Rik','Csda','Efha1','Mterfd2','Mterfd1', \
	'Efha2','Rps4y2','BC068281','Gm12060','9130206I24Rik','3110021A11Rik','0610007P14Rik','1110001J03Rik', \
	'2410089E03Rik','Ccdc94','Dnahc7a','D8Ertd82e','Gm5918','1810026J23Rik','Fam212a','Tmem48','AW549542', \
	'Gm5088','Ufd1l','Gm13247','Fam198a','Fam19a4','Fam198b','Gm6484','4930529M08Rik','BC016423', \
	'Fam132a','4930529F24Rik','9030617O03Rik','2210015D19Rik','AI480653','1110037F02Rik','Pet112','Mut', \
	'Slmo2','Ttc18','Carkd','Slmo1','Mkl2','2410016O06Rik','Gm17769','Gm17762','Mir143hg','B930003M22Rik', \
	'Tmem27','2810428I15Rik','Taf4a','Gm5607','Ppp2r4','BC048502','Acn9','Lepre1','2010012O05Rik','Mira', \
	'D630013N20Rik','BC005764','Sepw1','Dpcr1','Hiatl1','2410066E13Rik','Fam96a','Fam96b','Gm711','Lect1', \
	'Ppapdc1a','Pvrl1','Deb1','Gm19897','5930403L14Rik','1700026L06Rik','Fam21','Lphn3','Gm166','Park2', \
	'Lphn2','Fam175a','Fam175b','Ftsj2','Fam60a','AA987161','Sgol1','Mll3','8430410A17Rik','Mll1', \
	'Tmem55b','9430016H08Rik','Adc','Mll5','Nim1','8430408G22Rik','3110007F17Rik','Tmem110', \
	'5730422E09Rik','0610037L13Rik','1700016L04Rik','Gm3414','4833424O15Rik','4930524L23Rik','Ccrl1', \
	'0610009O20Rik','Ccdc79','4933403G14Rik','BC053749','1700019L03Rik','Wisp3','Nat6','1700045I19Rik', \
	'Ptplad2','1500002O10Rik','2410018M08Rik','Sdpr','Heatr2','Hnrpdl','2810417H13Rik','1810041L15Rik', \
	'Ppapdc2','Pvrl2','1700024F13Rik','Fam150a','Fam105b','4930444A02Rik','Fam150b','Mir5109','Ccbl1', \
	'Epb4.1l4a','Wdr8','Zfp259','Sfrs18','1700084C01Rik','1810013D10Rik','A030009H04Rik','A830080D01Rik', \
	'AI854517','C030030A07Rik','Ppap2a','AI314831','4921511H03Rik','Rn45s','Ftsjd2','2810407C02Rik', \
	'C630020P19Rik','Efcab4a','Tmem5','1110035M17Rik','LOC100861615','0610031J06Rik','Gm4827', \
	'6330403A02Rik','Bloc1s2a','Lrdd','2410137F16Rik','Ptplad1','Gpr44','3110001D03Rik','Glt25d2', \
	'Bet3l','4933411K20Rik','E130112N10Rik','Setd8','Gm7173','Fam58b','Spt1','Prosapip1','BC030336', \
	'4930538K18Rik','Atp5sl','BC048609','Darc','Ikbkap','Zcchc6','4930500J02Rik','Cpsf3l','Gm9199', \
	'Lrrc33','Smcr7l','Cml2','Gm9054','Vwa9','BC031361','Etohi1','1110018G07Rik','Cdk3-ps', \
	'9530091C08Rik','Gm216','2010111I01Rik','Fam19a5','BC022687','Sc4mol','Gm6194','Bre','Dpcd', \
	'4933422H20Rik','0610012H03Rik','2610017I09Rik','Ppapdc3','Eif2c4','A230056J06Rik','Atpbd4','Eif2c1', \
	'1110008F13Rik','Eif2c2','Fam69c','Fam69b','Fam69a','Tmem55a','Gm1661','1700112E06Rik', \
	'9430023L20Rik','Stra13','Gpr56'])
	#
	dOldNwNm = {}
	#
	for gnNm in sGnsMm:
		ext = "/xrefs/symbol/mus_musculus/%s?xref"%gnNm
		resp = requests.get(server+ext, headers={ "Content-Type" : \
		"application/json"})
		if not resp.ok:
			nwGnNm=None
		else:
			decoded = resp.json()
			if decoded:
				nwGnNm=str(decoded[0][u'id'])
			else:
				nwGnNm=None
		print gnNm,nwGnNm,dENSMBLtoGnNm.get(nwGnNm,None)
		dOldNwNm[gnNm]=dENSMBLtoGnNm.get(nwGnNm,None)
	#
	ooutFlOldNew = open(outFlOldNew,'w')
	ooutFlOldNew.write('Gene name to patch\tGene name in GRCm38_p6\n')
	ooutFlOldNew.write('\n'.join(['%s\t%s'%(gnNm,dOldNwNm[gnNm]) for gnNm \
	in sorted(dOldNwNm.keys())]))
	ooutFlOldNew.close()
	#NOTE: to avoid inconsistencies with already annotated genes, the
	#annotations of 1700012D01Rik, 2610301G19Rik, B3gnt1, F630111L10Rik, 
	#Gm3414, Gm4532, Pacrgl, Wbp7 1700016L04Rik, 1700024F13Rik, 1700034F02Rik, 
	#1700094D03Rik, 2210416O15Rik, 2310069G16Rik, 2610203C20Rik, 2610524H06Rik, 
	#2900008C10Rik, 9430076C15Rik, A230056J06Rik, A830010M20Rik, Abp1, Aim1, 
	#BC039771, Baat1, C130030K03Rik, D630013N20Rik, E130112N10Rik, Gm2518, 
	#Klhl17, Lsmd1, Nrp, Raet1a, Spt1 and Nov were forced to None



"""
Finally, the names of the genes are going to be changed from GRCh38_p12/GRCm38_p6 
to genecodeV28Comp/genecodeV18Comp in an effort to make them comparable with
our results. The next script corrects the enrichment files.
"""
######################
######################
if crrctIndxEnrchmntsFls:
	#p12/p6 versions to genecodeV28Comp/genecodeV18Comp files
	outFlCnvrt_mm10 = 'data/GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv'
	outFlCnvrt_hg38 = 'data/GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv'
	#Get ENSEMBL code for the p12/p6 annotations
	dhg38p12_gnCdV28 = dict([(l.split('\t')[0],l.split('\t')[1]) for l \
	in open(outFlCnvrt_hg38).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	dmm10p6_gnCdV18 =  dict([(l.split('\t')[0],l.split('\t')[1]) for l \
	in open(outFlCnvrt_mm10).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	sSelf_mm10 = set(dmm10p6_gnCdV18.values())
	sSelf_hg38 = set(dhg38p12_gnCdV28.values())
	#
	dSppGnsCnvrt = {'Human':dhg38p12_gnCdV28,'Mouse':dmm10p6_gnCdV18}
	dSppSself = {'Human':sSelf_hg38,'Mouse':sSelf_mm10}
	#
	inFlEnrchmnt = 'index.enrchmnts.txt'
	for l in open(inFlEnrchmnt,'r'):
		if l.strip() and l[0]!='#':
			sGnLst = set()
			l=l.splitlines()[0]
			infl = l.split('\t')[4]
			gnClmn = int(l.split('\t')[2])
			spp = l.split('\t')[5]
			dGnsCnvrt = dSppGnsCnvrt[spp]
			sGnsInSpp = dSppGns[spp]
			sGnsSelf = dSppSself[spp]
			dbTyp = l.split()[0]
			#	
			oInFl = open(infl,'w')
			fldrNm,flNm = os.path.split(l.split('\t')[4])
			inflDpt=os.path.join(fldrNm,'depot.d',flNm)
			assert os.path.exists(inflDpt)
			for lne in open(inflDpt,'r'):
				if lne.strip():
					if lne[0]!='#':
						lne = lne.splitlines()[0].split('\t')
						gnNm = lne[gnClmn]
						if dbTyp=='Self':#For self-DBs
							if gnNm in sGnsSelf:
								oInFl.write('%s\n'%'\t'.join(lne))
							else:
								sGnLst.add(gnNm)
								sGnsInSpp.add(gnNm)
						else:#For databases
							if dGnsCnvrt.has_key(gnNm):
								lne[gnClmn]=dGnsCnvrt[gnNm]#This converts
								#from GRCh38_p12/GRCm38_p6 to 
								#genecodeV28Comp/genecodeV18Comp
								oInFl.write('%s\n'%'\t'.join(lne))
							else:
								sGnLst.add(gnNm)
					else:
						oInFl.write(lne)
			oInFl.close()
			if sGnLst:
				print '\nGenes not found in file %s: %s'%(infl, \
				','.join(sorted(sGnLst)))
				print 'Genes not found in file %s and in p6/p12 annotations: %s'% \
				(infl,','.join(sorted(sGnLst.difference(sGnsInSpp))))
			#Results next
					


"""
Finally, the names of the genes are going to be changes from GRCh38_p12/GRCm38_p6 
to genecodeV28Comp/genecodeV18Comp in an effort to make them comparable with
our results. The next script corrects the correlation files.
"""
######################
######################
if crrctIndxCrrltnsFls:
	#p12/p6 versions to genecodeV28Comp/genecodeV18Comp files
	outFlCnvrt_mm10 = 'data/GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv'
	outFlCnvrt_hg38 = 'data/GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv'
	#Get ENSEMBL code for the p12/p6 annotations
	dhg38p12_gnCdV28 = dict([(l.split('\t')[0],l.split('\t')[1]) for l \
	in open(outFlCnvrt_hg38).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	dmm10p6_gnCdV18 =  dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(outFlCnvrt_mm10).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	sSelf_mm10 = set(dmm10p6_gnCdV18.values())
	sSelf_hg38 = set(dhg38p12_gnCdV28.values())
	#
	dSppGnsCnvrt={'Human':dhg38p12_gnCdV28,'Mouse':dmm10p6_gnCdV18}
	dSppSself = {'Human':sSelf_hg38,'Mouse':sSelf_mm10}
	#
	inFlCrrltn = 'index.crrltns.txt'
	for l in open(inFlCrrltn,'r'):
		if l.strip() and l[0]!='#':
			sGnLst = set()
			l=l.splitlines()[0]
			infl = l.split('\t')[3]
			gnClmn = 0
			spp = l.split('\t')[4]
			dGnsCnvrt = dSppGnsCnvrt[spp]
			sGnsInSpp = dSppGns[spp]
			sGnsSelf = dSppSself[spp]
			dbTyp = l.split()[0]
			#For databases	
			oInFl = open(infl,'w')
			fldrNm,flNm = os.path.split(l.split('\t')[3])
			inflDpt=os.path.join(fldrNm,'depot.d',flNm)
			# ~ print inflDpt
			assert os.path.exists(inflDpt)
			for lne in open(inflDpt,'r'):
				if lne.strip():
					if lne[0]!='#':
						lne = lne.splitlines()[0].split('\t')
						gnNm = lne[gnClmn]
						if dbTyp=='Self':#For self-DBs
							if gnNm in sGnsSelf:
								oInFl.write('%s\n'%'\t'.join(lne))
							else:
								sGnLst.add(gnNm)
								sGnsInSpp.add(gnNm)
						else:#For databases
							if dGnsCnvrt.has_key(gnNm):
								lne[gnClmn]=dGnsCnvrt[gnNm]#This converts
								#from GRCh38_p12/GRCm38_p6 to 
								#genecodeV28Comp/genecodeV18Comp
								oInFl.write('%s\n'%'\t'.join(lne))
							else:
								sGnLst.add(gnNm)
					else:
						oInFl.write(lne)
			oInFl.close()
			if sGnLst:
				print 'Genes not found in file %s: %s'%(infl, \
				','.join(sorted(sGnLst)))
				print 'Genes not found in file %s and in p6/p12 annotations: %s'% \
				(infl,','.join(sorted(sGnLst.difference(sGnsInSpp))))


"""
Next test if the differences in between the depot and the final files 
correspond to genes not present in GRCh38_p12/GRCm38_p6 but present
in genecodeV28Comp/genecodeV18Comp
"""
######################
######################
if tstIndxEnrchmntsFls:		
	#p12/p6 versions to genecodeV28Comp/genecodeV18Comp files
	outFlCnvrt_mm10 = 'data/GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv'
	outFlCnvrt_hg38 = 'data/GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv'
	#Get ENSEMBL code for the p12/p6 annotations
	dhg38p12_gnCdV28 = dict([(l.split('\t')[0],l.split('\t')[1]) for l \
	in open(outFlCnvrt_hg38).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	dmm10p6_gnCdV18 =  dict([(l.split('\t')[0],l.split('\t')[1]) for l \
	in open(outFlCnvrt_mm10).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	sSelf_mm10 = set(dmm10p6_gnCdV18.values())
	sSelf_hg38 = set(dhg38p12_gnCdV28.values())
	sSelf_mm10_p6 = set(dmm10p6_gnCdV18.keys())
	sSelf_hg38_p12 = set(dhg38p12_gnCdV28.keys())
	#
	dSppGnsCnvrt={'Human':dhg38p12_gnCdV28,'Mouse':dmm10p6_gnCdV18}
	dSppSself = {'Human':sSelf_hg38,'Mouse':sSelf_mm10}
	dSppSself_ps = {'Human':sSelf_hg38_p12,'Mouse':sSelf_mm10_p6}
	#
	inFlEnrchmnt = 'index.enrchmnts.txt'
	for l in open(inFlEnrchmnt,'r'):
		if l.strip() and l[0]!='#':
			sGnLst = set()
			l=l.splitlines()[0]
			infl = l.split('\t')[4]
			gnClmn = int(l.split('\t')[2])
			spp = l.split('\t')[5]
			dGnsCnvrt = dSppGnsCnvrt[spp]
			sGnsInSpp = dSppGns[spp]
			sGnsSelf = dSppSself[spp]
			sGnsSelf_ps = dSppSself_ps[spp]
			dbTyp = l.split()[0]
			#For final databases
			fnlNms=set([lne.splitlines()[0].split('\t')[gnClmn] for lne \
			in open(infl,'r') if lne.strip() and lne[0]!='#'])
			#For initial databases
			fldrNm,flNm = os.path.split(infl)
			inflDpt=os.path.join(fldrNm,'depot.d',flNm)
			if dbTyp=='Self':#For self-DBs1, where differences are due to depecrated names 	
				oriNms=set([lne.splitlines()[0].split('\t')[gnClmn] for \
				lne in open(inflDpt,'r') if lne.strip() and lne[0]!='#'])
			else:#For self-DBs1, where differences are due to names changes between GRCh38_p12/GRCm38_p6 and genecodeV28Comp/genecodeV18Comp codes.
				oriNms=set([dGnsCnvrt[lne.splitlines()[0].split('\t')[gnClmn]] \
				for lne in open(inflDpt,'r') if lne.strip() and \
				dGnsCnvrt.has_key(lne.splitlines()[0]. \
				split('\t')[gnClmn]) and lne[0]!='#'])
			#Test for differences
			print 'For file %s:'%flNm
			print 'Number of genes in the original file and the transformed file, respectively: %s %s'% \
			(len(oriNms),len(fnlNms))
			print 'Number of genes different in the original file and the transformed file: %s\n'% \
			len(oriNms.difference(fnlNms))
			try:
				assert not fnlNms.difference(oriNms)
			except:
				import pdb
				pdb.set_trace()



"""
Next test if the differences in between the depot and the final files 
correspond to genes not present in GRCh38_p12/GRCm38_p6 but present
in genecodeV28Comp/genecodeV18Comp
"""
######################
######################
if tstIndxCrrltnsFls:
	#p12/p6 versions to genecodeV28Comp/genecodeV18Comp files
	outFlCnvrt_mm10 = 'data/GRCm38_p6_to_genecodeV18Comp_ENSMBLgnNm.tsv'
	outFlCnvrt_hg38 = 'data/GRCh38_p12_to_genecodeV28Comp_ENSMBLgnNm.tsv'
	#Get ENSEMBL code for the p12/p6 annotations
	dhg38p12_gnCdV28 = dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(outFlCnvrt_hg38).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	dmm10p6_gnCdV18 =  dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
	open(outFlCnvrt_mm10).read().splitlines() if l.strip() and \
	l.split('\t')[0]!='Gene name' and l[0]!='#'])
	sSelf_mm10 = set(dmm10p6_gnCdV18.values())
	sSelf_hg38 = set(dhg38p12_gnCdV28.values())
	#
	dSppGnsCnvrt={'Human':dhg38p12_gnCdV28,'Mouse':dmm10p6_gnCdV18}
	dSppSself = {'Human':sSelf_hg38,'Mouse':sSelf_mm10}
	#
	inFlCrrltn = 'index.crrltns.txt'
	for l in open(inFlCrrltn,'r'):
		if l.strip() and l[0]!='#':
			sGnLst = set()
			l=l.splitlines()[0]
			infl = l.split('\t')[3]
			gnClmn = 0
			spp = l.split('\t')[4]
			dGnsCnvrt = dSppGnsCnvrt[spp]
			sGnsInSpp = dSppGns[spp]
			sGnsSelf = dSppSself[spp]
			dbTyp = l.split()[0]
			#For final databases
			fnlNms=set([lne.splitlines()[0].split('\t')[gnClmn] for lne \
			in open(infl,'r') if lne.strip() and lne[0]!='#'])
			#For initial databases
			fldrNm,flNm = os.path.split(infl)
			inflDpt=os.path.join(fldrNm,'depot.d',flNm)
			if dbTyp=='Self':#For self-DBs1, were differences are due to depecrated names 	
				oriNms=set([lne.splitlines()[0].split('\t')[gnClmn] for \
				lne in open(inflDpt,'r') if lne.strip() and lne[0]!='#'])
			else:#For self-DBs1, were differences are due to names changes between GRCh38_p12/GRCm38_p6 and genecodeV28Comp/genecodeV18Comp codes.
				oriNms=set([dGnsCnvrt[lne.splitlines()[0].split('\t') \
				[gnClmn]] for lne in open(inflDpt,'r') if lne.strip() \
				and dGnsCnvrt.has_key(lne.splitlines()[0]. \
				split('\t')[gnClmn]) and lne[0]!='#'])
			#Test for differences
			print 'For file %s:'%flNm
			print 'Number of genes in the original file and the transformed file, respectively: %s %s'% \
			(len(oriNms),len(fnlNms))
			print 'Number of genes different in the original file and the transformed file: %s\n'% \
			len(oriNms.difference(fnlNms))
			try:
				assert not fnlNms.difference(oriNms)
			except:
				import pdb
				pdb.set_trace()


