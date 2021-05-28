#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  counts_to_h5ad.py
#  
#  Copyright 2019 oscar <oscarbed@ki.se>
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


##################################
# Work on scanpy-based analysis  #                
##################################

import argparse,os

from numpy import array
from singlecell.scanpy_mod import mk_adata_frm_csv,addPAGODAtSNE

#Assert is running on python 3
import sys
assert sys.version_info.major>2#Assert is running on python 3

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
 
#input parameters  
parser.add_argument('-r', '--raw_cnts',default=None,help='csv file with raw counts')
parser.add_argument('-c', '--clstrs_PAGODAnrmlzd',default=None,help='tsv file with cells and clusters i.e. "cell1\tcluster1\ncell2\tcluster2\n..."')
parser.add_argument('-d', '--depot_prefix',default=None,help='prefix for depot file')
#Files to make dictionaries
parser.add_argument('-s', '--flSmplStg',default=None,help='tsv file with sample names and stages i.e. "Sample_name1\tStage1\nSample_name2\tStage2\n..."')
parser.add_argument('-O', '--flSmplOutcm',default=None,help='tsv file with sample names and outcomes i.e. "Sample_name1\tOutcome1\nSample_name2\tOutcome2\n..."')
#Optional PAGODA embeding results
parser.add_argument('-P','--pagoda_tSNEFl',default=None,help='PAGODA tSNE file')
parser.add_argument('-F','--flippingTsne',default=False,type=str2bool,help='If True will flip vertically the tSNE embedding')
parser.add_argument('-H','--flippingHrzntlTsne',default=False, type=str2bool,help='If True will flip horizontally the tSNE embedding')
parser.add_argument('-b','--w_bscs',default=True, type=str2bool,help='If True will write a h5ad file with the basics counts and PAGODA embedding')

args = parser.parse_args()

if args:
	#input paramenters
	raw_cnts = args.raw_cnts
	clstrs_PAGODAnrmlzd = args.clstrs_PAGODAnrmlzd
	depot_prefix = args.depot_prefix
	#Files to make dictionaries
	flSmplStg = args.flSmplStg
	flSmplOutcm = args.flSmplOutcm
	#Optional PAGODA embeding results
	pagoda_tSNEFl = args.pagoda_tSNEFl
	flippingTsne = args.flippingTsne
	flippingHrzntlTsne = args.flippingHrzntlTsne
	w_bscs = args.w_bscs
	

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

#For adrnlGlnd
dCll_clstrs_PAGODAnrmlzd = dict([(l.split()[0],l.split()[1]) for l in \
open(clstrs_PAGODAnrmlzd).read().splitlines()[1:] if l.strip()])

#Stages, expected outcome, cell cycle genes
dSmplStg = dict([l.splitlines()[0].split('\t') for l in open(flSmplStg,'r') \
if l.strip()])
dSmplOutcm = dict([l.splitlines()[0].split('\t') for l in open(flSmplOutcm,'r') \
if l.strip()])		


###########################
###########################
###########################
######   Switches   #######
###########################
###########################
###########################
#Write basic info filke
wrtBscs=0
if w_bscs:
	wrtBscs=1 # If True will write a h5ad file with the basics counts and PAGODA embedding



##########################
##########################
##########################
######   Execute   #######
##########################
##########################
##########################


####################
####################
#Make adata set
####################
####################
if wrtBscs:
	adata = mk_adata_frm_csv(raw_cnts)
	adata.obs['stage'] = \
	array([dSmplStg['.'.join(v.split('.')[:-1])] for v in \
	adata.obs_names])
	adata.obs['outcome'] = \
	array([dSmplOutcm['.'.join(v.split('.')[:-1])] for v in \
	adata.obs_names])
	adata.obs['samples'] = \
	array(['.'.join(v.split('.')[:-1]) for v in \
	adata.obs_names])
	adata.obs['PAGODA_hc'] = \
	array([dCll_clstrs_PAGODAnrmlzd.get(v,'absent_N') for v in \
	adata.obs_names])
	#Add PAGODA tSNE
	if pagoda_tSNEFl is not None:
		addPAGODAtSNE(adata,pagoda_tSNEFl,flippingTsne,addSfx='PAGODA_hc', \
		flippingHrzntlTsne=flippingHrzntlTsne)
	#Write basic shape
	if wrtBscs:
		outBscs = '%s_cnts.h5ad'%depot_prefix
		adata.write(outBscs)

