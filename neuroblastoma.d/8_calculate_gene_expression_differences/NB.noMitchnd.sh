#!/bin/bash -l

#Output folder parental. It has the results from 7_run_PAGODA and will include the output.
outFldrPrnt="NB.d"
#Number of clusters
lKmins="10"
#File 1 with sample\tgroup\tstage
smplsGrpFl="hg38_NB.tsv"
#Column in each smplsGrpFl with the sample name
smplsClmn=0
#Column in each smplsGrpFl with the group name
grpsClmn=2

python gnExprssnDffncs_mstr.py -o=$outFldrPrnt -k=$lKmins -g=$smplsGrpFl -s=$smplsClmn -G=$grpsClmn
