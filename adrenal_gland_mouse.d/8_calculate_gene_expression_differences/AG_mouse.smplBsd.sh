#!/bin/bash -l

#Output folder parental. It has the results from 7_run_PAGODA and will include the output.
outFldrPrnt="AG_mouse.d"
#Number of clusters
lKmins="19"
#File 1 with sample\tgroup\tstage
smplsGrpFl="smpl_AG_mouse.database.tsv"
#Column in each smplsGrpFl with the sample name
smplsClmn=0
#Column in each smplsGrpFl with the group name
grpsClmn=2

python gnExprssnDffncs_mstr.py -o=$outFldrPrnt -k=$lKmins -g=$smplsGrpFl -s=$smplsClmn -G=$grpsClmn
