#!/bin/bash -l

#The folllowing will execute in all databases
rn_inAllDB=TRUE

#Execute
python db_pop_marker_enrichment.del.py -D=$rn_inAllDB
