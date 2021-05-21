#!/bin/bash -l

mm10_fasta="mouse_index/STAR_indexes.d/chroms/chr1.fa mouse_index/STAR_indexes.d/chroms/chr2.fa mouse_index/STAR_indexes.d/chroms/chr3.fa mouse_index/STAR_indexes.d/chroms/chr4.fa mouse_index/STAR_indexes.d/chroms/chr5.fa mouse_index/STAR_indexes.d/chroms/chr6.fa mouse_index/STAR_indexes.d/chroms/chr7.fa mouse_index/STAR_indexes.d/chroms/chr8.fa mouse_index/STAR_indexes.d/chroms/chr9.fa mouse_index/STAR_indexes.d/chroms/chr10.fa mouse_index/STAR_indexes.d/chroms/chr11.fa mouse_index/STAR_indexes.d/chroms/chr12.fa mouse_index/STAR_indexes.d/chroms/chr13.fa mouse_index/STAR_indexes.d/chroms/chr14.fa mouse_index/STAR_indexes.d/chroms/chr15.fa mouse_index/STAR_indexes.d/chroms/chr16.fa mouse_index/STAR_indexes.d/chroms/chr17.fa mouse_index/STAR_indexes.d/chroms/chr18.fa mouse_index/STAR_indexes.d/chroms/chr19.fa mouse_index/STAR_indexes.d/chroms/chrM.fa mouse_index/STAR_indexes.d/chroms/chrX.fa mouse_index/STAR_indexes.d/chroms/chrY.fa mouse_index/STAR_indexes.d/chroms/chrYFP.fa mouse_index/STAR_indexes.d/chroms/chrICRE.fa mouse_index/STAR_indexes.d/chroms/ERCC-00002.fa mouse_index/STAR_indexes.d/chroms/ERCC-00003.fa mouse_index/STAR_indexes.d/chroms/ERCC-00004.fa mouse_index/STAR_indexes.d/chroms/ERCC-00009.fa mouse_index/STAR_indexes.d/chroms/ERCC-00012.fa mouse_index/STAR_indexes.d/chroms/ERCC-00013.fa mouse_index/STAR_indexes.d/chroms/ERCC-00014.fa mouse_index/STAR_indexes.d/chroms/ERCC-00016.fa mouse_index/STAR_indexes.d/chroms/ERCC-00017.fa mouse_index/STAR_indexes.d/chroms/ERCC-00019.fa mouse_index/STAR_indexes.d/chroms/ERCC-00022.fa mouse_index/STAR_indexes.d/chroms/ERCC-00024.fa mouse_index/STAR_indexes.d/chroms/ERCC-00025.fa mouse_index/STAR_indexes.d/chroms/ERCC-00028.fa mouse_index/STAR_indexes.d/chroms/ERCC-00031.fa mouse_index/STAR_indexes.d/chroms/ERCC-00033.fa mouse_index/STAR_indexes.d/chroms/ERCC-00034.fa mouse_index/STAR_indexes.d/chroms/ERCC-00035.fa mouse_index/STAR_indexes.d/chroms/ERCC-00039.fa mouse_index/STAR_indexes.d/chroms/ERCC-00040.fa mouse_index/STAR_indexes.d/chroms/ERCC-00041.fa mouse_index/STAR_indexes.d/chroms/ERCC-00042.fa mouse_index/STAR_indexes.d/chroms/ERCC-00043.fa mouse_index/STAR_indexes.d/chroms/ERCC-00044.fa mouse_index/STAR_indexes.d/chroms/ERCC-00046.fa mouse_index/STAR_indexes.d/chroms/ERCC-00048.fa mouse_index/STAR_indexes.d/chroms/ERCC-00051.fa mouse_index/STAR_indexes.d/chroms/ERCC-00053.fa mouse_index/STAR_indexes.d/chroms/ERCC-00054.fa mouse_index/STAR_indexes.d/chroms/ERCC-00057.fa mouse_index/STAR_indexes.d/chroms/ERCC-00058.fa mouse_index/STAR_indexes.d/chroms/ERCC-00059.fa mouse_index/STAR_indexes.d/chroms/ERCC-00060.fa mouse_index/STAR_indexes.d/chroms/ERCC-00061.fa mouse_index/STAR_indexes.d/chroms/ERCC-00062.fa mouse_index/STAR_indexes.d/chroms/ERCC-00067.fa mouse_index/STAR_indexes.d/chroms/ERCC-00069.fa mouse_index/STAR_indexes.d/chroms/ERCC-00071.fa mouse_index/STAR_indexes.d/chroms/ERCC-00073.fa mouse_index/STAR_indexes.d/chroms/ERCC-00074.fa mouse_index/STAR_indexes.d/chroms/ERCC-00075.fa mouse_index/STAR_indexes.d/chroms/ERCC-00076.fa mouse_index/STAR_indexes.d/chroms/ERCC-00077.fa mouse_index/STAR_indexes.d/chroms/ERCC-00078.fa mouse_index/STAR_indexes.d/chroms/ERCC-00079.fa mouse_index/STAR_indexes.d/chroms/ERCC-00081.fa mouse_index/STAR_indexes.d/chroms/ERCC-00083.fa mouse_index/STAR_indexes.d/chroms/ERCC-00084.fa mouse_index/STAR_indexes.d/chroms/ERCC-00085.fa mouse_index/STAR_indexes.d/chroms/ERCC-00086.fa mouse_index/STAR_indexes.d/chroms/ERCC-00092.fa mouse_index/STAR_indexes.d/chroms/ERCC-00095.fa mouse_index/STAR_indexes.d/chroms/ERCC-00096.fa mouse_index/STAR_indexes.d/chroms/ERCC-00097.fa mouse_index/STAR_indexes.d/chroms/ERCC-00098.fa mouse_index/STAR_indexes.d/chroms/ERCC-00099.fa mouse_index/STAR_indexes.d/chroms/ERCC-00104.fa mouse_index/STAR_indexes.d/chroms/ERCC-00108.fa mouse_index/STAR_indexes.d/chroms/ERCC-00109.fa mouse_index/STAR_indexes.d/chroms/ERCC-00111.fa mouse_index/STAR_indexes.d/chroms/ERCC-00112.fa mouse_index/STAR_indexes.d/chroms/ERCC-00113.fa mouse_index/STAR_indexes.d/chroms/ERCC-00116.fa mouse_index/STAR_indexes.d/chroms/ERCC-00117.fa mouse_index/STAR_indexes.d/chroms/ERCC-00120.fa mouse_index/STAR_indexes.d/chroms/ERCC-00123.fa mouse_index/STAR_indexes.d/chroms/ERCC-00126.fa mouse_index/STAR_indexes.d/chroms/ERCC-00130.fa mouse_index/STAR_indexes.d/chroms/ERCC-00131.fa mouse_index/STAR_indexes.d/chroms/ERCC-00134.fa mouse_index/STAR_indexes.d/chroms/ERCC-00136.fa mouse_index/STAR_indexes.d/chroms/ERCC-00137.fa mouse_index/STAR_indexes.d/chroms/ERCC-00138.fa mouse_index/STAR_indexes.d/chroms/ERCC-00142.fa mouse_index/STAR_indexes.d/chroms/ERCC-00143.fa mouse_index/STAR_indexes.d/chroms/ERCC-00144.fa mouse_index/STAR_indexes.d/chroms/ERCC-00145.fa mouse_index/STAR_indexes.d/chroms/ERCC-00147.fa mouse_index/STAR_indexes.d/chroms/ERCC-00148.fa mouse_index/STAR_indexes.d/chroms/ERCC-00150.fa mouse_index/STAR_indexes.d/chroms/ERCC-00154.fa mouse_index/STAR_indexes.d/chroms/ERCC-00156.fa mouse_index/STAR_indexes.d/chroms/ERCC-00157.fa mouse_index/STAR_indexes.d/chroms/ERCC-00158.fa mouse_index/STAR_indexes.d/chroms/ERCC-00160.fa mouse_index/STAR_indexes.d/chroms/ERCC-00162.fa mouse_index/STAR_indexes.d/chroms/ERCC-00163.fa mouse_index/STAR_indexes.d/chroms/ERCC-00164.fa mouse_index/STAR_indexes.d/chroms/ERCC-00165.fa mouse_index/STAR_indexes.d/chroms/ERCC-00168.fa mouse_index/STAR_indexes.d/chroms/ERCC-00170.fa mouse_index/STAR_indexes.d/chroms/ERCC-00171.fa"
gtf_mm10="mouse_index/mm10.genecodeV18Comp.ERCCYFPiCre.cfflinks.gnNms.biotyp.gtf"
genDir="mouse_index/STAR_indexes.d/mm10.genecodeV18Comp.ERCCYFPiCre.cfflinks.gnNms.biotyp.STAR.d"

STAR --runMode genomeGenerate --runThreadN 10 --genomeDir $genDir --genomeFastaFiles $mm10_fasta --sjdbGTFfile $gtf_mm10 --sjdbOverhang 40 --genomeChrBinNbits 16
