##########################
##########################
##########################
#####   BisqueRNA   ######
##########################
##########################
##########################

rtrnClstrs=function(outClstrs)
	#Method to read clustering and return clusters
	{
	dataset2=read.csv(outClstrs, header=T, sep="	")
	cell=dataset2[,1]
	cell=as.character(cell)
	dataset=data.frame(dataset2[,-1])
	rownames(dataset)=cell
	return(t(dataset))
	}

runBisqueRNAStndAln=function(cnts_inBlk_fl,cnts_SC_inFl,clstrs_SC_inFl, 
	outFlPrprtnsBisqueRNA,sep=',')
	{
	library(BisqueRNA)
	library(Biobase)
	#Load bulk
	cnts_blk_sampleID = rtrnCntsFrmFl(cnts_inBlk_fl)
	cnts_blk = cnts_blk_sampleID$gene
	#Load reference SC
	cnts_ref_sampleID = rtrnCntsFrmFl(cnts_SC_inFl,sep=sep)
	cnts_ref = cnts_ref_sampleID$gene
	ar_clstrsOri = rtrnClstrs(clstrs_SC_inFl)
	smplNames_pData_trngData = rtrnTrngDtStrctr(cnts_SC_inFl,clstrs_SC_inFl)
	smplNames = smplNames_pData_trngData$smplNames
	pData = smplNames_pData_trngData$pData
	trngData = smplNames_pData_trngData$trngData
	#
	bulk.eset <- Biobase::ExpressionSet(assayData = cnts_blk)
	sample.ids <- colnames(cnts_ref[,match(colnames(ar_clstrsOri),
	colnames(cnts_ref))])
	sc.pheno <- data.frame(check.names=F, check.rows=F,stringsAsFactors=F,
	row.names=sample.ids,SubjectName=c(smplNames),cellType=c(ar_clstrsOri))
	#
	sc.meta <- data.frame(labelDescription=c("SubjectName","cellType"),
	row.names=c("SubjectName","cellType"))
	sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno,varMetadata=sc.meta)
	#
	sc.eset <- Biobase::ExpressionSet(assayData=cnts_ref[,
	match(colnames(ar_clstrsOri),colnames(cnts_ref))],phenoData=sc.pdata)
	res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset,sc.eset, 
	markers=NULL,use.overlap=FALSE)
	#
	ref.based.estimates <- t(res$bulk.props)
	#Write matrices
	write.csv(ref.based.estimates,outFlPrprtnsBisqueRNA)
	}

rtrnCntsFrmFl=function(data_input, sep="	")
	{
	dataset2=read.csv(data_input,header=T,sep=sep,check.names=FALSE)
	genename=dataset2[,1]
	genename=as.character(genename)
	dataset=dataset2[,-1]
	rownames(dataset)=genename
	theColNames=colnames(dataset)	
	gene = as.matrix(sapply(dataset,as.numeric))
	colnames(gene)=theColNames
	rownames(gene)=rownames(dataset)
	#
	sampleID = as.character(lapply(strsplit(colnames(gene), split="-[0]"),head, n=1))
	#
	return(list('gene'=gene,'sampleID'=sampleID))
	}


rtrnTrngDtStrctr=function(cnts_SC_inFl,clstrs_SC_inFl,sep=',')
	{
	#Load reference SC
	cnts_ref_sampleID = rtrnCntsFrmFl(cnts_SC_inFl,sep=sep)
	cnts_ref = cnts_ref_sampleID$gene
	sampleID = cnts_ref_sampleID$sampleID
	ar_clstrsOri = rtrnClstrs(clstrs_SC_inFl)
	#assayData:
	smplNames = as.character(lapply(strsplit(as.character(colnames(ar_clstrsOri)), 
	split="\\.[A-Za-z]."),head, n=1))
	pData = data.frame(row.names=colnames(cnts_ref)[match(colnames(ar_clstrsOri), 
	colnames(cnts_ref))],sampleID=c(smplNames),cellType=c(ar_clstrsOri))
	trngData = ExpressionSet(cnts_ref[,match(colnames(ar_clstrsOri), 
	colnames(cnts_ref))], phenoData=AnnotatedDataFrame(data=pData))
	#
	return(list('smplNames'=smplNames,'pData'=pData,'trngData'=trngData))
	}


cnts_inBlk_fl="172_NB.allGns.noKIFIsfrms.counts.htsq.multi.csv"
cnts_SC_inFl="NB.noMitchnd.KEql10.cnts.csv"
clstrs_SC_inFl="NB.noMitchnd.KEql10.clstrsLbld.tsv"
outFlPrprtnsBisqueRNA="172_NB.BisqueRNA.prptns.csv"
runBisqueRNAStndAln(cnts_inBlk_fl,cnts_SC_inFl,clstrs_SC_inFl,outFlPrprtnsBisqueRNA)

#~ Decomposing into 10 cell types.
#~ Using 42215 genes in both bulk and single-cell expression.
#~ Converting single-cell counts to CPM and filtering zero variance genes.
#~ Filtered 1676 zero variance genes.
#~ Converting bulk counts to CPM and filtering unexpressed genes.
#~ Filtered 383 unexpressed genes.
#~ Generating single-cell based reference from 3212 cells.

