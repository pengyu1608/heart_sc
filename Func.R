require(Seurat)
require(dplyr)
require(Matrix)
require(ggplot2)
require(reshape2)
require(gridExtra)
require(tidyr)



doPCA<-function(x,title){

	require(ggplot2)
	x_select<-x[apply(x,1,var)>0.5,]

	df_PCA<-prcomp(t(x_select))
	df_out<-as.data.frame(df_PCA$x)

	df_out$group <- factor(sapply( strsplit(as.character(row.names(df_out)), "~"), "[[", 1))
	cell<- sapply(strsplit(rownames(df_out)[1],"~"),"[[",1)

	percentage <- round(df_PCA$sdev^2 / sum(df_PCA$sdev^2) * 100, 2)
	percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

	p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))+xlab(percentage[1]) + ylab(percentage[2])+ggtitle(title)
	print(p+geom_point())
	return(df_out)

}


library_normalize <-function (count_matrix,alignedRead){

	library_size <-median(alignedRead)+mad(alignedRead)
	adjust_factor<-alignedRead/library_size

	norm_count<-round(sweep(count_matrix,2,adjust_factor,FUN="/"))

	return(norm_count)

}


doCluster<-function(x,title,corMethod){

	x_select<-x[apply(x,1,var)>0.5,]
	cluster<-hclust(as.dist(1-cor(x_select,method=corMethod)))
	plot(as.dendrogram(cluster),cex.lab=0.1,cex=0.3,main=title)
	heatmap.2(x_select,distfun=function(x) as.dist(1-cor(t(x),method="s")),margins=c(5,5),keysize=0.8,density.info="none",trace="none",col=bluered(25),ColSideColors=CelColor,cexCol=0.3,labCol="",labRow="")

}

calTPM <-function(x){
    data<-apply(x,2,function(x){x/sum(x)*10**6})
    return(data)
}

calTPM2<-function(x,merge_info){

    	aligned_read<-matrix(nrow=1,merge_info$Aligned)
	aligned_read[aligned_read ==0]<- 1
	result<-sweep(x,2,aligned_read,FUN="/")*1e6
	return(result)
}

create_cell_umi_matrix <-function(umi_file,gene_info_file,prefix){

	umi_count<-read.csv(umi_file);
	merge_gene_info<-read.csv(gene_info_file,skip=8)
	merge_gene_info<-merge_gene_info[-nrow(merge_gene_info),]
	merge_gene_info[,2]<-factor(merge_gene_info[,2],prefix,levels=sort(unique(merge_gene_info[,2])))

	plate_id<-paste(merge_gene_info[,2],merge_gene_info[,6],merge_gene_info[,7],sep="~");
	
	cell_info<-data.frame(Barcode=merge_gene_info$Barcode,Type=merge_gene_info$Type,ID=paste(prefix,plate_id,sep="_"),Aligned=merge_gene_info$Aligned)
	
	umi_count_t<-t(umi_count)
	control_barcode_column<-which(merge_gene_info$Type %in% c("Pos_Ctrl","Neg_Ctrl"))
	
	if(length(control_barcode_column)){ 
	
		cell_umi<-umi_count_t[,-c(control_barcode_column)]
		cell_info<-cell_info[-c(control_barcode_column),]} else{
		
		cell_umi<-umi_count_t
	}

	merge_gene_info<-data.frame(Sample=prefix,merge_gene_info)
	
	colnames(cell_umi)<-cell_info$ID[match(colnames(cell_umi),cell_info$Barcode)]

	result<-list(umi=cell_umi,info=cell_info,ori_info=merge_gene_info)
	
	return(result)
	#return(cell_umi)
}



doTSNE<-function(data,merge_info,method="UMI",name,pc,resolution,min_gene=1000,marker=c()){

require(Seurat)
require(dplyr)
require(Matrix)

pdf_file<-paste(name,".pdf",sep="")
pdf(file=pdf_file)

if(method =="UMI"){ 
	all_cell_TPM<-calTPM2(data,merge_info)
	all_cell_log_TPM<-log2(all_cell_TPM+1)
	} else{
		all_cell_TPM<-calTPM(data)
		#all_cell_log_TPM<-log2(all_cell_TPM/10+1)
		all_cell_log_TPM<-log2(all_cell_TPM+1)
	}

all_cell_seurat <- new("seurat", raw.data = all_cell_log_TPM)
all_cell_seurat <- Setup(all_cell_seurat, min.cells = 10, min.genes = min_gene, do.logNormalize = F, project = "sc", names.delim ="~")

all_cell_seurat <- MeanVarPlot(all_cell_seurat ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, y.cutoff =1, do.contour = F)
all_cell_seurat <- PCA(all_cell_seurat, pc.genes = all_cell_seurat@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
all_cell_seurat <- ProjectPCA(all_cell_seurat)
#PCHeatmap(all_cell_seurat, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

#all_cell_seurat <- JackStraw(all_cell_seurat, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(all_cell_seurat, PCs = 1:16)

##use standard viariation explained by PCs to choose
print(PCElbowPlot(all_cell_seurat))

##Cluster the cells

#stash origin id

all_cell_seurat <- StashIdent(all_cell_seurat, save.name = "origin_id")

all_cell_seurat <- FindClusters(all_cell_seurat, pc.use = pc, resolution = resolution, print.output = 0, save.SNN = T,temp.file.location="/tmp/",algorithm=3)

##Run Non-linear dimensional reduction (tSNE)

all_cell_seurat <- RunTSNE(all_cell_seurat, dims.use = pc, do.fast = T)

#note that you can set do.label=T to help label individual clusters
#TSNEPlot(all_cell_seurat)

# find markers for every cluster compared to all remaining cells, report only the positive ones
all_cell_seurat.markers <- FindAllMarkers(all_cell_seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

##draw heatmap of markers 

all_cell_seurat.markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name

gene_marker<- if(length(marker)>0) marker else top20$gene

DoHeatmap(all_cell_seurat, genes.use = gene_marker, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = F,cex.row=0.5,density.info="none",keysize=0.8)
plot1<-TSNEPlot(all_cell_seurat, do.return = T, group.by = "origin_id", no.legend = F, do.label = F)
plot2<-TSNEPlot(all_cell_seurat, do.return = T, no.legend = F, do.label = T)
print(plot1)
print(plot2)
MultiPlotList(list(plot1, plot2), cols = 2)
dev.off()


result<-list("seurat"=all_cell_seurat,"markers"=all_cell_seurat.markers)
return(result)

#change cluster id 
#current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
#new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
#all_cell_seurat@ident <- plyr::mapvalues(all_cell_seurat@ident, from = current.cluster.ids, to = new.cluster.ids)
#TSNEPlot(all_cell_seurat, do.label = T, pt.size = 0.5)
}

output_TSNE<-function(data,name,marker=c()){

require(Seurat)
require(dplyr)
require(Matrix)

pdf_file<-paste(name,".pdf",sep="")
pdf(file=pdf_file)

all_cell_seurat<-data$seurat
all_cell_seurat.markers<-data$markers
all_cell_seurat.markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name

gene_marker<- if(length(marker)>0) marker else top20$gene

plot1<-TSNEPlot(all_cell_seurat, do.return = T, group.by = "origin_id", no.legend = F, do.label = F)
plot2<-TSNEPlot(all_cell_seurat, do.return = T, no.legend = F, do.label = T)
#print(plot1)
#print(plot2)
MultiPlotList(list(plot1, plot2), cols = 2)
DoHeatmap(all_cell_seurat, genes.use = gene_marker, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = F,cex.row=0.5,density.info="none",keysize=0.8)
dev.off()

}

#source("/bioware/bin/SingleCell/heatmap.3.R")


format_number<- function(x){
   
   if(x >1e-4){
		x<-round(x,3)
   } else{
        x<-sprintf("%.3e",x)
   }
   return(x)
}

draw_gene_plot<-function(new_marker_tpm,class_label,class_label2){

require(ggplot2)
require(gridExtra)

graph_object<-list()
graph_object2<-list()

for (i in 1:nrow(new_marker_tpm)){

    id=rownames(new_marker_tpm)[i]
	#id=paste(id,gene_alias[match(id,rownames(gene_alias),nomatch=0),2],sep=" ")
	tmp_df<-data.frame(value=new_marker_tpm[i,],class=class_label)
	tmp_df2<-data.frame(value=new_marker_tpm[i,],class=class_label2)
    #p<-ggplot(data=tmp_df,aes(x=class,y=value,fill=class))+scale_fill_manual(values=c("gold", "darkgreen"))+geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = boxplot.stats(tmp_df$value)$stats[c(1, 5)]*1.05)+ylab("logTPM")+ggtitle(id)+guides(fill=guide_legend(title=NULL))
	
	#p<-ggplot(data=tmp_df,aes(x=class,y=value,fill=class))+scale_fill_manual(values=c("gold", "darkgreen"))+geom_violin()+geom_boxplot(width=.1)+ylab("logTPM")+ggtitle(id)+guides(fill=guide_legend(title=NULL))
	
	#wilcox.pvalue<-wilcox.test(subset(tmp_df,class=="low",select=c("value"))[,1],subset(tmp_df,class=="high",select=c("value"))[,1])$p.value
	#wilcox.pvalue<-format_number(wilcox.pvalue)
	#id<-paste(id,wilcox.pvalue)
	p<-ggplot(data=tmp_df,aes(x=class,y=value,fill=class))+geom_violin(trim = T)+geom_boxplot(width=.1,outlier.shape = NA)+ylab("Expression")+ggtitle(id)+guides(fill=guide_legend(title=NULL))
	p2<-ggplot(data=tmp_df2,aes(x=class,y=value,fill=class))+geom_violin(trim = T)+geom_boxplot(width=.1,outlier.shape = NA)+ylab("Expression")+ggtitle(id)+guides(fill=guide_legend(title=NULL))
	
    graph_object[[i]]<-p	
	graph_object2[[i]]<-p2
    
	#p2<-p+theme_bw()+theme(plot.title = element_text(hjust = 0.5), axis.text.y  = element_text(size=12),axis.text.x = element_text(size=12),axis.title.x=element_text(size=12),axis.title.y = element_text(size=12),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())
	#print(p)	
	
	}
	
	result<-list(g1=graph_object,g2=graph_object2)
	return(result)
}

preTSNE<- function(input,cell_info,name,ifNorm=T,min.cell=10,min.gene=500,mito.perc=1,x.cutoff=3,y.cutoff=1,regress_factor=c(),doJackStraw=F){

	require(Seurat)
	require(dplyr)
	require(Matrix)
	require(ggplot2)
	require(reshape2)
	require(gridExtra)
	require(tidyr)
	
	pdf_file<-paste(name,"_preTSNE_QC.pdf",sep="")
	pdf(file=pdf_file)
	
	data <- CreateSeuratObject(raw.data = input, min.cells = min.cell, min.genes = min.gene, project = "sc", names.delim ="~")
		
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = data@data), value = TRUE,ignore.case =T)
	percent.mito <- Matrix::colSums(data@raw.data[mito.genes, ])/Matrix::colSums(data@raw.data)

	input_cell_info<-cell_info[colnames(input),]

	input_cell_info$newType<-input_cell_info$Type

	input_cell_info<- input_cell_info %>% separate(newType, c("condition","sample","group"), "_")

	#sample<-input_cell_info$sample
	#condition<-input_cell_info$condition
	#group<-input_cell_info$group

	meta_data<-data.frame(input_cell_info[,c("sample","condition","group")])

	data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")
	data <- AddMetaData(object = data, metadata = meta_data, col.name = c("sample","condition","group"))
	#data <- AddMetaData(object = data, metadata = condition, col.name = "condition")
	#data <- AddMetaData(object = data, metadata = group, col.name = "group")
	
	print(VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 2))
	par(mfrow = c(1, 2))
	GenePlot(object = data, gene1 = "nUMI", gene2 = "percent.mito")
	GenePlot(object = data, gene1 = "nUMI", gene2 = "nGene")

	#Filter by subset. E.g. filter by two factors nGene and percent.mito, define low and high threshold separately 

	data <- FilterCells(object = data, subset.names = c("nGene", "percent.mito"), 
		low.thresholds = c(min.gene, -Inf), high.thresholds = c(Inf, mito.perc))
		
	#normalize the data 
	#employ a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result

	par(mfrow = c(1,1))

	if(ifNorm){

		data <- NormalizeData(object = data, normalization.method = "LogNormalize", 
		scale.factor = 10000)
	}

	data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR, 
		x.low.cutoff = 0.0125, x.high.cutoff = x.cutoff, y.cutoff = y.cutoff)
    
	var.gene.num=length(data@var.genes)
	barplot(length(data@var.genes),main=paste(var.gene.num,"genes"))
	#3466

	#ScaleData now incorporates the functionality of the function formerly known as RegressOut (which regressed out given the effects of provided variables and then scaled the residuals)
	#scale and center the data. Setting center to TRUE will center the expression for each gene by subtracting the average expression for that gene. Setting scale to TRUE will scale the expression level for each gene by dividing the centered gene expression levels by their standard deviations 

	if(length(regress_factor)) {
		data <- ScaleData(object = data, vars.to.regress = regress_factor)
	} else{
		data <- ScaleData(object = data)
	}
	
	data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:5, 
		genes.print = 5)

	#Print genes correlates with each PC, to use all genes, set use.full=T
	
	#PrintPCA(object = data, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

	#Visulize the genes with highest correlation with PC 1 and 2

	#VizPCA(object = data, pcs.use = 1:2)
	#
	data <- ProjectPCA(object = data, do.print = FALSE)

	##print PC ElbowPlot to choose significant PCs.
	print(PCElbowPlot(data))

	if (doJackStraw){
		data <- JackStraw(object = data, num.replicate = 100, do.print = FALSE)
		print(JackStrawPlot(data, PCs = 1:16))
	}

	dev.off()
	return(data)

}

doTSNE_2 <-function(data,name,pc,resolution,marker=c(),k=30,algorithm=1){

	pdf_file<-paste(name,"_TSNE.pdf",sep="")
	pdf(file=pdf_file)
	
	require(Seurat)
	
	#find Clusters

	data <- StashIdent(data, save.name = "origin_id")
	data <- FindClusters(object = data, reduction.type = "pca", dims.use = pc, resolution = resolution, print.output = 0, save.SNN = TRUE,temp.file.location="/tmp/",k.param=k,algorithm=algorithm )

	#PrintFindClustersParams(object = data)

	data <- RunTSNE(object = data, dims.use = pc, do.fast = TRUE)

	#note that you can set do.label=T to help label individual clusters
	#TSNEPlot(data)

	# find markers for every cluster compared to all remaining cells, report only the positive ones
	#data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
	data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.2, thresh.use = 0.2, return.thresh =0.05)

	##draw heatmap of markers

	data.markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20
	# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name

	gene_marker<- if(length(marker)>0) marker else top20$gene
	#gene_marker<- top20$gene

	#print(DoHeatmap(data, genes.use = gene_marker, slim.col.label = TRUE, remove.key = TRUE,col.low = "#3399FF", col.mid = "white",col.high = "#993333"),cex.row=3,cex.col=4)
	print(DoHeatmap(data, genes.use = gene_marker, slim.col.label = TRUE, remove.key = TRUE,cex.row=3,cex.col=4))
	plot1<-TSNEPlot(data, do.return = T, group.by = "condition", no.legend = F, do.label = F)
	plot2<-TSNEPlot(data, do.return = T, no.legend = F, do.label = T)
	print(plot1)
	print(plot2)
	print(plot_grid(plot1, plot2))
	dev.off()
	
	result<-list("seurat"=data,"markers"=data.markers)
	return(result)
	
}


doTSNE_3 <-function(data,name,pc,resolution,marker=c()){

	pdf_file<-paste(name,"_TSNE.pdf",sep="")
	pdf(file=pdf_file)
	
	require(Seurat)
	
	#find Clusters

	data <- StashIdent(data, save.name = "origin_id")
	data <- FindClusters(object = data, reduction.type = "pca", dims.use = pc, resolution = resolution, print.output = 0, save.SNN = TRUE,temp.file.location="/tmp/",)

	#PrintFindClustersParams(object = data)

	data <- RunTSNE(object = data, dims.use = pc, do.fast = TRUE)

	#note that you can set do.label=T to help label individual clusters
	#TSNEPlot(data)

	# find markers for every cluster compared to all remaining cells, report only the positive ones
	#data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
	data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.2, thresh.use = 0.2, return.thresh =0.05)

	##draw heatmap of markers

	#data.markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20
	# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name

	#gene_marker<- if(length(marker)>0) marker else top20$gene
	#gene_marker<- top20$gene

	#print(DoHeatmap(data, genes.use = gene_marker, slim.col.label = T, remove.key = F))
	plot1<-TSNEPlot(data, do.return = T, group.by = "origin_id", no.legend = F, do.label = F)
	plot2<-TSNEPlot(data, do.return = T, no.legend = F, do.label = T)
	print(plot1)
	print(plot2)
	print(plot_grid(plot1, plot2))
	dev.off()
	
	result<-list("seurat"=data,"markers"=data.markers)
	return(result)
	
}


create_cell_umi_matrix_mod <-function(umi_file,gene_info_file,prefix,mod_info_matrix=c()){

        umi_count<-read.csv(umi_file);
        merge_gene_info<-read.csv(gene_info_file,skip=8)
        merge_gene_info<-merge_gene_info[-nrow(merge_gene_info),]

	#the control_barcode_column is not same in merge.report.csv and genematrix.csv !!! Don't use it!

	control_barcode_column<-which(merge_gene_info$Type %in% c("Pos_Ctrl","Neg_Ctrl"))

	if(length(control_barcode_column)){

		merge_gene_info<-merge_gene_info[-c(control_barcode_column),]
	}



	ori_cell_type<-merge_gene_info[,2]
	plate_info<-mod_info_matrix[mod_info_matrix$Plate ==prefix,]

	if(length(plate_info) ){

		oriID=as.character(plate_info$oriID)
		newID=as.character(plate_info$newID)
		merge_gene_info$Type<- plyr::mapvalues(x = merge_gene_info$Type, from = oriID, to = newID)
		
	}
	
	#
	plate_id <- prefix
	row_id <- merge_gene_info[,6]
	col_id <- merge_gene_info[,7]	

        cell_id<-paste("sc",plate_id,row_id,col_id,sep="_")
        #cell_id<-paste(merge_gene_info[,2],1:length(merge_gene_info[,2]),sep="~")

	cell_id<-gsub("-","_",cell_id)

	cell_id<-toupper(cell_id)

	merge_gene_info$mito.perc<- round(merge_gene_info$Mitochondria.Read.Pairs /merge_gene_info$Assigned,2)

        merge_gene_info<-data.frame(ID=cell_id,merge_gene_info)
	rownames(merge_gene_info)<-cell_id	

        umi_count_t<-t(umi_count)

	cell_umi<-umi_count_t

        #if(length(control_barcode_column)){

        #        cell_umi<-umi_count_t[,-c(control_barcode_column)]
                #cell_info<-cell_info[-c(control_barcode_column),]} else{
        #        merge_gene_info<-merge_gene_info[-c(control_barcode_column),]} else{

        #        cell_umi<-umi_count_t
        #}

        #merge_gene_info<-data.frame(Sample=prefix,merge_gene_info)

        colnames(cell_umi)<-merge_gene_info$ID[match(colnames(cell_umi),merge_gene_info$Barcode)]

	cell_umi<-cell_umi[,!is.na(colnames(cell_umi))]

        result<-list(umi=cell_umi,info=merge_gene_info)

        return(result)
        #return(cell_umi)
}

merge_plate_data<- function(mod_info_matrix=c()){
	
	sample_list<-list.files(pattern=".*_merged.report.csv")
	
	sample_data<-list()

	for (i in 1:length(sample_list)){
	
			info_file<-sample_list[i]
			matrix_file<-paste0(sub("_merged.report.*","",info_file),"_merged.genematrix.csv")
			plate_id<-sub("_.*","",info_file)
			plate_id<-sub("-.*","",plate_id)
			#print(c(matrix_file,plate_id,sample_id))

			#modify the name of cell types 
			if (! is.null(mod_info_matrix)){

				res<-create_cell_umi_matrix_mod(matrix_file,info_file,plate_id,mod_info_matrix)
			}else{

				res<-create_cell_umi_matrix_mod(matrix_file,info_file,plate_id)
			}
			sample_data[[plate_id]]<-res

	}


	total_gene<-vector()

	for(i in 1:length(sample_list)){

			total_gene<-c(total_gene,rownames(sample_data[[i]]$umi))
	}

	cell_info_combine<-data.frame()
	for(i in 1:length(sample_list)){

			cell_info_combine<-rbind(cell_info_combine,sample_data[[i]]$info)
	}


	total_gene_uniq<-sort(unique(total_gene))

	cell_umi_combine<-data.frame(row.names=total_gene_uniq)

	for(i in 1:length(sample_list)){

			cell_umi_combine<-data.frame(cell_umi_combine,sample_data[[i]]$umi[match(total_gene_uniq,rownames(sample_data[[i]]$umi)),])
	}

	cell_umi_combine[is.na(cell_umi_combine)]<-0
	
	result<-list(umi=cell_umi_combine,info=cell_info_combine)
	
	return(result)
}


create_cell_info_mod <-function(gene_info_file,prefix,info_matrix){

        merge_gene_info<-read.csv(gene_info_file,skip=8)
        merge_gene_info<-merge_gene_info[-nrow(merge_gene_info),]
        control_barcode_column<-which(merge_gene_info$Type %in% c("Pos_Ctrl","Neg_Ctrl"))

        ori_cell_type<-merge_gene_info[,2]
        plate_info<-info_matrix[info_matrix$Plate ==prefix,]

        #new_cell_type<-plate_info[match(ori_cell_type,plate_info[,"oriID"]),"newID"]

	if(nrow(plate_info)>0){
        	
		new_cell_type<-plate_info[match(ori_cell_type,plate_info[,"oriID"]),"newID"]
        	new_cell_type<-toupper(new_cell_type)
		merge_gene_info[,2]<-factor(new_cell_type,levels=sort(unique(new_cell_type)))

	}
	

	merge_gene_info$Plate=prefix;
	#merge_gene_info[,2]<-paste(merge_gene_info[,2],prefix,sep="~")

	merge_gene_info[,2]<-gsub("-","_",merge_gene_info[,2])

	merge_gene_info[,2]<-toupper(merge_gene_info[,2])

	#merge_gene_info[,2]<-paste(merge_gene_info[,2],1:length(merge_gene_info[,2]),sep="~")

        merge_gene_info$mito.perc<- round(merge_gene_info$Mitochondria.Read.Pairs /merge_gene_info$Assigned,2)

	row_id<-merge_gene_info[,6]
	col_id<-merge_gene_info[,7]

	merge_gene_info<-data.frame(ID=paste(merge_gene_info$Plate,row_id,col_id,sep="_"),merge_gene_info)

        return(merge_gene_info)
        #return(cell_umi)
}


draw_stack_perc<- function(x,column,mod_string=c(),orderBy=c()){

		require(ggplot2)
		

		data_df<-data.frame(cluster=x$data$ident,class=x$data[,as.character(column)],value=1)

		if(!is.null(mod_string)){

		  data_df$class=sub(mod_string,"",data_df$class)
		}

        cluster_sum<-aggregate(value~cluster, data_df,sum)
		group_sum<-aggregate(value~cluster+class, data_df,sum)
		group_percent<-data.frame(group_sum,percent=0)
		
		total_class_sum<-aggregate(value~class, data_df,sum)
		total_sum<-sum(data_df$value)
		
		cluster_total<-max(as.numeric(data_df$cluster))
		
		total_class_df<-data.frame(cluster=factor("total"),total_class_sum,percent=round(total_class_sum[,"value"]/total_sum,2)*100)
		
		for (i in names(table(group_percent$cluster))){
		
			group_percent[group_percent$cluster==i,"percent"]=round(group_percent[group_percent$cluster==i,"value"]/cluster_sum[cluster_sum$cluster==i,"value"]*100,2)
		
		}

		
		group_percent<-rbind(group_percent,total_class_df)
		
		#cluster by group count summary
		p1<-ggplot(group_sum,aes(x=factor(cluster),y=value,fill=class))+ geom_bar(stat="identity")+ylab("Cell Count")+xlab("Cluster")+geom_text(aes(label=value), vjust=-0.05)

		
		#percent summary 

		if(!is.null(orderBy)){ 
	
			if(orderBy %in% group_percent$class){

				#tmp_df<-subset(group_percent,class==as.character(orderBy))
				tmp_df<-subset(group_percent,class==orderBy & cluster !="total")

				tmp_df<-tmp_df[order(tmp_df$percent,decreasing=T),]

				tmp_df$cluster<-as.character(tmp_df$cluster)

				#group_percent$cluster<-factor(group_percent$cluster,levels=c(tmp_df$cluster,"total")) }
				group_percent$cluster<-factor(group_percent$cluster,levels=c(tmp_df$cluster,"total"))}

		}

			
		p2<-ggplot(group_percent,aes(x=factor(cluster),y=percent,label=paste(percent,"%"),fill=class))+ geom_bar(stat="identity")
		
		
		#class count summary

		p3<-ggplot(total_class_sum,aes(x=factor(class),y=value,fill=class))+ geom_bar(stat="identity")+ylab("Cell Count")+xlab("Cluster")+geom_text(aes(label=value), vjust=-0.3)


	        group_percent_pie<-subset(group_percent,cluster!="total")
		group_percent_pie$cluster=factor(group_percent_pie$cluster)

		if(column =="sample"){

			#class=paste0("N",sort(as.numeric(sub("N","",levels(group_percent_pie$class)))))

			#group_percent_pie$class=factor(group_percent_pie$class,levels=class)

		#labels=paste0("cluster ",1:length(unique(group_percent_pie$cluster)))
		#names(labels)<-1:length(unique(group_percent_pie$cluster))

			p4<-ggplot(data=group_percent_pie,aes(x="",y=percent,label=paste(percent,"%"),fill=class))+ geom_bar(width=1,stat="identity") + facet_wrap(~cluster,ncol=8)+coord_polar("y", start=0)+scale_x_discrete(expand = c(0, 0), drop = TRUE,labels=NULL)+scale_y_discrete(expand = c(0, 0), drop = TRUE,labels=NULL)+theme(axis.line=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),plot.title=element_text(size=14, face="bold"),strip.background=element_blank(),text=element_text(size=8))+scale_fill_discrete(name="person")
			print(p4)
		}

 		#p3<-ggplot(total_class_sum,aes(x=factor(class),y=value,fill=class))+ geom_bar(stat="identity")
		print(p1)
		print(p2)
		print(p3)
	 
		#ggplot(group_sum,aes(x=factor(cluster),y=percent,label=paste(percent,"%"),fill=class))+ geom_bar(stat="identity",position=position_dodge())
}


write_marker<-function(x,file,mod_string=c(),orderBy=c()){

	require(Seurat)
	require(ggplot2)
	require(dplyr)
	
	#x$markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> x_top20
	#x$markers %>% group_by(cluster) %>% top_n(200, avg_diff) -> x_top200
	x$markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> x_top20
	x$markers %>% group_by(cluster) %>% top_n(500, avg_logFC) -> x_top500
	x_top500<-subset(x_top500,p_val_adj<0.05)
	file_name=paste0(file,".csv")
	file_name=paste0(file,".csv")
	write.csv(x_top500,file=file_name)
	
	pdf_name=paste0(file,".pdf")
	pdf(file=pdf_name)
	print(DoHeatmap(x$seurat, genes.use = x_top20$gene, slim.col.label = T, remove.key = T,cex.row=3))
	plot_obj<-TSNEPlot(x$seurat, do.return = T, no.legend = F, do.label = T)	
	print(plot_obj)

	plot_obj$data$sample<-x$seurat@meta.data$sample
	plot_obj$data$condition<-x$seurat@meta.data$condition
	plot_obj$data$group<-x$seurat@meta.data$group


#	plot_obj$data$sample<-sub(".*CM~","",rownames(plot_obj$data))
#	plot_obj$data$sample<-sub("~.*","",plot_obj$data$sample)
#	plot_obj$data$condition<-sub("~.*","",rownames(plot_obj$data))
#	plot_obj$data$group<-sub("N_|HF_","",rownames(plot_obj$data))
#	plot_obj$data$group<-sub("_.*","",plot_obj$data$group)


	
	new_p<-ggplot(plot_obj$data,aes(x=tSNE_1,y=tSNE_2,color=group))+ geom_point()
	new_p1<-ggplot(plot_obj$data,aes(x=tSNE_1,y=tSNE_2,color=ident,pch=group))+ geom_point()
	new_p2<-ggplot(plot_obj$data,aes(x=tSNE_1,y=tSNE_2,color=sample,pch=group))+ geom_point()

        print(font_theme(new_p))
        print(font_theme(new_p1))
	print(font_theme(new_p2))

	if(!is.na(plot_obj$data$sample[1])){

	    draw_stack_perc(plot_obj,"sample",orderBy=orderBy)
	}

	draw_stack_perc(plot_obj,"condition",mod_string,orderBy=orderBy)

	 if(!is.na(plot_obj$data$group[1])){

            draw_stack_perc(plot_obj,"group",orderBy=orderBy)
        }


	tsne_cluster_correlation(x)

	dev.off()
	return(plot_obj)
}

tsne_cluster_correlation <- function(seurat_data){

	library(gplots)
	library(RColorBrewer)
	require(reshape2)

	cluster_vargene<-seurat_data$seurat@var.genes
	cluster_data_var<-seurat_data$seurat@data[cluster_vargene,]
	cluster_df_var<-data.frame(t(as.matrix(cluster_data_var)),ident=seurat_data$seurat@ident)
	cluster_df_var_melt<-melt(cluster_df_var,id="ident")

	cluster_var_avg<-aggregate(value~ident+ variable,cluster_df_var_melt,mean)
	cluster_var_avg_cast<-acast(cluster_var_avg, formula = variable~ident)

	cor_result<-cor(cluster_var_avg_cast,method="s")
	cor_result_pearson<-cor(cluster_var_avg_cast,method="p")

	
	#pdf(file=paste0(pdf_name,".pdf"))

	heatmap.2(cor_result_pearson,notecol="black",margins=c(5,5),keysize=0.9,density.info="none",trace="none",dendrogram="row",Colv="NA",col=brewer.pal(8,"GnBu"),cexRow =0.8,cexCol=0.8,main="Mean Pearson cor")

	heatmap.2(cor_result,notecol="black",margins=c(5,5),keysize=0.9,density.info="none",trace="none",dendrogram="row",Colv="NA",col=brewer.pal(8,"GnBu"),cexRow =0.8,cexCol=0.8,main="Mean spearman cor")

	#dev.off()
	
	result<-list(cor_sp=cor_result,cor_pearson=cor_result_pearson)
	return(result)

}


draw_heatmap_label <- function(tsne_res,marker_gene=c(),label_gene=c(),title,angle=0){

	require(Seurat)
	obj_plot<-DoHeatmap(tsne_res$seurat,genes.use=marker_gene,col.low = "#3399FF", col.mid = "white",col.high = "#993333",slim.col.label = TRUE, remove.key = F,cex.row=4,title=title)

	labelRow<- rev(unique(obj_plot$data$gene))
	labelRow<-as.character(labelRow)
	labelRow[!labelRow %in% label_gene] <-""
    
	#labelRow[labelRow !=""] <- paste0("-",labelRow[labelRow !=""])


	col.low = "#3399FF"
	col.mid = "white"
	col.high = "#993333"
	group.cex = 15
	group.spacing = 0.15
	cex.col = 8
	cex.row = 4

	key.direction <- "vertical"
	key.title.pos <- "left"

	

	heatmap <- ggplot(data = obj_plot$data, mapping = aes(x = cell,
			y = gene, fill = expression)) + geom_tile() + scale_fill_gradient2(low = col.low,
			mid = col.mid, high = col.high, name = "Expression",
			guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) +
			scale_y_discrete(position = "right", labels = labelRow) +
			theme(axis.line = element_blank(), axis.title.y = element_blank(),
				axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex,angle = -90),
				axis.text.y = element_text(size = cex.row,angle=angle,face="bold"), axis.ticks.x = element_blank(), axis.text.x = element_text(size = cex.col),axis.title.x = element_blank(),text=element_text(family="ArialMT"))

	switch <- "x"
	#heatmap <- heatmap + theme(axis.title.x = element_blank(),
	#			axis.text.x = element_blank(), axis.ticks.x = element_blank(),
	#			axis.line = element_blank(), axis.title.y = element_blank(),
	#			axis.ticks.y = element_blank(),strip.text.x = element_text(angle = -90))

	heatmap <- heatmap + facet_grid(facets = ~ident, drop = TRUE,
				space = "free", scales = "free", switch = switch,
				) + scale_x_discrete(expand = c(0, 0), drop = TRUE)
	panel.spacing <- unit(x = group.spacing, units = "lines")
				
	heatmap <- heatmap + theme(strip.background = element_blank(),
				panel.spacing = panel.spacing)

	heatmap <- heatmap + labs(title = title)			
	
	return(heatmap)
	
}

detachAllPackages <- function() {

  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")

  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]

  package.list <- setdiff(package.list,basic.packages)

  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

draw_heatmap_label2 <- function(tsne_res,marker_gene=c(),label_gene=c(),title,angle=0){

	require(Seurat)
	require(ggrepel)
	require(ggplot2)
	require(cowplot)
	
	theme_set(theme_gray())
	
	obj_plot<-DoHeatmap(tsne_res$seurat,genes.use=marker_gene,col.low = "#3399FF", col.mid = "white",col.high = "#993333",slim.col.label = TRUE, remove.key = F,cex.row=4,title=title)

	labelRow<- rev(unique(obj_plot$data$gene))
	labelRow<-as.character(labelRow)
	labelRow[!labelRow %in% label_gene] <-""

	col.low = "#3399FF"
	col.mid = "white"
	col.high = "#993333"
	group.cex = 15
	group.spacing = 0.15
	cex.col = 8
	cex.row = 6

	key.direction <- "vertical"
	key.title.pos <- "top"
	
	gene_num= length(labelRow)
	
	heatmap <- ggplot(data = obj_plot$data, mapping = aes(x = cell,
			y = gene, fill = expression)) + geom_tile() + scale_fill_gradient2(low = col.low,
			mid = col.mid, high = col.high, name = "Expression",
			guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) +
			scale_x_discrete(expand = c(0, 0), drop = TRUE,labels=NULL)+scale_y_discrete(expand = c(0, 0), drop = TRUE,labels=NULL)+
			theme(axis.line = element_blank(), axis.title= element_blank(),
				axis.ticks = element_blank(), strip.text.x = element_text(size = group.cex,angle=-90,family="ArialMT"),
				axis.text.x= element_text(),axis.text.y= element_text(),legend.justification = "left",
				plot.margin = margin(5.5, 0, 5.5, 5.5, "pt")) 
				
	switch <- "x"
	
	heatmap <- heatmap + facet_grid(facets = ~ident, drop = TRUE,
				space = "free", scales = "free", switch = switch
				)
				
	panel.spacing <- unit(x = group.spacing, units = "lines")
				
	heatmap <- heatmap + theme(strip.background = element_blank(),
				panel.spacing = panel.spacing)

	heatmap <- heatmap + labs(title = title)			
	
	# make the axis plot
	axis <- ggplot(data.frame(y = 1:gene_num,
							  gene = labelRow),
				   aes(x = 0, y = y, label = gene)) +
	  geom_text_repel(#min.segment.length =0.5,
					 min.segment.length = grid::unit(0, "pt"),
					 color = "grey30",  ## ggplot2 theme_grey() axis text
					 size = 0.8*11/.pt,  ## ggplot2 theme_grey() axis text
					 fontface = "bold.italic",family="ArialMT") +
	  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
						 breaks = NULL, labels = NULL, name = NULL) +
	  scale_y_continuous(limits = c(0.5, gene_num+0.5), expand = c(0, 0),
						 breaks = NULL, labels = NULL, name = NULL) +
	  theme(panel.background = element_blank(),
			plot.margin = margin(0, 0, 0, 0, "pt"))

	# align and combine
	aligned <- align_plots(heatmap + theme(legend.position = "none"), axis, align = "h", axis = "tb")
	aligned <- append(aligned, list(get_legend(heatmap)))
	plot_obj<- plot_grid(plotlist = aligned, nrow = 1, rel_widths = c(5, .5, .7))

	return(plot_obj)
	
}

font_theme<-function(p){

	p<- p+ theme(axis.line.x = element_line(size=1.5),axis.line.y = element_line(size=1.5),axis.text=element_text(size=14,family="ArialMT"),text=element_text(family="ArialMT",face="bold",size=14),legend.position = "bottom",legend.title = element_text(size=12, face="bold"),legend.text = element_text(size = 12, face = "bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())
	
	return(p)
}

font_theme2<-function(p){

	p<- p+ theme(axis.line.x = element_line(size=1.5),axis.line.y = element_line(size=1.5),axis.text=element_text(size=16,family="ArialMT"),text=element_text(family="ArialMT",face="bold",size=16),legend.position = "bottom",legend.title = element_text(size=14, face="bold"),legend.text = element_text(size = 14, face = "bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())
	
	return(p)
}

font_theme_slim<-function(p){

	p<- p+ theme(axis.line.x = element_line(size=1),axis.line.y = element_line(size=1),axis.text=element_text(size=12,family="ArialMT"),text=element_text(family="ArialMT",face="bold",size=12),legend.position = "bottom",legend.title = element_text(size=10, face="bold"),legend.text = element_text(size = 10, face = "bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())
	
	return(p)
}

create_mark_plot<-function(data,xaxis="tSNE_1",yaxis="tSNE_2",bygroup="",colors,title,legend_title,marker){

		p<-ggplot(data, aes(x=get(xaxis), y=get(yaxis), colour=get(bygroup))) + geom_point(size=0.5) + scale_colour_gradientn(colours=colors,name=legend_title)+labs(title=title,subtitle=paste(marker,collapse=","))+xlab("tSNE 1")+ylab("tSNE 2")+theme(plot.subtitle=element_text(size=6))
		p<-font_theme_slim(p)
		return(p)
}

capwords <- function(s, strict = FALSE) {
         cap <- function(s) paste(toupper(substring(s, 1, 1)),
                       {s <- substring(s, 2); if(strict) tolower(s) else s},
                                  sep = "", collapse = " " )
         sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
     }




	
draw_marker_tsne_plot<-function(data,species="hs"){
	
	library("RColorBrewer")
	library("gridExtra")
	
	cm_known_marker<-c("TNNT2","ACTN2","MYH6","MYH7","TNNI3","MYL2","MYL2")
	vsmc_known_marker<-c("MYH11","ACTA2","CNN1","TAGLN","MYOCD","CALD1","MYLK"); 
	ec_known_marker<-c("VWF","CDH5","THBD","TEK","ENG","PECAM1"); 
	fibroblast_known_marker<-c("FN1","VIM","COL1A1","COL1A2"); 
	macrophage_known_marker<-c("CD163","S100A8","CSF1R","C5AR1","CD74"); 
	t_cell_known_marker<-c("CD3D","CD69","CD96","PTPRC","CD2")
	
	if(species !="hs"){
	
		cm_known_marker<-sapply(cm_known_marker,capwords,strict=TRUE)
		vsmc_known_marker<-sapply(vsmc_known_marker,capwords,strict=TRUE)
		ec_known_marker<-sapply(ec_known_marker,capwords,strict=TRUE)
		fibroblast_known_marker<-sapply(fibroblast_known_marker,capwords,strict=TRUE)
		macrophage_known_marker<-sapply(macrophage_known_marker,capwords,strict=TRUE)
		t_cell_known_marker<-sapply(t_cell_known_marker,capwords,strict=TRUE)
	}

	cell_cm_exp_avg<-apply(data$seurat@data[cm_known_marker,],2,mean)
	cell_vsmc_exp_avg<-apply(data$seurat@data[vsmc_known_marker,],2,mean)
	cell_ec_exp_avg<-apply(data$seurat@data[ec_known_marker,],2,mean)
	cell_fibro_exp_avg<-apply(data$seurat@data[fibroblast_known_marker,],2,mean)
	cell_mp_exp_avg<-apply(data$seurat@data[macrophage_known_marker,],2,mean)
	cell_tcell_exp_avg<-apply(data$seurat@data[t_cell_known_marker,],2,mean)

	cell_marker_exp_avg<-cbind(cell_cm_exp_avg,cell_ec_exp_avg,cell_fibro_exp_avg,cell_vsmc_exp_avg,cell_mp_exp_avg,cell_tcell_exp_avg)
	colnames(cell_marker_exp_avg)<-c("cm","ec","fc","vsmc","mp","tcell")

	cell_marker_exp_avg<-data.frame(data$seurat@dr$tsne@cell.embeddings,cell_marker_exp_avg)
	
	exp_colors<-brewer.pal(5, "Purples")

	cm_mark_plot<- create_mark_plot(cell_marker_exp_avg,bygroup="cm",colors=exp_colors,title="CM marker",legend_title="logTPM",marker=cm_known_marker)

	ec_mark_plot<- create_mark_plot(cell_marker_exp_avg,bygroup="ec",colors=exp_colors,title="EC marker",legend_title="logTPM",marker=ec_known_marker)
	fibro_mark_plot<- create_mark_plot(cell_marker_exp_avg,bygroup="fc",colors=exp_colors,title="FC marker",legend_title="logTPM",marker=fibroblast_known_marker)

	vsmc_mark_plot<- create_mark_plot(cell_marker_exp_avg,bygroup="vsmc",colors=exp_colors,title="VSMC marker",legend_title="logTPM",marker=vsmc_known_marker)

	mp_mark_plot<- create_mark_plot(cell_marker_exp_avg,bygroup="mp",colors=exp_colors,title="Macrophage marker",legend_title="logTPM",marker=macrophage_known_marker)

	tcell_mark_plot<- create_mark_plot(cell_marker_exp_avg,bygroup="tcell",colors=exp_colors,title="T cell marker",legend_title="logTPM",marker=t_cell_known_marker)


	marker_plot_list<-list(cm_mark_plot,vsmc_mark_plot,ec_mark_plot,mp_mark_plot,fibro_mark_plot,tcell_mark_plot)

	ml<-marrangeGrob(marker_plot_list, nrow=2, ncol=3)
	
	return(ml)
	
}


draw_horizon_violin<-function(seurat_res,marker,scale="width",width=0.5){

	cluster_marker_exp<-seurat_res$seurat@data[marker,]
	cluster_ident<-seurat_res$seurat@ident

	select_df <- data.frame(t(as.matrix(cluster_marker_exp)),cluster=cluster_ident)
	select_df_melt<-melt(select_df,id="cluster")

	colnames(select_df_melt)<-c("cluster","gene","logTPM")
	
	p<-ggplot(select_df_melt,aes(x=cluster,y=logTPM,fill=cluster)) +geom_violin(scale=scale,width=width) +facet_grid(~gene) + coord_flip()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text.x = element_text(size=8, angle=75))+ guides(fill=FALSE)

	#p<-ggplot(select_df_melt,aes(x=cluster,y=logTPM,fill=cluster)) +geom_violin(scale=scale,width=width) +facet_grid(~gene) + coord_flip()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=4, fill=NA))
			
	p<-font_theme(p)
	
	return(p)

}
