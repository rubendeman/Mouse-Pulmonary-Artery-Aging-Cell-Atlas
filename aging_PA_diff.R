#This script requires R 4.4.1
setwd("/gpfs/gibbs/project/kaminski/rd796/hpaproj")
library('Seurat')
library('ggplot2')
library('cowplot')
library('patchwork')
library('stringr')
library('dplyr')
library('ggpubr')
library('tidyr')
library('ComplexHeatmap')

##############################
######## 1. Load Data ########
##############################
#hpa <- readRDS('mold.labeled_Seurat_object_6-2-2023.rds') #Old
#hpa <- readRDS('/gpfs/gibbs/project/kaminski/zc354/Azimuth_Analysis/Output/aging_PA_anno.Rds') #Old
#hpa <- readRDS('/gpfs/gibbs/project/kaminski/zc354/Artery_aging/data_share/merged_raw_3vs3_qc_anno_nov1.Rds') #Backup at RDM_AGING_PA.RData
hpa <- readRDS("/gpfs/gibbs/project/kaminski/zc354/Artery_aging/data_share/annotated_3vs4_artery_aging_afterqc_mar12.Rds")

#hpa@meta.data$cell.type<-hpa@meta.data$predicted.annotation.l1
hpa@meta.data$cell.type<-hpa@meta.data$cell_type_annotation_final
Idents(hpa) <- hpa@meta.data$cell.type
hpa<-subset(hpa,subset=cell.type %in% c('UnknownB','NeuronB','NeuronA','UnknownA','Chondrocyte'),invert=T)
age=grepl("2",hpa$orig.ident,fixed=T)
hpa@meta.data$age=factor(age,labels=c("Young","Aged"))

hpa<-JoinLayers(hpa)
DefaultAssay(hpa)<-'RNA'
hpa<-NormalizeData(hpa)
DimPlot(hpa,reduction = 'umap')

hpasub=subset(hpa,subset=cell.type %in% c('MacLung','MacPeriVasc','Myeloid'))
i1<-DimPlot(hpasub,reduction = 'umap',cols=c('#F8766D','#00A9FF','#00BFC4'))+xlim(c(-16,-8))

set.seed(2)
fill_df <- hpa@meta.data %>% select('cell.type') %>% unique() %>% mutate(color = sample(scales::hue_pal()(length(unique(hpa$cell.type)))))
#save(fill_df,file='aging_PA_colormat.RData')
a1<-DimPlot(hpa,reduction = 'umap',label=T,repel=T,group.by='cell.type',cols=fill_df$color,order=rev(fill_df$cell.type))+ggtitle(NULL)+theme(plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
a2<-DimPlot(hpa,reduction = 'umap',label=F,repel=T,group.by='age',cols=c('yellow','purple'))+ggtitle(NULL)+theme(plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
a3<-DimPlot(hpa,reduction = 'umap',label=F,repel=T,group.by='orig.ident')+ggtitle(NULL)+theme(plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
a1+a2+a3

#Marker Genes
cell.types=names(table(hpa$cell.type)[order(names(table(hpa$cell.type)))])
levels(hpa) <- c(rev(cell.types),levels(hpa)[!(levels(hpa) %in% cell.types)])
#marker.genes<-c('Top2a','Bank1','Jchain','Vwf','Fabp4','Calcrl','Dcn','Inmt','Mmrn1','Chil3','C1qa','Rbfox1','Retnla','S100a8','Acta2','Trbc2')
marker.genes<-c('Top2a','Bank1','Jchain','Vcam1','Vwf','Fabp4','Epcam','Inmt','Chil3','C1qa','Rbfox1','Retnla','S100a8','Hbb-bt','Acta2','Trbc2')
DotPlot(hpa, features = marker.genes, idents=cell.types)+ theme_minimal()+theme(axis.text.x = element_text(angle=90))+ RotatedAxis()+labs(y=NULL)

###################################
######## 2. Myeloid Subset ########
###################################
#Myeloid Subset
hpasub=subset(hpa,subset=cell.type %in% c('MacLung','MacPeriVasc','Myeloid'))

hpasub <- NormalizeData(hpasub)
hpasub <- ScaleData(hpasub)
hpasub <- FindVariableFeatures(hpasub)
hpasub <- RunPCA(hpasub,npcs = 10)
hpasub <- RunUMAP(hpasub,dims = 1:10)

i1<-DimPlot(hpasub,label=T,repel=T,group.by='cell.type',cols=c('#00BFC4','#F8766D','#00A9FF'))+ggtitle(NULL)
#i2<-DotPlot(hpasub,features=c('Chil3','Fabp5','C1qa','Pf4','Cd24a','Itgae'),group.by='cell.type')
i2<-DotPlot(hpasub,features=c('Chil3','Fabp5','C1qa','Pf4','Ccr2','Cd24a'),group.by='cell.type') #Ccr2 is myeloid marker
i1+i2

##########################################
######## 3. Add age/SenMayo score ########
##########################################
#SENMAYO
load('/gpfs/gibbs/project/kaminski/rd796/ageproj/senmayolist.RData') #SenMayo list
#cell=read.table('/gpfs/gibbs/project/kaminski/rd796/ageproj/CellAge.csv',sep=';',header=TRUE)[,2]
#fridman=read.csv('/gpfs/gibbs/project/kaminski/rd796/ageproj/FRIDMAN.csv')
#segura=read.csv('/gpfs/gibbs/project/kaminski/rd796/ageproj/SEGURA.csv')
#purcell=read.csv('/gpfs/gibbs/project/kaminski/rd796/ageproj/PURCELL.csv')
#consensus=c('C1qtnf1','C3','Ccnd2','Fer1l4','Hist1h1c','Hist1h2ac','Hist1h2bd','Pik3ip1','Slc1a1','Tcea2','Wdr63')
#csgene=c('Cdkn2b','Ccnd1','Hdac1','E2f1','E2f3','Ccne1','Jun','Pcna','Smad3','Mapk1','Sirt6','Tgfb2','Smad1','Cdkn2c','Akt1','Ccne2','Tp63','Prkdc','Rbl1','Mapk3')

#Add SenMayo score
#senmayolist=consensus #Or change to cell, fridman, purcell, or segura
senmayolist=str_to_title(unlist(senmayolist))
hpa <- ScaleData(object = hpa, features = (senmayolist)) #Need to have SenMayo genes as char vector
#Make sure all SenMayo genes in dataset
DefaultAssay(hpa)<-'RNA'
features.keep=unlist(senmayolist) %in% rownames(hpa)
senmayolist=senmayolist[features.keep]
#Now compute SenMayo score and add to meta-data
hpa <- AddModuleScore(hpa, features = list(senmayolist), slot='counts',name="SenMayo",ctrl=50) #*** Change name to reflect the senescence signature

#Let's look at SenMayo score distribution
Idents(hpa)=hpa@meta.data$cell.type
range(hpa@meta.data$SenMayo1)
#Density plot to determine threshold (i.e. High SenMayo cells are above this threshold)
plot(density(hpa@meta.data$SenMayo1),main='Senescence Score Distribution',xlab='Senescence Score',frame.plot=F,xlim=c(-0.2,0.6),xaxt="n",ylim=c(0,10))
thresh=0.3
abline(v=thresh,col="red")
axis(1,pos=0)

#Now we can add a binary SenMayo score
senbin <- rep("Low Sen Score",length(hpa@meta.data$SenMayo1))
senbin[hpa@meta.data$SenMayo1>thresh]="High Sen Score" #Cells above threshold are "High SenMayo cells"
hpa <- AddMetaData(object = hpa, metadata = senbin,col.name = 'senbin')

#s1=DimPlot(hpa,group.by='senbin',order='High Sen Score',cols=c('#4B88A2','#D3D4D9'))+ggtitle(NULL)
s1=DimPlot(hpa,group.by='senbin',cols=c('#4B88A2','#D3D4D9'))+ggtitle(NULL)
s2=DimPlot(hpa,label=T,repel=T,group.by='cell.type')+ggtitle(NULL)
s1+s2
table(hpa$senbin)

#################################
######## 4. YOUNG VS OLD ########
#################################
hpa$cell.type[hpa$cell.type %in% c('EndothelialA','EndothelialB','Endothelial')]<-'Endothelial'
cells.test=c('Fibroblast','SMC','MacPeriVasc','Myeloid','Endothelial')

#Violin Plot
Idents(hpa)<-hpa$cell.type
VlnPlot(hpa, features = "SenMayo1", group.by="cell.type",split.by = "age",idents=cells.test,pt.size = 0.1,cols=c('blue','lightblue'))+ggtitle(NULL)+stat_compare_means(aes(label = paste0(after_stat(p.signif))))+labs(x=NULL)
VlnPlot(hpa, features = "SenMayo1", group.by="cell.type",split.by = "orig.ident",idents=cells.test,pt.size = 0.1)+ggtitle(NULL)+stat_compare_means(aes(label = paste0(after_stat(p.signif))))+labs(x=NULL)
VlnPlot(hpa, features = "Cdkn1a", group.by="cell.type",split.by = "age",idents=cells.test,pt.size = 0)+ggtitle(NULL)+stat_compare_means(aes(label = paste0(after_stat(p.signif))))+labs(x=NULL)
VlnPlot(hpa, features = "Cdkn1a", group.by="cell.type",split.by = "orig.ident",idents=cells.test,pt.size = 0.1)+ggtitle(NULL)+stat_compare_means(aes(label = paste0(after_stat(p.signif))))+labs(x=NULL)

#Boxplots by cell type (cells)
rscale1=function(x){(x-min(x))/(max(x)-min(x))}
data=data.frame(Type=hpa@meta.data$cell.type, Age=hpa@meta.data$age, Score=rscale1(hpa@meta.data$SenMayo1), Gene=hpa@assays[["RNA"]]@layers$data[which(rownames(hpa)=="Cdkn1a"),])
data=data[data$Type %in% cells.test,]
ggplot(data, aes(x=Type, y=log10(Score), fill=Age)) + geom_boxplot(outlier.size=0.1) + theme_classic()+theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) + geom_point(position=position_jitter(),size=0.1)+scale_fill_manual(values=c('yellow','purple'))+
stat_compare_means(aes(label = paste0(after_stat(p.signif)))) + ylab("log SenMayo Score")

data %>% group_by(Type, Age) %>% summarise(score=median(Score))
ggplot(data, aes(x=Type, y=log10(Gene), fill=Age)) + geom_boxplot(outlier.size=0.1) + theme_classic()+theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) + geom_point(position=position_jitter(),size=0.1)+scale_fill_manual(values=c('blue','lightblue'))+
stat_compare_means(aes(label = paste0(after_stat(p.signif)))) + ylab("log SenMayo Score")

#Boxplots by cell type (samples)
tbl=data.frame(sid=hpa$orig.ident,ct=hpa$cell.type,sm=hpa$SenMayo1,Age=hpa$age, Gene=hpa@assays[["RNA"]]@layers$data[which(rownames(hpa)=="Cdkn1a"),])
tbl2=tbl %>% group_by(sid,ct,Age) %>% summarise(sm=mean(sm),gene=mean(Gene))
tbl2=tbl2[tbl2$ct %in% cells.test,]
ggplot(tbl2,aes(x=ct,y=sm,fill=Age))+geom_boxplot()+geom_point(position=position_jitterdodge())+ theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12))+scale_fill_manual(values=c('yellow','purple'))+
xlab("Type") + ylab("SenMayo Score")
ggsave('crime1.pdf', width = 7, height = 5)

ggplot(tbl2,aes(x=ct,y=gene,fill=Age))+geom_boxplot()+geom_point(position=position_jitterdodge())+ theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+scale_fill_manual(values=c('blue','lightblue'))+
stat_compare_means(aes(label = paste0(after_stat(p.signif))))

#Intelligently downsample to show cells equal per sample
cells.keep=NULL
set.seed(15)
for (k in cells.test){
    nametmp<-names(table(hpa$orig.ident))
    strtmp=lapply(nametmp,function(l){colnames(hpa)[hpa$cell.type==k&hpa$orig.ident==l]})
    mintmp=min(unlist(lapply(strtmp,length)))
    strtmp2=lapply(strtmp,function(l){try(sample(l,mintmp))})
    cells.keep=c(cells.keep,unlist(strtmp2))
}
hpa2=hpa[,cells.keep]
rscale1=function(x){(x-min(x))/(max(x)-min(x))}
data=data.frame(Type=hpa2@meta.data$cell.type, Age=hpa2@meta.data$age, Score=(hpa2@meta.data$SenMayo1), Gene=hpa2@assays[["RNA"]]@layers$data[which(rownames(hpa)=="Cdkn1a"),])
data=data[data$Type %in% cells.test,]
ggplot(data, aes(x=Type, y=(Score), fill=Age)) + geom_boxplot(outlier.size=0.1) + theme_classic()+theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) + geom_point(position=position_jitter(),size=0.1)+scale_fill_manual(values=c('blue','lightblue'))+
stat_compare_means(aes(label = paste0(after_stat(p.signif)))) + ylab("SenMayo Score")
data %>% group_by(Type,Age) %>% summarise(score=median(Score))

data2=data[data$Gene>0,]
ggplot(data2, aes(x=Type, y=Gene, fill=Age)) + geom_boxplot(outlier.size=0.1) + theme_classic()+theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) + geom_point(position=position_jitter(),size=0.1)+scale_fill_manual(values=c('blue','lightblue'))+
stat_compare_means(aes(label = paste0(after_stat(p.signif)))) + ylab("Gene Expression")

#P16 Percent
#Idents(hpa)<-hpa@meta.data$age
Idents(hpa)<-hpa@meta.data$orig.ident
a<-DotPlot(hpa,features='Cdkn1a')
a$data$id<-hpa@meta.data$age[match(a$data$id,hpa@meta.data$orig.ident)]
p1<-ggplot(a$data, aes(x=id, y=pct.exp)) + geom_boxplot() + labs(x='Age',y='% Cells Expressing p21')+theme_minimal()+ geom_point(position=position_jitter())
p2<-ggplot(a$data, aes(x=id, y=avg.exp)) + geom_boxplot() + labs(x='Age',y='Average Expression of p21')+theme_minimal()+ geom_point(position=position_jitter())
p1+p2

Idents(hpa)<-hpa$age
a<-DotPlot(hpa,features=c('Cdkn1a','Cdkn2a','Gdf15'),split.by='age')
a
View(a$data)

###################################################
######## 5. Differentially expressed genes ########
###################################################
cells.test=c('Fibroblast','SMC','MacLung','MacPeriVasc','Myeloid')
#Find markers (cell-type specific)
datalist = list() #Empty list for results
#Iterate through cell types
for (i in cells.test){
Idents(hpa)=hpa@meta.data$cell.type
tmp=subset(hpa, idents = i) #Subset dataset so we look at a single cell type
Idents(tmp)=tmp@meta.data$age #Make age the active Ident
agegenes <- FindMarkers(tmp, ident.1 = "Aged", ident.2 = "Young", verbose = FALSE)
datalist[[i]]=na.omit(agegenes) #Add to our list of data frames
}
#Put results in a data frame
names(datalist)=cells.test

#LOOK AT THE RESULTS
v=datalist[['Myeloid']]
View(v[!grepl('^Rp',rownames(v))&!grepl('^mt',rownames(v)),])

#GOST
#Make it into a loop
library('gprofiler2')
plist=list()
cells.test=c('Fibroblast','SMC','MacPeriVasc','Myeloid')
for (i in cells.test){
agemat=datalist[[i]]
genes.test.pos=rownames(agemat)[agemat$avg_log2FC>0&agemat$p_val<0.05]
genes.test.neg=rownames(agemat)[agemat$avg_log2FC<0&agemat$p_val<0.05]
genes.test.pos=genes.test.pos[!grepl('^Rp',genes.test.pos)&!grepl('^mt',genes.test.pos)][1:500]
genes.test.neg=genes.test.neg[!grepl('^Rp',genes.test.neg)&!grepl('^mt',genes.test.neg)][1:500]

pos_go<-gost(query=genes.test.pos, organism = "mmusculus")
neg_go<-gost(query=genes.test.neg, organism = "mmusculus")

kres=pos_go$result[pos_go$result$source %in% c('REAC'),]
kres=kres[order(kres$p_value,decreasing=F),]
kmat=kres[,c(3,11)]
if(nrow(kmat)>5){kmat<-kmat[1:5,]}
kmat$term_name<-factor(kmat$term_name,levels=kmat$term_name)
f1<-ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) +
    geom_col()+scale_fill_gradient2(low="gray",  high="yellow",breaks=c(mean(-log10(kmat$p_value))))+
    geom_text(aes(label = term_name,hjust=1)) +
    theme_classic()+theme(legend.position="none",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='REAC Term')+ggtitle(paste0(i,' (+)'))+theme(plot.title = element_text(hjust = 0.5))
kres=neg_go$result[neg_go$result$source %in% c('REAC'),]
kres=kres[order(kres$p_value,decreasing=F),]
kmat=kres[,c(3,11)]
if(nrow(kmat)>5){kmat<-kmat[1:5,]}
kmat$term_name<-factor(kmat$term_name,levels=kmat$term_name)
f2<-ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) +
    geom_col()+scale_fill_gradient2(low="gray",  high="purple",breaks=c(mean(-log10(kmat$p_value))))+
    geom_text(aes(label = term_name,hjust=1)) +
    theme_classic()+theme(legend.position="none",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='REAC Term')+ggtitle(paste0(i,' (-)'))+theme(plot.title = element_text(hjust = 0.5))
plist[[paste0(i,' (+)')]]=f1
plist[[paste0(i,' (-)')]]=f2
}
plot_grid(plotlist=plist,ncol=2,byrow=T)

#SENMAYO hmap by sample
hpa$cell.type[hpa$cell.type %in% c('EndothelialA','EndothelialB','Endothelial')]<-'Endothelial'
cells.test=c('Endothelial','Fibroblast', 'SMC', 'MacLung','MacPeriVasc','Myeloid')
smallobj=subset(hpa,subset=cell.type %in% cells.test)

Idents(hpa)<-paste(hpa$cell.type,hpa$age)
sm=senmayolist[senmayolist %in% rownames(hpa)]
agegenes <- FindMarkers(hpa, features=sm,ident.1 = "Fibroblast Aged", ident.2 = "Fibroblast Young",verbose = FALSE)
fib_genes=rownames(agegenes)[agegenes$avg_log2FC>0&agegenes$p_val<0.05]
fib_genes2=rownames(agegenes)[agegenes$avg_log2FC>0&agegenes$p_val_adj<0.05]
#fib_genes<-c("Il6st","Serpine2","Bmp6","Ctsb","Rps6ka5","Mmp3","Angpt1","Tnfrsf1b","Ptbp1","Egfr","Vegfc","Cxcl12")
agegenes <- FindMarkers(hpa, features=sm,ident.1 = "MacPeriVasc Aged", ident.2 = "MacPeriVasc Young",verbose = FALSE)
mac_genes=rownames(agegenes)[agegenes$avg_log2FC>0&agegenes$p_val<0.05]
mac_genes2=rownames(agegenes)[agegenes$avg_log2FC>0&agegenes$p_val_adj<0.05]
#mac_genes<-c("Igf1","Ptbp1","Lcp1","Ccl8","Il6st","Tnfrsf1a","Gmfg","Egfr","Cxcl12","Il15","Cxcl16","Bmp2")
agegenes <- FindMarkers(hpa, features=sm,ident.1 = "Myeloid Aged", ident.2 = "Myeloid Young",verbose = FALSE)
mye_genes=rownames(agegenes)[agegenes$avg_log2FC>0&agegenes$p_val<0.05]
mye_genes2=rownames(agegenes)[agegenes$avg_log2FC>0&agegenes$p_val_adj<0.05]


feats=sm
#feats=c('Gdf15','Cdkn1a')
rowsplit=NULL
meow<-AverageExpression(smallobj,assays='RNA',layer='data',feats,group.by = c('cell.type','age','orig.ident'))
meow<-as.matrix(meow$RNA)
#CREATE ANNOT
    ann <- data.frame(sub('_.*','',colnames(meow)),sub('.*_','',colnames(meow)))
    colnames(ann) <- c('Cell','Sample')
    ann=ann%>%mutate(Age=smallobj$age[match(ann$Sample,smallobj$orig.ident)])
    ann[is.na(ann)]<-'Young'
colAnn <- HeatmapAnnotation(df = ann)
meow=t(apply(meow,1,function(x){(x-min(x))/(max(x)-min(x))}))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
h1=Heatmap(na.omit(meow), name = "Expression",  column_split=data.frame(ann$Cell,factor(ann$'Age',levels=c('Young','Aged'))),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,
            cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = F,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45, row_split=rowsplit,
            show_column_names = FALSE, show_row_names=T,  use_raster = F)
draw(h1,merge_legend=T)

#hmap by cell
feats=sm
#NEWER sub
cells.test=c('Fibroblast','MacPeriVasc','Myeloid')
all=c(fib_genes2,mac_genes2,mye_genes2)
#dup=unique(all[duplicated(all)])
dup<-names(table(all))[which(table(all)==3)]
fib1=fib_genes[!fib_genes %in% dup][1:5]
mac1=mac_genes[!mac_genes %in% dup][1:5]
mye1=mye_genes[!(mye_genes %in% fib1 | mye_genes %in% mac1 | mye_genes %in% dup)][1:5]
#feats=c("Serpine2","Bmp6","Ctsb","Mmp3","Il6st","Ptbp1","Egfr","Cxcl12","Igf1","Lcp1","Ccl8","Il15")
feats=c(fib1,mac1,mye1,dup)

smallobj=subset(hpa,subset=cell.type %in% cells.test)
len=c(length(fib1),length(mac1),length(mye1),length(dup))
anno=data.frame(Cell=rep(c(cells.test,"All"),len))

#feats=c('Gdf15','Cdkn1a')
rowsplit=NULL
ord=order(paste(smallobj$cell.type,smallobj$age,smallobj$orig.ident))
ind=match(feats,rownames(smallobj))
meow=smallobj@assays$RNA@layers$data[ind,] %>% as.matrix() #check order
rownames(meow)<-rownames(smallobj@assays$RNA)[ind]
#meow=smallobj@assays$RNA@data[feats,ord] %>% as.matrix() #check order

#Col Annot
ann <- data.frame(smallobj$cell.type, smallobj$orig.ident, smallobj$age)
#ann <- data.frame(smallobj$cell.type[ord], smallobj$orig.ident[ord], smallobj$age[ord])
colnames(ann) <- c('Cell','Sample','Age')
            hcols=list(Cell=c(Fibroblast='#0BB702',MacPeriVasc='#F8766D',Myeloid='#00A9FF'),Age=c(Aged='purple',Young='yellow'))
#colAnn <- HeatmapAnnotation(df = ann)
colAnn <- HeatmapAnnotation(df = ann,col=hcols)

#Row Annot
            hcols=list(Cell=c(Fibroblast='#0BB702',MacPeriVasc='#F8766D',Myeloid='#00A9FF',All='black'))
HAleft <- rowAnnotation(df=anno,col=hcols,show_legend=T)
#HAleft=NULL

meow=t(apply(meow,1,function(x){(x-min(x))/(max(x)-min(x))}))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
h1=Heatmap(na.omit(meow), name = "Expression",  column_split=data.frame(ann$Cell,factor(ann$'Age',levels=c('Young','Aged'))),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,left_annotation=HAleft,
            cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = F,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45, row_split=rowsplit,
            show_column_names = FALSE, show_row_names=T,  use_raster = F)
pdf('crime2.pdf', width = 12, height = 6)       
draw(h1,merge_legend=T)
dev.off()
