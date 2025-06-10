setwd("/gpfs/gibbs/project/kaminski/rd796/hpaproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(NICHES)

##############################
######## 1. Load Data ########
##############################
# load from aging_PA_diff.R
hpa <- readRDS("/gpfs/gibbs/project/kaminski/zc354/Artery_aging/data_share/annotated_3vs4_artery_aging_afterqc_mar12.Rds")
DefaultAssay(hpa)<-"RNA"
hpa$cell.type<-hpa$cell_type_annotation_final
hpa<-subset(hpa,subset=cell.type %in% c('EndothelialA','EndothelialB','Endothelial','Fibroblast','Mesothelial','SMC','MacLung','MacPeriVasc','Myeloid'))
hpa$age<-"old"
hpa$age[hpa$orig.ident %in% c("3m","f3m","m3mo_nox")]<-"young"

### Downsample Method 1
all.cells.remove=list()
cells.remove=NULL
all.Types <- unique(hpa@meta.data$cell.type)
for(i in 1:length(all.Types)){
temp.meta <- hpa@meta.data %>%
filter(cell.type==all.Types[i])
nold=nrow(temp.meta[temp.meta$age=='old',])
nyoung=nrow(temp.meta[temp.meta$age=='young',])
if(nold>nyoung){
ind=(1:nrow(hpa@meta.data))[hpa@meta.data$cell.type==all.Types[i]&hpa@meta.data$age=='old']
cells.remove <- sample(ind, size=(nold-nyoung), replace=FALSE)
} else{
ind=(1:nrow(hpa@meta.data))[hpa@meta.data$cell.type==all.Types[i]&hpa@meta.data$age=='young']
cells.remove <- sample(ind, size=(nyoung-nold), replace=FALSE)
}
all.cells.remove[[i]]<-cells.remove
}
hpa<- hpa[,-unlist(all.cells.remove)]
all(rownames(hpa@meta.data) == colnames(hpa@assays$RNA))
table(paste(hpa$age,hpa$cell.type))

### Downsample Method 2
#Idents(hpa)<-hpa$age
#diffsub=max(table(hpa$age))-min(table(hpa$age))
#whichsub=names(which(table(hpa$age)==max(table(hpa$age))))
#ind=(1:nrow(hpa@meta.data))[hpa@meta.data$age==whichsub]
#cells.remove <- sample(ind, size=diffsub, replace=FALSE)
#hpa<- hpa[,-unlist(cells.remove)]

###########################
######## 3. NICHES ########
###########################
Idents(hpa)<-hpa$cell.type
cell.endo <- WhichCells(hpa, idents = c('EndothelialA','EndothelialB','Endothelial'))
cell.mes <- WhichCells(hpa, idents = c('Fibroblast','Mesothelial','SMC'))
cell.imm <- WhichCells(hpa, idents = c('MacLung','MacPeriVasc','Myeloid'))
hpa <- SetIdent(hpa, cells = cell.endo, value = 'Endothelial')
hpa <- SetIdent(hpa, cells = cell.mes, value = 'Mesenchyme')
hpa <- SetIdent(hpa, cells = cell.imm, value = 'Immune')
hpa <- AddMetaData(object = hpa, metadata = Idents(hpa), col.name = 'Class')

hpa[["RNA"]]<-JoinLayers(hpa[["RNA"]])
##########################################################
# Impute the dataset and save the output for later
# Pick which genes to impute
num.cells.per.feature <- 50
genes.to.impute <- rownames(hpa)[rowSums(hpa@assays$RNA@layers$counts>0)>num.cells.per.feature]
# impute using ALRA
options(warn = 1)
gc()
hpa <- SeuratWrappers::RunALRA(hpa, genes.use = genes.to.impute)
# Save imputed data for later
gc()
##########################################################

DefaultAssay(hpa)="RNA"
data.list <- SplitObject(hpa, split.by="orig.ident")
# Make sure that data is normalized
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x)
})

# Run NICHES on each system and store/name the outputs
scc.list <- list()
for(i in 1:length(data.list)){
    print(i)
    scc.list[[i]] <- RunNICHES(data.list[[i]],
                               LR.database="fantom5",
                               species="mouse",
                               assay="alra", #alra or RNA
                               cell_types = "cell.type",
                               meta.data.to.map = c('orig.ident','cell.type','age','Class'),
                               SystemToCell = F,
                               CellToCell = T)
}
names(scc.list) <- names(data.list)
# Merge outputs
temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToCell # Isolate CellToCell Signaling
  gc()
}
cell.to.cell <- merge(temp.list[[1]],temp.list[[2]])
cell.to.cell <- merge(cell.to.cell,temp.list[[3]])
cell.to.cell <- merge(cell.to.cell,temp.list[[4]])
cell.to.cell <- merge(cell.to.cell,temp.list[[5]])
cell.to.cell <- merge(cell.to.cell,temp.list[[6]])
cell.to.cell <- merge(cell.to.cell,temp.list[[7]])
cell.to.cell

# Format metadata for plotting
cell.to.cell$Sample <- factor(cell.to.cell$orig.ident.Sending)
cell.to.cell$Condition <- factor(cell.to.cell$age.Sending)

# Clean each data set
# CTC
VlnPlot(cell.to.cell,group.by='Sample',c('nFeature_CellToCell'),raster=F,pt.size = 0.1)
cell.to.cell <- subset(cell.to.cell,nFeature_CellToCell>40)
VlnPlot(cell.to.cell,group.by='Sample',c('nFeature_CellToCell'),raster=F,pt.size = 0.1)

##########################################################
#Integration
cell.to.cell$temp<-as.character(cell.to.cell$orig.ident=='m27mo')
data.list <- SplitObject(cell.to.cell, split.by="temp")
#data.list <- SplitObject(cell.to.cell, split.by="Sample")
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data.list)
immune.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
cell.to.cell <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(cell.to.cell) <- "integrated"
##########################################################

#UMAP
cell.to.cell <- JoinLayers(cell.to.cell)
cell.to.cell <- NormalizeData(cell.to.cell)
cell.to.cell <- ScaleData(cell.to.cell)
cell.to.cell <- FindVariableFeatures(cell.to.cell)
cell.to.cell <- RunPCA(cell.to.cell,npcs = 20)
cell.to.cell <- RunUMAP(cell.to.cell,dims = 1:20)

##########################################################
#save(cell.to.cell,file='aging_PA_celltocell.RData')
#save(cell.to.cell,file='aging_PA_celltocell_integrated.RData')
save(cell.to.cell,file='aging_PA_cc_3_17.RData')

#load('aging_PA_celltocell.RData')
#load('aging_PA_celltocell_integrated.RData')
load('aging_PA_cc_3_17.RData')
##########################################################

#Redo metadata
cell.to.cell$Condition<-as.character(cell.to.cell$Condition)
cell.to.cell$Condition[cell.to.cell$Condition=="old"]<-"Aged"
cell.to.cell$Condition[cell.to.cell$Condition=="young"]<-"Young"
cell.to.cell$Condition<-factor(cell.to.cell$Condition,levels=c("Young",'Aged'))

#Plot UMAPS with relevant metadata
load('aging_PA_colormat.RData')
c1<-DimPlot(cell.to.cell,group.by = 'Condition',raster=F,shuffle = T)+scale_color_manual(values=c('yellow','purple'))
c2<-DimPlot(cell.to.cell,group.by = 'Sample',raster=F,shuffle = T)#+scale_color_manual(values=c('blue','purple','yellow','yellow','red','yellow'))
c3<-DimPlot(cell.to.cell,group.by = 'Class.Sending',raster=F,shuffle = T)+scale_color_manual(values=c('#330C2F','#7B287D','#7067CF'))
c4<-DimPlot(cell.to.cell,group.by = 'Class.Receiving',raster=F,shuffle = T)+scale_color_manual(values=c('#330C2F','#7B287D','#7067CF'))
c5<-DimPlot(cell.to.cell,group.by = 'SendingType',cols=fill_df$color,order=rev(fill_df$cell.type),raster=F,shuffle = T)
c6<-DimPlot(cell.to.cell,group.by = 'ReceivingType',cols=fill_df$color,order=rev(fill_df$cell.type),raster=F,shuffle = T)
#plot_grid(c1,c2,c3,c4,c5,c6)
plot_grid(c1,c2,c5,c6)

#imm.mes <- subset(cell.to.cell,subset = Class.Sending == 'Immune' & Class.Receiving == 'Mesenchyme')
imm.mes <- subset(cell.to.cell,subset = Class.Sending == 'Immune' & ReceivingType %in% c('Fibroblast','SMC'))
#imm.mes <- subset(cell.to.cell,subset = Class.Sending == 'Immune')
imm.mes <- ScaleData(imm.mes)
imm.mes <- FindVariableFeatures(imm.mes)
imm.mes <- RunPCA(imm.mes,npcs = 20)
imm.mes <- RunUMAP(imm.mes,dims = 1:20)

d1<-DimPlot(imm.mes,group.by = 'Condition',raster=F,shuffle = T)+scale_color_manual(values=c('blue','lightblue'))
d2<-DimPlot(imm.mes,group.by = 'Sample',raster=F,shuffle = T)#+scale_color_manual(values=c('blue','purple','yellow','yellow','red','yellow'))
d3<-DimPlot(imm.mes,group.by = 'Class.Sending',raster=F,shuffle = T)
d4<-DimPlot(imm.mes,group.by = 'Class.Receiving',raster=F,shuffle = T)
d5<-DimPlot(imm.mes,group.by = 'SendingType',raster=F,shuffle = T)#+scale_color_manual(values=c('#00BBDB','#00C19C','#00A5FF'))#c('#99A800','#F8766D','#06A4FF'))
d6<-DimPlot(imm.mes,group.by = 'ReceivingType',raster=F,shuffle = T)#+scale_color_manual(values=c('#00BC56','#FB61D7','#00BFC4'))
plot_grid(d1,d2,d3,d4,d5,d6)

plot_grid(c3,c4,d5,d6,d1,d2,ncol=2,byrow=T)

#DGE the entire cluster
DefaultAssay(imm.mes)<-'CellToCell'
imm.mes<-JoinLayers(imm.mes)
Idents(imm.mes)<-paste(imm.mes$age.Sending,imm.mes$VectorType)
age.mrk1<-FindMarkers(imm.mes,ident.1='old MacPeriVasc—Fibroblast',ident.2='young MacPeriVasc—Fibroblast')
age.mrk2<-FindMarkers(imm.mes,ident.1='old MacPeriVasc—SMC',ident.2='young MacPeriVasc—SMC')
age.mrk3<-FindMarkers(imm.mes,ident.1='old Myeloid—Fibroblast',ident.2='young Myeloid—Fibroblast')
age.mrk4<-FindMarkers(imm.mes,ident.1='old Myeloid—SMC',ident.2='young Myeloid—SMC')
dlist=list(age.mrk1,age.mrk2,age.mrk3,age.mrk4)
dlist=lapply(dlist,function(x){x[x$p_val_adj<0.05&x$avg_log2FC>0,]})
#dlist=lapply(dlist,function(x){x[x$p_val<0.05&x$avg_log2FC>0,]}) #for upset
#dlist=lapply(dlist,function(x){x[x$p_val_adj<0.05,]})
dlist=lapply(dlist,function(x){cbind(mech=rownames(x),x)})
names(dlist)<-c('MacPeriVasc—Fibroblast','MacPeriVasc—SMC','Myeloid—Fibroblast','Myeloid—SMC')
dlist2=bind_rows(dlist,.id='id')
#write.csv(dlist2,file='agingpamrk.csv')
namelist=lapply(dlist,rownames)
names(namelist)<-c('MacPeriVasc—Fibroblast','MacPeriVasc—SMC','Myeloid—Fibroblast','Myeloid—SMC')

library(UpSetR)
u1<-upset(fromList(namelist),order.by='freq',sets=names(namelist),keep.order=T,empty.intersections="on",text.scale=1.5)
u1

namelistint=Reduce(intersect,namelist[c('MacPeriVasc—Fibroblast','MacPeriVasc—SMC')])
rnk=order(colMeans(rbind(match(namelistint, rownames(dlist[['MacPeriVasc—Fibroblast']])),match(namelistint, rownames(dlist[['MacPeriVasc—SMC']])))))
namelistint=namelistint[rnk]
#write.csv(namelistint,file='agingpamrk_overlap.csv')


#Endothelial-Mesenchymal
endo.mes <- subset(cell.to.cell,subset = Class.Sending == 'Endothelial' & ReceivingType %in% c('Fibroblast','SMC'))
DefaultAssay(endo.mes)<-'CellToCell'
endo.mes<-JoinLayers(endo.mes)
endo.mes@meta.data[,2:21]<-sapply(endo.mes@meta.data[,2:21],function(x){gsub("Endothelial.","Endothelial",x)})
Idents(endo.mes)<-paste(endo.mes$age.Sending,endo.mes$VectorType)
age.mrk<-FindMarkers(endo.mes,ident.1='old Endothelial—SMC',ident.2='young Endothelial—SMC')
namelistint<-rownames(age.mrk)[age.mrk$avg_log2FC>0&age.mrk$p_val_adj<0.05]
#write.csv(age.mrk,file='endo_smc_mrk.csv')

#GO
check=namelistint
golist <- unlist(strsplit(check,split = "—"))
noquote(golist)

library('gprofiler2')
go_out<-gost(query=golist, organism = "mmusculus")
#kres=go_out$result[go_out$result$source=='KEGG',]
kres=go_out$result[go_out$result$source=='KEGG'&go_out$result$term_size<200,]
kmat=kres[,c(3,11)]
if(nrow(kmat)>10){kmat<-kmat[1:10,]}
kmat$term_name<-factor(kmat$term_name,levels=rev(kmat$term_name))
ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) +
    geom_col()+scale_fill_gradient(low="lightblue",high="blue",limits=c(0,max(-log10(kmat$p_value))))+
    geom_text(aes(label = term_name),colour='white',hjust=1) +
    theme_classic()+theme(legend.position="right",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='KEGG Term')+ guides(fill=guide_legend(title="-log10(p-value)"))


#Violins
vector='Tgfb1—Itgb8'
#imm.mes$orig.ident=factor(imm.mes$orig.ident,levels=c('3m','3m1','m3mo','20m','24m','m27mo'))
v1<-VlnPlot(imm.mes,vector,split.by='Condition',group.by='orig.ident',cols=c('yellow','purple'))+labs(x='Signaling Archetype') # grouping by cluster
v1<-VlnPlot(imm.mes,vector,group.by='Condition',pt.size=0,cols=c('yellow','purple'))+labs(x='Signaling Archetype') # grouping by cluster
v2<-VlnPlot(cell.to.cell,vector,split.by='Condition',group.by='Class.Joint',pt.size=0,cols=c('yellow','purple'))+labs(x='Cell Class')+ggtitle(NULL) # grouping by class
v3<-VlnPlot(imm.mes,vector,split.by='Condition',group.by='cell.type.Joint',pt.size=0,cols=c('yellow','purple'))+labs(x='Cell Type')+ggtitle(NULL) # grouping by cell type
plot_grid(v1,v3,ncol=1)+theme(plot.margin = unit(c(0, 0, 0, 1), "cm"))

#Heatmap
load('aging_PA_colormat.RData')
tmp=unlist(fill_df$color)
names(tmp)<-fill_df$cell.type
tmp=c(tmp,old='blue',young='lightblue')
feats=sub.all$gene
ann=data.frame(Archetype=rep(unique(sub.all$cluster),each=nrow(sub.all)/length(unique(sub.all$cluster))))
h=CustomHeatmap(imm.mes,features=feats,rowanno=ann,selected.row.anotations = NULL,
secondary='Condition',tertiary = 'SendingType',quarternary = 'ReceivingType',secondary.cols=tmp[15:16],tertiary.cols=tmp,quarternary.cols=tmp,
labels=c('Signaling Archetype','Age','Sending Cell Type','Receiving Cell Type'),range.frac=0.5)
draw(h,merge_legend=T)

#Circuit Plot. Myeloid & Mesenchymal
split.lung <- SplitObject(hpa,split.by = 'age')
split.ctc <- SplitObject(cell.to.cell,split.by = 'Condition')
vector='Tgfb1—Itgb8'
nodes<-c('MacPeriVasc','Myeloid','MacLung','Fibroblast','SMC')
#cols<-fill_df$color[match(nodes,fill_df$cell.type)]
cols<-c("#0BB702", "#00BFC4", "#F8766D", "#00A9FF", "#00C19A")
c1=CircuitPlot(split.lung$Young,split.ctc$young,group.by='cell.type',global.node.list = nodes,cols.use=cols,feature=vector)+ggtitle(paste(vector,": Young"))
c2=CircuitPlot(split.lung$Aged,split.ctc$old,group.by='cell.type',global.node.list = nodes,cols.use=cols,feature=vector)+ggtitle(paste(vector,": Aged"))
c1+c2

#Circuit Plot. Endothelial & Mesenchymal
cell.to.cell@meta.data<-sapply(cell.to.cell@meta.data,function(x){gsub("Endothelial.","Endothelial",x)})
vector='Bmp6—Acvr2a'
nodes<-c('Endothelial','Fibroblast','SMC')
cols<-fill_df$color[match(nodes,fill_df$cell.type)]
c1=CircuitPlot(split.lung$Young,split.ctc$young,group.by='cell.type',global.node.list = nodes,cols.use=cols,feature=vector)+ggtitle(paste(vector,": Young"))
c2=CircuitPlot(split.lung$Aged,split.ctc$old,group.by='cell.type',global.node.list = nodes,cols.use=cols,feature=vector)+ggtitle(paste(vector,": Aged"))
c1+c2