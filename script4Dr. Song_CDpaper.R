#### code for Dr. Song's project


######################  Figure 4I: reanalysis of public data (Daniel et al) #############
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Matrix)

DAT <- Read10X(data.dir = '/fs/ess/PAS2579/NoJoon/GSE188666/out',gene.column=1)
DAT <- CreateSeuratObject(counts = DAT,project = "GSE18866",min.cells = 0, min.features = 0)

a <-sapply(strsplit(colnames(DAT),'_',1),'[',3)
DAT$time <- a


meta <- read.table('/fs/ess/PAS2579/NoJoon/GSE188666/out/GSE188666_scRNA_LCMV_metadata.tsv')

identical(colnames(DAT),rownames(meta))  # check if the order of cells are the same

DAT <- AddMetaData(DAT, meta)

Idents(DAT) <- DAT$ident
UMAP <- data.frame(UMAP_1=DAT$UMAP_1,UMAP_2=DAT$UMAP_2)
UMAP <- as.matrix(UMAP)
DAT[['umap']] <-CreateDimReducObject(embeddings = UMAP,key = 'UMAP_',assay = DefaultAssay(DAT))
DAT$clusters <- DAT$ident


DimPlot(DAT)

count <- DAT@assays$RNA@counts
UMAP <- DAT@reductions$umap@cell.embeddings

### NKRT signature
GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')   # alternative set2
GENE.count <- count[which(toupper(rownames(count))%in% toupper(GENE$gene)),]
GENE.count <- t(GENE.count)

a <- rep(1,ncol(GENE.count))
a <- as.matrix(a)
GENE.score <- GENE.count%*%a

summary(GENE.score[,1])  
hist(GENE.score[,1])
hist(log(GENE.score[,1]+1))


rownames(GENE.score) <- rownames(GENE.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(GENE.score)),]

df <- data.frame(UMAP.sub,GENE.score, log(GENE.score+1))

mid <- median(log(GENE.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(GENE.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()



## RAR

RAR.target <- read.csv('/fs/ess/PAS2579/NoJoon/set3.csv',header=T)

RAR.count <- count[which(rownames(count)%in% RAR.target$Set3_RAR),]
RAR.count <- t(RAR.count)

a <- rep(1,ncol(RAR.count))
a <- as.matrix(a)
RAR.score <- RAR.count%*%a

summary(RAR.score[,1])


rownames(RAR.score) <- rownames(RAR.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(RAR.score)),]

df <- data.frame(UMAP.sub,RAR.score)
mid <- median(log(RAR.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(RAR.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()



##### combine RAR and NKRT together

RAR <- RAR.target$Set3_RAR


TARGET.GENE <- c(GENE$gene,toupper(RAR))

TARGET.count <- count[which(toupper(rownames(count))%in% TARGET.GENE),]
TARGET.count <- t(TARGET.count)


a <- vector()

IND <- which(colnames(TARGET.count)%in% toupper(GENE$gene))

for(i in 1:ncol(TARGET.count)){
  if (i %in% IND){
    temp <- colnames(TARGET.count)[i]
    tt <- GENE[which(toupper(GENE$gene)==temp),]$exp
    a[i] <- ifelse(tt=='High',0.1,0.1)
  }else{
    a[i] <- 0.1
  }
}



a <- as.matrix(a)
TARGET.score <- TARGET.count %*% a

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(TARGET.score)),]
TARGET.score.mod <- log(TARGET.score[,1])

for (i in 1:length(TARGET.score.mod)){
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]>4.55,4.55,TARGET.score.mod[i])
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]<3.25,3.25,TARGET.score.mod[i])
}

hist(TARGET.score.mod)


df <- data.frame(UMAP.sub,TARGET.score,TARGET.score.mod)

mid <- median(TARGET.score.mod)
ggplot(df,aes(UMAP_1, UMAP_2,color=TARGET.score.mod))+geom_point(size=0.2) +
  scale_colour_gradient2(low = "royalblue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()



################## Figure 6 F,G,H (zheng et al., GSE221064)
## create seurat object
library(Seurat)
library(dplyr)
library(ggplot2)

dir.in <- '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/'

NAMES <- list.files(dir.in)

data.id <- grep('GSM',NAMES)
DATA.names <- NAMES[data.id]

for (i in 1:length(data.id)){
  data.dir <- paste0(dir.in, DATA.names[i])
  data <- Read10X(data.dir = data.dir)
  seurat.obj <- CreateSeuratObject(counts = data,project = DATA.names[i],min.features = 200,min.cells = 3)
  seurat.obj[['percent.mt']] <- PercentageFeatureSet(seurat.obj,pattern = '^mt-')
  
  saveRDS(seurat.obj,file = paste0(dir.in,'seurat/',DATA.names[i],'_seurat.rds'))
}


dir.in <- '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/seurat/'
FILES <-list.files(dir.in,pattern = 'rds')
FILES


i <-1
DATA1 <- readRDS(paste0(dir.in,FILES[i]))
DATA1 <- subset(DATA1,subset = percent.mt <20 & nCount_RNA <75000)

i <-2
DATA2 <- readRDS(paste0(dir.in,FILES[i]))
DATA2 <- subset(DATA2,subset = percent.mt <25 & nCount_RNA < 100000 )

i <-3
DATA3 <- readRDS(paste0(dir.in,FILES[i]))
DATA3 <- subset(DATA3,subset = percent.mt <20 & nCount_RNA < 100000 )

#### log normalization and combine

DATA1 <-NormalizeData(DATA1, normalization.method = 'LogNormalize',scale.factor = 10000)
DATA1 <-FindVariableFeatures(DATA1,selection.method='vst',nfeatures=2000)
DATA1 <-RenameCells(DATA1,add.cell.id = 'A7')  # add prefix to the colnames to avoid duplicated cell names during integration

DATA2 <-NormalizeData(DATA2, normalization.method = 'LogNormalize',scale.factor = 10000)
DATA2 <-FindVariableFeatures(DATA2,selection.method='vst',nfeatures=2000)
DATA2 <-RenameCells(DATA2,add.cell.id = 'A7.1')  # add prefix to the colnames to avoid duplicated cell names during integration

DATA3 <-NormalizeData(DATA3, normalization.method = 'LogNormalize',scale.factor = 10000)
DATA3 <-FindVariableFeatures(DATA3,selection.method='vst',nfeatures=2000)
DATA3 <-RenameCells(DATA3,add.cell.id = 'B7')  # add prefix to the colnames to avoid duplicated cell names during integration

DATA.list <- list(DATA1,DATA2,DATA3)

features <- SelectIntegrationFeatures(object.list = DATA.list, nfeatures = 3000)
DATA.anchors <- FindIntegrationAnchors(object.list = DATA.list,
                                       anchor.features = features)
AB.combined <- IntegrateData(anchorset = DATA.anchors)

saveRDS(AB.combined, file = '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/seurat/combined/AB.combined.rds')

AB.combined <- readRDS(file = '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/seurat/combined/AB.combined.rds')

AB.combined <- ScaleData(AB.combined)
AB.combined <- RunPCA(AB.combined, npcs = 30)
ElbowPlot(AB.combined)
AB.combined <- RunUMAP(AB.combined, reduction = "pca", dims = 1:30)

AB.combined <- FindNeighbors(AB.combined,dims = 1:30)
AB.combined <- FindClusters(AB.combined,resolution = 0.3)

# AB.combined <- FindClusters(AB.combined,resolution = 0.8)




###annotate the cells using scType

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
#tissue = "Heart" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = AB.combined[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(AB.combined@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(AB.combined@meta.data[AB.combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(AB.combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


AB.combined@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  AB.combined@meta.data$customclassif[AB.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


DimPlot(AB.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  


DefaultAssay(AB.combined) <- 'RNA'
Idents(AB.combined) <- AB.combined$customclassif

DATA.NK <- subset(AB.combined,subset = customclassif =='Natural killer  cells')

DATA.NK <- DATA.NK %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>%
  RunUMAP(dims = 1:20)

DATA.NK <- FindNeighbors(DATA.NK,dims = 1:20)

DATA.NK <- FindClusters(DATA.NK,resolution = 0.3) 

DimPlot(DATA.NK, label = T)

DATA.NK.sub <- subset(DATA.NK, subset = seurat_clusters %in% c(0,1,2,3,4))
Idents(DATA.NK.sub) <- DATA.NK.sub$seurat_clusters
DotPlot(DATA.NK.sub,features = c('Cd69','Eomes','Bcl2','Gzmb'))


#### heatmap 
GENES <- c('Il2ra','Il2rb','Il4r','Il10ra','Il17ra','Il12rb1','Il15ra','Ifngr2','Ifngr1','Ifnar1','Ifnar2')
#GENES <- c('Cd69','Itgam','Sell','Gzmb','Bcl2','Eomes','Tbx21')


averaged <- AverageExpression(DATA.NK.sub, assays = 'RNA',group.by = 'seurat_clusters')

PLOT <- averaged[[1]]
idx <- which(rownames(PLOT)%in% GENES)
FIND <- rownames(PLOT)[idx]
setdiff(GENES,FIND)
PLOT.sub <- as.matrix(PLOT[idx,])

temp <- apply(PLOT.sub,1,sum)
PLOT.sub <- PLOT.sub[-which(temp==0),]

pheatmap(PLOT.sub,scale='row')


################## Figure 7
#### Figure 7C: reanalysis of spatial data (zhang et al.)

# create seurat object

library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)


dir.in <- '/fs/ess/PAS2579/NoJoon/GSE238264_HCC_spatial/GSE238264_RAW/'

SAMPLES <- c('HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR')

for (i in 1:length(SAMPLES)){
  dir <- paste0(dir.in, SAMPLES[i])
  EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
  IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
  coord <- IMAGE@coordinates
  coord <- coord[order(match(rownames(coord),colnames(EXP))),]
  identical(rownames(coord),colnames(EXP))
  DATA <- CreateSeuratObject(counts = EXP,project = SAMPLES[i],assay = 'Spatial')
  DATA@images$image <- new(Class = 'SlideSeq', assay = 'Spatial', key = 'image_', coordinates = coord[, 2:3])
  saveRDS(DATA, file = paste0(dir,'/',SAMPLES[i],'_seurat.rds'))
}



## non-responder: HCC6NR
i <-6

dir <- paste0(dir.in, SAMPLES[i])
dir
EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
coord <- IMAGE@coordinates
coord <- coord[order(match(rownames(coord),colnames(EXP))),]

CD8.COUNT <- EXP[which(rownames(EXP) %in% c('CD3E','CD3D','CD247','CD8B','CD8A','LCK')),]
CD8.COUNT.sum <- colSums(CD8.COUNT)
summary(CD8.COUNT.sum)

NK.COUNT <- EXP[which(rownames(EXP) %in% c('KLRD1','KLRF1','GNLY','CD7','NCR1','NKG7','KLRK1')),]

NK.COUNT.sum <- colSums(NK.COUNT)
summary(NK.COUNT.sum)

# CD8+ spots
IDX2 <- rep(0,length(CD8.COUNT.sum))
IDX2[which(CD8.COUNT.sum>=1)] <-1
IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "red"))+
  ggtitle(paste0(SAMPLES[i],'_CD8','_cutoff=1'))+theme(legend.position="none")


# NK+ spots
IDX2 <- rep(0,length(NK.COUNT.sum))
IDX2[which(NK.COUNT.sum>=1)] <-1
IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "blue"))+
  ggtitle(paste0(SAMPLES[i],'_NK','_cutoff=1'))+theme(legend.position="none")


## responder:HCC2R

i <- 2

dir <- paste0(dir.in, SAMPLES[i])
dir
EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
coord <- IMAGE@coordinates
coord <- coord[order(match(rownames(coord),colnames(EXP))),]

CD8.COUNT <- EXP[which(rownames(EXP) %in% c('CD3E','CD3D','CD247','CD8B','CD8A','LCK')),]
CD8.COUNT.sum <- colSums(CD8.COUNT)
summary(CD8.COUNT.sum)

NK.COUNT <- EXP[which(rownames(EXP) %in% c('KLRD1','KLRF1','GNLY','CD7','NCR1','NKG7','KLRK1')),]

NK.COUNT.sum <- colSums(NK.COUNT)
summary(NK.COUNT.sum)

IDX2 <- rep(0,length(CD8.COUNT.sum))
IDX2[which(CD8.COUNT.sum>=1)] <-1


IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

# CD8+ spots
ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "red"))+
  ggtitle(paste0(SAMPLES[i],'_CD8','_cutoff=1'))+theme(legend.position="none")


# NK+ spots
IDX2 <- rep(0,length(NK.COUNT.sum))
IDX2[which(NK.COUNT.sum>=1)] <-1


IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "blue"))+
  ggtitle(paste0(SAMPLES[i],'_NK','_cutoff=1'))+theme(legend.position="none")


## ratio of CD8+ spots to NK+ spots across samples

RATIO <- vector()
for (i in 1:length(SAMPLES)){
  dir <- paste0(dir.in, SAMPLES[i])
  dir
  EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
  IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
  coord <- IMAGE@coordinates
  coord <- coord[order(match(rownames(coord),colnames(EXP))),]

  CD8.COUNT <- EXP[which(rownames(EXP) %in% c('CD3E','CD3D','CD247','CD8B','CD8A','LCK')),]
  CD8.COUNT.sum <- colSums(CD8.COUNT)
  summary(CD8.COUNT.sum)

  NK.COUNT <- EXP[which(rownames(EXP) %in% c('KLRD1','KLRF1','GNLY','CD7','NCR1','NKG7','KLRK1')),]

  NK.COUNT.sum <- colSums(NK.COUNT)
  RATIO[i] <- length(which(CD8.COUNT.sum>=1))/length(which(NK.COUNT.sum>=1))
  
}

DF <- data.frame(ratio = RATIO)
DF$group <- rep(c('R','NR'), c(4,3))
ggplot(DF, aes(x = group, y = ratio)) +  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1)) + theme_classic()


## Figure 7D: reanalysis of scRNA-seq data (zheng et al., pancancer)

CD8_integrated <- readRDS('/fs/ess/PAS2579/NoJoon/pancancer/int.CD8.S35.sce.merged.rds')
UMAP.harmony <- as.data.frame(reducedDims(CD8_integrated)[['harmony.umap']])
UMAP.harmony$group <- as.character(CD8_integrated$meta.cluster)

ggplot(UMAP.harmony, aes(harmony.umap_1, harmony.umap_2, color= group))+ geom_point() + theme_classic()

GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')  # alternative set2
IL2 <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IL2_pathway.csv')
IFNalpha <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IFNalpha_pathway.csv')
RAR.target <- read.csv('/fs/ess/PAS2579/NoJoon/set3.csv',header=T)    # IRIS3 RAR target


IL2 <- IL2$Symbol
IL2 <- IL2[nzchar(IL2)]

IFNalpha <- IFNalpha$Symbol
IFNalpha <- IFNalpha[nzchar(IFNalpha)]

## if IRIS3 RAR
RAR <- RAR.target$Set3_RAR


count <- CD8_integrated@assays$data$exprs

length(intersect(IL2,toupper(GENE$gene)))  # 0
length(intersect(IFNalpha,toupper(GENE$gene)))  # 0
length(intersect(toupper(RAR),toupper(GENE$gene)))  # 0


TARGET.GENE <- c(GENE$gene,IL2,IFNalpha,toupper(RAR))

TARGET.count <- count[which(toupper(rownames(count))%in% TARGET.GENE),]
TARGET.count <- t(TARGET.count)



a <- vector()

IND <- which(colnames(TARGET.count)%in% toupper(GENE$gene))

for(i in 1:ncol(TARGET.count)){
  if (i %in% IND){
    temp <- colnames(TARGET.count)[i]
    tt <- GENE[which(toupper(GENE$gene)==temp),]$exp
    a[i] <- ifelse(tt=='High',0.1,0.1)
  }else{
    a[i] <- 0.1
  }
}



a <- as.matrix(a)
TARGET.score <- TARGET.count %*% a




summary(TARGET.score[,1])  
hist(TARGET.score[,1])

TT <-rownames(TARGET.count)[-which(is.na(TARGET.score))]

TARGET.score <- TARGET.score[-which(is.na(TARGET.score))]
names(TARGET.score) <- TT
UMAP.sub <- UMAP.harmony[which(rownames(UMAP.harmony)%in% names(TARGET.score)),]


TARGET.score.mod <- TARGET.score

for (i in 1:length(TARGET.score.mod)){
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]>2,2,TARGET.score.mod[i])
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]<(-2),(-2),TARGET.score.mod[i])
}

hist(TARGET.score.mod)


df <- data.frame(UMAP.sub,TARGET.score,TARGET.score.mod)

mid <- median(TARGET.score.mod)
ggplot(df,aes( harmony.umap_1, harmony.umap_2,color=TARGET.score.mod))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()


## extract the signature from pancancer data
temp <- CD8_integrated@assays$data$exprs
CV <- apply(temp,1,sd)
INDS <- which(is.na(CV))

CD8_sub <- CD8_integrated[-INDS,]
stats <- modelGeneCV2(CD8_sub,assay.type = 'exprs')
HIGH <- getTopHVGs(stats,var.field = 'ratio')
FINAL <- union(HIGH,TARGET.GENE)
FINAL <- FINAL[which(FINAL%in% rownames(CD8_sub))]
CD8_sub <- CD8_sub[FINAL,]

left <- setdiff(1:nrow(df),INDEX)
library(scran)

RST <- list()
for (i in 10001:10020){
  set.seed(i)
  left.sub <- sample(left,1000)
  DF.sub <- DF[c(left.sub,INDEX),]
  #CD8.sub <- CD8_integrated[,which(colnames(CD8_integrated)%in% rownames(DF.sub))]
  CD8.sub <- CD8_sub[,which(colnames(CD8_sub)%in% rownames(DF.sub))]
  DF.sub <-DF.sub[order(match(rownames(DF.sub),colnames(CD8.sub))),]
  
  #LOW <- findMarkers(CD8.sub,groups = DF.sub$A,assay.type = 'exprs')
  #LOW <- findMarkers(CD8.sub,groups = as.factor(DF.sub$A),assay.type = 'exprs',test.type='t')
  LOW <-scoreMarkers(CD8.sub,groups = DF.sub$A,assay.type = 'exprs')
  LOW[[1]][order(LOW[[1]]$mean.AUC, decreasing=TRUE),1:4]   # this is the upregulated genes in group1
  RST[[i-10000]] <- LOW
  
}

saveRDS(RST,file = '/fs/ess/PAS2579/NoJoon/pancancer/low.rds')


low <-readRDS('/fs/ess/PAS2579/NoJoon/pancancer/low.rds')


MARKERS <- function(x){
  temp <- x[[1]]
  t2 <-temp[order(temp$mean.AUC, decreasing=TRUE),] 
  rownames(t2[1:50,])
}

low.marker <-lapply(low,MARKERS)
markers1 <- Reduce(intersect, low.marker)


c10 <- readRDS('/fs/ess/PAS2579/NoJoon/pancancer/high_C10.rds')

c10.marker <-lapply(c10,MARKERS)
markers.c10 <- Reduce(intersect, c10.marker)

c789 <- readRDS('/fs/ess/PAS2579/NoJoon/pancancer/high_C7-9.rds')

c789.marker <-lapply(c789,MARKERS)
markers.c789 <- Reduce(intersect, c789.marker)


####### pre
yost.pre.CD8 <- readRDS('/fs/ess/PAS2579/NoJoon/Yost/yost.pre.CD8.rds')
responder <- c('su001','su002','su003','su004','su009','su012')
noResponder <- c('su005','su006','su007','su008','su010')


yesCell <- colnames(yost.pre.CD8)[which(yost.pre.CD8$patient %in% responder)] 
NoCell <- colnames(yost.pre.CD8)[which(yost.pre.CD8$patient %in% noResponder)] 

count.data <- GetAssayData(yost.pre.CD8)

## low
low.count <- count.data[which(rownames(count.data)%in% markers1),]
low.count <- t(low.count)
a <- rep(1,ncol(low.count))
a <- as.matrix(a)
low.score <- low.count%*%a

summary(low.score[,1])  # there are ~700 NaN
hist(low.score[,1])

low.score.mod <- low.score[,1]
hist(low.score.mod)

for (i in 1:length(low.score.mod)){
  #low.score.mod[i] <- ifelse(low.score.mod[i]>45,45,low.score.mod[i])
  low.score.mod[i] <- ifelse(low.score.mod[i]<20,20,low.score.mod[i])
}

hist(low.score.mod)




UMAP <- yost.pre.CD8@reductions$umap@cell.embeddings

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(low.score)),]

df <- data.frame(UMAP.sub,low.score,low.score.mod)
RESPONSE <- rep(-1,nrow(df))
for (i in 1:nrow(df)){
  RESPONSE[i] <- ifelse(rownames(df)[i] %in% yesCell,1,0)
}

df$responder <- RESPONSE


ggplot(df, aes(x=as.factor(RESPONSE), y=low.score.mod)) + 
  geom_boxplot()+theme_classic()+ggtitle('Yost_pre_low')


##
c10.count <- count.data[which(rownames(count.data)%in% markers.c10),]
c10.count <- t(c10.count)
a <- rep(1,ncol(c10.count))
a <- as.matrix(a)
c10.score <- c10.count%*%a

summary(c10.score[,1])  # there are ~700 NaN
hist(c10.score[,1])

c10.score.mod <- c10.score[,1]
hist(c10.score.mod)

for (i in 1:length(c10.score.mod)){
  c10.score.mod[i] <- ifelse(c10.score.mod[i]>80,80,c10.score.mod[i])
  c10.score.mod[i] <- ifelse(c10.score.mod[i]<20,20,c10.score.mod[i])
}

hist(c10.score.mod)




UMAP <- yost.pre.CD8@reductions$umap@cell.embeddings

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(c10.score)),]

df <- data.frame(UMAP.sub,c10.score,c10.score.mod)
RESPONSE <- rep(-1,nrow(df))
for (i in 1:nrow(df)){
  RESPONSE[i] <- ifelse(rownames(df)[i] %in% yesCell,1,0)
}

df$responder <- RESPONSE

ggplot(df, aes(x=as.factor(RESPONSE), y=c10.score.mod)) + 
  geom_boxplot()+theme_classic()+ggtitle('Yost_pre_c10')

df.sub <- df[which(df$c10.score>=37&df$c10.score<=65),]

ggplot(df.sub, aes(x=as.factor(responder), y=c10.score.mod)) + 
  geom_boxplot()+theme_classic()+ggtitle('Yost_pre_c10')

t.test(df.sub[which(df.sub$responder==1),]$c10.score.mod,df.sub[which(df.sub$responder==0),]$c10.score.mod,var.equal = T)



############## Figure S11C: scRNA-seq DATA (Daniel et al)

library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)

DAT <- Read10X(data.dir = '/fs/ess/PAS2579/NoJoon/GSE188666/out',gene.column=1)
DAT <- CreateSeuratObject(counts = DAT,project = "GSE18866",min.cells = 0, min.features = 0)

a <-sapply(strsplit(colnames(DAT),'_',1),'[',3)
table(a)  # contains D8 and D21

DAT$time <- a


meta <- read.table('/fs/ess/PAS2579/NoJoon/GSE188666/out/GSE188666_scRNA_LCMV_metadata.tsv')
table(meta$tissue)   # liver, spleen and lung

identical(colnames(DAT),rownames(meta))  # check if the order of cells are the same

DAT <- AddMetaData(DAT, meta)

Idents(DAT) <- DAT$ident
UMAP <- data.frame(UMAP_1=DAT$UMAP_1,UMAP_2=DAT$UMAP_2)
UMAP <- as.matrix(UMAP)
DAT[['umap']] <-CreateDimReducObject(embeddings = UMAP,key = 'UMAP_',assay = DefaultAssay(DAT))
DAT$clusters <- DAT$ident


DimPlot(DAT,label=T)



count <- DAT@assays$RNA@counts
UMAP <- DAT@reductions$umap@cell.embeddings

### NKRT signature
GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')   # alternative set2
GENE.count <- count[which(toupper(rownames(count))%in% toupper(GENE$gene)),]
GENE.count <- t(GENE.count)

a <- rep(1,ncol(GENE.count))
a <- as.matrix(a)
GENE.score <- GENE.count%*%a

summary(GENE.score[,1])  
hist(GENE.score[,1])
hist(log(GENE.score[,1]+1))


rownames(GENE.score) <- rownames(GENE.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(GENE.score)),]

df <- data.frame(UMAP.sub,GENE.score, log(GENE.score+1))

mid <- median(log(GENE.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(GENE.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()



## RAR

RAR.target <- read.csv('/fs/ess/PAS2579/NoJoon/set3.csv',header=T)

RAR.count <- count[which(rownames(count)%in% RAR.target$Set3_RAR),]
RAR.count <- t(RAR.count)

a <- rep(1,ncol(RAR.count))
a <- as.matrix(a)
RAR.score <- RAR.count%*%a

summary(RAR.score[,1])


rownames(RAR.score) <- rownames(RAR.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(RAR.score)),]

df <- data.frame(UMAP.sub,RAR.score)
mid <- median(log(RAR.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(RAR.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()


## IL2 
IL2 <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IL2_pathway.csv')

IL2 <- IL2$Symbol
IL2 <- IL2[nzchar(IL2)]

IL2.count <- count[which(toupper(rownames(count))%in% IL2),]
IL2.count <- t(IL2.count)

a <- rep(1,ncol(IL2.count))
a <- as.matrix(a)
IL2.score <- IL2.count%*%a

summary(IL2.score[,1])


rownames(IL2.score) <- rownames(IL2.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(IL2.score)),]

df <- data.frame(UMAP.sub,IL2.score)
mid <- median(log(IL2.score+1))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(IL2.score+1)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()





## IFNalpha
IFNalpha <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IFNalpha_pathway.csv')
IFNalpha <- IFNalpha$Symbol
IFNalpha <- IFNalpha[nzchar(IFNalpha)]

IFN.count <- count[which(toupper(rownames(count))%in% IFNalpha),]
IFN.count <- t(IFN.count)

a <- rep(1,ncol(IFN.count))
a <- as.matrix(a)
IFN.score <- IFN.count%*%a

summary(IFN.score[,1])


rownames(IFN.score) <- rownames(IFN.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(IFN.score)),]

df <- data.frame(UMAP.sub,IFN.score)
mid <- median(log(IFN.score+1))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(IFN.score+1)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()


##### combine RAR and NKRT together


GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')  # alternative set2

TARGET.GENE <- c(GENE$gene,IL2,IFNalpha,toupper(RAR))

TARGET.count <- count[which(toupper(rownames(count))%in% TARGET.GENE),]
TARGET.count <- t(TARGET.count)


a <- vector()

IND <- which(colnames(TARGET.count)%in% toupper(GENE$gene))

for(i in 1:ncol(TARGET.count)){
  if (i %in% IND){
    temp <- colnames(TARGET.count)[i]
    tt <- GENE[which(toupper(GENE$gene)==temp),]$exp
    a[i] <- ifelse(tt=='High',0.1,0.1)
  }else{
    a[i] <- 0.1
  }
}



a <- as.matrix(a)
TARGET.score <- TARGET.count %*% a

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(TARGET.score)),]
TARGET.score.mod <- log(TARGET.score[,1])


df <- data.frame(UMAP.sub,TARGET.score,TARGET.score.mod)

mid <- median(TARGET.score.mod)
ggplot(df,aes(UMAP_1, UMAP_2,color=TARGET.score.mod))+geom_point(size=0.2) +
  scale_colour_gradient2(low = "royalblue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()

