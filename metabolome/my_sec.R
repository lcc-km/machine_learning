library(Rtsne)
library(umap)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)

# -- load dataset ####
write.path<- "/Users/lucc/Desktop/zyf/jak_stat/"
setwd(write.path)
sce <- readRDS(paste(write.path,"JAK_STAT_Norm_P7.P11.30.38.96_tsne1.rds",sep=""))
###
### 获得原始表达矩阵 ##
counts <- counts(sce)
#write.table(counts,'row_count.txt',quote = FALSE,sep = "\t")
#log
counts.LogNormalize <- log1p(t(t(counts) / colSums(counts)) * 10000)
# scale
scaled_data <- t(scale(t(counts.LogNormalize)))
####### pca
pca.results <- prcomp(x = scaled_data,rank. = 10)
pac.matrix <- pca.results$rotation

######## KNN  FindNeighbors
#1 使用10个PC
a=new(Class = RcppAnnoy::AnnoyEuclidean, f=10)
#2 求距离的输入数据为每个细胞的RNA表达量的前10个PC
pac.matrix
for(ii in seq(nrow(pac.matrix)) ){
  a$addItem(ii-1, pac.matrix[ii,])
}
#3 建立树的索引，默认建立50棵树
a$build(50) 
#4 对于每个细胞，求其前k个最近邻细胞
index=a #第3中建立的索引
query=pac.matrix #要搜索的数据，就是每个细胞的PC矩阵
k=20 #最近邻的几个点？
search.k=-1 #不限制搜索次数
include.distance=T #是否返回距离
nn <- matrix(data=rep(0, nrow(pac.matrix)**2),nrow=nrow(pac.matrix), ncol=nrow(pac.matrix))
dim(nn) 
for(x in 1:nrow(pac.matrix)){
  res <- index$getNNsByVectorList(pac.matrix[x, ], k, search.k, include.distance)
  res2=list(res$item + 1, res$distance) #C++下标是0-based，变为R的1-based
  nn[x, res2[[1]]] = 1
}
rownames(nn)=rownames(pac.matrix)
colnames(nn)=rownames(pac.matrix)

########### FindClusters
class(nn)
nn2=as(nn, "dgCMatrix")
# check: 检测大小
object.size(nn) 
object.size(nn2)
# 转为3列
df_nn=as.data.frame(summary(nn2))
# Louvain  算法
library(igraph)
edges <- df_nn
colnames(edges) <- c("from", "to", "weight")
#Create graph for the algorithms
g <- graph_from_data_frame(edges, directed = FALSE)
# Louvain
lc <- cluster_louvain(g)
head(membership(lc))
str(lc)
communities(lc)
plot(lc, g)
# make data.frame from Louvain result
lc_df=data.frame(cell=as.numeric(lc$names), cid=lc$membership)
lc_df=lc_df[order(lc_df$cell),]


#tsne
tsne_out = Rtsne(
  pac.matrix,
  dims = 2,
  theta = 0.1
)
tsne_zz <- tsne_out$Y
t.data_zz <- data.frame(lc_df,tsne_zz)
colnames(t.data_zz) <- c("cell", "state",  "tSNE1","tSNE2") 
ggplot(t.data_zz,aes(x=tSNE1,y=tSNE2))+
  geom_point(aes(color=factor(state)),size=2)+
  theme_classic()
#umap
umap_out <- umap(pac.matrix)
umap_zz <- umap_out$layout
umap.data_zz <- data.frame(lc_df,umap_zz)
colnames(umap.data_zz) <- c("cell", "state",  "umap1","umap2") 
ggplot(umap.data_zz,aes(x=umap1,y=umap2))+
  geom_point(aes(color=factor(state)),size=2)+
  theme_classic()

library(uwot)
result=umap(
	X = pac.matrix, #矩阵: 前dims个PC列
	#n_threads = nbrOfWorkers(), #线程数
	n_neighbors = as.integer(x = 30), #默认30
	n_components = as.integer(x = 2), #默认2
	metric = 'cosine', #这个参数
	n_epochs = NULL,
	learning_rate = 1,
	min_dist = 0.3,
	spread = 1,
	set_op_mix_ratio = 1,
	local_connectivity = 1,
	repulsion_strength = 1,
	negative_sample_rate = 5,
	a = NULL,
	b = NULL,
	fast_sgd = F, #uwot.sgd
	verbose = T,
	ret_model = F #return.model
)
umap_zz <- result
umap.data_zz2 <- data.frame(lc_df,umap_zz)
colnames(umap.data_zz2) <- c("cell", "state",  "umap1","umap2") 
ggplot(umap.data_zz2,aes(x=umap1,y=umap2))+
  geom_point(aes(color=factor(state)),size=2)+
  theme_classic()


