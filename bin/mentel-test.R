
#install.packages("vegan")
library(vegan)

lis<-list.files()

mat<-matrix(0,5074,25)
for (i in 1:25) {
  fst<-read.table(lis[i],header = T)
  mat[,i]<-fst$fst
}
write.csv(mat,file = "gro.csv",quote = F)

fst<-read.csv("gro.csv",header = T)
nfst<-fst[,-1]/(1-fst[,-1])
nfst$gro<-fst$gro

IBE<-read.csv("IBE.csv",header = T)
IBD<-read.csv("IBD.csv",header = T)
#IBD
#25
veg.dist <- vegdist(nfst[,c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBD[,-c(1,2)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBD.mentel.25.csv",quote = F,row.names = F)
#group1
veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 1)],c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBD[IBE$id[which(IBE$gro %in% 1)],-c(1,2)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBD.mentel.group1.csv",quote = F,row.names = F)
#group2
veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 2)],c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBD[IBE$id[which(IBE$gro %in% 2)],-c(1,2)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBD.mentel.group2.csv",quote = F,row.names = F)
#group3
veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 3)],c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBD[IBE$id[which(IBE$gro %in% 3)],-c(1,2)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBD.mentel.group3.csv",quote = F,row.names = F)
#IBE
veg.dist <- vegdist(nfst[IBE$id,c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBE[IBD$id,-c(1,2,3)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBE.mentel.group25.csv",quote = F,row.names = F)

#group1
veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 1)],c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBE[which(IBE$gro %in% 1),-c(1,2,3)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBE.mentel.group1.csv",quote = F,row.names = F)
#group2
veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 2)],c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBE[which(IBE$gro %in% 2),-c(1,2,3)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBE.mentel.group2.csv",quote = F,row.names = F)

#group2
veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 3)],c(1:5074)]) # Bray-Curtis
env.dist <- vegdist(IBE[which(IBE$gro %in% 3),-c(1,2,3)])
menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
result<-cbind(menr$statistic,menr$signif)
colnames(result)<-c("r","Pvalue")
write.csv(result,file = "IBE.mentel.group3.csv",quote = F,row.names = F)

save.image("mentel.RDdata")

#mantel.partial(env.dist,veg.dist, method = "pearson", permutations = 999)
#group3
result<-matrix(0,6,2)
for (i in 4:9){
	veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 3)],c(1:5074)]) # Bray-Curtis
	env.dist <- vegdist(IBE[which(IBE$gro %in% 3),i])
	menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
	result[i-3,]<-cbind(menr$statistic,menr$signif)
	colnames(result)<-c("r","Pvalue")
	write.csv(result,file = paste("IBE.mentel.group3","csv",sep="."),quote = F,row.names = F)
}


#group2
result<-matrix(0,6,2)
for (i in 4:9){
	veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 2)],c(1:5074)]) # Bray-Curtis
	env.dist <- vegdist(IBE[which(IBE$gro %in% 2),i])
	menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
	result[i-3,]<-cbind(menr$statistic,menr$signif)
	colnames(result)<-c("r","Pvalue")
	write.csv(result,file = paste("IBE.mentel.group2","csv",sep="."),quote = F,row.names = F)
}

#group1
result<-matrix(0,6,2)
for (i in 4:9){
	veg.dist <- vegdist(nfst[IBE$id[which(IBE$gro %in% 1)],c(1:5074)]) # Bray-Curtis
	env.dist <- vegdist(IBE[which(IBE$gro %in% 1),i])
	menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
	result[i-3,]<-cbind(menr$statistic,menr$signif)
	colnames(result)<-c("r","Pvalue")
	write.csv(result,file = paste("IBE.mentel.group1","csv",sep="."),quote = F,row.names = F)
}

#group1
result<-matrix(0,6,2)
for (i in 4:9){
	veg.dist <- vegdist(nfst[IBE$id,c(1:5074)]) # Bray-Curtis
	env.dist <- vegdist(IBE[,i])
	#env.dist <- vegdist(IBE[,i])
	menr<-mantel(env.dist,veg.dist, method="spearman", permutations=1000)
	result[i-3,]<-cbind(menr$statistic,menr$signif)
	colnames(result)<-c("r","Pvalue")
	write.csv(result,file = paste("IBE.mentel.group25","csv",sep="."),quote = F,row.names = F)
}


save.image("mentel.RDdata")



