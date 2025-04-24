data=read.csv("all.bins.abu",header = 1,row.names = 1,check.names = F,sep='\t')
drep_bins=read.csv("drep_bins_name.txt",header = F,check.names = F,sep='\t')
drep_bins$variable="TRUE"
rownames(drep_bins)=drep_bins$V1
data$variable=drep_bins[rownames(data),2]
data=na.omit(data)
data=data[,-25]
write.table(x = data,"all.drep_bins.abu.tsv",sep = '\t',row.names = T,quote = F)
group=read.table("../metadata.txt",header = 1,stringsAsFactors = T)

library(ggplot2)
library(ggsci)
relative_abundance=function(d){
  dta_sum=apply(d,2,function(x){x/sum(x)})
}
abu=as.data.frame(relative_abundance(data))

bacteria_taxonomy=read.table("prokaryotic_species_clusters_taxonomy.tsv",sep='\t',row.names = 1,header = 1)
phylum=strsplit(bacteria_taxonomy$classification,split = ";")
phylum=do.call(rbind,phylum)
bacteria_taxonomy$phylum=phylum[,2]
bacteria_taxonomy$class=phylum[,3]
##PCoA
library(dplyr)
library(vegan)
results5=t(abu)
results5=as.matrix(results5)
distance <- as.matrix(vegdist(results5, method= "bray",na.rm = T))
colnames(distance)=rownames(results5)
rownames(distance)=rownames(results5)
distance2=matrix(NA,nrow(distance),ncol(distance))
for(i in 1:nrow(distance)){
  for(j in 1:ncol(distance)){
    if ( i != j){
      distance2[i,j]=distance[i,j]
    }
  }
}

rownames(distance2)=rownames(distance)
colnames(distance2)=colnames(distance)
viral_distance=read.table("../viral_distance.tsv",sep='\t',check.names = F)
mantel(distance2,ydis = viral_distance,permutations = 9999,method = "spearman")

library(pheatmap)
pheatmap(distance2,cluster_rows = F,cluster_cols = F,clustering_method = "average",
         na_col = "lightblue",border_color = F,cellwidth = 10,cellheight = 8,
         color = colorRampPalette(c("white",
                                    "white",
                                    "white",
                                    "white",
                                    "#40CCA7"))(30))

library(reshape2)
bac_dist=melt(as.matrix(distance2))
vir_dist=melt(as.matrix(viral_distance))
bac_dist$vir_dist=vir_dist$value
ggplot(bac_dist,aes(value,vir_dist))+
  geom_point(alpha=0.8,shape=21,fill="green")+
  geom_smooth(formula = y~x,method = "lm",inherit.aes = T,color="red")+
  theme_bw()+labs(x="Bray distance among genome-centric \nmicrobial communities",
                  y="Bray distance among \nviral communities")

summary(lm(formula = value~vir_dist,data = bac_dist))

pcoa_v <- cmdscale(distance, k = (nrow(results5) - 1), eig = TRUE)
point_v <- data.frame(pcoa_v$point)
pcoa_eig_v <- (pcoa_v$eig)[1:2] / sum(pcoa_v$eig)
sample_site_v <- data.frame ({pcoa_v$point})[1:3]
sample_site_v$names <- rownames(sample_site_v)
names(sample_site_v)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')
rownames(group)=group$sample
sample_site_v$site = group[rownames(sample_site_v),2]
sample_site_v$stage=rep(c("early","late"),c(12,12))
write.table(sample_site_v,"viral_pcoa123.txt",sep = '\t',quote = F)

library(vegan)

dist = vegdist(t(abu), method = 'bray',na.rm = T)

p_value=anosim(x= dist,grouping = sample_site_v$stage,permutations = 9999)
summary(p_value)
p_value=anosim(x= dist,grouping = sample_site_v$site,permutations = 9999)
summary(p_value)
p_value_1=adonis2(formula = dist~sites,group,permutations = 9999)
p_value_1
p_value_2=adonis2(formula = dist~year,group,permutations = 9999)
p_value_2
p_value_3=adonis2(formula = dist~stage,sample_site_v,permutations = 9999)
p_value_3

dist_year=vegdist(group$year)
mantel(distance,dist_year)
p <- ggplot(sample_site_v, aes(PCoA1, PCoA2)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(aes(fill=site,shape=stage),size = 6, alpha = 0.8) + 
  scale_fill_d3() + 
  scale_shape_manual(values = c(21,24))+
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig_v[1], 2), '%'), 
       y = paste('PCoA axis2: ', round(100 * pcoa_eig_v[2], 2), '%'))+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = site),level = 0.95, show.legend = TRUE,inherit.aes = F)+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = site),level = 0.975, show.legend = TRUE,inherit.aes = F)+
  annotate('text', colour="#8766d")+
  guides(fill = guide_legend(override.aes = list(size = 5, shape = 21)))
p

##Candidatus_Methylomirabilis
tmp=strsplit(bacteria_taxonomy$classification,split = ";")
tmp=do.call(rbind,tmp)
bacteria_taxonomy$family=tmp[,5]
bacteria_taxonomy$order=tmp[,4]
abu$class=bacteria_taxonomy[rownames(abu),22]
abu$family=bacteria_taxonomy[rownames(abu),23]
abu$order=bacteria_taxonomy[rownames(abu),24]
Methylomirabilis=abu[grep("Methylomi",abu$class),]
apply(Methylomirabilis[,1:24], 2, sum)
library(reshape2)
Methylomirabilis$class=as.factor(Methylomirabilis$class)
Methylomirabilis=melt(as.matrix(Methylomirabilis[,1:24]))

library(ggplot2)
library(ggsci)
ggplot(Methylomirabilis,aes(x = Var2,y=value*100,fill=Var1))+
  geom_bar(stat = "identity",
           alpha=0.7,width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  scale_fill_aaas()+
  labs(y="Relative abundance (%)")

Nitrososphaeria=abu[grep("Nitrososphaeria",abu$class),]
apply(Nitrososphaeria[,1:24], 2, sum)
library(reshape2)
Nitrososphaeria=melt(as.matrix(Nitrososphaeria[,1:24]))
ggplot(Nitrososphaeria,aes(x = Var2,y=value*100,fill=Var1))+
  geom_bar(stat = "identity",
           alpha=0.7,width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  scale_fill_lancet()+
  labs(y="Relative abundance (%)")

Nitrospiria=abu[grep("Nitrospiria",abu$class),]
apply(Nitrospiria[,1:24], 2, sum)
library(reshape2)
library(ggplot2)
library(ggsci)
Nitrospiria=melt(as.matrix(Nitrospiria[,1:24]))
ggplot(Nitrospiria,aes(x = Var2,y=value*100,fill=Var1))+
  geom_bar(stat = "identity",
           alpha=0.7,width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  scale_fill_d3()+
  labs(y="Relative abundance (%)")

ggsave("Nitrospiria.pdf",width = 7.8,height = 4)

abu$phylum=bacteria_taxonomy[rownames(abu),21]
Desulfobacterota=abu[grep("Desulfobacterota",abu$phylum),]
Desulfobacterota=melt(as.matrix(Desulfobacterota[,c(1:24)]))
ggplot(Desulfobacterota,aes(x = Var2,y=value*100,fill=Var1))+
  geom_bar(stat = "identity",
           alpha=0.7,width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  scale_fill_d3(palette = "category20c")+
  labs(y="Relative abundance (%)")

Sulfur_oxider=abu[c(grep("Acidiferrobacterales",abu$order),grep("Sulfuricellaceae",abu$family)),]
apply(Sulfur_oxider[,1:24],2,sum)
library(reshape2)
Sulfur_oxider=melt(Sulfur_oxider[,c(1:24,29)])
ggplot(Sulfur_oxider,aes(x = variable,y=value*100,fill=order))+
  geom_bar(stat = "identity",
           alpha=0.7,width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  scale_fill_d3()+
  labs(y="Relative abundance (%)")

##share
abu_vc_2=abu
for (i in 1:ncol(abu_vc_2)){
  col=as.matrix(abu_vc_2[,i])
  rownames(col)=rownames(abu_vc_2)
  col=as.matrix(col[col[,1]>1,])
  colnames(col)=colnames(abu_vc_2)[i]
  nd=as.matrix(new.cbind(nd,col))
}

for(i in 1:nrow(nd)){
  for (j in 1:ncol(nd)){
    if(is.na(nd[i,j])==FALSE){
      nd[i,j]=1
    }
    else{
      nd[i,j]=0
    }
  }
}

share=matrix(NA,ncol(nd),ncol(nd))
for (i in 1:ncol(nd)){
  for (j in 1:ncol(nd)){
    print(j)
    sum=as.matrix(apply(nd[,c(i,j)],1,sum))
    shared=as.data.frame(nrow(nd[sum==2,]))
    if(nrow(shared)>0){
      a=sum(nd[,i])
      b=sum(nd[,j])
      viral_shared_content=((shared/a)+(shared/b))/2
      share[i,j]=viral_shared_content[1,1]
    }else{
      share[i,j]==0
    }
  }
}

colnames(share)=colnames(nd)
rownames(share)=colnames(nd)
share[is.na(share)]=0

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
share1=get_lower_tri(share)

group$stage=rep(x = c("early","later"),c(12,12))

library(reshape2)
share1=melt(share1)
share1=na.omit(share1)
share1=share1[share1$value!=1,]
share1=as.data.frame(share1)

nd=matrix(NA,nrow(share1),1)
for (i in 1:nrow(share1)){
  print(i)
  name=as.matrix(group[grep(share1[i,1],group$sample),4])
  name_1=as.matrix(group[grep(share1[i,2],group$sample),4])
  if ( name_1[1,1]==name[1,1]){
    if (name_1[1,1]=='early'){
     nd[i,1]='intra_early_stage' 
    }
    if (name_1[1,1]=="later"){
      nd[i,1]='intra_later_stage'
    }
  }else{
    nd[i,1]='inter_stage'
  }

}

share1$zone=nd
library(ggplot2)
library(ggpubr)
share1=as.data.frame(share1)
p=ggplot(share1,aes(zone,value))+
  geom_boxplot(fill="lightblue",outlier.size = 0,outlier.colour = 'white')+
  geom_jitter(width = 0.1,alpha=0.7,size=2)+
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4)
p
p+stat_compare_means(value ~ group, data = share1,
                     method = "wilcox.test",label="p.signif",paired = F,
                     comparisons =list(c("intra_early_stage","inter_stage"
                     ),c("intra_later_stage","intra_early_stage")))

share1=get_lower_tri(share)

library(pheatmap)
annotation_row=as.data.frame(group_non[,c(1,2)])
rownames(annotation_row)=group_non[,4]
colnames(annotation_row)=c("type","zone")

pheatmap(share,
         border_color = "grey",
         color = colorRampPalette(c("white","black"))(20),
         cluster_rows = F,cluster_cols = F,
         clustering_method = "average")

##alpha diversity
library(vegan)
chao=estimateR(t(data[,1:24]))[2,]%>%as.data.frame()
shannon=diversity(t(abu[,1:24]))%>%as.data.frame()
colnames(shannon)="shannon_index"
shannon$site=group[rownames(shannon),2]
shannon_aver=rowsum(shannon$shannon_index,shannon$site)/3
shannon_aver=as.data.frame(shannon_aver)
colnames(shannon_aver)="average"
shannon_aver$SD=tapply(shannon$shannon_index, shannon$site, sd)
group1=read.csv("../../huge_phage/group1.txt",sep = '\t',header = F,row.names = 1)
shannon_aver$year=group1[rownames(shannon_aver),1]
shannon_aver$site=rownames(shannon_aver)
library(ggplot2)
ggplot(shannon_aver,aes(x = year,y=average))+
  geom_point(stat = 'identity',size=4,shape=21,fill="lightblue")+
  geom_errorbar(aes(ymin=average-SD,ymax=average+SD))+
  theme_classic()+geom_line()+
  geom_text(aes(label=site),size=4,vjust=1.7)+
  labs(y="shannon index")
shannon$year=group[rownames(shannon),3]
model <- loess(average ~  year, data = shannon_aver)
summary(model)
sum_abu$year=group[rownames(sum_abu),3]
ggplot(data = shannon_aver,mapping = aes(x = year,y=average))+
  theme_classic()+
  labs(y="shannon index")+
  geom_smooth(method = 'loess', formula = y~x, span = 1,color="red")+
  geom_point(stat = 'identity',size=4,shape=21,fill="lightblue")+
  geom_errorbar(aes(ymin=average-SD,ymax=average+SD))+
  theme_classic()+geom_line()+
  geom_text(aes(label=site),size=4,vjust=1.7)+
  labs(y="shannon index")
ggplot()+
  theme_classic()+
  labs(y="Accumulative abundance")+
  geom_smooth(data =sum_abu,mapping = aes(x = year,y=sum_abu),
              method = 'loess', formula = y~x, span = 1,color="red")

####kraken
kraken=read.csv("all.final_report.profile",header = 1,row.names = 1,check.names = F,sep='\t')
kraken_order=kraken[grep("o_",rownames(kraken)),]
kraken_order=kraken_order[-grep("f_",rownames(kraken_order)),]
kraken_order=kraken_order[-grep("g_",rownames(kraken_order)),]
kraken_order=kraken_order[-grep("s_",rownames(kraken_order)),]
kraken_order=relative_abundance(kraken_order)
kraken_SOB=rbind(
  kraken_order[grep(pattern = "Acidiferrobacterales",rownames(kraken_order)),]
)
kraken_SOB=kraken_SOB*100
library(ggplot2)
library(reshape2)
library(ggsci)
kraken_SOB=melt(kraken_SOB)
ggplot(kraken_SOB,aes(Var2,y = value))+
  geom_bar(stat = 'identity',fill="#1F77B4")+
  theme_classic()
kraken_phylum=kraken[grep("p_",rownames(kraken)),]
kraken_phylum=kraken_phylum[-grep("o_",rownames(kraken_phylum)),]
kraken_phylum=kraken_phylum[-grep("c_",rownames(kraken_phylum)),]
kraken_phylum=kraken_phylum[-grep("g_",rownames(kraken_phylum)),]
kraken_phylum=kraken_phylum[-grep("s_",rownames(kraken_phylum)),]

kraken_phylum=relative_abundance(kraken_phylum)
sum=apply(kraken_phylum,1,sum)
kraken_phylum=kraken_phylum[order(sum,decreasing = T),]
kraken_phylum_top10=kraken_phylum[1:10,]
kraken_phylum_top10=melt(kraken_phylum_top10)
ggplot(kraken_phylum_top10,aes(Var2,value*100,fill=Var1))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_d3()+theme_bw()+
  theme(axis.text.x.bottom = element_text(angle = 60,hjust = 1))+
  labs(y="Relative abundance (%)")
  
  
