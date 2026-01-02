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

summary(lm(formula = value~vir_dist,data = bac_dist))

pcoa_v <- cmdscale(distance, k = (nrow(results5) - 1), eig = TRUE)
point_v <- data.frame(pcoa_v$point)
pcoa_eig_v <- (pcoa_v$eig)[1:2] / sum(pcoa_v$eig)
sample_site_v <- data.frame ({pcoa_v$point})[1:3]
sample_site_v$names <- rownames(sample_site_v)
names(sample_site_v)[1:3] <- c('PCoA1', 'PCoA2','PCoA3')
rownames(group)=group$sample
sample_site_v$site = group[rownames(sample_site_v),2]
sample_site_v$stage=rep(c("early","late"),c(9,15))
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

##Candidatus_Desulfobacterota
tmp=strsplit(bacteria_taxonomy$classification,split = ";")
tmp=do.call(rbind,tmp)
bacteria_taxonomy$family=tmp[,5]
bacteria_taxonomy$order=tmp[,4]
abu$class=bacteria_taxonomy[rownames(abu),22]
abu$family=bacteria_taxonomy[rownames(abu),23]
abu$order=bacteria_taxonomy[rownames(abu),24]

Desulfobacterota=abu[grep("Desulfobacterota",abu$phylum),]
Desulfobacterota=melt(as.matrix(Desulfobacterota[,c(1:24)]))
ggplot(Desulfobacterota,aes(x = Var2,y=value*100,fill=Var1))+
  geom_bar(stat = "identity",
           alpha=0.7,width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  scale_fill_d3(palette = "category20c")+
  labs(y="Relative abundance (%)")

##SOB
SOB=read.table("SOB.txt",sep='\t',header = 1,row.names = 1)
SOB=as.data.frame(SOB)
SOB[,2:29]=abu[rownames(SOB),1:28]
library(reshape2)
Sulfur_oxider=melt(SOB[,c(2:29)])
ggplot(Sulfur_oxider,aes(x = variable,y=value*100,fill=order))+
  geom_bar(stat = "identity",
           alpha=0.7,width = 0.8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  scale_fill_d3()+
  labs(y="Relative abundance (%)")
env=read.table("../../202505批次/sulfate.txt",header = 1,row.names = 1,sep='\t',check.names = F)
SOB_sum=data.frame(SOB=apply(SOB[,2:25],2,sum))
env$SOB_sum=SOB_sum[rownames(env),1]

summary(lm(formula = SOB_sum~Sulfite_2025,env))
summary(lm(formula = SOB_sum~Sulfate_2021,env))
summary(lm(formula = SOB_sum~Sulfate_202505,env))

write.table(shannon,"shannon_index.tsv",sep="\t",quote = F)
env$stage=rep(c("S1-S3","S4-S8"),c(9,15))
library(ggplot2)
ggplot(env)+
  geom_point(mapping = aes(SOB_sum,Sulfate_2021,fill = stage),shape=21,size=2)+
  geom_point(aes(SOB_sum,Sulfate_202505,fill = stage),shape=22,size=2)+
  geom_smooth(aes(SOB_sum,Sulfate_2021),method = "lm")+
  geom_smooth(aes(SOB_sum,Sulfate_202505),method = "lm")+
  labs(y="Concerntration of Sulfate (mg/kg)",x="Relative abundance of rDSRB")+
  theme_bw()
ggplot(env)+
  geom_point(mapping = aes(SOB_sum,Sulfite_2025,fill = stage),shape=21,size=2)+
  geom_smooth(aes(SOB_sum,Sulfite_2025),method = "lm")+
  labs(y="Concerntration of Sulfate (mg/kg)",x="Relative abundance of rDSRB")+
  theme_bw()
  

