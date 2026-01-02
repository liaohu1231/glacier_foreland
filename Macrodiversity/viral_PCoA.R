data=read.csv("all.GFVD.rpkm.tsv",header = 1,row.names = 1,check.names = F,sep='\t')
group=read.table("metadata.txt",header = 1,stringsAsFactors = T,sep='\t')

library(ggplot2)
library(ggsci)
relative_abundance=function(d){
  dta_sum=apply(d,2,function(x){x/sum(x)})
}
abu=relative_abundance(data)
abu=as.data.frame(abu)
quality=read.table("quality_summary.tsv",sep='\t',header = T,row.names = 1)
quality[grep("Yes",quality$provirus),7]="proviruses"
quality$checkv_quality=factor(quality$checkv_quality,levels = rev(c("Complete","High-quality",
                                                                "Medium-quality","proviruses",
                                                                "Low-quality","Not-determined")))
ggplot(quality,aes(checkv_quality,contig_length))+
  stat_boxplot(geom = "errorbar",
               width=0.3)+
  geom_boxplot(fill="#008B00",outlier.colour = "gray",
               outlier.shape = 21,outlier.fill = "white",outlier.size = 1)+
  theme_light()+ scale_y_log10()+ coord_flip()+labs(y="length (bp)",x="")+
  theme(axis.text = element_text(color = "black"),
                                                      axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
                                                      legend.key.size = unit(5, 'mm'),
                                                      legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'white', linetype = 2, size = 0), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'))
  
table(quality$checkv_quality)

###virsort
virsort=read.table("glacial-final-viral-score0.9.tsv",sep = '\t',header = 1,row.names = 1)
genetic=as.data.frame(table(virsort$max_score_group))
genetic$Var2=as.factor("total")
ggplot(genetic, aes(x=Var2, y=Freq,fill=Var1)) +
  geom_bar(stat="identity",width=1,color="black")+
  coord_polar("y",start=0)+
  theme_void()

library(ggpubr)
ggplot(virsort,aes(max_score_group,length))+
  stat_boxplot(geom = "errorbar",
               width=0.3)+
  geom_boxplot(fill="#008B00",outlier.colour = "gray",
               outlier.shape = 21,outlier.fill = "white",outlier.size = 1)+
  theme_light()+ scale_y_log10()+ 
  labs(y="length (bp)",x="")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'white', linetype = 2, size = 0), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'))+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ssDNA","dsDNAphage")))

###two tools
viral_number=data.frame(tools=as.factor(c("vibrant","virsort2","merged")),number=c(282829,354839,500296))
viral_number$tools=factor(viral_number$tools,levels = viral_number$tools)
ggplot(viral_number,aes(tools,number))+
  geom_bar(stat = "identity",width = 0.8,fill="#008B00")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'white', linetype = 2, size = 0), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'))

##PCoA
library(dplyr)
library(vegan)
results5=t(abu)
results5=as.matrix(results5)
library(vegan)
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

write.table(x = distance2,"viral_distance.tsv",quote = F,row.names = T,sep='\t')
library(pheatmap)
pheatmap(distance2,cluster_rows = F,cluster_cols = F,clustering_method = "average",
         na_col = "lightblue",border_color = F,cellwidth = 10,cellheight = 8,
         color = colorRampPalette(c("white","white","white","white",
                                    "white","white","white","white","white",
                                    "white","white","white","white","white","#40CCA7"))(30))

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
p_value=anosim(x= dist,grouping = sample_site_v$site,permutations = 9999)
summary(p_value)
p_value_1=adonis2(formula = dist~sites,group,permutations = 9999)
p_value_1
p_value_2=adonis2(formula = dist~year,group,permutations = 9999)
p_value_2
p_value_3=adonis2(formula = dist~stage,sample_site_v,permutations = 9999)
p_value_3
distance_early=vegdist(results5[1:12,], method= "bray",na.rm = T)
distance_later=vegdist(results5[13:24,], method= "bray",na.rm = T)
dist_env_early=vegdist(group$Sulfate_21_25[1:12])
dist_env_later=vegdist(group$Sulfate_21_25[13:24])
mantel(distance_early,dist_env_early)
mantel(distance_later,dist_env_later)
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

#lifestyle
lifestyle=read.table("GTVD_phatyp_prediction.tsv",sep='\t',header = 1,row.names = 1)
lifestyle[lifestyle$PhaTYPScore==0,2]="unclassified"
abu$lifestyle=lifestyle[rownames(abu),2]
abu$length=lifestyle[rownames(abu),1]
abu_lifestyle=aggregate.data.frame(abu[,1:24],list(abu$lifestyle),sum)
library(reshape2)
abu_lifestyle=melt(abu_lifestyle)
colnames(abu_lifestyle)[1]="group"
library(ggplot2)
library(ggsci)
ggplot(abu_lifestyle,aes(variable,value,fill=group))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_d3()+
  labs(y="Relative abundance",x="Samples")

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
group1=read.csv("../huge_phage/group1.txt",sep = '\t',header = F,row.names = 1)
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


###geographic range
data_rel=abu
for (i in 1:nrow(data_rel)){
  for(j in 1:ncol(data_rel)){
    if (data_rel[i,j]>0){
      data_rel[i,j]=1
    }
  }
}

data_rel=as.data.frame(t(data_rel))
rownames(group)=group$sample
group
data_rel$site=group[rownames(data_rel),2]
group$stage=rep(c("early","late"),c(12,12))
data_rel$stage=group[rownames(data_rel),4]

data_sites=as.data.frame(matrix(NA,nrow = 8,ncol = ncol(data_rel)-2))
colnames(data_sites)=rownames(abu)
for (i in 1:nrow(data)) {
  print(i)
  data_sites[,i]=tapply(data_rel[,i],INDEX = data_rel$site,sum)
  for(j in 1:nrow(data_sites)){
    if(data_sites[j,i]>0){
      data_sites[j,i]=1
    }
  }
}
tmp=as.data.frame(tapply(data_rel[,1],INDEX = data_rel$site,sum))
rownames(data_sites)=rownames(tmp)
data_stages=as.data.frame(matrix(NA,2,ncol = ncol(data_rel)-2))
colnames(data_stages)=rownames(abu)

for (i in 1:nrow(data)){
  print(i)
  data_stages[,i]=tapply(data_rel[,i],INDEX = data_rel$stage,sum)
  for(j in 1:nrow(data_stages)){
    if(data_stages[j,i]>0){
      data_stages[j,i]=1
    }
  }
}
tmp=as.data.frame(tapply(data_rel[,1],INDEX = data_rel$stage,sum))
rownames(data_stages)=rownames(tmp)

vOTUs_geo=data.frame(geogra_range=rep(c(NA),c(nrow(data))))
rownames(vOTUs_geo)=rownames(data)
for (i in 1:nrow(data)){
  if (sum(data_rel[,i])==1 && data_stages[1,i]==1){
    vOTUs_geo[i,]='sample_specific (early)'
  }
  else if(sum(data_rel[,i])==1 && data_stages[2,i]==1){
    vOTUs_geo[i,]='sample_specific (later)'
  }else{
    if(sum(data_sites[,i])==1 && data_stages[1,i]==1){
      vOTUs_geo[i,]='site_specific (early)'
    }
    else if(sum(data_sites[,i])==1 && data_stages[2,i]==1){
      vOTUs_geo[i,]='site_specific (later)'
    }else{
      if(sum(data_stages[,i])==2){
        vOTUs_geo[i,]="inter_stage"
      }else{
        if(data_stages[1,i]==1){
          vOTUs_geo[i,]='intra_early_stage'
        }
        if(data_stages[2,i]==1){
          vOTUs_geo[i,]='intra_later_stage'
        }
      }
    }
  }
}
write.table(vOTUs_geo,"GFVD_vOTUs_geographic_range_s12_12.tsv",sep='\t',quote = F,row.names = T)

###viral_host
hosts=read.table("host/viral_host_network.tsv",sep='\t',header = 1)
nrow(as.data.frame(table(hosts$host)))
Sulfuricaulis=hosts[grep("Sulfuricaulis",hosts$genus),]
Methylobacterium=hosts[grep("Methylobacterium",hosts$genus),]
Rokubacteriales=hosts[grep("Rokubacteriales",hosts$order),]
rownames(Rokubacteriales)=Rokubacteriales$viruses
Rokubacteriales[,9:32]=abu[rownames(Rokubacteriales),1:24]
abu_Rokubacteriales=data.frame(sum_abu=apply(Rokubacteriales[,9:32],2,sum))
abu_Rokubacteriales$samples=rownames(abu_Rokubacteriales)
abu_Rokubacteriales$stage=rep(c("S1-S3","S4-S8"),c(9,15))
library(ggplot2)
library(ggpubr)
ggplot(abu_Rokubacteriales,aes(stage,sum_abu*100))+
  geom_boxplot(fill="#FF7F0E")+theme_bw()+
  geom_point()+
  labs(y="Relative abundance (%)")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("S1-S3","S4-S8")))
ggsave(filename = "Rokubacteriales.pdf",device = 'pdf',width = 3,height = 4.5)

##SOB_viruses
SOB_v=read.table("host/SOB_viruses_host_uniq.tsv",sep='\t')
SOB_v=SOB_v[!duplicated(SOB_v$V1),]
SOB_v[,6:29]=abu[SOB_v$V1,]
SOB_v=SOB_v[,-c(2:4)]
SOB_v_sum=data.frame(sum=apply(SOB_v[,3:26],2,sum))
SOB_v_sum$sum=SOB_v_sum$sum*100
SOB_v_sum$site=group[rownames(SOB_v_sum),2]
SOB_v_ms=data.frame(mean=tapply(SOB_v_sum$sum,SOB_v_sum$site,mean),
                    sd=tapply(SOB_v_sum$sum, SOB_v_sum$site, sd))
SOB_v_ms$site=rownames(SOB_v_ms)
SOB_v_ms$group="group"
sum(abu$`1-1`)
library(reshape2)
library(ggplot2)
ggplot(SOB_v_ms,aes(site,mean))+
  geom_point(stat = 'identity',shape=21,size=5,fill="#FF9999")+
  geom_line(stat = 'identity',aes(group=group),size=1)+
  theme_classic()+
  labs(y="relative_abundance (%)\n of rDSRB-associated viruses")
 
SOB_v$lifestyle=lifestyle[SOB_v$V1,2]
SOB_lifestyle=aggregate.data.frame(SOB_v[,3:26],list(SOB_v$lifestyle),sum)
rownames(SOB_lifestyle)=SOB_lifestyle$Group.1
SOB_lifestyle=relative_abundance(SOB_lifestyle[,2:25])
SOB_lifestyle=t(SOB_lifestyle)
SOB_lifestyle=as.data.frame(SOB_lifestyle)
SOB_lifestyle$site=group[rownames(SOB_lifestyle),2]
SOB_temperate_mean=data.frame(mean=tapply(SOB_lifestyle$temperate,SOB_lifestyle$site,mean))
SOB_temperate_mean$group="group"
SOB_temperate_mean$site=rownames(SOB_temperate_mean)
library(ggsci)
ggplot(SOB_temperate_mean,aes(site,mean*100))+
  geom_point(stat = 'identity',shape=21,size=5,fill="#FF9999")+
  geom_line(stat = 'identity',aes(group=group),size=1)+
  theme_classic()+
  labs(y="relative_abundance (%) \nof temperate viruses")
mean(SOB_temperate_mean$mean[3:8])
