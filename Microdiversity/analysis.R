scaffold_nucl=read.csv("all.scaffold.nucl_diversity.tsv",sep='\t',
                       header = 1,row.names = 1,check.names = F)

nucl_freq=data.frame(freq=apply(scaffold_nucl,1,function(x){sum(x>0)}))
for (i in 1:nrow(scaffold_nucl)){
  for (j in 1:ncol(scaffold_nucl)){
    if ( scaffold_nucl[i,j]==0 ){
      scaffold_nucl[i,j]=NA
    }
  }
}
lifestyle=read.table("../GTVD_phatyp_prediction.tsv",sep='\t',header = 1,row.names = 1)
lifestyle[lifestyle$PhaTYPScore==0,2]="unclassified"
scaffold_nucl$lifestyle=lifestyle[rownames(scaffold_nucl),2]

average_scal_nucl=data.frame(total_average=NA)
for (i in 1:ncol(scaffold_nucl)){
  average_scal_nucl[i,1]=median(scaffold_nucl[,i],na.rm = T)
}

rownames(average_scal_nucl)=colnames(scaffold_nucl)

length=read.table("GTVD_quality_summary.filter10k_95-80_1500.length",
                  row.names = 1)
quality=read.csv("../quality_summary.tsv",sep='\t',header = 1,row.names = 1)
quality=as.data.frame(quality)
scaffold_nucl$completeness=quality[as.character(rownames(scaffold_nucl)),9]
scaffold_nucl$completeness=as.numeric(scaffold_nucl$completeness)
scaffold_nucl$quality=quality[as.character(rownames(scaffold_nucl)),7]
for(i in 1:nrow(scaffold_nucl)){
  if ( is.na(scaffold_nucl$quality[i]) == TRUE ){
    scaffold_nucl$quality[i]="proviruses"
  }
}
scaffold_nucl_com=scaffold_nucl[grep(pattern = "Complete",scaffold_nucl$quality),]
scaffold_nucl_high=scaffold_nucl[grep(pattern = "High",scaffold_nucl$quality),]
scaffold_nucl_pro=scaffold_nucl[grep(pattern = "proviruses",scaffold_nucl$quality),]
scaffold_nucl_length=rbind(scaffold_nucl_com,scaffold_nucl_high)
scaffold_nucl_length$length=length[as.character(rownames(scaffold_nucl_length)),1]
scaffold_nucl_length$length_cate=NA
for(i in 1:nrow(scaffold_nucl_length)){
  if(scaffold_nucl_length$length[i]<20000 && scaffold_nucl_length$length[i]>1500){
    scaffold_nucl_length$length_cate[i]="1.5-20kbp"
  }
  if(scaffold_nucl_length$length[i]<60000 && scaffold_nucl_length$length[i]>20000){
    scaffold_nucl_length$length_cate[i]="20-60kbp"
  }
  if(scaffold_nucl_length$length[i]>60000 && scaffold_nucl_length$length[i]<120000){
    scaffold_nucl_length$length_cate[i]="60-120kbp"
  }
  if(scaffold_nucl_length$length[i]>120000 && scaffold_nucl_length$length[i]<200000){
    scaffold_nucl_length$length_cate[i]="120-200kbp"
  }
  if(scaffold_nucl_length$length[i]>200000){
    scaffold_nucl_length$length_cate[i]="jumbophage"
  }
}
scaffold_nucl_length$length_cate=as.factor(scaffold_nucl_length$length_cate)

scaffold_nucl_length[,30:32]=quality[rownames(scaffold_nucl_length),4:6]
scaffold_nucl_length$host_percent=as.numeric(scaffold_nucl_length$host_genes)/as.numeric(scaffold_nucl_length$gene_count)
scaffold_nucl_length=scaffold_nucl_length[scaffold_nucl_length$host_percent<0.1,]
average_nucl=as.data.frame(matrix(data = NA,nrow = 5,ncol = 24))
phage_type=c("1.5-20kbp","20-60kbp","60-120kbp","120-200kbp","jumbophage")
rownames(average_nucl)=phage_type
for(i in 1:24){
  print(i)
  tmp=scaffold_nucl_length[,c(i,28)]
  tmp=na.omit(tmp)
  for (j in 1:nrow(average_nucl)){
    a=tmp[grep(rownames(average_nucl)[j],tmp$length_cate),]
    a_1=matrix(NA,1000,1)
    for (l in 1:1000){
      if(nrow(a)>20){
        a_1[l,]=mean(a[sample(1:nrow(a),size=20,replace=F),1])
      }else{
        a_1[l,]=mean(a[,1])
      }
    }
   average_nucl[j,i]=mean(a_1[,1])
  }
}
table(scaffold_nucl_length$length_cate)
colnames(average_nucl)=colnames(scaffold_nucl_length)[1:24]

average_pro_nucl=as.data.frame(matrix(data = NA,nrow = 1,ncol = 24))
phage_type=c("prophage")
rownames(average_pro_nucl)=phage_type
for(i in 1:24){
  print(i)
  tmp=scaffold_nucl_pro[,c(i,26)]
  tmp=na.omit(tmp)
  a_1=matrix(NA,1000,1)
  for (l in 1:1000){
    if(nrow(tmp)>50){
      a_1[l,1]=mean(tmp[sample(1:nrow(tmp),size=50,replace=F),1])
    }else{
      a_1[l,1]=mean(tmp[,1])
    }
  }
  average_pro_nucl[,i]=mean(a_1[,1])
}
colnames(average_pro_nucl)=colnames(scaffold_nucl_length)[1:24]
average_nucl=rbind(average_nucl,average_pro_nucl)
average_nucl=as.matrix(average_nucl)

##lifestyle
average_nucl_lifestyle=as.data.frame(matrix(data = NA,nrow = 3,ncol = 24))
phage_type=c("virulent","temperate","unclassified")
rownames(average_nucl_lifestyle)=phage_type
for(i in 1:24){
  print(i)
  tmp=scaffold_nucl[,c(i,25)]
  #tmp=scaffold_nucl_length[grep("20-60kbp",scaffold_nucl_length$length_cate),c(i,25)]
  tmp=na.omit(tmp)
  for (j in 1:nrow(average_nucl_lifestyle)){
    a=tmp[grep(rownames(average_nucl_lifestyle)[j],tmp$lifestyle),]
    a_1=matrix(NA,1000,1)
    for (l in 1:1000){
      if(nrow(a)>20){
        a_1[l,]=mean(a[sample(1:nrow(a),size=20,replace=F),1])
      }else{
        a_1[l,]=mean(a[,1])
      }
    }
    average_nucl_lifestyle[j,i]=mean(a_1[,1])
  }
}

colnames(average_nucl_lifestyle)=colnames(scaffold_nucl_length)[1:24]
library(reshape2)
average_nucl=melt(as.matrix(average_nucl_lifestyle))
average_nucl=na.omit(average_nucl)
metadata$stage2=rep(c("early","late"),c(12,12))
average_nucl$stage2=metadata[average_nucl$Var2,4]
library(ggpubr)
library(ggplot2)
library(ggsci)
ggplot(average_nucl,aes(stage2,value,
                        fill=stage2))+
  geom_boxplot(color="#087984",outlier.size = 0.7)+
  theme_bw()+
  facet_wrap(~Var1)+
  scale_fill_d3(palette = "category20")+
  theme(axis.text.x = element_text(angle = 30,hjust = 0.75,face = "bold"))+
  labs(y="Microdiversity",x="viral types")+
  stat_compare_means(method = "wilcox.test",comparisons = 
                       list(c("early","late")
                       ))
  

  
ggplot(average_nucl[grep("early",average_nucl$stage2),],
       aes(stage2,value,fill=stage2))+
  geom_boxplot(color="#087984",outlier.size = 0.7)+
  theme_bw()+
  scale_fill_d3(palette = "category20")+
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 0.75,face = "bold"))+
  stat_compare_means(method = "wilcox.test",comparisons = 
                     list(c("virulent","temperate")
                     ))+
  labs(y="log10(Microdiversity)",x="viral types")

ggplot(average_nucl[grep("pro",average_nucl$Var1),],
       aes(stage,log10(value),fill=stage))+
  geom_boxplot(color="#2BEF0F",outlier.size = 0.7)+
  theme_bw()+
  scale_fill_d3(palette = "category20")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,face = "bold"))+
  stat_compare_means(method = "wilcox.test",comparisons = 
                       list(c("putative_lytic","lysogenic")
                       ))+
  labs(y="log10(Microdiversity)",x="viral types")



ggplot(average_nucl[grep("late",average_nucl$stage2),],aes(Var1,log10(value),fill=stage2))+
  geom_boxplot(color="#087984",outlier.size = 0.7)+
  theme_bw()+
  scale_fill_d3(palette = "category20")+
  theme(axis.text.x = element_text(angle = 30,hjust = 0.75,face = "bold"))+
stat_compare_means(method = "wilcox.test",comparisons = 
                     list(c("putative_lytic","lysogenic")))+
  labs(y="log10(Microdiversity)",x="viral types")

library(rstatix)
average_nucl$log10=log10(average_nucl$value)
stat.test <- average_nucl %>%
  wilcox_test(
     log10~ Var1,alternative = "greater",
    p.adjust.method = "bonferroni"
  )

metadata=read.table("../../huge_phage/metadata.txt",sep='\t',header = 1,row.names = 1)

###bacterial
bacteria=read.table("all_bacterial_genome.nucl_diversity.tsv",
                    sep = '\t',header = 1,row.names = 1,
                    check.names = F)
ave_bacteria=as.data.frame(matrix(data = NA,nrow = 1,ncol = 24))
rownames(ave_bacteria)="bacteria"
for(i in 1:24){
  print(i)
  tmp=data.frame(value=bacteria[,c(i)])
  tmp=na.omit(tmp)
  a_1=matrix(NA,1000,1)
  for (l in 1:1000){
    if(nrow(tmp)>50){
      a_1[l,1]=mean(tmp[sample(1:nrow(tmp),size=50,replace=F),1])
    }else{
      a_1[l,1]=mean(tmp[,1])
    }
  }
  ave_bacteria[,i]=mean(a_1[,1])
}
colnames(ave_bacteria)=colnames(bacteria)

ave_viruses=as.data.frame(matrix(data = NA,nrow = 1,ncol = 24))
rownames(ave_viruses)="viruses"
for(i in 1:24){
  print(i)
  tmp=data.frame(value=scaffold_nucl[,c(i)])
  tmp=na.omit(tmp)
  a_1=matrix(NA,1000,1)
  for (l in 1:1000){
    if(nrow(tmp)>50){
      a_1[l,1]=mean(tmp[sample(1:nrow(tmp),size=50,replace=F),1])
    }else{
      a_1[l,1]=mean(tmp[,1])
    }
  }
  ave_viruses[,i]=mean(a_1[,1])
}
colnames(ave_viruses)=colnames(scaffold_nucl)[1:24]

ave_bacteria=rbind(ave_bacteria,ave_viruses)

ave_bacteria=as.matrix(ave_bacteria)
ave_bacteria=melt(ave_bacteria)
ave_bacteria$stage=metadata[ave_bacteria$Var2,3]
library(ggsci)
ggplot(ave_bacteria,aes(x=value,group=Var1,fill=Var1))+
  #geom_boxplot(color="#2BEF0F",outlier.size = 0.7)+
  geom_density(alpha=0.6)+
  theme_bw()+
  scale_fill_d3(palette = "category20b")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,face = "bold"))+
  #stat_compare_means(method = "wilcox.test",comparisons = 
   #                    list(c("bacteria","viruses")
    #                   ))+
  labs(y="Density of microdiversity",x="Microdiversity")
median(ave_bacteria[grep("bacteria",ave_bacteria$Var1),3])
sd(ave_bacteria[grep("bacteria",ave_bacteria$Var1),3])
median(ave_bacteria[grep("viruses",ave_bacteria$Var1),3])
sd(ave_bacteria[grep("viruses",ave_bacteria$Var1),3])

bacteria_taxonomy=read.table("../bacteria/prokaryotic_species_clusters_taxonomy.tsv",sep="\t")
class=strsplit(bacteria_taxonomy$classification,split = ";")
class=do.call(rbind,class)
bacteria_taxonomy$class=class[,3]
bacteria_taxonomy$phylum=class[,2]
bacteria$class=bacteria_taxonomy[rownames(bacteria),21]
bacteria$phylum=bacteria_taxonomy[rownames(bacteria),24]
table_class=as.data.frame(table(bacteria$class))
table_class=table_class[order(table_class$Freq,decreasing = T),]
bacteria$class=factor(bacteria$class,levels = as.factor(table_class$Var1))
library(reshape2)
bacteria_class=melt(bacteria)
bacteria_class$stage=metadata[bacteria_class$variable,4]
bacteria_class=na.omit(bacteria_class)
bacteria_class_top=matrix(data = NA,nrow = 0,ncol = ncol(bacteria_class))
for (i in 1:6) {
  patt=table_class[i,1]
  tmp=bacteria_class[grep(patt,bacteria_class$class),]
  bacteria_class_top=rbind(bacteria_class_top,tmp)
}
library(ggplot2)
library(ggsci)
library(ggpubr)
ggplot(bacteria_class_top,aes(x = class,y = value,fill=stage))+
  geom_boxplot(width=0.6,outlier.size = 0.7,outlier.shape = 21)+
  coord_flip()+
  theme_classic()+
  scale_fill_lancet()+
  theme(axis.text.x.bottom = element_text(angle = 40,hjust = 1))+
  stat_compare_means(inherit.aes = T,method = "wilcox.test",
                     comparisons = list(c("c__Gammaproteobacteria","c__Alphaproteobacteria"),
                                        c("c__Bacteroidia","c__Gemmatimonadetes"),
                                        c("c__Gemmatimonadetes","c__Acidobacteriae"),
                                        c("c__Acidobacteriae","c__Actinomycetia"),
                                        c("c__Gammaproteobacteria","c__Actinomycetia")))+
  labs(y="Microdiversity")

Gamma=bacteria[grep("c__Gammaproteobacteria",bacteria$class),]
order=strsplit(bacteria_taxonomy$classification,split = ";")
order=do.call(rbind,order)
bacteria_taxonomy$order=order[,4]
bacteria_taxonomy$family=order[,5]
Gamma$order=bacteria_taxonomy[rownames(Gamma),22]
Gamma$family=bacteria_taxonomy[rownames(Gamma),23]
Gamma=melt(Gamma)
Gamma$stage=metadata[Gamma$variable,4]
Gamma=na.omit(Gamma)
Gamma$sites=metadata[Gamma$variable,1]
ggplot(Gamma[grep("Acid",Gamma$order),],aes(x = sites,y = value,fill=sites))+
  geom_boxplot(width=0.6,outlier.size = 0.7,outlier.shape = 21)+
  theme_classic()+
  scale_fill_lancet()+
  theme(axis.text.x.bottom = element_text(angle = 40,hjust = 1))+
  labs(y="Microdiversity")+
  stat_compare_means(comparisons = list(c("S1","S2")))

Roku=bacteria[grep("Methylom",bacteria$class),]
Roku=melt(Roku)
Roku=na.omit(Roku)
Roku$site=metadata[Roku$variable,1]
ggplot(Roku,aes(x = site,y = value,fill=site))+
  geom_boxplot(width=0.6,outlier.size = 0.7,outlier.shape = 21)+
  theme_classic()+
  scale_fill_d3()+
  theme(axis.text.x.bottom = element_text(angle = 40,hjust = 1))+
  labs(y="Microdiversity")+
  stat_compare_means(comparisons = list(c("S4","S6")))
###quality filtering
length=read.table("GTVD_quality_summary.filter10k_95-80_1500.length",
                  row.names = 1)
length[,2:14]=quality[rownames(length),]

length$host_percent=length$host_genes/length$gene_count
length2=as.data.frame(matrix(NA,nrow = 0,ncol = ncol(length)))

for(i in 1:nrow(length)){
  if(is.na(length$host_percent[i])==FALSE){
    if(length$host_percent[i]<0.1){
      length2=rbind(length2,length[i,])
    }
  }else{
    length2=rbind(length2,length[i,])
  }
}

write.table(length2[,1],file = "GTVD_contigs_filtering.tsv",row.names = T,quote = F,sep = '\t')
jumbophage_length=length2[length2$V2>200000,]
write.table(jumbophage_length,file = "GTVD_jumbophage.tsv",row.names = T,quote = F,sep = '\t')


###conANI
compare=read.table("compare_comparisonsTable.tsv",sep='\t',header = 1)
compare=compare[compare$percent_genome_compared>0.5,]
metadata$stage=rep(c("early","later"),c(9,15))

compare$site1=metadata[as.character(compare$name1),1]
compare$site2=metadata[as.character(compare$name2),1]

for (i in 1:nrow(compare)){
  if (compare$site1[i]==compare$site2[i]){
    compare$group1[i]="intra_site"
  }
  if(compare$site1[i]!=compare$site2[i]){
    compare$group1[i]="inter_site"
  }
}

library(ggplot2)
library(ggpubr)
ggplot(compare,aes(x = group1,y = conANI))+
  geom_boxplot(fill="lightblue",outlier.size = 0.5,
               outlier.colour = "blue",outlier.fill = "blue",outlier.alpha = 0.7)+
  theme_bw()+
  stat_compare_means(method = "wilcox.test",
                    comparisons = list(c("inter_site","intra_site")))
ggsave("compare_site_conANI.pdf",device = "pdf",width = 2.2,height = 2.7)

compare$stage1=metadata[as.character(compare$name1),3]
compare$stage2=metadata[as.character(compare$name2),3]
for (i in 1:nrow(compare)){
  if (compare$stage1[i]==compare$stage2[i]){
    if(compare$site1[i]==compare$site2[i]){
      compare$group2[i]="intra_site"
    }else{
      if(compare$stage2[i]=="later"){
        compare$group2[i]="intra_later"
      }
      if(compare$stage1[i]=='early'){
        compare$group2[i]='intra_early'
      }
    }

  }
  if(compare$stage1[i]!=compare$stage2[i]){
    compare$group2[i]="inter_stage"
  }
}
compare$year1=metadata[compare$name1,2]
compare$year2=metadata[compare$name2,2]
compare$year_dist=1+abs(compare$year1-compare$year2)
#compare$normal_conANI=compare$conANI/compare$year_dist
ggplot(compare,aes(x = group2,y =conANI))+
  geom_boxplot(fill="lightblue",outlier.size = 0.5,
               outlier.colour = "blue",outlier.fill = "blue",outlier.alpha = 0.7)+
  theme_bw()+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("intra_site","intra_later"),c("intra_early","intra_later"),
                                        c("inter_stage","intra_early")))
ggsave("compare_stage_conANI.pdf",device = "pdf",width = 3.2,height = 2.7)
compare$length=length2[compare$scaffold,1]
compare=na.omit(compare)
for (i in 1:nrow(compare)){
  if(compare$length[i]>100000){
    compare$phages[i]="huge_phage"
  }
  if(compare$length[i]<200000){
    compare$phages[i]="normal_phages"
  }
}

ggplot(compare[grep(pattern = "inter_site",compare$group1),],aes(phages,conANI))+
  geom_boxplot(fill="lightblue",outlier.size = 0.5,
               outlier.colour = "blue",outlier.fill = "blue",outlier.alpha = 0.7)+
  theme_bw()+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("huge_phage","normal_phages")))
ggsave("compare_phages_conANI.pdf",device = "pdf",width = 2.2,height = 2.7)
compare$group2=as.factor(compare$group2)
library(dplyr)
library(ggplot2)
library(ggsci)
lm(formula = conANI~year_dist,data = compare[grep("intra_later",compare$group2),])%>%summary()
ggplot(compare,aes(year_dist,conANI))+
  geom_point(aes(fill=group2,color=group2),
             shape=21,stat = 'identity',size=1,alpha=0.7)+
  #geom_boxplot(aes(group=year_dist,color=group2),outlier.colour = "white",outlier.alpha = 0,outlier.size = 0)+
  geom_smooth(aes(group=group2,color=group2),formula = y~x,method = "lm")+
  theme_bw()+
  #geom_smooth(data = compare,aes(year_dist,conANI,group=phages,color=phages),
   #           method = 'lm',formula = y~x,inherit.aes = F)+
  scale_fill_aaas()+
  scale_color_d3()
ggsave("compare_lm_dist_year_conANI.pdf",device ="pdf",width = 5.5,height = 3)
ggplot(compare[grep("huge",compare$phages),],aes(year_dist,conANI))+
  geom_point(aes(fill=group2,color=group2),
             shape=21,stat = 'identity',size=1,alpha=0.7)+
  #geom_boxplot(aes(group=year_dist,color=group2),outlier.colour = "white",outlier.alpha = 0,outlier.size = 0)+
  geom_smooth(aes(group=group2,color=group2),formula = y~x,method = "lm")+
  theme_bw()+
  #geom_smooth(data = compare,aes(year_dist,conANI,group=phages,color=phages),
  #           method = 'lm',formula = y~x,inherit.aes = F)+
  scale_fill_aaas()+
  scale_color_d3()

###gene
viruses_gene=read.table("all.viruses_genes_nucl1_5.tsv",
                        header = 1,row.names = 1,check.names = F)
bacterial_gene=read.table("all.bacterial_gene_info_filter17.tsv",sep='\t',header = 1,row.names = 1,check.names = F)
library(dplyr)
viruses_gene_1=na.omit(select(viruses_gene,'1-1'))
viruses_gene_1$`1-1`=as.numeric(viruses_gene_1$`1-1`)
library(reshape2)
viruses_gene=melt(viruses_gene)
viruses_gene$site=metadata[viruses_gene$variable,1]
viruses_gene$stage=metadata[viruses_gene$variable,3]
library(ggplot2)
library(ggsci)
ggplot(viruses_gene,aes(x=log10(value),fill=stage))+
  geom_density(lwd = 0.5,
               linetype = 2,
               alpha=0.6,
               adjust=1.75)+
  scale_fill_d3(palette = "category20")+
  theme_classic()
###bacterial_genes
bacterial_gene=read.table("all.MAGs_gene_nucl_diversity.tsv",sep='\t',header = 1,row.names = 1,check.names = F)
bacteria_function=read.table("../functions/SOB_contigs.txt",sep='\t',header = 1,row.names = 1)
bacterial_gene$KO=bacteria_function[rownames(bacterial_gene),1]
library(reshape2)
Dsr=melt(bacterial_gene)
Dsr=na.omit(Dsr)
Dsr=Dsr[-grep("K00471",Dsr$KO),]
Dsr=Dsr[-grep("K01011",Dsr$KO),]
Dsr=Dsr[-grep("K23077",Dsr$KO),]
library(ggplot2)
library(ggpubr)
ggplot(Dsr,aes(KO,value))+
  geom_boxplot(fill="lightblue",outlier.size = 0.8)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 40,hjust = 1))+
  labs(y="Microdiversity")+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("K11180","K08358"),
                                      c("K11181","K18500")))


#pNpS
bacterial_pNpS=read.table("all.MAGs_pNpS.tsv",sep='\t',header = 1,row.names = 1,check.names = F)
viral_pNpS=read.table("all_GFVD_genes.pNpS.tsv",sep='\t',header = 1,row.names = 1,check.names = F)

average_pNpS=as.data.frame(matrix(data = NA,nrow = 2,ncol = 24))
phage_type=c("bacteria","viruses")
rownames(average_pNpS)=phage_type
for(i in 1:24){
  print(i)
  tmp=as.data.frame(bacterial_pNpS[,i])
  tmp=na.omit(tmp)
  a_1=matrix(NA,1000,1)
  for (l in 1:1000){
    if(nrow(tmp)>100){
      a_1[l,1]=mean(tmp[sample(1:nrow(tmp),size=100,replace=F),1])
    }else{
      a_1[l,1]=mean(tmp[,1])
    }
  }
  average_pNpS[1,i]=mean(a_1[,1])
}
for(i in 1:24){
  print(i)
  tmp=as.data.frame(viral_pNpS[,i])
  tmp=na.omit(tmp)
  a_1=matrix(NA,1000,1)
  for (l in 1:1000){
    if(nrow(tmp)>100){
      a_1[l,1]=mean(tmp[sample(1:nrow(tmp),size=100,replace=F),1])
    }else{
      a_1[l,1]=mean(tmp[,1])
    }
  }
  average_pNpS[2,i]=mean(a_1[,1])
}
colnames(average_pNpS)=colnames(viral_pNpS)
library(reshape2)
average_pNpS_melt=melt(as.matrix(average_pNpS))
average_pNpS_melt$stage=metadata[average_pNpS_melt$Var2,4]

library(ggplot2)
library(ggsci)
library(ggpubr)
ggplot(average_pNpS_melt,aes(Var1,value))+
  geom_boxplot(aes(fill=Var1,color=Var1),alpha=0.4,outlier.shape=21)+
  facet_wrap(~stage,ncol=3) +
  labs(y="pNpS")+theme_bw()+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("bacteria","viruses")))

ggplot(average_pNpS_melt,aes(Var2,value))+
  geom_bar(aes(fill=Var1,color=Var1),stat = 'identity',alpha=0.4)+
  labs(y="pNpS")+theme_bw()
  
