SOB_vt=read.table("all_SOB_viruses_gene_trans.mean2",sep='\t',header = T,check.names = F,row.names = 1)
sum=apply(SOB_vt, 1, sum)
SOB_vt=SOB_vt[sum>0,]
sum_col=data.frame(sum_sample=apply(SOB_vt,2,sum))
metadata=read.table("../metadata.txt",sep = '\t',header = 1,check.names = F)
rownames(metadata)=metadata$transcriptomes_2025_sample
sum_col$site=metadata[rownames(sum_col),2]
sum_col_mean=data.frame(mean=tapply(sum_col$sum_sample,sum_col$site,mean),
                        sd=tapply(sum_col$sum_sample, sum_col$site, sd))
sum_col_mean$site=rownames(sum_col_mean)
sum_col_mean$group="group"
library(ggplot2)
ggplot(sum_col_mean,aes(site,mean))+
  geom_point(stat = 'identity',shape=21,size=5,fill="#FF9999")+
  geom_line(stat = 'identity',aes(group=group),size=1)+
  theme_classic()+
  labs(y="Accumulative mean \ncoverage")
ggsave(filename = "all_SOB_viruses_gene_trans.pdf",device = "pdf",width = 6,height = 2.2)
SOB_v_ann=read.table("SOB_viruses_best_exec.tsv",sep = '\t',header = F,row.names = 2)
SOB_vt$anno=SOB_v_ann[rownames(SOB_vt),6]
SOB_vt_ann=na.omit(SOB_vt)
write.table(x = SOB_vt_ann,"Table_S9.tsv",sep='\t',quote = F,row.names = T)
