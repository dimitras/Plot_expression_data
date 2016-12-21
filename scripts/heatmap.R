library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)

ipadata <- read.csv("~/Dropbox/lihong/IPA/WT-KOcomparison0.6.24.pathway.edited.for.heatmap.filtered.csv")
ipadata %>% 
  arrange(WT0hKO0h) %>% 
  # arrange(desc(WT0hKO0h)) %>% #order descending
  rename(`0h`=WT0hKO0h, `6h`=WT6hKO6h, `24h`=WT24hKO24h) %>%
  mutate(Name=factor(Name,levels=as.character(Name))) %>% 
  gather(TP,Pvalue,`0h`:`24h`) %>%
  mutate(TP=factor(TP,levels=c("0h","6h","24h"))) %>% 
  ggplot(aes(x = TP, y = Name, fill=Pvalue)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high="darkblue") +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("~/Dropbox/lihong/IPA/WT-KOcomparison0.6.24.pathway.edited.for.heatmap.filtered.tiff", height = 25, width = 15, units ="cm")
ggsave("~/Dropbox/lihong/IPA/WT-KOcomparison0.6.24.pathway.edited.for.heatmap.filtered.eps", height = 25, width = 15, units ="cm")




ipadata <- read.csv("~/Dropbox/lihong/IPA/WT-KOcomparison0.6.24.disease.edited.for.heatmap.filtered.csv")
ipadata %>% 
  arrange(WT0hKO0h) %>% 
  rename(`0h`=WT0hKO0h, `6h`=WT6hKO6h, `24h`=WT24hKO24h) %>%
  mutate(Name=factor(Name,levels=as.character(Name))) %>% 
  gather(TP,Pvalue,`0h`:`24h`) %>%
  mutate(TP=factor(TP,levels=c("0h","6h","24h"))) %>% 
  ggplot(aes(x = TP, y = Name, fill=Pvalue)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgreen", high="darkgreen") +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("~/Dropbox/lihong/IPA/WT-KOcomparison0.6.24.disease.edited.for.heatmap.filtered.tiff", height = 30, width = 10, units ="cm")
ggsave("~/Dropbox/lihong/IPA/WT-KOcomparison0.6.24.disease.edited.for.heatmap.filtered.eps", height = 30, width = 10, units ="cm")




genedata <- read.csv("~/Dropbox/lihong/results/WT0h-KO0h_vs_WT6h-KO6h_vs_WT24h-KO24h_common_genes.xlsx.csv")
gene_data.for_hclust = genedata %>% select(id, starts_with("WT")) %>% column_to_rownames("id") %>% as.matrix
# gene_data.for_hclust = genedata %>% column_to_rownames("id") %>% as.matrix
genedata.hclust_results = hclust(dist(gene_data.for_hclust))
# genedata.hclust_results = hclust(as.dist((1 - cor(t(gene_data.for_hclust)))))

fc.genedata = genedata %>%
  gather(Sampleid,Count,WT1_0h:KO4_24h) %>%
  mutate(TP=gsub(".*([[:digit:]]_[[:digit:]]+h)","\\1",Sampleid),
         TP=factor(TP,levels=c("1_0h", "2_0h","3_0h","4_0h" ,"1_6h", "2_6h","3_6h","4_6h","1_24h", "2_24h","3_24h","4_24h")),
         gt=gsub("^(WT|KO).*","\\1",Sampleid),
         gt=factor(gt,levels=c("WT","KO")),
         logcount = log10(Count),
         id = factor(id, levels=genedata.hclust_results$labels[genedata.hclust_results$order])
  ) %>%
  ggplot(aes(x = TP, y = id, fill=logcount)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high="darkblue") +
  facet_grid(.~gt) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        strip.background = element_blank())
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_abundance.tiff", height = 40, width = 15, units ="cm")
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_abundance.eps", height = 30, width = 15, units ="cm")






# heatmap with zscore 
genedata <- read.csv("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes.csv")
gene_data.for_hclust = genedata %>% column_to_rownames("gene_symbol") %>% select(WT1_0h:KO4_24h) %>% as.matrix
genedata.hclust_results = hclust(as.dist((1 - cor(t(gene_data.for_hclust)))))
# genedata.hclust_results = hclust(dist(gene_data.for_hclust))
# gene_data.for_hclust = genedata %>% select(id, starts_with("WT")) %>% column_to_rownames("id") %>% as.matrix

genedata %>%
  select(id, gene_symbol, WT1_0h:KO4_24h) %>% 
  gather(Sampleid,Count,WT1_0h:KO4_24h) %>%
  mutate(TP=gsub(".*([[:digit:]]_[[:digit:]]+h)","\\1",Sampleid), 
         TP=factor(TP,levels=c("1_0h", "2_0h","3_0h","4_0h" ,"1_6h", "2_6h","3_6h","4_6h","1_24h", "2_24h","3_24h","4_24h")),
         gt=gsub("^(WT|KO).*","\\1",Sampleid),
         gt=factor(gt,levels=c("WT","KO")),
         gene_symbol = factor(gene_symbol, levels=genedata.hclust_results$labels[genedata.hclust_results$order])
  ) %>%

  group_by(gene_symbol) %>%
  mutate(rmean = mean(Count, na.rm=TRUE),
         rSD = sd(Count, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(zscore = (Count - rmean) / rSD) %>% 
         # zscore = ifelse(zscore > 2, 2, ifelse(zscore < -2, -2, zscore))) %>%
    
  ggplot(aes(x = TP, y = gene_symbol , fill=zscore)) +
  geom_tile() +
  # scale_fill_gradient(low = "lightblue", high="darkblue") +
  scale_fill_gradientn(colours = c("red", "orange", "yellow", "green", "blue"), values = rescale(c(3,2,1,0,-1))) +
  facet_grid(.~gt) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 15),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.grid = element_blank(),
        strip.background = element_blank())
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_zscore_corhclust_ncolor.tiff", height = 55, width = 25, units ="cm")
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_zscore_corhclust_ncolor.eps", height = 55, width = 25, units ="cm")







genedata <- read.csv("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_full.csv")

# gene_data.for_hclust = genedata %>% select(id, starts_with("WT")) %>% column_to_rownames("id") %>% as.matrix
# genedata.hclust_results = hclust(dist(gene_data.for_hclust))

fc.genedata = genedata %>%
  gather(Sampleid,Count,WT1_0h:KO4_24h) %>%
  mutate(TP=gsub(".*([[:digit:]]_[[:digit:]]+h)","\\1",Sampleid), 
         TP=factor(TP,levels=c("1_0h", "2_0h","3_0h","4_0h" ,"1_6h", "2_6h","3_6h","4_6h","1_24h", "2_24h","3_24h","4_24h")),
         TP.norep=gsub(".*_([[:digit:]]+h)","\\1",Sampleid), 
         TP.norep=factor(TP.norep,levels=c("0h", "6h", "24h")),
         gt=gsub("^(WT|KO).*","\\1",Sampleid),
         gt=factor(gt,levels=c("WT","KO")),
         id = factor(id, levels=genedata.hclust_results$labels[genedata.hclust_results$order])
         ) %>%
  group_by(id, gt, TP.norep) %>%
  summarize(repmean = mean(Count, na.rm=TRUE)) %>%
  ungroup() %>% spread(gt, repmean) %>% 
  mutate(log2fc = log2((KO+0.01)/(WT+0.01)),
         log2fc = ifelse(log2fc>4, 4, log2fc)
         )
  # # cluster by logFC
  gene_data.for_hclust = fc.genedata %>% select(-WT, -KO) %>% spread(TP.norep, log2fc) %>% data.frame %>% rownames_to_column("X") %>%  column_to_rownames("id") %>% select(-X) %>% as.matrix
  # genedata.hclust_results = hclust(as.dist((1 - cor(t(gene_data.for_hclust)))))
  genedata.hclust_results = hclust(dist(gene_data.for_hclust))

fc.genedata  %>% mutate(id = factor(id, levels=genedata.hclust_results$labels[genedata.hclust_results$order]))  %>% ggplot(aes(x = TP.norep, y = id, fill=log2fc)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high="darkblue", name = "log2FC\nKO/WT", labels = c("2", "3", ">4"), breaks = c(2,3,4)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        strip.background = element_blank())
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_logFC.pdf")
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_logFC.tiff", height = 40, width = 15, units ="cm")
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_logFC.eps", height = 30, width = 15, units ="cm")





# Heatmap for the union
uniondata <- read.csv("~/Dropbox/lihong/results/uniongenes_fullexpression.csv") #union3DEcomparisons_fc1.5pv0.1_union.csv
# uniondata.for_hclust = uniondata %>% column_to_rownames("id") %>% select(WT1_0h.C1:KO4_0h.C1, KO1_6h.C2:WT4_6h.C2, WT1_24h.C3:KO4_24h.C3) %>%  as.matrix
uniondata.for_hclust = uniondata %>% column_to_rownames("id") %>% select(WT1_0h:KO4_24h) %>%  as.matrix
# uniondata.hclust_results = hclust(dist(uniondata.for_hclust))
uniondata.hclust_results = hclust(as.dist((1 - cor(t(uniondata.for_hclust)))))

uniondata %>%
  # select(id, WT1_0h.C1:KO4_0h.C1, KO1_6h.C2:WT4_6h.C2, WT1_24h.C3:KO4_24h.C3) %>% 
  # gather(Sampleid,Count,WT1_0h.C1:KO4_24h.C3) %>%
  # mutate(TP=gsub(".*([[:digit:]]_[[:digit:]]+h).C[[:digit:]]","\\1",Sampleid), 
  select(id, WT1_0h:KO4_24h) %>% 
  gather(Sampleid,Count,WT1_0h:KO4_24h) %>%
  mutate(TP=gsub(".*([[:digit:]]_[[:digit:]]+h)","\\1",Sampleid), 
         TP=factor(TP,levels=c("1_0h", "2_0h","3_0h","4_0h" ,"1_6h", "2_6h","3_6h","4_6h","1_24h", "2_24h","3_24h","4_24h")),
         gt=gsub("^(WT|KO).*","\\1",Sampleid),
         gt=factor(gt,levels=c("WT","KO")),
         id = factor(id, levels=uniondata.hclust_results$labels[uniondata.hclust_results$order])
  ) %>%
  
  group_by(id) %>%
  mutate(rmean = mean(Count, na.rm=TRUE),
         rSD = sd(Count, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(zscore = (Count - rmean) / rSD) %>% 
  # zscore = ifelse(zscore > 2, 2, ifelse(zscore < -2, -2, zscore))) %>%
  
  ggplot(aes(x = TP, y = id, fill=zscore)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred" ,"red", "orange", "yellow", "green", "blue", "purple"), values = rescale(c(4,3,2,1,0,-1,-2))) +
  # scale_fill_manual(breaks=c("\[Inf,3)", "\[3,2)", "\[2,1)", "\[1,0)", "\[0,-1)", "\[-1,-Inf)"), values = c("red", "orange", "yellow", "white", "green", "blue")) + #ncolor
  # scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0) + #3color
  # scale_fill_gradient(low = "lightblue", high = "darkblue") + #2color
  facet_grid(.~gt) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 15),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.grid = element_blank(),
        strip.background = element_blank())
# ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.union_zscore.tiff", height = 40, width = 15, units ="cm")
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.union_fc1.25pv0.1_zscore_corhclust_ncolor.tiff", height = 45, width = 20, units ="cm")
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.union_fc1.25pv0.1_zscore_corhclust_ncolor.eps", height = 45, width = 20, units ="cm")








# FC in the input file ????????

genedata <- read.csv("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes.csv")
row.names(genedata) = genedata$id
gene_data.for_hclust = genedata %>% as.matrix
genedata.hclust_results = hclust(dist(gene_data.for_hclust))

fc.genedata = genedata %>%
  gather(FCperTP,logFC,logFC_0h:logFC_24h) %>%
  mutate(log2fc = logFC,
         log2fc = ifelse(logFC>4, 4, log2fc),
         id = factor(id, levels=genedata.hclust_results$labels[genedata.hclust_results$order])
  )
# # cluster by logFC
gene_data.for_hclust = fc.genedata %>% spread(FCperTP, log2fc) %>% data.frame %>% rownames_to_column("X") %>%  column_to_rownames("id") %>% select(-X) %>% as.matrix
genedata.hclust_results = hclust(dist(gene_data.for_hclust))

fc.genedata  %>% mutate(id = factor(id, levels=genedata.hclust_results$labels[genedata.hclust_results$order]))  %>% ggplot(aes(x = FCperTP, y = id, fill=log2fc)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high="darkblue", name = "log2FC\nKO/WT", labels = c("2", "3", ">4"), breaks = c(2,3,4)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        strip.background = element_blank())
ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_log2FC.tiff", height = 40, width = 10, units ="cm")

ggsave("~/Dropbox/lihong/results/WT-KOcomparison0.6.24.common_genes_logFC.eps", height = 30, width = 15, units ="cm")
