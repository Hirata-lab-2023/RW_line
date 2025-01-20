library(ggplot2)
library(openxlsx)
library(patchwork)
library(ggsignif)
library(dplyr)
library(ggbreak)
library(stringr)
library(UpSetR)
library(reshape2)

########################################################################
#SNPs_per
########################################################################

a = read.table("vcf/merge/merged.vcf", sep = "\t", header = F)
all = a
strain_name=list("strain" = c(rep("AB", 3), 
                              rep("TU", 3),
                              rep("TL", 3),
                              rep("WIK", 3),
                              rep("SAT", 3),
                              rep("NHGRI-1", 3),
                              rep("PET", 3),
                              rep("RW", 3),
                              rep("M-AB", 3),
                              rep("*AB", 3),
                              rep("IM", 3),
                              rep("India", 3)))
hetero = NULL
homo = NULL
for (i in 1:36) {
  #hetero
  b = a[grep("0/1",a[,i+9]),]
  hetero_count = c(hetero_count,nrow(b))
  c = (nrow(b)/1345101831)*100
  d = c(as.data.frame(strain_name$strain[i]),c)
  e = as.data.frame(d)
  colnames(e)[1] = "name"
  colnames(e)[2] = "hetero_genome"
  hetero = rbind(hetero,e)
  #homo
  b = a[grep("1/1",a[,i+9]),]
  homo_count = c(homo_count,nrow(b))
  c = (nrow(b)/1345101831)*100
  d = c(as.data.frame(strain_name$strain[i]),c)
  e = as.data.frame(d)
  colnames(e)[1] = "name"
  colnames(e)[2] = "homo_genome"
  homo = rbind(homo,e)
  rm(b)
  rm(c)
  rm(d)
  rm(e)
  cat(i)
}
strain = unique(strain_name$strain)

#hetero
hetero$name = factor(hetero$name,levels = strain)

hetero_mean  = group_by(hetero, name) %>% 
  summarise_all(funs(mean))

sd  = group_by(hetero, name) %>% 
  summarise_all(funs(sd))

se = sd[,2]/sqrt(3)

hetero_mean = cbind (hetero_mean, sd$hetero_genome, se)
hetero_mean = as.data.frame(hetero_mean)
colnames(hetero_mean) = c("name", "mean", "sd", "se")

hetero_mean$name = factor(hetero_mean$name,levels = strain)

hetero_g = ggplot(hetero_mean, aes(x = name, y = mean))+
  geom_bar(stat = "identity", fill ="white", color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),width = .1)+
  ggtitle("Hetero%SNPs")+
  xlab("Labstoks")+
  ylab("Hetero%SNPs")+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13))

hetero_g = hetero_g + 
  geom_jitter(data = hetero, aes(x = name, y = hetero_genome), size = 1) +
  theme_bw()+
  theme(legend.position = "none")

pdf("hetero_per.pdf", height = 5, width = 8)
hetero_g
dev.off()

#homo
homo$name = factor(homo$name,levels = strain)

homo_mean  = group_by(homo, name) %>% 
  summarise_all(funs(mean))

sd  = group_by(homo, name) %>% 
  summarise_all(funs(sd))

se = sd[,2]/sqrt(3)

homo_mean = cbind (homo_mean, sd$homo_genome, se)
homo_mean = as.data.frame(homo_mean)
colnames(homo_mean) = c("name", "mean", "sd", "se")

homo_mean$name = factor(homo_mean$name,levels = strain)

homo_g = ggplot(homo_mean, aes(x = name, y = mean))+
  geom_bar(stat = "identity", fill ="white", color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),width = .1)+
  ggtitle("Homo%SNPs")+
  xlab("Labstoks")+
  ylab("homo%SNPs")+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13))

homo_g = homo_g + 
  geom_jitter(data = homo, aes(x = name, y = homo_genome), size = 1) +
  theme_bw()+
  theme(legend.position = "none")

pdf("homo_per.pdf", height = 5, width = 8)
homo_g
dev.off()

############################################################################################
# Anotation
# wget https://ftp.ensembl.org/pub/release-109/gtf/danio_rerio/Danio_rerio.GRCz11.109.gtf.gz
############################################################################################

vcf.ano2genename = function(snpeff){
  tmp = strsplit(snpeff$V8,"\\|")
  tmp = unlist(lapply(tmp, FUN = function(x){
    return(x[4])
  }))
  return(unique(tmp))
}
vcf.ano2geneID = function(snpeff){
  tmp = strsplit(snpeff$V8,"\\|")
  tmp = unlist(lapply(tmp, FUN = function(x){
    return(x[5])
  }))
  return(unique(tmp))
}
vcf.transcriptID = function(snpeff){
  tmp = strsplit(snpeff$V8,"\\|")
  tmp = unlist(lapply(tmp, FUN = function(x){
    return(x[7])
  }))
  return(unique(tmp))
}

gtf <- read.table("Danio_rerio.GRCz11.109.gtf.gz", header = FALSE, sep = "\t")
gene_all = NULL
gene = gtf[gtf$V3=="exon",]
for (i in 1:length(gene$V1)) {
  gene_4 = NULL
  gene_1 = t(as.data.frame(unlist(strsplit(gene$V9[i], "; "))))
  gene_2 = t(as.data.frame(unlist(strsplit(gene_1[,6]," "))))
  gene_2 = gene_2[,2]
  gene_3 = t(as.data.frame(unlist(strsplit(gene_1[,3]," "))))
  gene_3 = gene_3[,2]
  gene_4 = cbind(gene_2, gene_3) 
  gene_all = rbind(gene_all,gene_4)
  cat(i,"\n")
}

gene_all = data.frame(gene_all)

library(dplyr)
gene_all_new = gene_all %>% distinct(gene_2,gene_3, .keep_all = T)
a = table(unlist(gene_all_new[,1]))
a = as.data.frame(a)
name = c("Gene_name", "Splicing_variant")
colnames(a) = name
Splicing_variant_name = a



strain_list = list(c("AB",
                     "TU",
                     "TL",
                     "WIK",
                     "SAT", 
                     "NHGRI",
                     "PET",
                     "RW",
                     "*AB",
                     "India",
                     "M-AB",
                     "IM",
                     "OK"))

#strain
i=1
ALL_homo_high = list()
for (i in 1:length(strain_list[[1]])){
  c=list()
  ano_high_path = sprintf("vcf/%s",strain_list[[1]][i]) 
  ano_high_list = list(list.files(ano_high_path, pattern="HIGH.vcf.gz", full.names=T))
  for (t in 1:length(ano_high_list[[1]])) {
    a = read.table((ano_high_list[[1]][t]), sep = "\t", header = F)
    b = a[grep("1/1",a$V10),]
    b = list(b)
    c[sprintf("%s_%d",as.character(strain_list[[1]][i]),t)] = b
  }
  ALL_homo_high[[as.character(strain_list[[1]][i])]] = c 
}

####
i=1
t=1
z=1
All_homo_high_list = list()
for (z in 1:length(strain_list[[1]])) {
  out_all = list()
  for (t in 1:3) {
    out = NULL
    a = as.data.frame(ALL_homo_high[[z]][t])
    for (gn in 1:length(a[,1])) {
      b = a[gn,]
      p = cbind(b[,1],b[,2])
      b = b[,8]
      q = unlist(strsplit(b,"\\|"))
      i = unlist(grep("HIGH",q))
      g = NULL
      for (n in c(i)) {
        c = b
        e = unlist(lapply(strsplit(c,"\\|"),FUN = function(x){return(x[c(n-1, n+1, n+4, n+7, n+8, n+9, n+10)])}))
        e = matrix(e,nrow = 1)
        d = Splicing_variant_name[Splicing_variant_name[,1]==e[2],2]
        if (length(d) == 0) {
          d =0
        }
        f = NULL
        f = cbind(p,e,d)
        g = rbind(g,f)
      }
      u = as.data.frame(table(g[,4]))
      for (o in 1:length(g[,1])) {
        for (i in 1:length(u[,1])) {
          x = g[g[o,4]==u[i,1],10]
          if (length(x) == 0) {
            x =0
          }
          if (x == u[,2]) {
            out = rbind(out,g[o,])
          }
        } 
      }
    }
    out = list(as.data.frame(out))
    out_all[sprintf("%s_%d",as.character(strain_list[[1]][z]),t)] = out
  }
  All_homo_high_list[[as.character(strain_list[[1]][z])]] = out_all
}

i=1
All_homo_high_gene_unique_list = list()
for (i in 1:length(strain_list[[1]])) {
  out = list()
  for (t in 1:length(All_homo_high_list[[i]])) {
    g = All_homo_high_list[[i]][[t]][[4]]
    g = unique(g)
    g = g[grep("si",g,invert = T)]
    g = g[grep("zgc",g,invert = T)]
    gnl = g[grep("[A-Z]",g,invert = T)]
    gnl = list(gnl)
    out[sprintf("%s_%d",as.character(strain_list[[1]][i]),t)] = gnl
  }
  All_homo_high_gene_unique_list[[as.character(strain_list[[1]][i])]] = out
}

All_homo_high_gene_unique_list[[1]]

library(sets)
All_homo_high_all_unique_list = list()
for (i in 1:length(strain_list[[1]])) {
  unique <- All_homo_high_gene_unique_list[[i]][[1]]
  for (q in 2:length(All_homo_high_gene_unique_list[[i]])) {
    unique <- set_intersection(unique, All_homo_high_gene_unique_list[[i]][[q]])
    unique = as.character(unique)
  }
  a = as.data.frame(All_homo_high_list[[i]][[1]])
  c = NULL
  for (q in 1:length(unique)) {
    b = a[a[,4]==unique[q],]
    c = rbind(c,b)
  }
  c = list(c)
  All_homo_high_all_unique_list[as.character(strain_list[[1]][i])] = c
}

gene_strain_list=list()
for(i in 1:length(strain_list[[1]])) {
  gene_strain_list[as.character(strain_list[[1]][i])] = All_homo_high_all_unique_list[[i]][4]
}
gene_strain_list
upset(fromList(gene_strain_list),nsets = 100, nintersects = 100,order.by = "freq")


gene_uniq_strain_list <- list()
group_names <- names(gene_strain_list)
for (group_name in group_names) {
  unique_strings <- gene_strain_list[[group_name]]
  for (other_group_name in setdiff(group_names, group_name)) {
    unique_strings <- setdiff(unique_strings, gene_strain_list[[other_group_name]])
  }
  gene_uniq_strain_list[[group_name]] <- unique_strings
}

print(gene_uniq_strain_list)

library(ggplot2)

gene_strain_count <- sapply(group_names, function(group) length(unique(gene_strain_list[[group]])))
gene_uniq_strain_count <- sapply(group_names, function(group) length(gene_uniq_strain_list[[group]]))
data <- data.frame(
  Group = factor(group_names, levels = group_names),
  ALL = c(gene_strain_count),
  Uniq = c(gene_uniq_strain_count)
)

library(ggplot2)
g = ggplot(NULL) +
  geom_bar(data = data,
           aes(x = Group, y = ALL),
           stat = "identity", position = "dodge", 
           width = 0.7, 
           fill = "white",
           color = "black", 
           size = 0.7) + 
  ylim(c(0,200))
g = g +
  geom_bar(data = data, 
           aes(x = Group, y = Uniq, fill = Group),
           stat = "identity", 
           position = "dodge", 
           width = 0.7, 
           color = "black", 
           size = 0.7)+ 
  labs(title = "Count of Strain Lists",
       x = "Group",
       y = "Count") +
  theme_bw()

pdf("desp_strain_all.pdf", width = 8, height = 6)
g
dev.off()

library(openxlsx)

wb <- createWorkbook()
for (i in 1:length(strain_list[[1]])) {
  addWorksheet(wb, strain_list[[1]][i])
  writeData(wb, sheet = i, rowNames = F, gene_strain_list[i])
}
saveWorkbook(wb, "n3_desp_gene.xlsx", overwrite = TRUE)


wb2 <- createWorkbook()
for (i in 1:length(strain_list[[1]])) {
  addWorksheet(wb2, strain_list[[1]][i])
  writeData(wb2, sheet = i, rowNames = F, gene_uniq_strain_list[i])
}
saveWorkbook(wb2, "strain_uni_desp_gene.xlsx", overwrite = TRUE)

i=1
t=1
t=2
gene_uniq_all_strain_list = NULL
for (i in 1:length(All_homo_high_list)) {
  a = as.data.frame(All_homo_high_list[[i]][1])
  c=NULL
  for (t in 1:length(gene_uniq_strain_list[[i]])) {
    b = a[a[,4] == gene_uniq_strain_list[[i]][t],]
    c = rbind(c, b)
  }
  gene_uniq_all_strain_list[[i]] <- c
}


wb <- createWorkbook()
for (i in 1:length(strain_list[[1]])) {
  addWorksheet(wb, strain_list[[1]][i])
  writeData(wb, sheet = i, rowNames = T, gene_uniq_all_strain_list[i])
}
saveWorkbook(wb, "gene_uniq_all_strain_list.xlsx", overwrite = TRUE)

for (i in 1:length(strain_list[[1]])) {
  a=as.data.frame(gene_uniq_all_strain_list[i])
  b=a[,5]
  write(x = b, file = sprintf("%s_%s",as.character(strain_list[[1]][i]), "gene_id.txt"), sep = "¥t")
}