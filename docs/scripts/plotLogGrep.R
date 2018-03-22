# Make the input file from a bunch of runs
# for i in annofilt_prokka_outputs/Lys*/annofilt.log 
#   do 
#     line=`cat ${i} | grep "Total from"`
#     echo $i $line >> logs.logs
# done
data <- read.table("~/results/2018-03-20-annofilt/logs.logs" ,skip = 6)
file <- as.character(gsub(".*/(.*?)/annofilt.log", "\\1", data[, 1]))
contigs <- as.numeric(data[, 9])
total <- as.numeric(data[, 11])
kept <- as.numeric(data[, 13])
lost <- as.numeric(data[, 15])
datacols <- c(9, 11, 13, 15)

wide <- data.frame(file, contigs, total, kept, lost)
wide <- wide %>%
  group_by(file) %>%
  mutate(genes_searched = contigs * 2,
         genes_searched_kept = genes_searched - lost,
         percent_rejected = 1 - (genes_searched_kept / genes_searched))

wide$pool <- as.factor(ifelse(grepl("Lys\\d*", wide$file), "Soil", "Enteric"))

library(ggplot2)
library(reshape2)
library(dplyr)
tall <- melt(
  wide, id.vars = c("file", "pool"), 
  variable.name = "category",
    value.name = "count")

ggplot(tall[tall$category %in% c("total", "kept"),], aes(fill=category)) + 
  geom_histogram(aes(x=count), alpha=.5 )
ggplot(tall[, ], aes(fill=category)) + 
  geom_density(aes(x=count), alpha=.5 ) +
  facet_wrap(~category, scales="free")
f <- tall[tall$category %in% c("genes_searched", "genes_searched_kept", "percent_rejected"), ]
summary(f)

ggplot(wide, aes(fill=pool,
                 x=as.numeric(as.factor(reorder(file, genes_searched))))) + 
  geom_area(aes(y=genes_searched_kept)) +
  geom_line(aes(y=genes_searched), color="black") +
  facet_wrap(~pool, scales = "free") +
  scale_fill_manual(values=c("#d2797966", "#b3ff6666"), guide=F) +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17)) +
  scale_x_continuous( expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  labs(x="Genomes", y="Count")

ggplot(wide, aes(fill=pool,
                 x=as.numeric(as.factor(reorder(file, genes_searched))))) + 
  geom_area(aes(y=genes_searched), fill="grey80") +
  # geom_line(aes(y=genes_searched), color="grey20", size=1) +
  geom_area(aes(y=genes_searched_kept)) +
  facet_wrap(~pool, scales = "free") +
  scale_fill_manual(values=c("#d27979", "#b3ff66"), guide=F) +
  theme_minimal()+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17),
        strip.text = element_text(size=17)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_log10(expand = c(0, 0)) +
  labs(x="Genomes", y="Count") 



# plotting the percent genes that are rejected by annofilt
ggplot(f[f$category=="percent_rejected", ],
       aes(x=pool, y=count * 100, group=pool, fill=pool)) + 
  geom_boxplot(lwd=1, color="black", outlier.shape = 2) +
  geom_jitter(width=.1, height=.1, alpha=.2) +
  scale_fill_manual(values=c("#d2797966", "#b3ff6666"), guide=F) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(limit=c(0, 100), expand = c(0, 0)) +
  theme_minimal()+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17)) +
  labs(
    # title="Percent genes rejected by annofilt",
     y="Percent Genes Rejected", 
     x="Collection")
  

