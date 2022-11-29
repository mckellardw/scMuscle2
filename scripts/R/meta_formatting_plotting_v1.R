setwd("C:/Users/dwm269/github/scMuscle2")
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggsci)

source("C:/Users/dwm269/github/DWM_utils/plotting_utils/scThemes.R")
scTheme <- scThemes()

ffq = read.csv("metadata_ffq/merged_metadata.csv")
meta = read.csv("scMuscle2_metadata_v1-0.csv")
spatial = read.csv("spatial_metadata_v1-0.csv")

# clean up...
meta = meta[meta$tissue != "",]

meta$tissue[meta$tissue=="embryo"] <- "whole"

# formatted for CZI
# meta = read.csv("scMuscle2_metadata_v1-0_forCZI.csv")

#
#     format meta ------------

meta$species <- factor(
  meta$species,
  levels = table(meta$species) %>% 
    sort(decreasing = T) %>%
    names()
)

meta$tissue <- factor(
  meta$tissue,
  levels = table(meta$tissue) %>% 
    sort(decreasing = T) %>%
    names()
)


#
# Add to meta from ffq  --------------
rows_to_parse = 1:nrow(meta)
rows_to_parse = rows_to_parse[!is.null(meta$SRR.accession) & !meta$SRR.accession%in%c("#TODO","NOT_PUBLIC") & meta$SRR.accession!=" "]

tmp.meta.col = vector(mode = "character",length=nrow(meta))
for(i in rows_to_parse){
  srrs = stringr::str_split(string = meta[i,'SRR.accession'], pattern="; ") %>% unlist()
  print(i)
    
  tmp = ffq$experiment.instrument[ffq$accession %in% srrs]  
    
  if(length(tmp)>1){
    tmp = unique(tmp) %>% paste(sep="; ")
  }else if(length(tmp)<1){
    tmp = " "
  }
  print(tmp)
  
  tmp.meta.col[i] <- tmp
  
}

meta$sequencing.instrument <- tmp.meta.col




table(meta$species[!meta$chemistry %in% c("CelSeq2","Drop-seq", "ddSeq","5prime") ])

# Plotting -------------------



# Species bar plot
ggplot(
  meta,
  aes(
    x=species
  )
)+
  geom_bar(fill="black")+
  # scale_y_log10()+
  scale_y_continuous(expand = expansion(c(0.005,0)))+
  scTheme$bar +
  theme(
    axis.text.x = element_text(face = "italic"),
    axis.title.x = element_blank()
  )+
  labs(
    y="Number of datasets"
  )


# Species pie chart
p.species =
  ggplot(
  meta,
  aes(
    y=1,
    fill=species
  )
)+
  geom_bar(
    # width=0.1, color="white"
  )+
  coord_polar()+
  scTheme$pie +
  theme(
    legend.position="right",
    legend.text = element_text(face="italic"),
    legend.title = element_text(face="bold")
  )+
  labs(
    fill="Species"
  )+
  scale_fill_manual(
    values=pal_npg()(length(unique(meta$species)))%>%rev()
  )


# Tissue pie chart
p.tissue =
  ggplot(
  meta,
  aes(
    y=1,
    fill=tissue
  )
)+
  geom_bar(
    # width=1, color="black"
  )+
  coord_polar()+
  scTheme$pie +
  theme(
    legend.position="right",
    legend.title = element_text(face="bold")
  )+
  labs(
    fill="Tissue source"
  )+
  scale_fill_manual(
    values=rainbow(n=unique(meta$tissue)%>%length(),end = 0.75)%>% 
      shades::brightness(values = 0.95) %>%
      shades::saturation(values = 0.7)
  )


# patch & save
p.species+p.tissue

ggsave(
  filename = "C:\\Users\\dwm269\\Box\\DWM\\Proposals\\2022\\czi_resubmission\\figures\\pies_v2.pdf",
  device="pdf",
  units="in",
  height=3,
  width=18
)









