---
title: "R Notebook"
output: html_notebook
---

# 1. Introduction

The purpose of this Markdown is to document variant processing and QC of the 450k exome release of the UK Biobank data and subsequent testing of T2D. The majority of the data files that are presented in this file contain protected UKBB participant data, and thus cannot be included in this repository. Please access our returned results via the UK Biobank and application 9905 to get this data.

## 1A. Setup

```{r setup}

library(data.table)
library(tidyverse)
library(patchwork)
library(lemon)
library(readxl)
library(broom)
library(svglite)
library(lubridate)

theme <- theme(panel.background=element_rect(fill="white"),line=element_line(size=1,colour="black",lineend="round"),axis.line=element_line(size=1),text=element_text(size=16,face="bold",colour="black"),axis.text=element_text(colour="black"),axis.ticks=element_line(size=1,colour="black"),axis.ticks.length=unit(.1,"cm"),strip.background=element_rect(fill="white"),axis.text.x=element_text(angle=45,hjust=1),legend.position="blank",panel.grid.major=element_line(colour="grey",size=0.5),legend.key=element_blank())

## Default theme w/legend
theme.legend <- theme + theme(legend.position="right")

```

# 2. Variant Processing and Filtering

To perform variant processing and quality control, I used the pipeline documented here: <https://github.com/mrcepid-rap#vcf-filtering-and-rare-variant-burden-testing>. Please see this GitHub repository for more details as to how this pipeline works. All work for this project is being completed on the UK Biobank Research Access Platform (RAP) project MRC - Variant Filtering 450k (project-G6BJF50JJv8p4PjGB9yy7YQ2). Analysis will proceed in 6 steps (after '-' is the applet being used on the RAP):

1.  Splitting BCFs - [mrcepid-bcfsplitter](http://github.com/mrcepid-rap/mrcepid-bcfsplitter)
2.  Filtering and VEP annotation - [mrcepid-filterbcf](https://github.com/mrcepid-rap/mrcepid-filterbcf)
3.  CADD annotation and statistics output - [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-annotatecadd)
4.  Collapsing variants according to pre-defined filters - [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants)
5.  Merging collapsed variants into resource bundles for association testing - [mrcepid-mergecollapsevariants](https://github.com/mrcepid-rap/mrcepid-mergecollapsevariants)
6.  Association testing - [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting)

There is also one supplementary step that prepares genetic resource files and sample lists for association testing:

-   Building of GRMs and sample inclusion/exclusion lists - [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms)

# 3. Analysis Functions

## 3A. Loading Transcripts

This is just to enable easy plotting of Manhattan plots. This code section also calculates the "mean" position of each gene genome-wide

```{r load transcripts}

transcripts <- fread("data_files/transcripts.tsv.gz")
setnames(transcripts,"#chrom","chrom")

# Get chromosome locations for the plot
mean.chr.pos <- transcripts[,mean(manh.pos), by = chrom]
chrom.to.omit <- c("17","19","21")
mean.chr.pos[,labels:=if_else(chrom %in% chrom.to.omit, "", chrom)]
setkey(mean.chr.pos,V1)
mean.chr.pos[,chrom:=factor(chrom, levels = mean.chr.pos[,chrom])]

```

## 3B. Functions for Loading Data

This is just a master function for loading BOLT results into R.

```{r loading data}

load.and.plot.data <- function(file.names = c(), p.val.col, AC.col, tool.name, marker.file = NULL, ymin = 0, ymax = 15) {
  
  # Read all data files and cat
  result.table <- data.table()
  for (file in file.names) {
    curr.tab <- fread(file)
    result.table <- rbind(result.table,curr.tab)
  }
  
  # Set gene and p. value column to something standard:
  setnames(result.table, p.val.col, "p.value.selected")
  result.table[,log.p:=-log10(p.value.selected)] # Add -log10 p value
  result.table[,chrom:=factor(chrom, levels = mean.chr.pos[,chrom])] 
  
  masks <- result.table[MASK != "",unique(MASK)]
  mafs <- result.table[MAF != "",unique(MAF)]
  
  # Make standard plots
  plots = list()
  
  for (mask in masks) {
    for (maf in mafs) {
      name = paste(mask,maf,sep="-")
      manh.plot <- plot.manh(result.table[get(AC.col) > 10 & MASK == mask & MAF == maf], "log.p",ymin=ymin,ymax=ymax)
      qq.plot <- plot.qq(result.table[get(AC.col) > 10 & MASK == mask & MAF == maf], "p.value.selected",ymin=ymin,ymax=ymax)
      comb.plot <- manh.plot + qq.plot + plot_layout(ncol = 2, nrow = 1, widths = c(3,1.3)) + plot_annotation(title = paste(tool.name, mask, maf, sep = " - "), theme = theme(plot.title=element_text(size=18,face="bold",colour="black")))
      
      plots[[name]] = list('manh.plot' = manh.plot,
                           'qq.plot' = qq.plot,
                           'comb.plot' = comb.plot)
    }
  }
  
  if (!is.null(marker.file)) {
    return(list('gene.table' = result.table, 
                'marker.file' = marker.file,
                'plots' = plots,
                'masks' = masks))
  } else {
    return(list('gene.table' = result.table, 
                'marker.file' = NULL,
                'plots' = plots,
                'masks' = masks))
  }
  
}

```

## 3C. Functions for Plotting

These two functions make manhattan and qq plots, respectively.

```{r plotting functions}

plot.manh <- function(stats, p.var, ymin = 0, ymax = 15) {
 
  manh.plot <- ggplot(stats, aes(manh.pos, get(p.var), colour = chrom)) + 
    geom_point() +
    geom_hline(yintercept = -log10(1.6e-6), colour = "red", linetype = 2) +
    geom_text(inherit.aes = F, data = stats[get(p.var)>-log10(1.4e-6)], aes(manh.pos, get(p.var), label = SYMBOL), position = position_nudge(0.01,0.2),hjust=0,angle=45) +
    scale_x_continuous(name = "Chromosome", label = mean.chr.pos[chrom != "Y",labels], breaks = mean.chr.pos[chrom != "Y",V1]) +
    scale_y_continuous(name = expression(bold(-log[10](italic(p)))), limits = c(ymin,ymax)) +
    scale_colour_manual(values = c(rep(c("#53878D","#7AC6CC"),11),"#53878D")) +
    coord_capped_cart(bottom="both",left="both") + # Comes from the "lemon" package
    theme + theme(panel.grid.major = element_blank())

  return(manh.plot)
}

plot.qq <- function(stats, p.var, ymin = 0, ymax = 15) {
  
  ## QQplot
  qqplot.data <- data.table(observed = stats[,get(p.var)],
                            SYMBOL = stats[,SYMBOL])
  setkey(qqplot.data,observed)
  qqplot.data <- qqplot.data[!is.na(observed)]
  qqplot.data[,observed:=-log10(observed)]
  qqplot.data[,expected:=-log10(ppoints(nrow(qqplot.data)))]
  
  qq.plot <- ggplot(qqplot.data, aes(expected, observed)) +
    geom_point() +
    geom_text(inherit.aes = F, data = qqplot.data[observed > -log10(1.4e-6)], aes(expected, observed, label = SYMBOL), position = position_nudge(-0.04,0),hjust=1) +
    geom_abline(slope = 1, intercept = 0,colour="red") +
    scale_x_continuous(name = expression(bold(Expected~-log[10](italic(p)))), limits = c(0,5)) +
    scale_y_continuous(name = "", limits = c(ymin,ymax)) +
    theme + theme(panel.grid.major = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
  return(qq.plot)
  
}

tabix.header <- c("varID","CHROM","POS","REF","ALT","ogVarID","FILTER","AF","F_MISSING","AN","AC","MANE","ENST","ENSG","BIOTYPE","SYMBOL","CSQ","gnomAD_AF","CADD","REVEL","SIFT","POLYPHEN","LOFTEE","PARSED_CSQ","MULTI","INDEL","MINOR","MAJOR","MAF","MAC","AA","AApos","BOLT_MAF","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM","p.value.selected","BOLT_AC")

Sys.setenv(PATH="/Users/eg15/Applications/htslib-1.10.2/bin/:/Users/eg15/Applications/bedtools2/bin/:/bin/:/usr/bin/")

plot.gene <- function(mask, maf, symbol, table, marker.file, query) {
  gene.info <- table[SYMBOL == symbol & MASK == mask & MAF == maf]
  system(paste("tabix", marker.file, gene.info[,coord],"> /tmp/parsed.tsv", sep = " "))
  variants <- fread("/tmp/parsed.tsv")
  system("rm /tmp/parsed.tsv")
  setnames(variants, names(variants), tabix.header)
  
  # Make sure the right gene is included
  variants <- variants[ENST == gene.info[,ENST]]
  
  # Make sure the right variants are include
  variants <- variants[eval(parse(text=query))]
  variants <- variants[BOLT_AC > 0]

  return(list(plot.gene.model(variants, 
                              gene.info[,ENST],
                              gene.info[,SYMBOL],
                              gene.info[,MASK],
                              gene.info[,MAF],
                              "BOLT"), 
              variants))
}

```

# 4. Phenotype Associations

## 4A. T2D

### ExWAS Results

```{r T2D, fig.height=6, fig.width=15}
bolt.T2D.ret <- load.and.plot.data(file.names = c("ukbb_data/T2D_ExWAS/T2D.bolt.genes.BOLT.stats.tsv.gz"),
                                   p.val.col="P_BOLT_LMM_INF",
                                   tool.name = "BOLT",
                                   AC.col = "AC",
                                   marker.file = "ukbb_data/T2D_ExWAS/T2D.bolt.markers.BOLT.stats.tsv.gz",
                                   ymax = 60)

for (mask in names(bolt.T2D.ret$plots)[grepl("MAF_01", names(bolt.T2D.ret$plots))]) {
  print(bolt.T2D.ret$plots[[mask]]$comb.plot)
}

```

### Additional Analyses

#### Protein Domain Associations

Here we are reading in individual-level variant information to determine protein domain-specific associations.

```{r protein domains}

pheno_covar <- fread("ukbb_data/pheno_tables/T2D.phenotypes_covariates.formatted.tsv")
IGF1R.carriers <- fread("ukbb_data/T2D_ExWAS/indv_gene_results/MISS_REVEL0_7.IGF1R.T2D.carriers_formated.tsv")
IGF1R.variants <- fread("ukbb_data/T2D_ExWAS/indv_gene_results/MISS_REVEL0_7.IGF1R.T2D.variant_table.tsv")
IGF1R.carriers <- merge(IGF1R.carriers, IGF1R.variants[,c("varID","AA","AApos")],by="varID")

pheno_covar[,has.var:=if_else(IID%in%IGF1R.carriers[,IID],1,0)]
pheno_covar[,has.var_kinase:=if_else(IID%in%IGF1R.carriers[AApos>=999 & AApos <=1274,IID],1,0)]
pheno_covar[,has.var_not_kinase:=if_else(IID%in%IGF1R.carriers[AApos <999 | AApos >1274,IID],1,0)]
pheno_covar[,has.var_cytoplasmic:=if_else(IID%in%IGF1R.carriers[AApos >= 960,IID],1,0)]
pheno_covar[,has.var_extracellular:=if_else(IID%in%IGF1R.carriers[AApos >= 741 & AApos <= 935,IID],1,0)]

model <- glm(T2D_for_GWAS ~ has.var + age + age_squared + sex + wes_batch + PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_covar, family="binomial")
model <- data.table(tidy(model))
domain.table <- model[term == "has.var"]

model <- glm(T2D_for_GWAS ~ has.var_kinase + age + age_squared + sex + wes_batch + PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_covar, family="binomial")
model <- data.table(tidy(model))
domain.table <- rbind(domain.table, model[term == "has.var_kinase"])

model <- glm(T2D_for_GWAS ~ has.var_not_kinase + age + age_squared + sex + wes_batch + PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_covar, family="binomial")
model <- data.table(tidy(model))
domain.table <- rbind(domain.table, model[term == "has.var_not_kinase"])

model <- glm(T2D_for_GWAS ~ has.var_cytoplasmic + age + age_squared + sex + wes_batch + PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_covar, family="binomial")
model <- data.table(tidy(model))
domain.table <- rbind(domain.table, model[term == "has.var_cytoplasmic"])

model <- glm(T2D_for_GWAS ~ has.var_extracellular + age + age_squared + sex + wes_batch + PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_covar, family="binomial")
model <- data.table(tidy(model))
domain.table <- rbind(domain.table, model[term == "has.var_extracellular"])

domain.table[,OR:=exp(estimate)]
domain.table[,ci.upper:=exp(estimate+(1.96*std.error))]
domain.table[,ci.lower:=exp(estimate-(1.96*std.error))]

ztest <- (domain.table[term == 'has.var_kinase', estimate]) - (domain.table[term == 'has.var_not_kinase', estimate])/((domain.table[term == 'has.var_kinase', std.error])^2 + (domain.table[term == 'has.var_kinase', std.error])^2)^(1/2) 
pnorm(ztest) * 2
```

#### Excluding T2D Cases

Purpose here is to exclude T2D cases and see if "controls" still have elevated glucose/hba1c for novel genes

```{r T2D Case Exclusion}

# Load carriers
IGF1R.carriers <- fread("ukbb_data/T2D_ExWAS/indv_gene_results/MISS_REVEL0_7.IGF1R.T2D.carriers_formated.tsv")
MLXIPL.carriers <- fread("ukbb_data/T2D_ExWAS/indv_gene_results/MISS_REVEL0_7.MLXIPL.T2D.carriers_formated.tsv")
ZEB2.carriers <- fread("ukbb_data/T2D_ExWAS/indv_gene_results/MISS_REVEL0_5.ZEB2.T2D.carriers_formated.tsv")
TNRC6B.carriers <- fread("ukbb_data/T2D_ExWAS/indv_gene_results/HC_PTV.TNRC6B.T2D.carriers_formated.tsv")

# Get T2D status
pheno_t2d <- fread("ukbb_data/pheno_tables/T2D.phenotypes_covariates.formatted.tsv")

# hba1c first
pheno_hba1c <- fread("ukbb_data/pheno_tables/HBA1C.phenotypes_covariates.formatted.tsv")
pheno_hba1c <- merge(pheno_hba1c, pheno_t2d[,c("IID","T2D_for_GWAS")], by = "IID")

pheno_hba1c[,has.IGF1R:=if_else(IID %in% IGF1R.carriers[,IID], 1, 0)]
pheno_hba1c[,has.MLXIPL:=if_else(IID %in% MLXIPL.carriers[,IID], 1, 0)]
pheno_hba1c[,has.ZEB2:=if_else(IID %in% ZEB2.carriers[,IID], 1, 0)]
pheno_hba1c[,has.TNRC6B:=if_else(IID %in% TNRC6B.carriers[,IID], 1, 0)]

# glucose second
pheno_glucose <- fread("ukbb_data/pheno_tables/GLUCOSE.phenotypes_covariates.formatted.tsv")
pheno_glucose <- merge(pheno_glucose, pheno_t2d[,c("IID","T2D_for_GWAS")], by = "IID")

pheno_glucose[,has.IGF1R:=if_else(IID %in% IGF1R.carriers[,IID], 1, 0)]
pheno_glucose[,has.MLXIPL:=if_else(IID %in% MLXIPL.carriers[,IID], 1, 0)]
pheno_glucose[,has.ZEB2:=if_else(IID %in% ZEB2.carriers[,IID], 1, 0)]
pheno_glucose[,has.TNRC6B:=if_else(IID %in% TNRC6B.carriers[,IID], 1, 0)]

# This function does glucose/hba1c level test w & w/o T2D carriers
run_t2d_test <- function(pheno.table, y.var, gene, control_var, model.name) {
  gene.var <- paste("has",gene,sep=".")
  covariates <- c(gene.var, control_var, "age", "age_squared", "sex", "wes_batch", paste0("PC",c(1:10)))
  formatted.formula <- as.formula(paste(y.var,paste(covariates, collapse ="+"),sep="~"))
  if (model.name == "all") {
    model <- glm(formatted.formula, data = pheno.table, family="gaussian")
  } else {
    model <- glm(formatted.formula, data = pheno.table[T2D_for_GWAS==0])
  }
  model <- data.table(tidy(model))
  model[,gene:=gene]
  model[,model:=model.name]
  model[,y.var:=y.var]
  return(model[term == gene.var])
}

model.table <- rbind(run_t2d_test(pheno_hba1c, "hba1c_qc", "MLXIPL", "hba1c_medications", "all"),
                     run_t2d_test(pheno_hba1c, "hba1c_qc", "IGF1R", "hba1c_medications", "all"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "MLXIPL", "glucose_medications", "all"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "IGF1R", "glucose_medications", "all"),
                     
                     run_t2d_test(pheno_hba1c, "hba1c_qc", "ZEB2", "hba1c_medications", "all"),
                     run_t2d_test(pheno_hba1c, "hba1c_qc", "TNRC6B", "hba1c_medications", "all"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "ZEB2", "glucose_medications", "all"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "TNRC6B", "glucose_medications", "all"),
                     
                     run_t2d_test(pheno_hba1c, "hba1c_qc", "MLXIPL", "hba1c_medications", "no.T2D"),
                     run_t2d_test(pheno_hba1c, "hba1c_qc", "IGF1R", "hba1c_medications", "no.T2D"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "MLXIPL", "glucose_medications", "no.T2D"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "IGF1R", "glucose_medications", "no.T2D"),
                     
                     run_t2d_test(pheno_hba1c, "hba1c_qc", "ZEB2", "hba1c_medications", "no.T2D"),
                     run_t2d_test(pheno_hba1c, "hba1c_qc", "TNRC6B", "hba1c_medications", "no.T2D"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "ZEB2", "glucose_medications", "no.T2D"),
                     run_t2d_test(pheno_glucose, "glucose_qc", "TNRC6B", "glucose_medications", "no.T2D"))
```

#### REVEL

Just gets the REVEL score percentiles so we can mention in the paper

```{r REVEL, fig.height=5, fig.width=12}

REVEL.scores <- fread("ukbb_data/misc/REVEL_scores.txt",col.names = c("REVEL"))
REVEL.scores[,dummy:=1]
REVEL.scores[,REVEL.bin:=cut(REVEL,breaks=seq(0,1,by=0.02))]

revel.bins <- REVEL.scores[,sum(dummy),by=REVEL.bin]
setkey(revel.bins,REVEL.bin)
revel.bins[,REVEL.bin:=if_else(is.na(REVEL.bin), "NA", as.character(REVEL.bin))]
revel.bins[,REVEL.bin:=factor(REVEL.bin, levels = revel.bins[,REVEL.bin])]

revel.bins.dens <- data.table()
tot.miss <- revel.bins[,sum(V1)]
for (i in c(1:nrow(revel.bins))) {
  curr.row <- revel.bins[i]
  prop <- curr.row[,V1] / tot.miss
  min <- as.double(str_remove(str_split(curr.row[,REVEL.bin], ",", simplify = T)[1], "\\("))
  max <- as.double(str_remove(str_split(curr.row[,REVEL.bin], ",", simplify = T)[2], "\\]"))
  cum.dens <- revel.bins[i:nrow(revel.bins),sum(V1)] / tot.miss
  revel.bins.dens <- rbind(revel.bins.dens,
                           data.table(REVEL.bin = curr.row[,REVEL.bin],
                                      min = min,
                                      max = max,
                                      tot = curr.row[,V1],
                                      prop = prop,
                                      cum.dens = cum.dens))
}
revel.bins.dens[,min:=if_else(is.na(min), -0.05, min)]
revel.bins.dens[,max:=if_else(is.na(max), 0, max)]


prop.plot <- ggplot(revel.bins.dens, aes(REVEL.bin, prop)) + 
  geom_col() + 
  scale_x_discrete(breaks = c("NA", levels(revel.bins.dens[,REVEL.bin]))) + 
  scale_y_continuous(name = "Proportion of Missense Variants") + 
  theme

dens.plot <- ggplot(revel.bins.dens, aes(min, cum.dens)) + 
  geom_vline(data = revel.bins.dens[min == 0.7 | min == 0.5 | min == 0.2], aes(xintercept = min), colour = "red", linetype = 2) + 
  geom_text(data = revel.bins.dens[min == 0.7 | min == 0.5 | min == 0.2], aes(x = min, y = 0.9, label = sprintf("%0.2f%%", cum.dens * 100)), hjust = -0.2, colour = "red") + 
  geom_line(size = 2) + 
  scale_x_continuous(name = "bin max REVEL score", breaks = c(-0.05, 0, 0.25, 0.5, 0.75, 1), labels = c("NA", 0, 0.25,0.5,0.75,1)) +
  scale_y_continuous(name = "Proportion of Variants in\nBin and Above") +
  theme

prop.plot + dens.plot + plot_layout(ncol = 2, nrow = 1)

```

# 5. Figures

Code block just sets constants used for figure making

```{r Setting Constants}

# valid masks
valid.masks <- c("HC_PTV","MISS_REVEL0_5","MISS_REVEL0_7","SYN")

# Abstract
# p. value
p.val.thresh <- 0.05 / nrow(bolt.T2D.ret$gene.table[MASK %in% valid.masks & AC > 30])

```

## 5A. Main Text

### Figure 1

This figure is pulled into Inkscape/Adobe Illustrator for further artistic editing for the final figure. Also calculate lambda for each model for the paper.

```{r Figure 1, fig.height=5, fig.width=10}

plot.table <- bolt.T2D.ret$gene.table[AC > 30 & MASK %in% valid.masks & MASK != "SYN"]
shape.table <- data.table(crossing(MAF=c("MAF_01","AC_1"), MASK = c("HC_PTV", "MISS_REVEL0_5", "MISS_REVEL0_7")))
shape.table[,mask.shape:=c(0, 1, 2, 15, 16, 17)]
plot.table <- merge(plot.table,shape.table,by=c("MASK","MAF"))

manh.plot <- ggplot(plot.table, aes(manh.pos, log.p, colour = chrom)) + 
  geom_point(aes(size = if_else(log.p>-log10(p.val.thresh), 1, 0.5), shape = mask.shape)) +
  geom_hline(yintercept = -log10(p.val.thresh), colour = "red", linetype = 2) +
  scale_x_continuous(name = "Chromosome", label = mean.chr.pos[chrom!="Y",labels], breaks = mean.chr.pos[chrom!="Y",V1]) +
  scale_y_continuous(name = expression(bold(-log[10](italic(p)))), limits = c(0,60)) +
  scale_colour_manual(values = rep(c("#53878D","#7AC6CC"),12)) +
  coord_capped_cart(bottom="both",left="both") + # Comes from the "lemon" package
  scale_shape_identity() +
  scale_size_continuous(range=c(0,2)) +
  guides(colour="none",size="none") +
  theme.legend + theme(text = element_text(size = 12), panel.grid.major = element_blank(), legend.position=c(0.20,0.80))

## QQplot
figure.QQ <- function(mask, show.x = F, show.text = F) {
  
  qqplot.data <- data.table(observed = bolt.T2D.ret$gene.table[AC > 30 & MASK == mask & MAF == "MAF_01" & !is.na(log.p), p.value.selected],
                          SYMBOL = bolt.T2D.ret$gene.table[AC > 30 & MASK == mask & MAF == "MAF_01" & !is.na(log.p),SYMBOL])
  setkey(qqplot.data,observed)
  qqplot.data <- qqplot.data[!is.na(observed)]
  qqplot.data[,chisq:=qchisq(1-observed,1)]
  lambda <- qqplot.data[,median(chisq) / qchisq(0.5,1)]
  print(paste0("Lambda for mask ", mask, " : ", lambda))
  qqplot.data[,observed:=-log10(observed)]
  qqplot.data[,expected:=-log10(ppoints(nrow(qqplot.data)))]
  
  qqplot <- ggplot(qqplot.data, aes(expected, observed)) +
    geom_point(size = 0.5) +
    geom_abline(slope = 1, intercept = 0,colour="red") +
    scale_x_continuous(name = if_else(show.x, expression(bold(Expected~-log[10](italic(p)))), expression("")), limits = c(0,5)) +
    scale_y_continuous(name = "", limits = c(0,60)) +
    theme + theme(text = element_text(size = 10), panel.grid.major = element_blank(), axis.ticks.y = element_blank(), axis.title = element_text(size = 10))
  if (show.text) {
    qqplot <- qqplot + geom_text(inherit.aes = F, data = qqplot.data[observed > -log10(p.val.thresh)], aes(expected, observed, label = SYMBOL), position = position_nudge(-0.04,0),hjust=1,size = 3)
  }
  return(qqplot)
  
}

qqplot.PTV <- figure.QQ("HC_PTV", F)
qqplot.MISS05 <- figure.QQ("MISS_REVEL0_5", F)
qqplot.MISS07 <- figure.QQ("MISS_REVEL0_7", F)
qqplot.SYN <- figure.QQ("SYN", T)

right.plot <- qqplot.PTV + qqplot.MISS05 + qqplot.MISS07 + qqplot.SYN + plot_layout(ncol = 1)
figure1 <- manh.plot + right.plot + plot_layout(ncol = 2, widths = c(0.85,0.15))
figure1
ggsave("Figures/Figure1_pre.png",figure1,width = 10, height = 4, dpi = 1200)

manh.plot <- manh.plot + geom_text(inherit.aes = F, data = bolt.T2D.ret$gene.table[AC > 30 & MASK != "SYN" & log.p>-log10(1.4e-6)], aes(manh.pos, log.p, label = SYMBOL), position = position_nudge(0.01,0.2),hjust=0,angle=45,size=4)
qqplot.PTV <- figure.QQ("HC_PTV", F, T)
qqplot.MISS05 <- figure.QQ("MISS_REVEL0_5", F, T)
qqplot.MISS07 <- figure.QQ("MISS_REVEL0_7", F, T)
qqplot.SYN <- figure.QQ("SYN", T, T)

right.plot <- qqplot.PTV + qqplot.MISS05 + qqplot.MISS07 + qqplot.SYN + plot_layout(ncol = 1)
fake.figure1 <- manh.plot + right.plot + plot_layout(ncol = 2, widths = c(0.85,0.15))
fake.figure1

```

### Figure 2

#### Run Individual Genes

Just an example of the command used with the MRCEpid-RAP project to generate summary stats for all genes. This files are downloaded from DNANexus and placed into the ukbb_data folder (which is not included in this repo due to data protection rules)

```{bash indvidual genes}

dx run mrcepid-runassociationtesting --destination t2d_genes/ -imode=extract -iassociation_tarballs=file-GFkfYj0JJv8fQgGG4z458P60 -iis_binary=true -ioutput_prefix=T2D -igene_ids=ENST00000403799,ENST00000257555,ENST00000678049,ENST00000454349,ENST00000650285,ENST00000313375,ENST00000627532 -iinclusion_list=file-GBVYYPQJ57y5ZV904gzFb873 -iphenofile=file-G5PFv00JXk80Qfx04p9X56X5 -isex=2

```

#### Plot OR vs MAF

```{r fig.height=4, fig.width=6}

load.model <- function(gene) {
  
  # select best model:
  min.result <- gene.table[ENST == gene & p.value.selected == gene.table[ENST==gene,min(p.value.selected)]]
  
  table <- fread(paste0("ukbb_data/T2D_ExWAS/indv_gene_results/",min.result[,MASK], ".", min.result[,SYMBOL], ".T2D.association_stats.tsv"))
  
  return(list(table[,n_car],
              table[,n_model],
              table[,effect],
              table[,std_err]))
  
}

genes[,c("n_car","n_model","effect","std_err"):=load.model(ENST), by=1:nrow(genes)]
T2D.table <- copy(genes)

T2D.table[,OR:=exp(effect)]
T2D.table[,ci.upper:=exp(effect+(1.96*std_err))]
T2D.table[,ci.lower:=exp(effect-(1.96*std_err))]
T2D.table[,ci.log.upper:=effect + (1.96 * std_err)]
T2D.table[,ci.log.lower:=effect - (1.96 * std_err)]
T2D.table[,cMAF:=n_car / (n_model * 2)]
T2D.table <- merge(T2D.table, shape.table, by =c("MAF","MASK"))

GCK.label <- paste0('GCK - OR=',
                    sprintf("%0.2f",T2D.table[SYMBOL == "GCK",OR]), ' [',
                    sprintf("%0.2f",T2D.table[SYMBOL == "GCK",ci.lower]),'-',
                    sprintf("%0.2f",T2D.table[SYMBOL == "GCK",ci.upper]), ']')

figure.2 <- ggplot(T2D.table[SYMBOL != "GCK"], aes(cMAF, OR, shape = mask.shape)) + 
  geom_point(size = 4, colour = "#53878D") + 
  geom_text(aes(label = SYMBOL), vjust = -0.25, hjust = -0.25) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), colour = "#53878D", width = 0) + 
  annotate('segment', x = T2D.table[SYMBOL == "GCK",cMAF], xend = T2D.table[SYMBOL == "GCK",cMAF], y=27.5, yend=30, size = 1, arrow=arrow(length = unit(0.25,"cm"))) +
  annotate('text', x = T2D.table[SYMBOL == "GCK",cMAF], y=28.75, label = GCK.label, hjust = -0.05) +
  scale_x_log10(name = "cMAF", limits=c(3e-5,1e-3), breaks = c(3e-5,1e-4,3e-4,1e-3), labels = c("0.003%", "0.01%", "0.03%","0.1%")) +
  scale_y_continuous(name = "Type 2 Diabetes Risk (OR)", limits = c(1,30)) +
  scale_shape_identity() +
  scale_colour_discrete(guide = guide_legend(direction = "vertical", title = "Variant Class"), labels=c("HC PTVs", "Missense REVEL ≥ 0.5", "Missense REVEL ≥ 0.7")) + 
  theme.legend + theme(legend.position=c(0.75,0.75),legend.text = element_text(size=10,face="bold",colour="black"))

figure.2

ggsave("Figures/Figure2.svg", figure.2, height = 4, width = 6)

```

### Figure 3

From Katherine Kentistou (Locus Zoom plot)

## 5B. Supplementary Figures

### Supplementary Figure 1

```{r Supp Fig 1, fig.height=4, fig.width=8}

plot.table <- bolt.T2D.ret$gene.table[AC > 30 & MASK %in% valid.masks & MASK == "SYN"]
manh.plot <- ggplot(plot.table[MASK == "SYN"], aes(manh.pos, log.p, colour = chrom)) + 
  geom_point() +
  geom_hline(yintercept = -log10(p.val.thresh), colour = "red", linetype = 2) +
  scale_x_continuous(name = "Chromosome", label = mean.chr.pos[chrom!="Y",labels], breaks = mean.chr.pos[chrom!="Y",V1]) +
  scale_y_continuous(name = expression(bold(-log[10](italic(p)))), limits = c(0,60)) +
  scale_colour_manual(values = rep(c("#53878D","#7AC6CC"),12)) +
  coord_capped_cart(bottom="both",left="both") + # Comes from the "lemon" package
  scale_shape_identity() +
  scale_size_continuous(range=c(0,2)) +
  guides(colour="none",size="none") +
  theme.legend + theme(text = element_text(size = 12), panel.grid.major = element_blank(), legend.position=c(0.20,0.80))

manh.plot
ggsave("Figures/SuppFig1.png", manh.plot, width = 8, height = 4, dpi = 450)
```

### Supplementary Figure 2

```{r Supp Fig 2, fig.height=6, fig.width=8}

staar.T2D.ret <- load.and.plot.data(file.names = c("ukbb_data/T2D_ExWAS/T2D.staar.genes.STAAR.stats.tsv.gz"),
                                    p.val.col="staar.O.p",
                                    tool.name = "STAAR",
                                    AC.col = "cMAC",
                                    ymax = 60)

glm.T2D.ret <- load.and.plot.data(file.names = c("ukbb_data/T2D_ExWAS/T2D.glm.genes.glm.stats.tsv.gz"),
                                    p.val.col="p_val_full",
                                    tool.name = "GLM",
                                    AC.col = "n_car",
                                    ymax = 60)


query.data <- function(symbol, mask, maf) {
  
  staar.p <- staar.T2D.ret$gene.table[SYMBOL == symbol & MASK == mask & MAF == maf,p.value.selected]
  glm.p <- glm.T2D.ret$gene.table[SYMBOL == symbol & MASK == mask & MAF == maf]
  return(list(staar.p,glm.p[,p.value.selected],glm.p[,effect],glm.p[,std_err]))
  
}

T2D.table[,c("staar.p","glm.p"):=query.data(SYMBOL, MASK, MAF),by=1:nrow(T2D.table)]

suppfigure2.table <- T2D.table[,c("SYMBOL","MASK","MAF","p.value.selected","staar.p","glm.p")]
setnames(suppfigure2.table,"p.value.selected","bolt.p")
suppfigure2.table <- as.data.table(pivot_longer(suppfigure2.table,-c(SYMBOL,MASK,MAF)))

suppfigure2.table[,mask.text:=if_else(MASK == "HC_PTV", "HC PTV", if_else(MASK == "MISS_REVEL0_5", "Missense REVEL ≥ 0.5", if_else(MASK == "MISS_REVEL0_7", "Missense REVEL ≥ 0.7", "NA")))]
suppfigure2.table[,maf.text:=if_else(MAF == "MAF_01", "MAF < 0.1%", "Singletons")]

suppfigure2.table[,axis.name:=paste0(SYMBOL, "\n(",mask.text, " ", maf.text, ")")]

supp.fig2 <- ggplot(suppfigure2.table, aes(axis.name, -log10(value), group = name, colour = name)) + 
  geom_hline(yintercept = -log10(p.val.thresh), size = 1, linetype = 2, colour = "red") +
  geom_point(position=position_dodge(0.5), size = 2) +
  coord_flex_cart(bottom = brackets_horisontal(tick.length = 0)) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = expression(bold(-log[10](italic(p))))) +
  scale_colour_manual(values = c("#53878D","#7AC6CC", "#E8B587"), labels=c("BOLT","GLM","STAAR"), guide=guide_legend(title = "Tool")) +
  theme.legend + theme(panel.grid.major.x=element_blank())
supp.fig2

ggsave(filename = "Figures/SuppFig2.png", supp.fig2, width = 8, height = 6, dpi = 450)

```

### Supplementary Figure 3

```{r Supp Fig 3, fig.height=6, fig.width=9}

bolt.height.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/height.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)
bolt.bmi.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/BMI.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)
bolt.glucose.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/glucose.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)
bolt.hba1c.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/hba1c.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)
bolt.ldl.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/ldl.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)
bolt.hdl.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/hdl.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)
bolt.triglycerides.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/triglycerides.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)
bolt.whr_bmi.ret <- load.and.plot.data(file.names = c("ukbb_data/phewas/whr_bmi.bolt.genes.BOLT.stats.tsv.gz"),p.val.col="P_BOLT_LMM_INF",tool.name = "BOLT",AC.col = "AC",ymax = 75)

bolt.height.ret$gene.table[,pheno:="height"]
bolt.bmi.ret$gene.table[,pheno:="bmi"]
bolt.glucose.ret$gene.table[,pheno:="glucose"]
bolt.hba1c.ret$gene.table[,pheno:="hba1c"]
bolt.ldl.ret$gene.table[,pheno:="ldl"]
bolt.hdl.ret$gene.table[,pheno:="hdl"]
bolt.triglycerides.ret$gene.table[,pheno:="trigs"]
bolt.whr_bmi.ret$gene.table[,pheno:="whr_bmi"]

search.tables <- function(gene, maf, mask) {
  
  ret.table <- rbind(bolt.height.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask],
        bolt.bmi.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask],
        bolt.glucose.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask],
        bolt.hba1c.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask],
        bolt.ldl.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask],
        bolt.hdl.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask],
        bolt.triglycerides.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask],
        bolt.whr_bmi.ret$gene.table[SYMBOL == gene & MAF == maf & MASK == mask])
  
  return(ret.table)
}

phewas_table <- rbind(search.tables("GCK","MAF_01","HC_PTV")[,c("SYMBOL", "MASK","MAF","pheno","BETA","SE","p.value.selected","log.p")],
      search.tables("GIGYF1","MAF_01","HC_PTV")[,c("SYMBOL", "MASK","MAF","pheno","BETA","SE","p.value.selected","log.p")],
      search.tables("HNF1A","MAF_01","HC_PTV")[,c("SYMBOL", "MASK","MAF","pheno","BETA","SE","p.value.selected","log.p")],
      search.tables("IGF1R","MAF_01","MISS_REVEL0_7")[,c("SYMBOL", "MASK","MAF","pheno","BETA","SE","p.value.selected","log.p")],
      search.tables("MLXIPL","MAF_01","MISS_REVEL0_7")[,c("SYMBOL", "MASK","MAF","pheno","BETA","SE","p.value.selected","log.p")],
      search.tables("TNRC6B","AC_1","HC_PTV")[,c("SYMBOL", "MASK","MAF","pheno","BETA","SE","p.value.selected","log.p")],
      search.tables("ZEB2","AC_1","MISS_REVEL0_5")[,c("SYMBOL", "MASK","MAF","pheno","BETA","SE","p.value.selected","log.p")])

phewas_table[,mask.text:=if_else(MASK == "HC_PTV", "HC PTV", if_else(MASK == "MISS_REVEL0_5", "Missense REVEL ≥ 0.5", if_else(MASK == "MISS_REVEL0_7", "Missense REVEL ≥ 0.7", "NA")))]
phewas_table[,maf.text:=if_else(MAF == "MAF_01", "MAF < 0.1%", "Singletons")]
phewas_table[,axis.name:=paste0(SYMBOL, "\n(",mask.text, " ", maf.text, ")")]
phewas_table[,log.p:=if_else(log.p>-log10(p.val.thresh), -log10(p.val.thresh), log.p)]
phewas_table[,ci.lower:=BETA-(1.96*SE)]
phewas_table[,ci.upper:=BETA+(1.96*SE)]

supp.fig.3 <- ggplot(phewas_table, aes(axis.name, log.p, group = pheno, colour = pheno)) +
  geom_hline(yintercept=-log10(p.val.thresh),linetype=2, colour = "red") +
  geom_point(position=position_dodge(0.5)) +
  coord_flex_cart(bottom = brackets_horisontal(tick.length = 0)) +
  scale_color_discrete(labels = c("BMI", "c. Glucose", "c. HbA1c", "c. HDL", "Height", "c. LDL", "c. Trig.", "WHR adj. BMI"),guide=guide_legend(title="Phenotype")) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = expression(bold(-log[10](italic(p))))) +
  theme.legend + theme(panel.grid.major.x=element_blank())
supp.fig.3

ggsave("Figures/SuppFig3.png",height=6, width=9, dpi=450)

```

### Supplementary Figure 4

```{r Supp Fig 4, fig.height=5, fig.width=8}

IGF1_funnel <- fread("data_files/IGF1_height_funnel.txt")
IGF1_funnel[,P.Adult_height:=as.double(P.Adult_height)]
IGF1_funnel[,ci.upper:=Beta.Adult_height + (1.96 * SE.Adult_height)]
IGF1_funnel[,ci.lower:=Beta.Adult_height - (1.96 * SE.Adult_height)]

setkey(IGF1_funnel,"Beta.Adult_height")
IGF1_funnel[,SNP:=factor(SNP, levels=IGF1_funnel[,SNP])]
IGF1_funnel[,pt.colour:=if_else(P.Adult_height < 0.05 & Beta.Adult_height > 0, "blue",
                                if_else(P.Adult_height < 0.05 & Beta.Adult_height < 0, "red","black"))]
IGF1_funnel[,adj.beta:=if_else(Beta.Adult_height > 1, 1, Beta.Adult_height)]

supp.fig.4 <- ggplot(IGF1_funnel, aes(SNP, adj.beta,colour=pt.colour)) + 
  geom_point(size=0.5) +
  geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width = 0) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "β Adult Height", limits = c(-1,1)) +
  scale_colour_identity() +
  coord_capped_cart(left='both') +
  theme + theme(panel.grid.major.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())
  
supp.fig.4
ggsave("Figures/SuppFig4.png", supp.fig.4, height = 5, width = 8, dpi = 450)

```

## 5C. Supplementary Tables

Only Supp Table 1 is automated. Other tables are manually formatted in MS Excel.

### Supplementary Table 1

#### Define Initial Table

```{r Table 1}

supp.table1 <- bolt.T2D.ret$gene.table[AC > 30 & MASK %in% valid.masks & MASK != "SYN" & p.value.selected < p.val.thresh]
supp.table1[,c("staar.p","glm.p","logOR","std.err"):=query.data(SYMBOL, MASK, MAF),by=1:nrow(supp.table1)]
supp.table1 <- supp.table1[,c("ENST","SYMBOL","MASK","MAF","p.value.selected","staar.p","glm.p","AC","logOR","std.err")]

```

#### T2D carrier prop

```{r carrier prop}

# Get prop of cases per gene:
get_T2D_prop <- function(MASK, SYMBOL, MAF){

  carriers <- fread(paste0("ukbb_data/T2D_ExWAS/indv_gene_results/T2D.",MASK,"-",MAF,".", SYMBOL, ".carriers_formatted.tsv"))
  
  carrier_pheno_covar <- fread("ukbb_data/pheno_tables/T2D.phenotypes_covariates.formatted.tsv")
  carrier_pheno_covar[,has.var:=IID %in% carriers[,IID]]
  
  num <- nrow(carrier_pheno_covar[has.var == T & T2D_for_GWAS == 1])
  den <- nrow(carrier_pheno_covar[has.var == T & T2D_for_GWAS == 0])
  prop <- (num / (num + den)) * 100
  
  return(prop)
  
}

supp.table1[,"Prop. T2D":=get_T2D_prop(MASK, SYMBOL, MAF), by = 1:nrow(supp.table1)]
```

#### Number of Variants / Gene

```{r variants per gene}

# Get number of variants per gene:
get_var_count <- function(MASK, SYMBOL, MAF){

  gene_variants <- fread(paste0("ukbb_data/T2D_ExWAS/indv_gene_results/T2D.",MASK,"-",MAF,".", SYMBOL, ".variant_table.tsv"))
  gene_variants <- gene_variants[MAF_tested != 0]
  return(nrow(gene_variants))
  
}

supp.table1[,"Num. Variants":=get_var_count(MASK, SYMBOL, MAF), by = 1:nrow(supp.table1)]

```

#### Formatting Table

```{r Supp Table 1}

# Convert to human readable...
supp.table1[,MASK:=if_else(MASK == "HC_PTV", "HC PTV", if_else(MASK == "MISS_REVEL0_5", "Missense REVEL ≥ 0.5", if_else(MASK == "MISS_REVEL0_7", "Missense REVEL ≥ 0.7", "NA")))]
supp.table1[,MAF:=if_else(MAF == "MAF_01", "MAF < 0.1%", "Singletons")]

setnames(supp.table1,names(supp.table1),c("ENST", "Symbol", "Mask", "MAF Cutoff", "p. value BOLT", "p. value STAAR", "p. value Logistic Model", "Num. Carriers", "log OR", "standard err.", "Prop. T2D", "Num. Variants"))

fwrite(supp.table1, "Figures/SuppTab1.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

```

## 5D. Numbers Catalogue

```{r Numbers Catalogue}
paste0("p. value threshold: ", sprintf("%0.1e", p.val.thresh))

# Formatted N/p. value/ORs
T2D.table[,paste0(SYMBOL, " (N=",  n_car, "; OR ", sprintf("%0.1f", OR), " [", sprintf("%0.1f", ci.lower), "-", sprintf("%0.1f", ci.upper),"]; p=", sprintf("%0.1e", p.value.selected), ")")]

# IGF1R Phenos:
phewas_table[SYMBOL == "IGF1R" & pheno == "height",paste0(SYMBOL, " for ", pheno, sprintf(" %0.1f", BETA), " [", sprintf("%0.1f", ci.lower), "-", sprintf("%0.1f", ci.upper),"]; p=", sprintf("%0.1e", p.value.selected), ")")]

# TNRC6B carriers
TNRC6B.pheno <- fread("ukbb_data/pheno_tables/T2D.phenotypes_covariates.formatted.tsv")
TNRC6B.carriers <- fread("ukbb_data/T2D_ExWAS/indv_gene_results/HC_PTV.TNRC6B.T2D.carriers_formated.tsv")
TNRC6B.pheno[,has_var:=IID %in% TNRC6B.carriers[,IID]]
TNRC6B.inc <- sprintf("%0.1f", prop.table(table(TNRC6B.pheno[,c("has_var","T2D_for_GWAS")]),margin = 1) * 100)
paste0("Incidence of Diabetes in non-carriers: ", TNRC6B.inc[3])
paste0("Incidence of Diabetes in TNRC6B carriers: ", TNRC6B.inc[4])

# TNRC6B exon:
paste0("Number of individuals with PEXT exon variants in TNRC6B: ", length(TNRC6B.carriers[(POS > 40264529 & POS < 40267175), unique(IID)]))
TNRC6B.pheno[,has.exon.var:=if_else(IID %in% TNRC6B.carriers[(POS > 40264529 & POS < 40267175) == F, IID], 1, 0)]
TNRC6B.exclude_exon <- data.table(tidy(glm(T2D_for_GWAS ~ has.exon.var + sex + wes_batch + age + age_squared + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = TNRC6B.pheno, family = "binomial")))
paste0("Exclude variants in PEXT exon of TNRC6B: ", TNRC6B.exclude_exon[term == "has.exon.var", sprintf("%0.1e", p.value)])

# IGF1R protein domains
domain.table[,paste0(term, " (OR=", sprintf("%0.1f", OR), " [", sprintf("%0.1f", ci.lower), "-", sprintf("%0.1f", ci.upper),"]; p=", sprintf("%0.1e", p.value), ")")]
```
