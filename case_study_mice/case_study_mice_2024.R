#' ---
#' title: "Identifying cis-driver genes for insulin and islet-specific eQTL hotspot on chr2 of  Mus musculus"
#' author: "Qi Zhang"
#' date: "Feb 18, 2024"
#' ---
#'

#'
#' In this demo, we analyze the mouse f2 cross data for diabetes study from Dr. Alan Attie's lab at University of Wisconsin Madison. There are many publications involved (e.g., Wang et al 2011, Tu et al 2012, and Tian et al 2015). The data analyzed in this demo is provided in folder data_mouse.
#'
#' We have mentioned that there is an islet specific eQTL hotspot and a QTL peak for insulin level colocalized on chr2. We decide to focus on a region from 120Mbps to 170Mbps on chr2, which include the center of the eQTL hotspot and the QTL peak based on the previous literature. We define an eQTL as cis- if the distance between the transcription starting site of the gene being regulated and the genetic marker is less than 10Mbps. This leads to 193 cis-regulated genes as exposures, and 4561 trans-regulated genes as the mediators. There are 2057 genetic markers across the whole genome, and we used them as confounders to remove the effect of the population structure.
#'
#' The dataset includes three independent measures of the blood insulin level at sacrifice (named ``10wk",``10wkRep", and ``rbm", respectively), which provides us a unique opportunity to evaluate statistical methods based the reproducibility of the analysis results across independent measurements of identical biological signal.  We conduct three high dimensional mediation analysis using the same confounders, exposures and mediators, and each of the three measurements of insulin as the outcome. We evaluate the performance of the proposed method MedDiC based on the reproducibility across the three analysis, and the overlap between the significant exposures and the known transcription factors.
#'
#' ## 1. Data preprocessing
#'
#'
#' #### 1.1 Creating the merged list of known transcription factors for evaluation.
#'
#' We first need to construct a merged list of known transcription factors based on three diverse sources.
#'
#' The three lists of of known TFs for mammals are downloaded from their publications, and are in folder data_mouse. We will use gene symbols for matching.
fantom5 = read.csv("fantom5_Transcription_Factors_mm9.csv",header = T)
gifford = read.table("mouse_ensemble_tfs_from_lambertetal_isyes_unique_gifford.txt",header = F)
animaltfdb = read.table("Mus_musculus_TF_animaltfdb3.txt",header = T,sep='\t')
animaltfdb[1,]
fantom5[1,]
gifford[1,]
#' The dataset from Gifford lab did not come with gene symbols. biomaRt package was used to retrieve the gene symbols.
#' The users will need to install R pckage biomaRt from Bioconductor.
#'
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")
#require(biomaRt)
#mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
#names(mouse)
#dtmouse = useDataset(dataset = 'mmusculus_gene_ensembl',mart = mouse)
#dtgifford = getBM(attributes =  c('mgi_symbol','ensembl_gene_id','gene_biotype'),filters = 'ensembl_gene_id',values = gifford$V1,mart = mouse)
#dim(dtgifford)
#dtgifford[1,]
#' Once we have the gene symbols from all three sources, we can pool them. There are 1929 unique TFs
#tf_known = setdiff(union(union(tolower(fantom5$Symbol),tolower(animaltfdb$Symbol)),tolower(dtgifford$mgi_symbol)),"")
#length(tf_known)
#save(tf_known,file='tf_known.rdata')
#' We provide this file for the user's convenience.

#'
#'## 1.2 Creating the datasets for the three mediation analysis.
#'
#' The locus 120Mbs-170Mbs on chr2 is selected based on the QTL peak and eQTL hot spot reported in the literature.
loci2 = c(120,170)
#' The eQTL results were downloaded with lod score cutoff 3 from http://diabetes.wisc.edu/f2.php on 06/06/2018. It was also provided in folder data_mouse
loci.eqtl = read.csv('eQTL_islet3.csv',sep='\t')
#' We further filter the eQTL to only chr2 markers, and LOD>=5
loci2.eqtl = subset(loci.eqtl,peak.chr==2&peak.pos.Mb>loci2[1]&peak.pos.Mb<loci2[2]&peak.score>=5)
chr2.eqtl = subset(loci2.eqtl,trans.chr==2)
dim(chr2.eqtl)
#' Next, we separate cis and trans eqtls.
dis.chr2 = with(chr2.eqtl,abs(peak.pos.Mb-trans.pos.Mb))
#' The cis-regulated genes are genes on chr2 whose distance from the associated markers is less than 10Mbs. The other genes are labeled as trans-regulated genes. We will use the probe IDs of these genes to subset the transcriptomics data.
cis.eqtl = subset(loci2.eqtl,trans.chr==2&(abs(peak.pos.Mb-trans.pos.Mb)<10))
dim(cis.eqtl)
cis.eqtl[1,]
trans.eqtl = subset(loci2.eqtl,!(trans.chr==2&(abs(peak.pos.Mb-trans.pos.Mb)<10)))
dim(trans.eqtl)
trans.eqtl[1,]
save(trans.eqtl,cis.eqtl,file='locus_info.rdata')

#' We read the genotype marker data.
gta = read.csv('genotype.csv',header=T,row.names = 1)
#' We read the gene expression data from a .csv file.
expr = read.csv('islet_expression.csv',header=T,row.names = 1)
expr[1:3,1:5]
colnames(expr) = gsub('X','',colnames(expr))
#' We read the three measures of insulin levels.
insulin_10wk = read.csv('insulin_10wk.csv',header=T,row.names = 1)
insulin_10wkRep = read.csv('insulin_10wkRep.csv',header=T,row.names = 1)
insulin_rbm = read.csv('insulin_rbm.csv',header=T,row.names = 1)
#' We will use sex of the mice as an additional confounder.
sex=read.csv('sex.csv',header = T,row.names = 1)
#' The numbers of mouse in different data types are different, and they can be matched using mouse ID (the row names of the data). There are also some missing values in each data type. We will analyze the subjects that the phenotype is not missing.
dim(insulin_10wk)
dim(insulin_10wkRep)
dim(insulin_rbm)
dim(gta)
dim(expr)
gta[1:3,1:5]
expr[1:3,1:5]
head(insulin_10wk)
mean(is.na(gta))
mean(is.na(expr))

#' The following is a utility function that create the dataset for the analysis of each phenotype. Since there are not many missing values in genotype and gene expressions, we impute using the mean of the non-missing values for each gene.

imputeMean <- function(x){
  x[is.na(x)] = mean(x,na.rm=T)
  return(x)
}

getMedData <- function(pheno_name,pheno,gta,expr,sex,cis,trans){
  y=pheno[,1]
  names(y) = rownames(pheno)
  y = y[!is.na(y)]
  mouse.common = intersect(names(y),intersect(rownames(gta),rownames(expr)))
  y=y[mouse.common]
  sex=sex[mouse.common,1]
  gta=apply(gta,2,imputeMean)
  gm = gta[mouse.common,]
  Z=expr[mouse.common,cis]
  M=expr[mouse.common,trans]
  Z = apply(Z,2, imputeMean)
  M = apply(M,2, imputeMean)
  dt4med = list(y=as.matrix(y),X=cbind(sex,gm),Z=Z,M=M)
  save(dt4med,file=paste0(pheno_name,'_input.rdata'))
  return(sapply(dt4med,dim))
}

getMedData('10wk',insulin_10wk,gta,expr,sex,cis=as.character(cis.eqtl$a.gene.id),trans = as.character(trans.eqtl$a.gene.id))


getMedData('10wkRep',insulin_10wkRep,gta,expr,sex,cis=as.character(cis.eqtl$a.gene.id),trans = as.character(trans.eqtl$a.gene.id))


getMedData('rbm',insulin_rbm,gta,expr,sex,cis=as.character(cis.eqtl$a.gene.id),trans = as.character(trans.eqtl$a.gene.id))

#'
#'## 2. Data analysis
#'
#' Next, we start the real data analysis for the three related phenotypes.
#'

library(Rcpp)
library(RcppArmadillo)
sourceCpp('..\\meddic.cpp')
source('..\\meddic.R')


#' This is a utility function that converts the matrix of the adjusted p-value to a tag for the effect type. e.g., an exposure with significant indirect effect and total effect but not direct effect has a tag "tot_indirect"
padj2tag =function(padj,pthr=0.1){
  ptot = padj$totPadj
  pdir = padj$directPadj
  pind = padj$indirectPadj
  out = rep('',length(ptot))
  out = ifelse(ptot<=pthr,paste(out,'tot',sep='_'),out)
  out = ifelse(pdir<=pthr,paste(out,'direct',sep='_'),out)
  out = ifelse(pind<=pthr,paste(out,'indirect',sep='_'),out)
  return(out)
}

ci2center = matrix(c(0.5,0.5,0,0,0,0,0,0,0.5,0.5,0,0,0,0,0,0,0.5,0.5),6,3)

padj.thr=0.1

load('locus_info.rdata')

#' The three mediation analysis in a for loop
for(pheno in c('10wk','10wkRep','rbm')){
load(paste0(pheno,'_input.rdata'))

# Calculate the score matrix for the marginal outcome model (without mediators), and save it in a file.
X = unique(dt4med$X,MARGIN = 2)
dim(X)

(fscore = paste0(pheno,'_score.rdata'))
# (1) prepare the data for meddic analysis
dt.meddic <- meddic_data_prep(dt4med$y,X,dt4med$Z,dt4med$M,thr.lse=0.2,no.penalty.exposure=0,no.penalty.confounder=0)
# (2) calculate the scores
if(!file.exists(fscore)){
  score.meddic <- meddic_get_score(dt.meddic)
  # The score calculation for the debiased lasso is the most time consuming step, especially when p and q are large. We recommend to   save the results in a file.
  save(score.meddic,file=fscore)
}else{
  load(fscore)
}
# (3)Estimating the direct, indirect and total effects
resmeddic <- meddic_get_results(dt.meddic,score.meddic)
# (4) Perform inference: calculating the raw p-values without adjusting multiple testing and 95% CIs.
inference_out <- inference.elem(resmeddic,conf.level = 0.95)
resmeddic$results_inference_alpha05 <- inference_out
save(resmeddic, file=paste0(pheno,'_meddic.rdata'))


# Create the data.frame for the significant exposures
pval.scl.asymp = resmeddic$results_inference_alpha05$pval.asymp
ci.scl.asymp = resmeddic$results_inference_alpha05$ci.asymp
est.scl.asymp = as.matrix(ci.scl.asymp)%*%ci2center
padj.scl.asymp = data.frame(sapply(pval.scl.asymp,p.adjust,method='fdr'))
colnames(pval.scl.asymp) = c('totPraw','directPraw','indirectPraw')
colnames(est.scl.asymp) = c('totEst','directEst','indirectEst')
colnames(padj.scl.asymp) = c('totPadj','directPadj','indirectPadj')
idsig = c(1,which(apply(padj.scl.asymp,1,min)<=padj.thr))
results.scl = cbind(est.scl.asymp,pval.scl.asymp,padj.scl.asymp)
if(length(idsig)>1){
  idtag = padj2tag(padj.scl.asymp[idsig,],pthr=padj.thr)
  loc.meddic= cbind(idsig,idtag,est.scl.asymp[idsig,],padj.scl.asymp[idsig,],pval.scl.asymp[idsig,],ci.scl.asymp[idsig,])[-1,]
  loc.meddic = data.frame(geneid=cis.eqtl$a.gene.id[loc.meddic[,'idsig']],genesymbol=cis.eqtl$symbol[loc.meddic[,'idsig']],loc.meddic)
}

write.csv(loc.meddic,row.names = F,quote = F,file=paste0(pheno,'_meddic_loci.csv'))

}

#'
#'## 3 Evaluations
#'
#'
#'
#'#### 3.1 Evaluation based on the reproducibility across the three analysis
#'

pheno.all = c('10wk','10wkRep','rbm')

#' This is a function that extracts the Gene IDs for each effect type
getGeneID = function(pheno){
    loc.meddic=read.csv(paste0(pheno,'_meddic_loci.csv'),header=T)

    output = list(indirect=subset(loc.meddic,grepl('_indirect',idtag))$geneid,direct=subset(loc.meddic,grepl('_direct',idtag))$geneid,tot=subset(loc.meddic,grepl('_tot',idtag))$geneid)
    return(output)
  }

geneID.all = lapply(pheno.all, getGeneID)
names(geneID.all)=pheno.all

#' This is the function that performs the overlaps between each pair of phenotypes.
getPairwise = function(nm){
    id.all = lapply(geneID.all,getElement,name=nm)
    ng = sapply(id.all,length)
    ins12 = intersect(id.all[[1]],id.all[[2]])
    ins13 = intersect(id.all[[1]],id.all[[3]])
    ins23 = intersect(id.all[[2]],id.all[[3]])
    output = c(ng,length(ins12),length(ins13),length(ins23))
    names(output) = c(pheno.all,paste(pheno.all[1],pheno.all[2],sep='_vs_'),paste(pheno.all[1],pheno.all[3],sep='_vs_'),paste(pheno.all[2],pheno.all[3],sep='_vs_'))
    return(output)
  }
#' Here is the overlap results.
t(sapply(names(geneID.all[[1]]),getPairwise))

#'
#'#### 3.2 Evaluation based on their overlaps with the known transcription factors.
#'
#' Now start comparison by checking the proportion of transcription factors.
load('tf_known.rdata')
#' We intersect the model outputs with this list to find the overlaps. We expect the genes with indirect effects include more TFs, and not so much for the genes with no indirect effects but only direct or total effects.

for(pheno in pheno.all){
  print(pheno)
  loc.meddic=read.csv(paste0(pheno,'_meddic_loci.csv'),header=T)

print(intersect(tf_known,tolower(with(loc.meddic,genesymbol[grepl('indirect',idtag)]))))

print(intersect(tf_known,tolower(with(loc.meddic,genesymbol[!grepl('indirect',idtag)]))))

}



