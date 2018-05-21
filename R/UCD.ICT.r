# three ICT datasets
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

library(ROCR)
library(caTools)
library(data.table)

qnorm.array <- function(mat)     ############# correct one
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}
rank.array <- function(mat)     ############# correct one
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average")/length(mat);
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}
setwd("/cbcb/project2-scratch/jooslee/Ayelet/pyrimidine/github/")

# (1) Van Allen dataset
load("data/Van.Allen.RData")
mut=prob$mut
patient=prob$patient
mRNA=prob$mRNA

# rank-normalization
mRNA.rank = lapply(1:ncol(mRNA), function(tt) rank.array(mRNA[,tt]))
mRNA.rank = t(do.call(rbind, mRNA.rank))

# UCD-score
UC.genes=c("ASL","ASS1","CPS1","OTC","SLC25A13","SLC25A15")
iuc=match(UC.genes,prob$genes)
q=c(-1,-1,1,-1,1,-1)
score.M=NULL
for (i in seq(length(UC.genes))){
	score.M=rbind(score.M,prob$mRNA.rank[iuc[i],]*q[i])
}
u.score=apply(score.M,2,sum,na.rm=T)
patient$u.score=u.score

icad=match("CAD",prob$genes)
cad=patient$cad=prob$mRNA.rank[icad,]

ratio=matrix(NA,4,length(prob$samples))
for (i in seq(length(prob$samples))){
	ix=which(mut$patient==prob$samples[i])
	bb=mut[ix,]
	ratio[1,i]= sum(bb$mja %in% c("A","G") & bb$mna %in% c("A","G"),na.rm=T)
	ratio[2,i]=	sum(bb$mja %in% c("A","G") & bb$mna %in% c("C","T"),na.rm=T)
	ratio[3,i]=	sum(bb$mja %in% c("C","T") & bb$mna %in% c("A","G"),na.rm=T)
	ratio[4,i]=	sum(bb$mja %in% c("C","T") & bb$mna %in% c("C","T"),na.rm=T)
	print(i)
}
patient$mut.bias=(ratio[2,]-ratio[3,])/apply(ratio,2,sum,na.rm=T)
patient$mut.load=apply(ratio,2,sum,na.rm=T)

pval=cbind(patient$u.score,patient$mut.load,patient$mut.bias)	
flag=1*(patient$RECIST %in% c("PR","CR"))
	ix=which(patient$RECIST %in% c("PR","CR","PD"))
	flag=flag[ix]
	pval=pval[ix,]

auc=prr=rep(NA,dim(pval)[2])

for (i in seq(dim(pval)[2])){
	pred <- prediction( pval[,i], flag)
	perf <- performance(pred,"tpr","fpr")
	auc.perf = performance(pred, "tpr","fpr",measure = "auc")
	auc[i]=auc.perf@y.values[[1]]
	perf1 <- performance(pred, "prec", "rec")
	rec=perf1@x.values[[1]]
	prec=perf1@y.values[[1]]
	prr[i] <- trapz(rec[2:length(rec)], prec[2:length(prec)])
}

wilcox.test(pval[flag==1,1],pval[flag==0,1],alternative="greater") #UCD: 0.04095
wilcox.test(pval[flag==1,2],pval[flag==0,2],alternative="greater") #mutational load: 0.1247
wilcox.test(pval[flag==1,3],pval[flag==0,3],alternative="greater") #PTMB: 0.5088


# (2) Hugo dataset
patient2=fread("data/GSE78220.phenotype.txt",header=T)
genes=fread("data/GSE78220.genes.txt",header=F)$V1
mRNA=fread("data/GSE78220_PatientFPKM.csv",header=T)

samples0=do.call(rbind,strsplit(colnames(mRNA),"[.]"))[,1]
samples=intersect(samples0,patient2$"Patient ID")
ix=match(samples,patient2$"Patient ID")
patient2=patient2[ix,]
iy=match(samples,samples0)
mRNA=mRNA[,iy,with=F]
mRNA=as.matrix(mRNA);class(mRNA) <- "numeric"
mRNA.rank = lapply(1:ncol(mRNA), function(tt) rank.array(mRNA[,tt]))
mRNA.rank = t(do.call(rbind, mRNA.rank))

UC.genes=c("ASL","ASS1","CPS1","OTC","SLC25A13","SLC25A15")
iuc=match(UC.genes,genes)
q2=c(-1,-1,1,-1,1,-1)

score.M=NULL
for (i in seq(length(UC.genes))){
	score.M=rbind(score.M,mRNA.rank[iuc[i],]*q2[i])
}
u.score=apply(score.M,2,sum,na.rm=T)

ir2=which(patient2$irRECIST %in% c("Complete Response","Partial Response"))
ir0=which(patient2$irRECIST %in% c("Progressive Disease"))
wilcox.test(u.score[ir2],u.score[ir0],alternative="greater") #0.0284

load("data/Hugo.mutation.RData")
pval=cbind(mut$bias,mut$ml,mut$score)
flag=(mut$response=="response")*1

na.inx=which(!is.na(pval[,1]))
	pval=pval[na.inx,]
	flag=flag[na.inx]

auc=prr=rep(NA,dim(pval)[2])
for (i in seq(dim(pval)[2])){
	pred <- prediction( pval[,i], flag)
	perf <- performance(pred,"tpr","fpr")
	auc.perf = performance(pred, "tpr","fpr",measure = "auc")
	auc[i]=auc.perf@y.values[[1]]
	perf1 <- performance(pred, "prec", "rec")
	rec=perf1@x.values[[1]]
	prec=perf1@y.values[[1]]
	prr[i] <- trapz(rec[2:length(rec)], prec[2:length(prec)])
}

# (3) Roh et al (Wargo) data
mut=fread("data/Roh.mutation.call.csv")
upatient=unique(unlist(mut[,1,with=F]))
samples=fread("data/Roh.samples.csv",header=T)
ix=match(upatient,samples$sample)
samples=samples[ix,]
patient2=samples

mut$mja=mut$ref_allele
mut$mna=mut$alt_allele

na.inx=which(mut$aaannotation!="")
na.pc=mut$aaannotation[na.inx]
pc=do.call(rbind,strsplit(na.pc,"[.]"))[,2]
mut$p.mja=mut$p.mna=rep(NA,nrow(mut))
mut$p.mja[na.inx]=substr(pc,1,1)
mut$p.mna[na.inx]=substr(pc,nchar(pc),nchar(pc))
mut$non.syn=(mut$p.mja!=mut$p.mna)*1
mut=mut[mut$non.syn==1,]

ratio=matrix(NA,4,length(upatient))
for (i in seq(length(upatient))){
	ix=which(unlist(mut[,1,with=F])==upatient[i])
	bb=mut[ix,]
	ratio[1,i]= sum(bb$mja %in% c("A","G") & bb$mna %in% c("A","G"),na.rm=T)
	ratio[2,i]=	sum(bb$mja %in% c("A","G") & bb$mna %in% c("C","T"),na.rm=T)
	ratio[3,i]=	sum(bb$mja %in% c("C","T") & bb$mna %in% c("A","G"),na.rm=T)
	ratio[4,i]=	sum(bb$mja %in% c("C","T") & bb$mna %in% c("C","T"),na.rm=T)
	print(i)
}
patient2$mut.bias=(ratio[2,]-ratio[3,])/apply(ratio,2,sum,na.rm=T)
patient2$mut.load=apply(ratio,2,sum,na.rm=T)

pval=cbind(patient2$mut.load,patient2$mut.bias)
flag=1*(samples$aPD1_response=="response")
#	pval=pval[!is.na(samples$aPD1_response) & samples$time=="postCTLA4_prePD1",]
#	flag=flag[!is.na(samples$aPD1_response) & samples$time=="postCTLA4_prePD1"]
	pval=pval[!is.na(samples$aPD1_response),]
	flag=flag[!is.na(samples$aPD1_response)]

pval=do.call(cbind,lapply(seq(dim(pval)[2]),function(x) {rank(pval[,x],na.last="keep")/length(pval[,x])}))
auc=prr=rep(NA,dim(pval)[2])

for (i in seq(dim(pval)[2])){
	pred <- prediction( pval[,i], flag)
	perf <- performance(pred,"tpr","fpr")
	auc.perf = performance(pred, "tpr","fpr",measure = "auc")
	auc[i]=auc.perf@y.values[[1]]
	perf1 <- performance(pred, "prec", "rec")
	rec=perf1@x.values[[1]]
	prec=perf1@y.values[[1]]
	prr[i] <- trapz(rec[2:length(rec)], prec[2:length(prec)])
}
