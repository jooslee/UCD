# peptidomics analysis
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

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

library(data.table)
library(parallel)
library(Peptides)

setwd("/cbcb/project2-scratch/jooslee/Ayelet/pyrimidine/github/")
hep=fread("data/Hepg2_shOTC.pairedVCF.peptides.txt")
lox=fread("data/Lox_citrin_OE.pairedVCF.peptides.txt")
u2o=fread("data/U2os_shASS.pairedVCF.peptides.txt")

hep$ucd=hep$"Intensity Hepg2- shOTC"/sum(hep$"Intensity")
hep$ctr=hep$"Intensity Hepg2- Ev"/sum(hep$"Intensity")
hep$ucd.count="Experiment Hepg2- shOTC"
hep$ctr.count="Experiment Hepg2- Ev"

lox$ucd=lox$"Intensity Lox- citrin OE"/sum(lox$Intensity)
lox$ctr=lox$"Intensity Lox-Ev"/sum(lox$Intensity)
lox$ucd.count=lox$"Experiment Lox- citrin OE"
lox$ctr.count=lox$"Experiment Lox-Ev"

u2o$ucd=u2o$"Intensity U2os- shASS"/sum(u2o$Intensity)
u2o$ctr=u2o$"Intensity U2os- Ev"/sum(u2o$Intensity)
u2o$ucd.count=u2o$"Experiment U2os- shASS"
u2o$ctr.count=u2o$"Experiment U2os- Ev"

unp1=hep$ctr[!is.na(hep$ucd.count) & !is.na(hep$ctr.count)]
unp2=lox$ctr[!is.na(lox$ucd.count) & !is.na(lox$ctr.count)]
unp3=u2o$ctr[!is.na(u2o$ucd.count) & !is.na(u2o$ctr.count)]
per1=hep$ucd[!is.na(hep$ucd.count) & !is.na(hep$ctr.count)]
per2=lox$ucd[!is.na(lox$ucd.count) & !is.na(lox$ctr.count)]
per3=u2o$ucd[!is.na(u2o$ucd.count) & !is.na(u2o$ctr.count)]

dat=data.frame(x=c(unp1,unp2,unp3,per1,per2,per3),
			cond=c(rep("UC-WT",length(c(unp1,unp2,unp3))),rep("UCD",length(c(per1,per2,per3)))))
dat=dat[!is.na(dat$x),]
dat$x=log10(dat$x)
dat=dat[is.finite(dat$x),]
wilcox.test(dat$x[dat$cond=="UCD"],dat$x[dat$cond=="UC-WT"],alternative="greater") #0.001001

# hydrophobicity
hp.score2 <- function(peptides){
	pp=toupper(peptides)
	return(hydrophobicity(pp,"Janin"))
}

hep$hp.score0=do.call(c,lapply(hep$Sequence,hp.score))
lox$hp.score0=do.call(c,lapply(lox$Sequence,hp.score))
u2o$hp.score0=do.call(c,lapply(u2o$Sequence,hp.score))

hep$hp.score=do.call(c,lapply(hep$Sequence,hp.score2))
lox$hp.score=do.call(c,lapply(lox$Sequence,hp.score2))
u2o$hp.score=do.call(c,lapply(u2o$Sequence,hp.score2))

per1=(hep$hp.score*hep$ucd)[!is.na(hep$ucd.count)]
unp1=(hep$hp.score*hep$ctr)[!is.na(hep$ctr.count)]
per2=(lox$hp.score*lox$ucd)[!is.na(lox$ucd.count)]
unp2=(lox$hp.score*lox$ctr)[!is.na(lox$ctr.count)]
per3=(u2o$hp.score*u2o$ucd)[!is.na(u2o$ucd.count)]
unp3=(u2o$hp.score*u2o$ctr)[!is.na(u2o$ctr.count)]

per1=hep$hp.score[hep$ucd!=0]
unp1=hep$hp.score[hep$ctr!=0]
per2=lox$hp.score[lox$ucd!=0]
unp2=lox$hp.score[lox$ctr!=0]
per3=u2o$hp.score[u2o$ucd!=0]
unp3=u2o$hp.score[u2o$ctr!=0]

dat=data.frame(x=c(unp1,unp2,unp3,per1,per2,per3),
			cond=c(rep("UC-WT",length(c(unp1,unp2,unp3))),rep("UCD",length(c(per1,per2,per3)))))
dat=dat[!is.na(dat$x),]
wilcox.test(dat$x[dat$cond=="UCD"],dat$x[dat$cond=="UC-WT"],alternative="greater")
#0.0006999

# neo-antigen
diff0=c(per1-unp1,per2-unp2,per3-unp3)

AA=fread("data/neoantigen.txt",header=F)
hp.score3 <- function(peptides){
	pp=toupper(peptides)
	return(hydrophobicity(pp,"Janin"))
}

hps2=do.call(c,lapply(AA$V2,hp.score3))
hps3=do.call(c,lapply(AA$V3,hp.score3))
wilcox.test(hps2[c(1:15)],hps3[c(1:15)],paired=T,alternative="less")$p.value #0.03456403

