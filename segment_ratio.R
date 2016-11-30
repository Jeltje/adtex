# ----------------------------------------------------------------------#
# Copyright (c) 2013, Kaushalya Amarasinghe.
#
# Note: this code has been altered slightly to fit in the Adtex docker container
# by Jeltje van Baren
#
# Original copyright notice:
#
# > Source License <
# This file is part of ADTEx.
#
#    ADTEx is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ADTEx is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ADTEx.  If not, see <http://www.gnu.org/licenses/>.
#
# 
#-----------------------------------------------------------------------#

require("wmtsa")
require("splus2R")

options <- commandArgs(trailingOnly = T)
sourceFunc = options[1]
control = options[2]
tumor = options[3]
thresh = as.numeric(options[4])
outputLoc = options[5]
bafIn = as.logical(options[6])
baf = options[7]
ploidyIn = as.logical(options[8])
outputRatio = options[9]
outputSnpSegs = options[10]
chrom=unlist(strsplit(options[11],","))

options(scipen=999) # removing scientific notation. fixed notation will be preferred
                    #unless it is more than scipen digits wider
source(sourceFunc)

control<-read.delim(control,header=F)
tumor<-read.delim(tumor,header=F)

ch<-sapply(chrom,function(x)(ifelse(sum(control[,1]==x & control[,4]>thresh)>100,1,0)))
ch<-ch[sapply(1:length(ch),function(x)(which(chrom[x]==names(ch))))]
chrom<-names(ch)[ch>0]
ch<-sapply(chrom,function(x)(ifelse((
  sum(tumor[,1]==x & tumor[,4]<thresh)/sum(tumor[,1]==x))==1,0,1)))
ch<-ch[sapply(1:length(ch),function(x)(which(chrom[x]==names(ch))))]
chrom<-names(ch)[ch>0]

result<-matrix(data=NA,ncol=8,nrow=sum(control[,4]>thresh & control[,1]%in%chrom))
colnames(result)<-c("chr","exon_start","exon_end","tumor_DOC","control_DOC",
                    "ratio_before_smoothing","ratio_after_smoothing",
                    "rationormalized_after_smoothing")

result<-data.frame(result)
f<-(control[,4]>thresh & control[,1]%in%chrom)
result[,1:3]<-control[f,1:3]
result[,5]<-control[f,4]
result[,4]<-tumor[f,4]
control<-control[f,]
tumor<-tumor[f,]
control[,5]<-control[,4]/mean(control[,4])
tumor[,5]<-tumor[,4]/mean(tumor[,4])
ratio<-tumor[,5]/control[,5]
df<-data.frame(cbind(chr=result[,1],control=control[,4],ratio=ratio))
df$chr<-result[,1]
df<-smoothData(df,chrom)
result$ratio_before_smoothing<-df$ratio
peak<-selectPeak(df$ratio_after_smoothing)
result$rationormalized_after_smoothing<-df$ratio_after_smoothing/peak
result$ratio_after_smoothing<-df$ratio_after_smoothing
write.table(result,file=outputRatio,quote=F,sep="\t",
            row.names=F)

if(bafIn & !ploidyIn){
  require("DNAcopy")
    
  baf<-read.delim(baf)
  snp<-data.frame(baf[,1])
  snp[,2]<-baf$SNP_loc-1
  snp[,3]<-baf$SNP_loc
  snp[,4]<-baf$tumor_BAF
  snp[,5]<-baf$control_BAF
  
  s<-snp[,4]
  med_o<-median(s)
  s[s>0.5]<-1-s[s>0.5]
  med_m<-median(s)
  if(med_o>0.42 & med_m>0.4){
    s2<-sapply(s,function(x)(min(0.5,((med_o-med_m)*x/med_o)+x)))
  }else{
    s2<-s
  }
  CNA.object<-CNA(s2,snp[,1],snp[,2],data.type="logratio",sampleid="tumour")
  smoothed.CNA.object<-smooth.CNA(CNA.object)
  segment.data<-segment(smoothed.CNA.object)
  write.table(segment.data$output[,2:6],file=outputSnpSegs,
              quote=F,sep="\t",row.names=F,col.names=F)
}


