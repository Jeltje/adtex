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

require("DNAcopy",quietly=T)

options <- commandArgs(trailingOnly = T)
sourceFunc = options[1]
ratio_data = options[2]
ploidy = as.numeric(options[3])
thresh = as.numeric(options[4])
bafin=options[5]
baf = options[6]
plot = options[7]
chrom=unlist(strsplit(options[8],","))
#totReads=options[9]
#print(c(control,tumor,ploidy,thresh,outputLoc,plot,chrom))

write(paste("Analysing tumour sample with base ploidy : ",as.character(ploidy), sep=""), stderr())
source(sourceFunc)
ratio_data<-read.delim(ratio_data)
result<-matrix(data=NA,ncol=12,nrow=nrow(ratio_data))
colnames(result)<-c("chr","exon_start","exon_end","cnv","tumor_DOC","control_DOC",
                    "ratio_before_smoothing","ratio_after_smoothing",
                    "rationormalized_after_smoothing","CNV_start",
                    "CNV_end","seg_mean")
result<-data.frame(result)
result[,1:3]<-ratio_data[,1:3]
result[,5:9]<-ratio_data[,4:8]

for (i in chrom){
  #print(paste("Analysing chromosome ",i," :",sep=""))
  f<-result$chr==i
  if(i==chrom[1]){
    fit.cnv<-predict.cnv(result$rationormalized_after_smoothing[f],d=ploidy,
                         models=list(model1=c(0:4)))
    mu<-fit.cnv$mu
    model<-fit.cnv$model
    id<-which(model==ploidy)
    model<-c(ploidy-1,ploidy,ploidy+1)
    mu<-mu[(id-1):(id+1)]
  }
  else{
    fit.cnv<-predict.cnv(result$rationormalized_after_smoothing[f],d=ploidy,
                         predict=list(mu=mu,model=model),
                         models=list(model1=c(0:4)))
    mu1<-fit.cnv$mu
    model1<-fit.cnv$model
    id<-which(model1==ploidy)
    model<-c(ploidy-1,ploidy,ploidy+1)
    mu<-(mu1[(id-1):(id+1)]+mu)/2
    #print(mu,model)
  }
  result$cnv[f]<-fit.cnv$cnv
  #result$ratio_before_smoothing[f]<-tumour/normal
  #result$ratio_after_smoothing[f]<-data$ratio*mean(tumour/normal)/mean(data$ratio)
  segments<-combine.segments(result$rationormalized_after_smoothing[f],fit.cnv$cnv)
  result$seg_mean[f]<-segments[,1]
  result$CNV_start[f]<-result$exon_start[f][c(segments[,2])]
  result$CNV_end[f]<-result$exon_end[f][c(segments[,3])]
  result$cnv[f]<-segments[,4]
}

state.mean<-tapply(result$rationormalized_after_smoothing,result$cnv,mean)
temp<-tapply(result$rationormalized_after_smoothing,result$cnv,sd)
temp<-temp[!is.na(temp)&!temp>1]
sd.thresh<-ifelse(length(temp)<3,mean(temp),mean(temp[1:3]))


for(i in chrom){
  f<-result$chr==i
  path<-result$cnv[f]
  obs<-result$rationormalized_after_smoothing[f]
  changes<-which(path[-length(path)]!=path[-1])+1
  if(length(changes)>1){
    temp<-F
    if(length(changes)>10 & length(table(path))<4){
      CNA.object<-CNA(obs,result$chr[f],result$exon_start[f],
                      data.type="logratio",sampleid="tumor")
      smoothed.CNA.object<-smooth.CNA(CNA.object)
      segment.data<-segment(smoothed.CNA.object,verbose=0,undo.splits="sdundo",
                            undo.SD=1)
      if(nrow(segment.data$output)==1){
        result$seg_mean[f]<-mean(obs)
        result$CNV_start[f]<-result$exon_start[f][1]
        result$CNV_end[f]<-result$exon_end[f][length(obs)]
        k<-as.numeric(names(which.min(abs(mean(obs)-state.mean))))
        result$cnv[f]<-k
        temp<-T
      }
    }
    if(!temp){
      range.start<-c(1,changes)
      range.end<-c(changes-1,length(obs))
      range.mean<-sapply(1:length(range.end),
                         function(x)(mean(obs[range.start[x]:range.end[x]])))
      range.cnv<-path[range.start]
      mean.diff<-abs(range.mean-state.mean[as.character(range.cnv)])
      for(i in 1:length(range.cnv)){
        if(mean.diff[i]>sd.thresh){
          k<-as.numeric(names(which.min(abs(range.mean[i]-state.mean))))
          path[range.start[i]:range.end[i]]<-as.numeric(k)
        }
      }
      segments<-cnv.segments(obs,path)
      result$seg_mean[f]<-segments[,1]
      result$CNV_start[f]<-result$exon_start[f][c(segments[,2])]
      result$CNV_end[f]<-result$exon_end[f][c(segments[,3])]
      result$cnv[f]<-path
    }
  }
}

if(bafin){
  baf<-read.delim(baf)
  state.mean<-tapply(result$rationormalized_after_smoothing,result$cnv,mean)
  baf$mirrored_BAF<-baf$tumor_BAF
  baf$mirrored_BAF[baf$tumor_BAF>0.5]<-1-baf$tumor_BAF[baf$tumor_BAF>0.5]
  
  for(i in chrom){
    f<-result$chr==i
    path<-result$cnv[f]
    obs<-result$rationormalized_after_smoothing[f]
    changes<-which(path[-length(path)]!=path[-1])+1
    if(length(changes)>10 & length(table(path))<4 & sum(path!=ploidy)>100){
      CNA.object<-CNA(baf$mirrored_BAF[baf$chrom==i],baf$chrom[baf$chrom==i],
                      baf$SNP_loc[baf$chrom==i],
                      data.type="logratio",sampleid="tumor")
      smoothed.CNA.object<-smooth.CNA(CNA.object)
      segment.data<-segment(smoothed.CNA.object,verbose=0,undo.splits="sdundo",
                            undo.SD=1)
      if(nrow(segment.data$output)<6){
        for(j in 1:nrow(segment.data$output)){
          kk<-which(segment.data$output$loc.start[j]<=result$exon_start[f] & segment.data$output$loc.end[j]>=result$exon_end[f])
          seg.m<-mean(obs[kk])
          if(length(kk)>0){
            k<-as.numeric(names(which.min(abs(seg.m-state.mean))))
            if(k==1 & segment.data$output$seg.mean[j]>0.3){
              k<-as.numeric(names(which.min(abs(seg.m-state.mean[-k]))))
            }
            path[kk]<-k
          }
        }
      }
    }
    segments<-cnv.segments(obs,path)
    result$seg_mean[f]<-segments[,1]
    result$CNV_start[f]<-result$exon_start[f][c(segments[,2])]
    result$CNV_end[f]<-result$exon_end[f][c(segments[,3])]
    result$cnv[f]<-path
  }
}


colnames(result)<-c("chr","exon_start","exon_end","cnv","tumor_DOC","control_DOC",
                    "ratio_before_smoothing","ratio_after_smoothing",
                    "rationormalized_after_smoothing","CNV_start",
                    "CNV_end","seg_mean")
result<-subset(result,select=c(-ratio_after_smoothing,-ratio_before_smoothing))

write.table(result,
            quote=F,sep="\t",row.names=F)
