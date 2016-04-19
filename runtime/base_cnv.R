# ----------------------------------------------------------------------#
# Note: this code has been altered slightly to fit in the Adtex docker container
# by Jeltje van Baren
# 
# Original copyright notice:
#
# Copyright (c) 2013, Kaushalya Amarasinghe.
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

options <- commandArgs(trailingOnly = T)
outF = options[1]
outdir = options[2]


alpha<-seq(from=0.1,to=0.9,by=0.05)

fit_score<-function(outF,i){
  cnv<-read.delim(paste(outF,"/cnv",i,"_baf.txt",sep=""),header=F)
  l<-as.numeric(cnv[,3])-as.numeric(cnv[,2])
  g<-paste(cnv[,1],cnv[,4],cnv[,5],sep="")
  ll<-round((cnv[,4]/tapply(l,g,sum)[g])*l)
  ll[ll==0]<-1
  ll[is.na(ll)]<-1
  t<-median(rep(cnv[,5],ll[]))
  t<-ifelse(t>0.41,t,0.45)
  cnv[cnv[,5]>t,5]<-0.5
  w<-rep(1,nrow(cnv))
  w[cnv[,5]>t]<-0.05
  err<-rep(100000,length(alpha))
  b_err<-rep(100000,length(alpha))
  thresh<-sum(0.3^2*ll*w)
  #print(thresh)
  
  for(i in 1:length(alpha)){
    a<-alpha[i]
    NB.hat<-(((2*a+(1-a)*cnv[,9])*cnv[,5]) - a)/(1-a)
    if(all(round(NB.hat)>=0) & all(round(NB.hat)<=(cnv[,9]%/%2))){
      b<-(a+(1-a)*round(NB.hat))/(2*a+(1-a)*cnv[,9])
      r<-(2*a+(1-a)*cnv[,9])/i
      #err[i]<-sum((b-cnv[,10])^2*cnv[,9]*w)
      #print(sum((r-cnv[,10])^2))
      #print(table(round(NB.hat)))
      err[i]<-sum((round(NB.hat)-NB.hat)^2*ll*w)
      b_err[i]<-sum((b-cnv[,5])^2*ll*w)
      #print((-100/thresh)*err[i]+100)
      #print(mean((-100/thresh)*(round(NB.hat)-NB.hat)^2*ll*w+100))
      err[i]<-ifelse(err[i]>thresh,100000,err[i])
    }
  }
  return(list(err=err,b_err=b_err))
}

e<-matrix(NA,ncol=length(alpha),nrow=3)
b_e<-matrix(NA,ncol=length(alpha),nrow=3)
for(i in 2:4){
  res<-fit_score(outF,i)
  e[(i-1),]<-res$err
  b_e[(i-1),]<-res$b_err
}

j<-which.min(e)
k<-which.min(b_e)
j_row<-row(e==e[j])[j]+1
k_row<-row(b_e==b_e[k])[k]+1
if(j==k | (j_row+k_row)%%2!=0 | e[j]/e[k]<0.7){
  print(paste("estimated base ploidy =",j_row))
  p=j_row
}else{
  print(paste("estimated base ploidy =",k_row))
  p=k_row
}
writeLines(as.character(p),con=paste(outdir,"/ploidy",sep=""))
