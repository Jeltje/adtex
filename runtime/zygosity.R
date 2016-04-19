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

options <- commandArgs(trailingOnly = T)
baf_file = options[1]
outputDir = options[2]
inputCNV = options[3]
thresh = as.numeric(options[4])
plotstate = options[5]
chrom = unlist(strsplit(options[6],","))


start.probs<-c(10,10,40,10,10,10,20,10,10,10,10,10,10,10,10)/190

pt_fun<-function(d){
  L<-2e6
  return(1 - 0.5*(1-exp(-d/(2*L))))
}


tran.matrix<-function(d,cn,nstates=15){
  #state transition matrix
  #d - distance between two observations in bp
  #cn - c(previous_cn,current_cn)
  #adapted from APOLLOH paper
  #same zygosity states i and j or i=j have A(i,j)= pt 
  #otherwise A(i,j)=(1-pt)/(nstates-1)
  A <- matrix(NA,ncol=nstates,nrow=nstates)
  colnames(A)<-rownames(A)<-c("LOH","LOH","HET","LOH","ASCNA","LOH","HET",
                              "ASCNA","LOH","ASCNA","ASCNA","LOH","ASCNA","HET",
                              "ASCNA")
  CN<-matrix(0,ncol=nstates,nrow=nstates) # copy number specific transitions
  prev_cn<-cn[1] # cn of previous position
  cur_cn<-cn[2] # cn of current position
  cn_s<-c(1,2,2,3,3,4,4,4,5,5,5,6,6,6,6) # copy number states
  CN[which(cn_s==prev_cn),which(cn_s==cur_cn)]<-1
  genos<-c(1,2,2,2,2,3,3,3,3,3,3,4,4,4,4)
  pt<-pt_fun(d)
  for (i in 1:nstates){
    z<-rownames(A)[i]
    k<-which(colnames(A)==z)
    A[i,k]<-pt
    if((genos[i]-1)>0){
      A[i,-k]<-(1-pt)/(genos[i]-1)
    }
    else{
      A[i,-k]<-(1-pt)
    }
    A[i,]<-A[i,]*CN[i,]
    t<-sum(A[i,])
    if(t>0){
      A[i,]<-A[i,]/t
    }
  }
  return(A)
}

forward_back<-function(Pi,B,d,cn,nstates=15){
  #Pi - initial state distribution for each chrom
  #B - emission probabilities (# of observations x # of states)
  #d - distance between current and previous SNP in bp for all SNPs (# of obs x 1)
  #cn - copy number of each SNP (# of obs x 1)
  #nstates - number of states
  #chr - chr of each SNP - this is for initiate the starting probabilities
  
  n<-nrow(B)
  m<-nstates
  lalpha<-lbeta<-matrix(NA,m,n)
    
  temp<-Pi*B[1,]
  sumtemp<-sum(temp)
  lscale<-log(sumtemp)
  temp<-temp/sumtemp
  lalpha[,1]<-log(temp)+lscale
  
  for(i in 2:n){
    A<-tran.matrix(d[i],c(cn[i-1],cn[i]))
    temp<-temp%*%A*B[i,]
    sumtemp<-sum(temp)
    lscale<-lscale+log(sumtemp)
    temp<-temp/sumtemp
    lalpha[,i]<-log(temp)+lscale
  }
  lbeta[,n]<-rep(0,m)
  temp<-rep(1/m,m)
  lscale<-log(m)
  
  for(i in (n-1):1){
    A<-tran.matrix(d[i+1],c(cn[i],cn[i+1]))
    temp<-A%*%(B[i+1,]*temp)
    lbeta[,i]<-log(temp)+lscale
    sumtemp<-sum(temp)
    temp<-temp/sumtemp
    lscale<-lscale+log(sumtemp)
  }
  list(la=lalpha,lb=lbeta)
}

estimate.param<-function(baf,m,pi,maxiter=1000,tol=1e-6,alpha,copy,d,chr){
  # baf - mirrored BAF at each SNP
  # m - number of states
  # pi - initial state distribution
  # maxiter - maximum iterations
  # tol - tolerance
  # alpha - initial guess for contamination level
  #copy - copy number at each SNP
  #d - distance between two SNPs
  #chr - chr of each SNP - this is for initiate the starting probabilities
  
  min_tol=1e7
  min_iter=NA
  alpha<-alpha
  
  # theoretical tumor BAFs in each state
  uT<-c(0,0,0.5,0,1/3,0,0.5,0.25,0,0.2,0.4,0,1/6,0.5,2/3)
  cT<-c(1,2,2,3,3,4,4,4,5,5,5,6,6,6,6)
  
  #means of normalized BAF with contamination
  mu_b<-(alpha+((1-alpha)*cT*uT))/((2*alpha) + ((1-alpha)*cT))
  sigma_b<-rep(0.1,m)
  
  sigma.b.next<-sigma_b
  alpha.next<-alpha
  mu_b.next<-mu_b
  mu_het<-0.5
  
  n<-length(baf) # number of observations
  pi.mat<-matrix(pi,ncol=m,nrow=sum(table(chr)[]>1),byrow=T)
  pi.next <-pi.mat
  chrom<-names(table(chr))[table(chr)[]>1]
  changes<-match(chrom,chr)
  chrom<-chrom[order(changes)]
  changes<-changes[order(changes)]
  
   
  for (iter in 1:maxiter){
    lallprobs<-matrix(data=NA,ncol=m,nrow=n)
    allprobs<-matrix(data=NA,ncol=m,nrow=n)
    for (i in 1:m){
      lallprobs[,i]<-dnorm(baf,mean=mu_b[i],sd=sigma_b[i],log=T)
      allprobs[,i]<-dnorm(baf,mean=mu_b[i],sd=sigma_b[i])
    }
    
    la<-lb<-matrix(NA,ncol=n,nrow=m)
    llk<-rep(0,n)
    
    for(j in 1:length(chrom)){
      f<-chr==chrom[j]
      fb<-forward_back(Pi=pi.mat[j,],B=allprobs[f,],d=d[f],cn=copy[f],nstates=m)
      la[,f]<-fb$la
      lb[,f]<-fb$lb
      c<-max(fb$la[,length(d[f])])
      llk_l<-c+log(sum(exp(fb$la[,length(d[f])]-c))) 
      llk[f]<-llk_l
    }
    
    copy.cT<-2*alpha+(1-alpha)*cT
    
    mu_b<-(alpha+((1-alpha)*cT*uT))/((2*alpha) + ((1-alpha)*cT))
    mu_b[c(3,7,14)]<-mu_het
    alpha.calc1<-rep(0,m)
    alpha.calc2<-rep(0,m)
    mu_het1<-0
    mu_het2<-0
    
    for(j in 1:m){
      if(j!=3 & j!=7 & j!=14){
        alpha.calc1[j]<-sum(exp(la[j,]+lb[j,]-llk)*(baf-(cT[j]*uT[j]/copy.cT[j])))*(1-cT[j]*uT[j])/copy.cT[j]
        alpha.calc2[j]<-sum(exp(la[j,]+lb[j,]-llk))*((1-cT[j]*uT[j])/copy.cT[j])^2
      }
      else{
        mu_het1<-sum(exp(la[j,]+lb[j,]-llk)*baf)+mu_het1
        mu_het2<-sum(exp(la[j,]+lb[j,]-llk))+mu_het2
      }
    }
    mu_het<-mu_het1/mu_het2
    alpha.next<-sum(alpha.calc1,na.rm=T)/sum(alpha.calc2,na.rm=T)
    alpha.next<-ifelse(alpha.next>0.7,0.7,ifelse(alpha.next<0,0,alpha.next))
    
    for(j in 1:length(changes)){
      pi.next[j,]<-t(exp(la[,changes[j]]+lb[,changes[j]]-llk[changes[j]]))
      pi.next[j,]<-pi.next[j,]/ifelse(sum(pi.next[j,])>0,sum(pi.next[j,]),1)
    }
    
    crit<-sum(sum(abs(pi.mat-pi.next))+sum(abs(alpha-alpha.next)))
    if(crit<min_tol){
      min_tol<-crit
      min_iter<-iter
      min_param<-list(mu_b=mu_b,alpha=alpha,sigma_b=sigma_b,pi=pi.mat)
    }
    
    if(crit<tol){
      return(list(mu_b=mu_b,alpha=alpha,sigma_b=sigma_b,pi=pi.mat))
    }
    
    pi.mat<-pi.next
    alpha<-alpha.next
  }
  print(paste("No convergence after",maxiter,"iterations at",tol,"tolerance level"))
  print(paste("Selecting parameters at",min_tol,"tolerance level"))
  return(min_param)
}

predict_geno<-function(baf,mu_b,sigma_b,pi,d,cn){
  # baf - normalized BAF observations
  # mu_b - mean of bafs at each state
  # sigma_b - sd of bafs at each state
  # pi - initial state distribution
  # d - distance between two SNPs
  # cn - copy number of each SNP
  
  n<-length(baf)
  m<-length(pi)
  B<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:m){
    B[,i] <- dnorm(baf,mean=mu_b[i],sd=sigma_b[i])
  }
  
  xi<-matrix(0,n,m)
  temp<-pi*B[1,]
  xi[1,]<-temp/sum(temp)
  for(i in 2:n){
    A<-tran.matrix(d[i],c(cn[i-1],cn[i]))
    temp<-apply(xi[i-1,]*A,2,max)*B[i,]
    xi[i,]<-temp/sum(temp)
  }
  iv<-numeric(n)
  iv[n]<-which.max(xi[n,])
  for(i in (n-1):1){
    A<-tran.matrix(d[i+1],c(cn[i],cn[i+1]))
    iv[i]<-which.max(A[,iv[i+1]]*xi[i,])
  }
  iv
}

snp_copy<-function(baf,cnv_res,chrom){
  cnv<-read.delim(cnv_res)
  
  baf$cn<-NA
  baf<-baf[baf$chrom%in%chrom,]
  k<-names(table(baf$chrom))[which(table(baf$chrom)[]==1)]
  if(length(k)>0){
    chrom<-chrom[-which(chrom==k)]
  }
  
  for(i in chrom){
    d1<-baf[baf$chrom==i,]
    d2<-cnv[cnv$chr==i,]
    seg<-unique(d2$CNV_start)
    seg_i<-match(seg,d2$CNV_start)
    for(j in 1:length(seg)){
      k<-which(d1$SNP_loc>=seg[j] & d1$SNP_loc<=d2$CNV_end[seg_i[j]])
      if(length(k)>0){
        if(d2$cnv[seg_i[j]]==0){
          baf$cn[baf$chrom==i][k]<-1
        }
        else{
          baf$cn[baf$chrom==i][k]<-d2$cnv[seg_i[j]]
        }
      }
    }
  }
  baf<-baf[!is.na(baf$cn),]
  return(baf)
}

plot.result<-function(bafcn,cnv_res,dir,chr){
  print("Plotting results...")
  for(i in chr){
    png(filename=paste(dir,"/",i,"_zygosity.png",sep=""),height=750,width=600,
        bg="white")
    par(mfrow=c(2,1),oma=c(0,0,0,2))
    f<-cnv_res$chr==i
    ch.len<-cnv_res$exon_end[f][sum(f)]
    col_mat<-col2rgb(cnv_res$cnv[f]+1)
    x<-(cnv_res$exon_start[f]+cnv_res$exon_end[f])/2/ch.len
    x.at<-seq(0,1,0.1)
    x.lb<-paste(round(x.at*ch.len/1000000,0),"M",sep="")
    plot(x,cnv_res$rationormalized_after_smoothing[f],
         col=rgb(col_mat[1,],col_mat[2,],col_mat[3,],50,maxColorValue=255),
         xlab="",ylab="DOC ratio",xlim=c(-0.05,1.05),
         main=paste("Chromosome",i,sep=""),pch=20,cex=0.5,ylim=c(0,3),axes=F)
    axis(2)
    axis(1,at=x.at,labels=x.lb,las=2,cex.axis=0.8)
    par(xpd=NA)#writing outside of the figure margins
    legend(x=1.06,y=3,legend=as.numeric(names(table(cnv_res$cnv[f]))),
           col=as.integer(names(table(cnv_res$cnv[f])))+1, pch=20)
    par(xpd=FALSE)#default behaviour
    ff<-bafcn$chrom==i
    x<-(bafcn$SNP_loc[ff])/ch.len
    col_mat<-col2rgb(bafcn$z[ff]+1)
    plot(x,bafcn$tumor_BAF[ff],xlab="",ylab="BAF",
         col=rgb(col_mat[1,],col_mat[2,],col_mat[3,],100,maxColorValue=255),
         pch=20,cex=0.6,ylim=c(0,1),xlim=c(-0.05,1.05),axes=F)
    axis(2)
    par(xpd=NA)#writing outside of the figure margins
    legend(x=1.06,y=1,legend=c("LOH","HET","ASCNA"),col=c(2:4), pch=20)
    par(xpd=FALSE)#default behaviour
    dev.off()
  }
  png(filename=paste(dir,"result.png",sep=""),height=750,width=1024,bg="white")
  par(mfrow=c(2,1))
  chrom<-names(table(cnv_res$chr))[table(cnv_res$chr)[]>0]
  changes<-match(chrom,cnv_res$chr)
  chrom<-chrom[order(changes)]
  changes<-changes[order(changes)]
  mt<-match(chr,chrom)
  limits<-match(chrom,cnv_res$chr)
  limits<-c(limits[-1]-1,nrow(cnv_res))
  limits<-limits[mt]
  tot<-sum(as.numeric(cnv_res$exon_end[limits]))
  end_limit<-cumsum(as.numeric(cnv_res$exon_end[limits]))
  plot(c(0,tot),c(0,max(cnv_res$rationormalized_after_smoothing)),type="n",xlab="",
       ylab="DOC ratio",xaxt="n",main="Ratio profile",ylim=c(0,3))
  end_limit2<-c(0,end_limit[-length(end_limit)])[match(chrom,chr)]
  len<-rep(end_limit2,table(cnv_res$chr)[chrom][])
  x<-(cnv_res$exon_start+cnv_res$exon_end)/2 + len
  points(x,cnv_res$rationormalized_after_smoothing,pch=20,cex=0.8,
         col=rgb(0,0,255,50,maxColorValue=255))
  abline(v=0,lty=1,col="lightgrey")
  text(end_limit[1]/2,0.5,chr[1], pos = 1, cex = 1)
  for (j in 2:length(end_limit)) {
    vpos = end_limit[j-1];
    tpos = (end_limit[j-1]+end_limit[j])/2
    text(tpos,0.5,chr[j], pos = 1, cex = 1)
    abline(v=vpos,lty=1,col="lightgrey")
  }
  
  
  plot(c(0,tot),c(0,max(bafcn$tumor_BAF)),type="n",xlab="",
       ylab="B allele frequency",xaxt="n",main="BAF profile",ylim=c(0,1))
  len<-rep(c(0,end_limit[-length(end_limit)]),table(bafcn$chrom)[chr][])
  x<-(bafcn$SNP_loc+ len)
  points(x,bafcn$tumor_BAF,col=rgb(0,0,255,50,maxColorValue=255),pch=20)
  abline(v=0,lty=1,col="lightgrey")
  text(end_limit[1]/2,0.05,chr[1], pos = 1, cex = 1)
  for (j in 2:length(end_limit)) {
    vpos = end_limit[j-1];
    tpos = (end_limit[j-1]+end_limit[j])/2
    text(tpos,0.05,chr[j], pos = 1, cex = 1)
    abline(v=vpos,lty=1,col="lightgrey")
  }
  dev.off()
}

run_hmm<-function(baf,cnv_res,dest,chrom,alpha=0.1){
  bafcn<-snp_copy(baf,cnv_res,chrom)
  params<-estimate.param(baf=bafcn$mirrored_BAF,m=15,pi=start.probs,maxiter=1000,
                         tol=1e-6,alpha=alpha,copy=bafcn$cn,
                         d=c(0,diff(bafcn$SNP_loc)+1),chr=bafcn$chrom)
  predict<-rep(0,nrow(bafcn))
  chr<-names(table(bafcn$chrom))[table(bafcn$chrom)[]>1]
  changes<-match(chr,bafcn$chrom)
  chr<-chr[order(changes)]
  for(i in 1:length(chr)){
    f<-bafcn$chrom==chr[i]
    predict[f]<-predict_geno(baf=bafcn$mirrored_BAF[f],mu_b=params$mu_b,
                             sigma_b=params$sigma_b,pi=params$pi[i,],
                             d=c(0,diff(bafcn$SNP_loc[f])+1),
                             cn=bafcn$cn[f])
  }
  bafcn$z<-predict
  bafcn$z[which(predict%in% c(1,2,4,6,9,12))]<-1
  bafcn$z[which(predict%in% c(3,7,14))]<-2
  bafcn$z[which(predict%in% c(5,8,10,11,13,15))]<-3
  cnv_res<-read.delim(cnv_res)
  if(plotstate){
    plot.result(baf=bafcn,cnv_res=cnv_res,dir=dest,chr=chr)
  }
  z<-c("LOH","HET","ASCNA")
  bafcn$zygosity<-z[bafcn$z]
  write.table(bafcn,file=paste(dest,"/zygosity.res",sep=""),quote=F,row.names=F,
              sep="\t")
  writeLines(as.character(params[c("alpha")]),con=paste(dest,"/contamination",sep=""))
}

baf_file<-read.delim(baf_file)
filter<-baf_file$control_doc>thresh
baf_file<-baf_file[filter,]
filter<-(baf_file$control_BAF<0.3)|(baf_file$control_BAF>0.7)
baf_file<-baf_file[!filter,]

baf_file$mirrored_BAF<-baf_file$tumor_BAF
baf_file$mirrored_BAF[baf_file$tumor_BAF>0.5]<-1-baf_file$tumor_BAF[baf_file$tumor_BAF>0.5]
bafextract=paste(outputDir,"/extracted.baf",sep="")
write.table(baf_file,file=bafextract,quote=F,sep="\t",row.names=F)

run_hmm(baf=baf_file,cnv_res=inputCNV,dest=outputDir,chrom=chrom)
unlink(bafextract)
