# ----------------------------------------------------------------------#
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

require("ifultools")

forward_back<-function(Pi,A,B,nstates){
  #Pi - initial state distribution
  #A - transition matrix
  #B - emission probabilities (# of observations x # ofstates)
  #nstates - number of states
  
  n<-nrow(B)
  m<-nstates
  lalpha<-lbeta<-matrix(NA,m,n)
  
  temp<-Pi*B[1,]
  sumtemp<-sum(temp)
  lscale<-log(sumtemp)
  temp<-temp/sumtemp
  lalpha[,1]<-log(temp)+lscale
  for(i in 2:n){
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
    temp<-A%*%(B[i+1,]*temp)
    lbeta[,i]<-log(temp)+lscale
    sumtemp<-sum(temp)
    temp<-temp/sumtemp
    lscale<-lscale+log(sumtemp)
  }
  list(la=lalpha,lb=lbeta)
}

EM<-function(x,m,A,pi,maxiter=1000,tol=1e-6,mu,sigma,d){
  # x - observations
  # m - number of states
  # A - transition probabilities
  # pi - initial state distribution
  # maxiter - maximum iterations
  # tol - tolerance
  # mu - initial guesses for mean of Gaussian emission
  # sigma - initial guesses for sd of Gaussian emission
  # d - ploidy
  
  min_tol=1e7
  min_iter=NA
  
  thresh<-c(0.2,0.15,0.1)
  #the above levels are for ploidy 4 only
  thresh.l<-c(0.1,0.7,0.9,1.15,1.4) # lower cutoff of mu
  thresh.u<-c(0.65,0.85,1.1,1.3,1.8) # upper cutoff of mu
  
  mu.next <-mu
  mu.theory<-mu
  sigma.next <-sigma
  A.next <-A
  pi.next <-pi
  n<-length(x) # number of observations
  
  for (iter in 1:maxiter){
    m<-length(mu)
    lallprobs<-matrix(data=NA,ncol=m,nrow=n)
    allprobs<-matrix(data=NA,ncol=m,nrow=n)
    for (i in 1:m){
      lallprobs[,i] <- dnorm(x,mean=mu[i],sd=sigma[i],log=T)
      allprobs[,i] <- dnorm(x,mean=mu[i],sd=sigma[i])
    }
    for(i in 1:n){
      if(sum(allprobs[i,])<=1e-100){
        x[i]<-ifelse(i>1&i<n,sum(x[i-2]+x[i+3])/6,ifelse(i==n,x[i-1],x[i+1]))
        for (j in 1:m){
          lallprobs[i,j] <- dnorm(x[i],mean=mu[j],sd=sigma[j],log=T)
          allprobs[i,j] <- dnorm(x[i],mean=mu[j],sd=sigma[j])
        }
      }
    }
    fb<-forward_back(Pi=pi,A=A,B=allprobs,nstates=m)
    la<-fb$la
    lb<-fb$lb
    c<-max(la[,n])
    llk <- c + log(sum(exp(la[,n]-c)))
    for(j in 1:m){
      for(k in 1:m){
        A.next[j,k]<-A[j,k]*sum(exp(la[j,1:(n-1)]+lallprobs[2:n,k]+lb[k,2:n]-llk))
      }
      if(sum(A.next[j,])==0){
        A.next[j,]<-A[j,]
      }
      mu.next[j]<-sum(exp(la[j,]+lb[j,]-llk)*x)/sum(exp(la[j,]+lb[j,]-llk))
      if(!is.na(mu.next[j])){
        if(d>=4){
          f<-(mu.next[j]>thresh.u[j] | mu.next[j]<thresh.l[j])
        }
        else{
          f<-abs(mu.next[j]-mu.theory[j])>thresh[d-1]
        }
        if(f){
          mu.next[j]<-ifelse(j!=m,mu[j],
                             ifelse(mu.next[j]>mu.theory[j],mu.next[j],
                                    mu[j]))
        }
      }
      else{
        mu.next[j]<-mu[j]
      }
    }
    zeros<-which(apply(A.next,1,sum)==0)
    A.next<-A.next/apply(A.next,1,sum)
    pi.next<-exp(la[,1]+lb[,1]-llk)
    pi.next<-pi.next/sum(pi.next)
    if(length(zeros)>0){
      A.next<-A.next[-zeros,-zeros]
      A<-A[-zeros,-zeros]
      pi.next<-pi.next[-zeros]
      pi<-pi[-zeros]
      mu.next<-mu.next[-zeros]
      mu<-mu[-zeros]
      sigma.next<-sigma.next[-zeros]
      sigma<-sigma[-zeros]
    }
    crit<-sum(abs(mu-mu.next))+sum(abs(A-A.next))+sum(abs(pi-pi.next))
    if(crit<min_tol){
      min_tol<-crit
      min_iter<-iter
      np<- m*m+m-1
      AIC<--2*(llk-np)
      BIC<--2*llk+np*log(n)
      min_param<-list(mu=mu,sigma=sigma,A=A,pi=pi,mllk=-llk,AIC=AIC,BIC=BIC,x=x)
    }
    if(crit<tol){
      np<- m*m+m-1
      AIC<--2*(llk-np)
      BIC<--2*llk+np*log(n)
      return(list(mu=mu,sigma=sigma,A=A,pi=pi,mllk=-llk,AIC=AIC,BIC=BIC,x=x))
    }
    mu<-mu.next
    sigma<-sigma.next
    A<-A.next
    pi<-pi.next
  }
  print(paste("No convergence after",maxiter,"iterations at",tol,"tolerance level"))
  print(paste("Selecting parameters at",min_tol,"tolerance level"))
  return(min_param)
}

viterbi<-function(x,mu,sigma,pi,A){
  #x -observations
  #mu - mean of gaussian emission
  #sigma - sd of gaussian emission
  #pi - initial state distribution
  #A - transition matrix
  
  n<-length(x)
  m<-nrow(A)
  B<-matrix(NA,nrow=n,ncol=m)
  for (i in 1:m){
    B[,i] <- dnorm(x,mean=mu[i],sd=sigma[i])
  }
  #To remove NAs if divided by zero
  for(i in 1:n){
    if(sum(B[i,])==0){
      B[i,]<-B[(i-1),]
    }
  }
  xi<-matrix(0,n,m)
  temp<-pi*B[1,]
  xi[1,]<-temp/sum(temp)
  for(i in 2:n){
    temp<-apply(xi[i-1,]*A,2,max)*B[i,]
    xi[i,]<-temp/sum(temp)
  }
  iv<-numeric(n)
  iv[n]<-which.max(xi[n,])
  for(i in (n-1):1){
    iv[i]<-which.max(A[,iv[i+1]]*xi[i,])
  }
  iv
}


predict.cnv<-function(x,d=2,goto.cnv=1e-3,goto.normal=1/20,goto.prob=1e-10,
                      models=list(model1=c(0:4),model2=c(0:3),model3=c(1:4),
                                  model4=c(1:3)),predict=list()){
  #x - ratios
  #d - normal state
  #goto.cnv - initial probability assumption to produce transition from normal 
  #           to cnv
  #goto.normal- initial probability assumption to produce transition from cnv 
  #           to normal
  #goto.prob - initial probability assumption for other transitions
  #models - list of models to be fitted
  #x<-madWins(x,2.5,5)$ywin
  BIC<-100000
  model.no<-0
  cont<-0
  mean<-1
  for(model in models){
    model.no=model.no+1
    S<-(model+d-2)
    normal.state <- which(S==d)
    nstates<-length(S)
    mu <- sapply(1:nstates,function(j)(mean*(cont + (1-cont)*(S[j]/d + ifelse(S[j]==0,0.1,0)))))
    sigma<-rep(sd(x),nstates)
    #if(sigma[1]>0.9){
      #sigma<-rep(0.2,nstates)
    #}
    pi<-rep(1,nstates)/nstates
    A<-matrix(NA,ncol=nstates,nrow=nstates)
    for(k in 1:ncol(A)){
      A[k,k]<-0.5
      A[k,-k]<-0.5/(nstates-1)
    }
    if(model.no==1){
      fitted.model<-list(mu=mu,sigma=sigma,A=A,pi=pi,model=S)
    }
    kk<-EM(x=x,m=length(S),A=A,pi=pi,maxiter=1000,tol=0.01,mu=mu,sigma=sigma,d=d)
    if(kk$BIC<BIC){#selecting best model with minimum BIC (i.e. min. information loss)
      if(all(diff(kk$mu)>=0)){ 
        BIC<-kk$BIC
        fitted.model<-kk
        fitted.model$model<-S
        fitted.model$mu.init<-mu
        if(fitted.model$sigma[1]>1){
            fitted.model$sigma<-rep(0.5,nstates)
          }
        if(abs(kk$mu[normal.state]-mean)>0.1){
          fitted.model$mu[normal.state]<-mean
        }
      }
    }
  }
  #print(fitted.model$mu)
  #print(fitted.model$model)
  kk.path<-viterbi(x=fitted.model$x,mu=fitted.model$mu,sigma=fitted.model$sigma,
                   pi=fitted.model$pi,A=fitted.model$A)
  cnv<-sapply(kk.path,function(j)(fitted.model$model[j]))
  fitted.model$cnv<-cnv
  return(fitted.model)
}


generate.data<-function(control,tumor,exome){
  data<-exome
  data$control<-control
  data$tumor<-tumor
  data$width<-exome$End-exome$Start
  d1<-matrix(c(data$control,data$tumor),ncol=2,nrow=nrow(data))
  med<-apply(d1,1,median)
  data$bg<-med
  cutoff<-median(data$control)
  f<-(data$control<cutoff)
  data$ratio<-data$tumor/data$control
  f<-data$control<quantile(data$control,p=0.5)
  data$ratio[f]<-wavShrink(data$ratio[f], wavelet="haar",thresh.fun="adaptive",shrink.fun="soft",
                           thresh.scale=1, xform="dwt",
                           n.level=ilogb(length(data$ratio[f]),base=2))
  return(data)
}

cnv.segments<-function(obs,path){
  changes<-which(path[-length(path)]!=path[-1])+1
  range.start<-c(1,changes)
  range.end<-c(changes-1,length(obs))
  range.mean<-sapply(1:length(range.end),
                     function(x)(mean(obs[range.start[x]:range.end[x]])))
  segments<-matrix(data=NA,ncol=3,nrow=length(obs))
  for(i in 1:length(range.end)){
    segments[range.start[i]:range.end[i],1]<-range.mean[i]
    segments[range.start[i]:range.end[i],2]<-range.start[i]
    segments[range.start[i]:range.end[i],3]<-range.end[i]
  }
  return(segments)
}

combine.segments<-function(obs,path){
  end<-T
  while(end){
    changes<-which(path[-length(path)]!=path[-1])+1
    tcount<-as.numeric(names(which.max(tapply(obs,path,length))))
    range.start<-c(1,changes)
    range.end<-c(changes-1,length(obs))
    range.mean<-sapply(1:length(range.end),
                       function(x)(mean(obs[range.start[x]:range.end[x]])))
    len<-range.end-range.start+1
    f<-which(len<5)
    if(length(f)>0){
      #print(f)
      table(path)
      for(j in f){
        if(path[range.start[j]]!=tcount){
          if(range.start[j]>2 & range.end[j]<length(obs)){
            s<-ifelse(which.min(c(abs(range.mean[j]-range.mean[j-1]),abs(range.mean[j]-range.mean[j+1])))==1,
                      path[range.start[j-1]],path[range.start[j+1]])
          }
          else if(range.start[j]==1){
            s<-path[range.start[j+1]]
          }
          else{
            s<-path[range.start[j-1]]
          }
          path[range.start[j]:range.end[j]]<-s
        }
        else{
          end<-F
        }
      }
    }
    else{
      end<-F
    }
  }
  segments<-matrix(data=NA,ncol=4,nrow=length(obs))
  segments[,1:3]<-cnv.segments(obs,path)
  segments[,4]<-path
  return(segments)
}

smoothData<-function(data,chrom=c(1:22)){
  data$ratio_after_smoothing<-data$ratio
  for(i in chrom){
    f<-data$chr==i
    df<-data[f,]
    f1<-df$control<median(df$control)
    df$ratio[f1]<-wavShrink(df$ratio[f1], wavelet="haar",thresh.fun="adaptive",shrink.fun="soft",
                            thresh.scale=1, xform="dwt",
                            n.level=ilogb(length(df$ratio[f1]),base=2))
    data$ratio_after_smoothing[f]<-df$ratio
  }
  return(data)
}

selectPeak<-function(ratio){
  dens<-density(ratio)
  all.peaks<-peaks(dens$y)
  peak.value<-dens$y[all.peaks]
  peak.max<-peak.value[which.max(peak.value)]
  sig.peaks<-(peak.value/peak.max) > 0.70
  peak.ratio<-dens$x[all.peaks][sig.peaks]
  main.peak<-peak.ratio[which.min(abs(peak.ratio-1))]
  #print(main.peak)
  return(main.peak)
}


####### Following codes are from ASCAT to remove outliers##############

madWins <- function(x,tau,k){
  xhat <- medianFilter(x,k)
  d <- x-xhat
  SD <- mad(d)
  z <- tau*SD
  xwin <- xhat + psi(d, z)
  outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
  return(list(ywin=xwin,sdev=SD,outliers=outliers))
}

psi <- function(x,z){
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}


#########################################################################
# Function to calculate running median for a given a window size
#########################################################################

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  
  runMedian <- runmed(x,k=filtWidth,endrule="median")
  
  return(runMedian)
  
}

######### Above codes are from ASCAT to remove outliers############

######################################################################
## Similar method to ASCAT ###########################################
######################################################################

BAFseg_obj<-function(baf_seg,cnv_res){
  chr<-c(1:22)
  kk<-which(baf_seg$data$chrom %in% chr)
  if(length(kk)>0){
    snp<-matrix(0,nrow=length(kk),ncol=6)
    snp<-data.frame(snp)
    colnames(snp)<-c("chrom","loc","normalized_BAF","seg_mean","cnv","doc_ratio")
    snp$cnv<-rep(NA,length(kk))
  }
  id=1
  for(i in chr){
    k<-which(baf_seg$data$chrom==i)
    id2<-length(k)
    if(id2>0){
      t<-id+id2-1
      snp$chrom[id:t]<-baf_seg$data$chrom[k]
      snp$loc[id:t]<-baf_seg$data$maploc[k]
      snp$normalized_BAF[id:t]<-baf_seg$data[,3][k]
      k2<-which(baf_seg$output$chrom==i)
      snp$seg_mean[id:t]<-populate_mean(snp[id:t,],baf_seg$output[k2,])
      id=id+id2
    }
  }
  k<-NULL
  for(c in chr){
    d1<-cnv_res[cnv_res$chr==c,]
    f<-snp$chrom==c
    d2<-snp[f,]
    for(j in 1:length(d2[,1])){
      k<-which(d2$loc[j]>=d1$CNV_start & d2$loc[j]<=d1$CNV_end)
      if(length(k)>0){
        snp$cnv[f][j]<-d1$cnv[k][1]
        snp$doc_ratio[f][j]<-d1$seg_mean[k][1]
      }
    }
  }
  snp<-snp[!is.na(snp$cnv),]
  return(snp)
}

##this function will add seg mean for each snp position
populate_mean<-function(snp,bafseg){
  if(length(snp$chrom)==sum(bafseg$num.mark)){
    seg_mean<-rep(-1,length(snp$chrom))
  }
  id=1
  for(i in 1:length(bafseg$ID)){
    id2<-bafseg$num.mark[i]
    seg_mean[id:(id+id2-1)]<-rep(bafseg$seg.mean[i],id2)
    id=id+id2
  }
  return(seg_mean)
}

makedistance_matrix<-function(BAFsegobj){
  snp<-BAFsegobj
  alpha<-seq(0.1,0.9,0.01)
  d<-rep(NA,length(alpha))
  
  for(i in 1:length(alpha)){
    a<-alpha[i]
    NB_hat<-((2*a+(1-a)*snp$cnv)*snp$seg_mean - a)/(1-a)
    d[i]<-sum((NB_hat - round(NB_hat))^2*ifelse(snp$seg_mean==0.5,0.05,1))
  }
  return(d)
}
