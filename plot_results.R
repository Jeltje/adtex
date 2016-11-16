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
outputLoc=options[1]
chrom=unlist(strsplit(options[2],","))

result<-read.delim(paste(outputLoc,"/cnv.result",sep=""))
print("Plotting the results :")
for(i in chrom){
  png(filename=paste(outputLoc,"/",i,".png",sep=""),height=600,width=600,
      bg="white")
  f<-result$chr==i
  col_mat<-col2rgb(result$cnv[f]+1)
  plot(result$exon_start[f],result$rationormalized_after_smoothing[f],
       col=rgb(col_mat[1,],col_mat[2,],col_mat[3,],50,maxColorValue=255),
       main=paste("Chromosome",i),xlab="Exon location",pch=20,cex=0.8,
       ylab="DOC normalized ratio",ylim=c(0,3))
  legend("topright",legend=as.numeric(names(table(result$cnv[f]))),
         col=as.integer(names(table(result$cnv[f])))+1, pch=20)
  dev.off()
}

