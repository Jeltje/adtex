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
outF = options[1]

read_file<-function(outF,i){
  cnv<-read.delim(paste(outF,"/cnv.result",i,sep=""))
  chr<-cnv[,1]
  cnv<-matrix(c(cnv$chr,cnv$CNV_start,cnv$CNV_end,cnv$cnv,cnv$seg_mean),ncol=5)
  cnv<-as.data.frame(cnv)
  cnv[,1]<-chr
  cnv<-cnv[!duplicated(cnv[,1:3]),]
  write.table(cnv,file=paste(outF,"/cnv",i,sep=""),quote=F,sep="\t",row.names=F,
              col.names=F)
}

for(i in 2:4){
  read_file(outF,i)
}
