#!/usr/bin/python2.7

# ----------------------------------------------------------------------#
# Note from Jeltje van Baren, 2015: slightly altered original code to avoid
# copying large input files, and removed concatenation of png files into a pdf
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

import os
import argparse
import sys
import subprocess
import shlex
from multiprocessing import Process, Manager
from getMeanCoverage import *


#absolute script path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))

class Options:

	def __init__(self):
		self.parser = argparse.ArgumentParser("Aberration Detection in Tumour EXome")
		self.parser.add_argument('-n','--normal', help = 'Matched normal sample in BAM format or bed formatted coverage [REQUIRED], to generate bed formatted coverage please see documentation',dest='control')
		self.parser.add_argument('-t','--tumor', help = 'Tumor sample in BAM format or bed format DOC coverage [REQUIRED], to generate bed formatted coverage please see documentation',dest='tumor')
		self.parser.add_argument('-b','--bed',help = 'BED format of the targeted regions [REQUIRED]', dest='bed')
		self.parser.add_argument('-o','--out',help="Output folder path name to store the output of analysis [REQUIRED]",
			action="store", dest="outFolder") 
		self.parser.add_argument("--DOC",
			help="If specified, matched normal and tumor inputs  will be in BED formatted coverage [False]",
			action="store_true", dest = "doc",default="False")
		self.parser.add_argument("--ploidy",help="Most common ploidy in the tumour sample [2]", dest="ploidy",action="store")
		self.parser.add_argument("--estimatePloidy", help="If provided, --baf must be specified to estimate base ploidy [FALSE]",dest="p_est", action="store_true",default="False")
		self.parser.add_argument("--minReadDepth", 
			help="The threshold for minimum read depth for each exon [10]",
			action="store", dest="minReadDepth", default=10) 
		self.parser.add_argument("-p", "--plot", 
			help="Plots each chromosome with CNV estimates [False]", 
			action="store_true", dest="plot", default="False")
		self.parser.add_argument("--baf", help="File containing B allele frequencies at heterozygous loci of the normal [optional]", dest="baf",action="store")
		
		args = self.parser.parse_args()
		if args.control:
			self.control = args.control
		else:
			self.parser.error("Matched normal sample not supplied")
		if args.tumor:
			self.tumor = args.tumor
		else:
			self.parser.error("Tumor sample not supplied")
		if args.bed:
			self.bed = args.bed
		else:
			self.parser.error('Targeted regions file not supplied')
		if args.baf:
			self.baf = args.baf
			self.bafin = "True"
		else:
			self.bafin = "False"
			self.baf = "False"
		if args.outFolder:
			self.outFolder = str(args.outFolder)
		else:
			self.parser.error('Output folder path not supplied')
		if str(args.doc)=="True":
			#self.totReads = args.doc
			self.docInput = "True"
		else:
			self.docInput = "False"
		if args.ploidy:
			self.ploidy = args.ploidy
			self.p_est="False"
			if str(args.p_est)=="True":
				print "Ploidy provided. Estimation step won't be extecuted"
			self.ploidyIn = "True"
		else:
			self.ploidyIn = "False"
		if str(args.p_est)=="True" and str(self.ploidyIn)=="False":
			if args.baf:
				self.p_est="True"
			else:
				self.parser.error("--baf must be provided if base ploidy estimation is used")
		elif str(self.ploidyIn)=="False":
			self.ploidy=2
			self.p_est = "False"			
		if args.minReadDepth:
			self.minRead = args.minReadDepth
		if args.plot:
			self.plot = args.plot


def splitBam(inF,outFolder,chroms):
	try:
		os.mkdir(outFolder + "/chr/")
	except:
		print "Folder exists"
		 					
	for c in chroms:
		outputfile = outFolder + "/chr/" + c+".bam"
		subprocess.call("samtools view -bh %s %s > %s" %(inF,c,outputfile),shell=True)
		

def splitBed(inF,outFolder,chroms):
	try:
		os.mkdir(outFolder + "/chr/")
	except:
		print "Folder exists"
	inF = open(inF,"r")
	check = "1"
	outputfile = outFolder + "/chr/" + str(check)+".bed"
	outfile = open(outputfile,"w")
	
	for row in inF:
		cols = row.split()
		chr = cols[0]
		if (chr != check):
			outfile.close()
			check = chr
			outputfile = outFolder + "/chr/" + str(check)+".bed"
			outfile = open(outputfile,"w")
		outfile.write(row.rstrip()+"\n")

	outfile.close()
		
def getCoverage(outF,bedF,chroms):
	outFile = outF+"/coverage.txt"
	iOutFile = open(outFile, "w")
	for c in chroms:
		inFile = outF + "/chr/" +c+".bam"
		targets = bedF + "/chr/"+c+".bed"
		args = shlex.split("coverageBed -abam %s -d -b %s" %(inFile, targets))
		output = subprocess.Popen(args, stdout = subprocess.PIPE).communicate()[0]
		iOutFile.write(output)
	iOutFile.close()
	subprocess.call("sort  -V -k1 -k2n -k3n -k4n %s > %s" %(outFile,outFile+".sorted"),shell=True)
	
def getCoveragefromBAM(outF,bedF,inF):
	outFile = outF+"/coverage.txt"
	targets = bedF
	inFile = inF
	
	#to remove duplicates
	#outbam = outF+"/noduplicates.bam"
	#subprocess.call("samtools view -F 0x400 %s > %s" %(inFile,outbam),shell=True)
	#subprocess.call("coverageBed -abam %s -d -b %s > %s" %(outbam,targets,outFile),shell=True)
	
	subprocess.call("coverageBed -abam %s -d -b %s > %s" %(inFile,targets,outFile),shell=True)
	subprocess.call("sort -V -k1 -k2n -k3n -k4n %s > %s" %(outFile,outFile+".sorted"),shell=True)


def sortFile(inF,fileN):
	inFile = inF+fileN
	subprocess.call("sort  -V -k1 -k2n -k3n -k4n %s > %s" %(inFile,inFile+".sorted"),shell=True)

def getTotReads(inF,folder):
	subprocess.call("samtools view %s | wc -l > %s/tot_reads.txt" %(inF,folder),shell=True)

def analyseCNV(params,ratio_data,outF,chroms):
	print "Analysing CNV..."
	rScriptName = os.path.join(scriptPath, "cnv_analyse.R")
	rFunctionsPath = os.path.join(scriptPath, "RFunction.R")
	
	def runCNV(par,ratios,outLoc,chrom,rScr,rFunctions,p):
		args = shlex.split("Rscript %s %s %s %s %s %s %s %s %s %s"
		%(rScr,rFunctions, ratios,p,par.minRead,outLoc,par.bafin,par.baf,par.plot,chrom))
		rscr = subprocess.call(args)
		
	if(str(params.p_est)=="False"):
		args = shlex.split("Rscript %s %s %s %s %s %s %s %s %s %s"
		%(rScriptName,rFunctionsPath, ratio_data,params.ploidy,params.minRead,outF,params.bafin,params.baf,params.plot,chroms))
		rscr = subprocess.call(args)
		args = shlex.split("mv %s %s" %(outF+"/temp/cnv.result"+str(params.ploidy),outF+"/cnv.result"))
		rscr = subprocess.call(args)
	else:
		cnv2 = Process(target= runCNV, args=(params,ratio_data,outF,chroms,rScriptName,rFunctionsPath,2))
		cnv3 = Process(target= runCNV, args=(params,ratio_data,outF,chroms,rScriptName,rFunctionsPath,3))
		cnv4 = Process(target= runCNV, args=(params,ratio_data,outF,chroms,rScriptName,rFunctionsPath,4))
		cnv2.start()
		cnv3.start()
		cnv4.start()
		cnv2.join()
	    	cnv3.join()
	    	cnv4.join()
		

def segmentRatio(params,cCoverage,tCoverage,outF,chroms):
	rScriptName = os.path.join(scriptPath,"segment_ratio.R")
	rFunctionsPath = os.path.join(scriptPath,"RFunction.R")
	args = shlex.split("Rscript %s %s %s %s %s %s %s %s %s %s" %(rScriptName,rFunctionsPath,cCoverage,tCoverage,params.minRead,outF,params.bafin,params.baf,params.ploidyIn,chroms))
	rscr = subprocess.call(args)

def zygosity(params,outF,chroms):
	print "Predicting Zygosity states"
	rScriptName = os.path.join(scriptPath,"zygosity.R")
	args = shlex.split("Rscript %s %s %s %s %s %s" %(rScriptName,params.baf,outF,params.minRead,params.plot,chroms))
	rscr = subprocess.call(args)

def getChroms(inF,outF):
	outFile=outF+"/targets.sorted"
	subprocess.call("sort  -V -k1 -k2n -k3n -k4n %s > %s" %(inF,outFile),shell=True)
	infile=open(outFile)
	

	chr=[]
	chr_prev=0

	for line in infile:
		chr_current=line.rstrip().split("\t")[0]
		if (chr_current!=chr_prev and chr_current[0]!="G"):
			chr_prev=chr_current
			chr.append(chr_current)
	
	chr=",".join(chr)
	return(chr)
	
def main():

	subprocess.call("date",shell=True)
	
	options = Options()
	
	control = options.control
	tumor = options.tumor
	targets = options.bed
	outF = options.outFolder
	docInput = options.docInput
	bafIn = options.bafin
	
	print "Creating output folder"
	
	if outF[len(outF)-1] == "/":
		outF = outF[:len(outF)-1]
	try:
		os.mkdir(outF)
	except:
		print "cannot create folder '%s'" %outF
		print "if folder already exist, please specify other folder"
		sys.exit(1)
		
	os.mkdir(outF+"/temp")
	os.mkdir(outF+"/temp/control")
	os.mkdir(outF+"/temp/tumor")
	
	chroms = getChroms(targets,outF)
	targets = outF+"/targets.sorted"
		
	if (str(docInput)=="True"):
		print "Generating mean coverage files..."
                os.symlink(os.path.abspath(control), os.path.join(outF, 'temp', 'control', 'coverage.txt'))
                os.symlink(os.path.abspath(tumor), os.path.join(outF, 'temp', 'tumor', 'coverage.txt'))
		ctrSort=Process(target= sortFile, args=(outF+"/temp/control","/coverage.txt"))
		tmrSort=Process(target= sortFile, args=(outF+"/temp/tumor","/coverage.txt"))
		ctrSort.start()
		tmrSort.start()
		ctrSort.join()
		tmrSort.join()
	   	
	else:
		print "Creating coverage files"
	
                # Following codes are for generating coverage for the whole bam at once
		ctrDOC = Process(target= getCoveragefromBAM, args=(outF+"/temp/control",targets,control))
		tmrDOC = Process(target= getCoveragefromBAM, args=(outF+"/temp/tumor",targets,tumor))
		ctrDOC.start()
		tmrDOC.start()
		ctrDOC.join()
		tmrDOC.join()
	    	
    	
    	ctrDOC = Process(target= getMeanCoverage, args=(outF+"/temp/control/coverage.txt.sorted",outF+"/control.coverage"))
	tmrDOC = Process(target= getMeanCoverage, args=(outF+"/temp/tumor/coverage.txt.sorted",outF+"/tumor.coverage"))
	ctrDOC.start()
	tmrDOC.start()
	ctrDOC.join()
    	tmrDOC.join()
    	
    	subprocess.call("rm -rf %s" %(outF+"/temp"),shell=True)
    	
    	os.mkdir(outF+"/temp")
    	  	
    	segmentRatio(options,outF+"/control.coverage",outF+"/tumor.coverage",outF,chroms)
    	
    	chroms=open(outF+"/chrom").readline()
    	
    	analyseCNV(options,outF+"/ratio.data",outF,chroms)
    	
    	if(str(options.p_est)=="True"):
   		print "Estimating base ploidy..."
    		rScriptName = os.path.join(scriptPath,"extract_cnv.R")
		args = shlex.split("Rscript %s %s" %(rScriptName,outF+"/temp"))
		rscr = subprocess.call(args)
		
		subprocess.call("intersectBed -a %s -b %s -wb > %s" %(outF+"/temp/snp_segments",outF+"/temp/cnv2",outF+"/temp/cnv2_baf.txt"),shell=True)
		subprocess.call("intersectBed -a %s -b %s -wb > %s" %(outF+"/temp/snp_segments",outF+"/temp/cnv3",outF+"/temp/cnv3_baf.txt"),shell=True)
    		subprocess.call("intersectBed -a %s -b %s -wb > %s" %(outF+"/temp/snp_segments",outF+"/temp/cnv4",outF+"/temp/cnv4_baf.txt"),shell=True)
    		
    		rScriptName = os.path.join(scriptPath,"base_cnv.R")
		args = shlex.split("Rscript %s %s" %(rScriptName,outF+"/temp"))
		rscr = subprocess.call(args)
    		
    		ploidy=open(outF+"/temp/ploidy").readline()
    		   		
    		args = shlex.split("mv %s %s" %(outF+"/temp/cnv.result"+str(ploidy),outF+"/cnv.result"))
		rscr = subprocess.call(args)
    	
    	subprocess.call("rm -rf %s" %(outF+"/temp"),shell=True)
    	
    	if str(options.plot)=="True":
		rScriptName = os.path.join(scriptPath,"plot_results.R")
		args = shlex.split("Rscript %s %s %s" %(rScriptName,outF,chroms))
		rscr = subprocess.call(args)
    		
    	if(str(bafIn) =="True"):
    		subprocess.call("mkdir %s" %(outF+"/zygosity"),shell=True)
    		zygosity(options,outF,chroms)
    		
    	subprocess.call("rm %s %s %s" %(outF+"/chrom",outF+"/ratio.data",outF+"/targets.sorted"),shell=True)
    	
    	subprocess.call("date",shell=True)
	

if __name__ == "__main__":
	main()
