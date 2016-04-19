#!/usr/bin/python2.7 -u

# ----------------------------------------------------------------------#
# Note from Jeltje van Baren, 2015: 
# reworked the code to make more robust.
# Create temporary output on the docker container instead of in the working directory
# Reimplement DNAcopy to take in log scores instead of ratios (code doesn't work properly on those)
# Catch runtime errors
# Altered original code to avoid
# copying large input files, and removed concatenation of png files into a pdf
# Also removed option to input bam files, since the bedtools call fills up the drive
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
import shutil
#import shlex   # remove me
#from multiprocessing import Process, Manager  # remove me


#absolute script path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))

class Options:
    """
        Note that True and False are set to strings - this is done to pass them to R programs in the pipeline
    """
    def __init__(self):
        self.parser = argparse.ArgumentParser('Aberration Detection in Tumour EXome')
        self.parser.add_argument('-n','--normal', help = 'Matched normal sample in bed formatted coverage [REQUIRED], to generate bed formatted coverage please see documentation',dest='control')
        self.parser.add_argument('-t','--tumor', help = 'Tumor sample in bed format DOC coverage [REQUIRED], to generate bed formatted coverage please see documentation',dest='tumor')
        self.parser.add_argument('-b','--bed',help = 'BED format of the targeted regions [REQUIRED]', dest='bed')
        self.parser.add_argument('-o','--out',help='Output folder name to store temporary files, folder will be created as a subfolder of the working directory if it does not exist [REQUIRED]', action='store', dest='outFolder') 
        self.parser.add_argument('--keeptemp',help='do not delete temporary output folder adtexout', dest='keeptemp',action='store_true')
        self.parser.add_argument('--ploidy',help='Most common ploidy in the tumour sample [2]', dest='ploidy',action='store')
        self.parser.add_argument('--estimatePloidy', help='If provided, --baf must be specified to estimate base ploidy [FALSE]',dest='p_est', action='store_true',default='False')
        self.parser.add_argument('--minReadDepth', 
            help='The threshold for minimum read depth for each exon [10]',
            action='store', dest='minReadDepth', default=10) 
        self.parser.add_argument('-p', '--plot', 
            help='Plots each chromosome with CNV estimates [False]', 
            action='store_true', dest='plot', default='False')
        self.parser.add_argument('--baf', help='File containing B allele frequencies at heterozygous loci of the normal [optional]', dest='baf',action='store')
        
        args = self.parser.parse_args()
        if args.control:
            self.control = args.control
        else:
            self.parser.error('Matched normal sample not supplied')
        if args.tumor:
            self.tumor = args.tumor
        else:
            self.parser.error('Tumor sample not supplied')
        if args.bed:
            self.bed = args.bed
        else:
            self.parser.error('Targeted regions file not supplied')
        if args.baf:
            self.baf = args.baf
            self.bafin = 'True'
        else:
            self.bafin = 'False'
            self.baf = 'False'
        if args.outFolder:
            self.outFolder = str(args.outFolder)
        else:
            self.parser.error('Output folder path not supplied')
        if args.keeptemp:
            self.keeptemp = True
        if args.ploidy:
            self.ploidy = args.ploidy
            self.p_est='False'
            if str(args.p_est)=='True':
                print "Ploidy provided. Estimation step won't be executed"
            self.ploidyIn = 'True'
        else:
            self.ploidyIn = 'False'
        if str(args.p_est)=='True' and str(self.ploidyIn)=='False':
            if args.baf:
                self.p_est='True'
            else:
                self.parser.error('--baf must be provided if base ploidy estimation is used')
        elif str(self.ploidyIn)=='False':
            self.ploidy=2
            self.p_est = 'False'            
        if args.minReadDepth:
            self.minRead = args.minReadDepth
        if args.plot:
            self.plot = 'True'

# Convenience functions used in the pipeline
def mkdir_p(path):
    """
    Allow overwriting an existing directory
    """
    if not os.path.isdir(path):
        try:
            os.makedirs(path)
        except:
            raise

def sortFile(infile,outfile):
    """
    Sort input coverage file by chromosome and start, end positions
    """
    with open(outfile, 'w') as o:
        subprocess.check_call(['sort', '-V', '-k1', '-k2n', '-k3n', '-k4n', infile], stdout=o)

def getChroms(infile):
    """
    Returns a list of chromosome IDs found in the first column of the input file, in the same order
    """
    chr=[]
    chr_prev=0

    with open (infile, 'r') as f:
        for line in f:
            chr_current=line.rstrip().split("\t")[0]
            if (chr_current!=chr_prev and chr_current[0]!="G"):
                chr_prev=chr_current
                chr.append(chr_current)
    
    return ",".join(chr)
    

def getMeanCoverage(infile,outfile):
    """
    Calculate mean coverage in input file, outputs a file with chrom, startpos, endpos, and mean coverage
    """
    with open(infile, 'r') as inF, open (outfile, 'w') as outF:

        tot = 0
        count = 0

        for line in inF:
            l = line.rstrip().split('\t')
            ll = len(l)
        
            tot = tot + float(l[ll-1])
            if (float(l[ll-1]) <10):
                count = count + 1
        
            if (int(l[2]) - int(l[1]) == int(l[ll-2])):
                #if count > 10:
                if count > -1:
                    chr = str(l[0])
                    start  = str(l[1])
                    end = str(l[2])
                    mean = round(float(tot)/int(l[ll-2]),0)
                    outF.write(chr+"\t"+start+"\t"+end+"\t"+str(mean)+"\n")
                tot = 0
                count =0
                

def segmentRatio(params, cCoverage, tCoverage, outF, chroms, ratioOutfile, snpSegfile):
    """
    Calls segment_ratio.R
    if bafin is positive but no ploidy is asked, the program outputs a second file snp_segments, 
    which is then called later from ADTEx.py
    chroms is a comma separated list of chromosomes
    """
    rScriptName = os.path.join(scriptPath,'segment_ratio.R')
    rFunctionsPath = os.path.join(scriptPath,'RFunction.R')
    subprocess.check_call(['Rscript', rScriptName, rFunctionsPath, cCoverage, 
       tCoverage, str(params.minRead), outF, params.bafin, params.baf, 
       params.ploidyIn, ratioOutfile, snpSegfile, chroms])


def analyseCNV(params, ratio_data, outdir, tmpdir,  chroms):
    """
    Calls cnv_analyse.R to estimate copynumber. When ploidy estimation is required, it runs the script three times, once
    for ploidy 2, 3, and 4. Output (in outdir) is called cnv.result or (in tmpdir) cnv.result2, cnv.result3, cnv.result4
    """
    rScriptName = os.path.join(scriptPath, 'cnv_analyse.R')
    rFunctionsPath = os.path.join(scriptPath, 'RFunction.R')
    outfiles = []
    # if we don't want to estimate ploidy
    if(str(params.p_est)=="False"):
        outfile = os.path.join(outdir, "cnv.result")
        outfiles.append(outfile)
        with open(outfile, 'w') as o:
            subprocess.check_call(['Rscript', rScriptName, rFunctionsPath, 
                   ratio_data, str(params.ploidy), str(params.minRead),  
                   params.bafin, params.baf, params.plot, chroms], stdout=o)
    else:
        # if we DO want to estimate ploidy, run CNV three times, once for each ploidy and output files as temporaries
	# in tmpdir
        for ploidy in ['2', '3', '4']:
            outfile = os.path.join(tmpdir, "cnv.result" + ploidy)
            outfiles.append(outfile)
            with open(outfile, 'w') as o:
                subprocess.check_call(['Rscript', rScriptName, rFunctionsPath, 
                    ratio_data, ploidy, str(params.minRead),  
                    params.bafin, params.baf, params.plot, chroms], stdout=o)
    return outfiles

def zygosity(params, outdir, cnv, chroms):
    """
    Calls zygosity.R, which outputs two files: contamination (with a contamimation number) and zygosity.res, which is
    a table with zygosity numbers. If plots are requested, these are put in outdir too.
    """
    rScriptName = os.path.join(scriptPath,"zygosity.R")
    subprocess.check_call(['Rscript', rScriptName, params.baf, outdir, cnv, str(params.minRead), params.plot, chroms])

def estimatePloidy(tmpdir, workdir, snpSegfile):
    """
    Runs extract_cnv.R, bedtools intersect, and base_cnv.R.
    extract_cnv.R expects cnv.result<ploidy> and outputs cnv<ploidy>
    bedtools then intersects the cnv<ploidy> file with the snpSegfile created in segmentRatio
    base_cnv uses the intersectfiles to determine the correct ploidy, which it writes to a file named ploidy
    The corresponding file is then moved to the working directory
    """
    rScriptName = os.path.join(scriptPath,"extract_cnv.R")
    subprocess.check_call(['Rscript', rScriptName, tmpdir])
    for i in ["2", "3", "4"]:
        cnvfile = os.path.join(tmpdir, 'cnv' + i)
        outfile = os.path.join(tmpdir, 'cnv' + i + "_baf.txt")
        with open(outfile, 'w') as o:
            subprocess.check_call(['bedtools', 'intersect', 
                '-a', snpSegfile,
                '-b', cnvfile,
                '-wb'
            ], stdout = o)
        
    rScriptName = os.path.join(scriptPath,"base_cnv.R")
    subprocess.check_call(['Rscript', rScriptName, tmpdir, workdir])

    # now move the cnv results with the selected ploidy to the output file
    ploidy=open(os.path.join(workdir, "ploidy")).readline().strip()
    shutil.copyfile(os.path.join(tmpdir, "cnv.result" + ploidy), os.path.join(workdir, "cnv.result"))         

def main():

    subprocess.call("date")
    
    options = Options()
    
    control = options.control
    tumor = options.tumor
    targets = options.bed
    bafIn = options.bafin

    # These folders are local to the docker container
    workdir = os.path.join('/data', outFolder)
    tmpdir = os.path.join(workdir, 'tmp')
    mkdir_p(tmpdir)
    
    chromstring= getChroms(targets)
        
    print 'Sorting input files...'
    sortedcontrol = os.path.join(tmpdir, 'control.coverage.sorted')
    sortedtumor   = os.path.join(tmpdir, 'tumor.coverage.sorted')
    sortFile(os.path.abspath(control), sortedcontrol)
    sortFile(os.path.abspath(tumor), sortedtumor)

    print 'Generating mean coverage files...'
    coveredcontrol = os.path.join(tmpdir, 'control.coverage')
    coveredtumor = os.path.join(tmpdir, 'tumor.coverage')
    getMeanCoverage(sortedcontrol, coveredcontrol)
    getMeanCoverage(sortedtumor, coveredtumor)

    print 'Creating segment ratios...'
    ratiofile = os.path.join(tmpdir, 'ratio.data')
    snpSegfile = os.path.join(tmpdir, 'snp_segments')
    segmentRatio(options, coveredcontrol, coveredtumor, workdir, chromstring, ratiofile, snpSegfile)
    
    print 'Analyzing CNV...'
    cnvFiles = analyseCNV(options, ratiofile, workdir, tmpdir, chromstring)

    if(str(options.p_est)=="True"):

        print "Estimating base ploidy..."
        estimatePloidy(tmpdir, workdir, snpSegfile)
    	
    if str(options.plot)=="True":
        #Runs plot_results.R, which expects a file named cnv.results in the input directory
	#Creates one plot per chromosome.

        rScriptName = os.path.join(scriptPath,"plot_results.R")
        subprocess.check_call(['Rscript', rScriptName, workdir, chromstring])

    # if baf is provided but ploidy does not need to be estimated
    if(str(bafIn) =="True"):
	print "Predicting Zygosity states"
        zygosity(options, workdir, cnvFiles[0], chromstring) 

    # Delete temporary directory	
    if not options.keeptemp:
        shutil.rmtree(tmpdir)
    subprocess.call("date",shell=True)
	

if __name__ == "__main__":
	main()

