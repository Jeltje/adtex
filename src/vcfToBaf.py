#!/usr/bin/env python

import sys, os, re, argparse

parser = argparse.ArgumentParser(description="Create BAF file as needed by ADTEx zygosity from input VCF")
parser.add_argument('-i', '--input', required=True, type=str, help="Extract B allele frequencies from input VCF")
parser.add_argument('-c', '--control', required=True, type=str, help="ID of control sample as used in the VCF header")
parser.add_argument('-t', '--tumor', required=True, type=str, help="ID of tumor sample as used in the VCF header")

# From the MuTect header:
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=BQ,Number=A,Type=Float,Description="Average base quality for reads supporting alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# So:
# The B allele frequency is present under the FA tag in the MuTect VCF file. This number is based on reads covering
# the ref and alt alleles under the AD tag, with the depth of coverage under DP (this number is the same as
# ref+alt)
# We want to print chrom, pos, control FA, tumor FA, control DP, tumor DP

def getColumn(header, id):
    """Determine which columns holds the sample"""
    try:
        tpos = header.index(id)
    except ValueError:
        print >>sys.stderr, "ERROR: cannot find {} in VCF header".format([id])
        sys.exit()
    return tpos

def getBaf(formatstring, valuestring):
    """
    Extract location of FA and DP tags in the fcol string, such as GT:AD:BQ:DP:FA:SS
    And return the corresponding values from the valuestring
    """
    ffields = formatstring.split(':')
    try:
        afloc = ffields.index('FA')
        dploc = ffields.index('DP')
    except ValueError:
        print >>sys.stderr, 'ERROR: cannot find DP or FA in the FORMAT column'
        sys.exit()
    sfields = valuestring.split(':')
    return sfields[afloc], sfields[dploc]
    
# Main
args = parser.parse_args()
print 'chrom\tSNP_loc\tcontrol_BAF\ttumor_BAF\tcontrol_doc\ttumor_doc'
with open(args.input, 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue

        fields = line.strip().split('\t')
        if line.startswith('#'):
            # determine tumor and control column
            tpos = getColumn(fields, args.tumor)
            cpos = getColumn(fields, args.control)
            # this should be column 8 but check anyway
            fcol = getColumn(fields, 'FORMAT')
            continue


        # now we should be in the main body
        # do NOT get the somatic mutations
        if fields[6] == "PASS":
            continue
    
        cbaf, cdoc = getBaf(fields[fcol], fields[cpos])
        # if the allele frequency is 0 or 1 in the control, don't bother
        if(cbaf == '1.00' or cbaf == '0.00'):
            continue

        tbaf, tdoc = getBaf(fields[fcol], fields[tpos])
        print '{}\t{}\t{}\t{}\t{}\t{}'.format(fields[0], fields[1], cbaf, 
                  tbaf, cdoc, tdoc)





