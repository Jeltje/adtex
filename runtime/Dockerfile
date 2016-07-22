FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
	r-base \
	vim

# DNAcopy and wmtsa package (move to R library location)
# Note: cannot use ADD because this untars the file
COPY DNAcopy_1.44.0.tar.gz /tmp/
RUN R CMD INSTALL /tmp/DNAcopy_1.44.0.tar.gz 
COPY splus2R_1.2-0.tar.gz /tmp/
RUN R CMD INSTALL /tmp/splus2R_1.2-0.tar.gz
COPY MASS_7.3-36.tar.gz /tmp/
RUN R CMD INSTALL /tmp/MASS_7.3-36.tar.gz
COPY ifultools_2.0-1.tar.gz /tmp/
RUN R CMD INSTALL /tmp/ifultools_2.0-1.tar.gz
COPY wmtsa_2.0-0.tar.gz /tmp/
RUN R CMD INSTALL /tmp/wmtsa_2.0-0.tar.gz

# bedtools
ADD exec/bedtools /usr/local/bin/

# Adtex suite 
ADD exec/RFunction.R /usr/local/bin/

# Slightly modified versions of Adtex suite code
ADD base_cnv.R /usr/local/bin/
ADD cnv_analyse.R /usr/local/bin/
ADD extract_cnv.R /usr/local/bin/
ADD plot_results.R /usr/local/bin/
ADD segment_ratio.R /usr/local/bin/
ADD zygosity.R /usr/local/bin/

# Heavily modified Adtex run script and post processing step
ADD ADTEx.py /usr/local/bin/
ADD basicDNAcopy.R /usr/local/bin/

# Boilerplate
RUN mkdir /opt/adtex
ADD wrapper.sh /opt/adtex/

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/adtex/wrapper.sh"]
CMD ["--help"]

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* 


