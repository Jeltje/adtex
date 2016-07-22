FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
	wget \
	g++ \
	make \
        dpkg-dev \
	zlib1g-dev  

# python3 is the default
RUN ln -s /usr/bin/python3 /usr/bin/python

# Create a new source directory
WORKDIR /home
RUN mkdir adtex-src
RUN mkdir exec
WORKDIR /home/adtex-src

# get adtex 
RUN wget http://downloads.sourceforge.net/project/adtex/ADTEx.v2.0/ADTEx.v.2.0.tar.gz
RUN tar xf ADTEx.v.2.0.tar.gz
RUN mv ADTEx.v.2.0/* /home/exec

# updates made some packages incompatible. Use older versions from jeltje/adtex git repo
# wmtsa dependencies. Note: use archived version of MASS to be compatible with R 3.0.2
#RUN wget -O /home/exec/splus2R_1.2-0.tar.gz https://cran.r-project.org/src/contrib/splus2R_1.2-0.tar.gz
#RUN wget -O /home/exec/MASS_7.3-36.tar.gz https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-36.tar.gz
#RUN wget -O /home/exec/ifultools_2.0-1.tar.gz https://cran.r-project.org/src/contrib/ifultools_2.0-1.tar.gz
# DNAcopy and wmtsa package 
# NOTE: this is NOT a stable way of doing this, bioconductor updates this package regularly
#RUN wget -O /home/exec/DNAcopy_1.44.0.tar.gz http://www.bioconductor.org/packages/release/bioc/src/contrib/DNAcopy_1.44.0.tar.gz
#RUN wget -O /home/exec/wmtsa_2.0-0.tar.gz https://cran.r-project.org/src/contrib/wmtsa_2.0-0.tar.gz

# bedtools - use the last bedtools, not any of the bedtools2 (numbering v2.18.0 and up) because the output
# format is different
RUN wget https://github.com/arq5x/bedtools/archive/v2.17.0.tar.gz && \
    tar xf v2.17.0.tar.gz
WORKDIR /home/adtex-src/bedtools-2.17.0
RUN make 
RUN mv bin/bedtools /home/exec

# change ownership of outputfiles so they can be deleted
RUN chmod -R 777 /home/exec

# Set working directory
WORKDIR /home
 
# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /home/adtex-src

