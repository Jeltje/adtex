FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
	wget \
	g++ \
	make \
        dpkg-dev \
	zlib1g-dev \ 
	r-base \
	vim


WORKDIR /home

# bedtools - use the last bedtools, not any of the bedtools2 (numbering v2.18.0 and up) because the output
# format is different
RUN wget https://github.com/arq5x/bedtools/archive/v2.17.0.tar.gz && \
    tar xf v2.17.0.tar.gz
WORKDIR /home/bedtools-2.17.0
RUN make 
RUN mv bin/bedtools /usr/local/bin/

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

# The code below is based on http://downloads.sourceforge.net/project/adtex/ADTEx.v2.0/ADTEx.v.2.0.tar.gz
# These are slightly modified versions of Adtex suite code
ADD RFunction.R /usr/local/bin/
ADD base_cnv.R /usr/local/bin/
ADD cnv_analyse.R /usr/local/bin/
ADD extract_cnv.R /usr/local/bin/
ADD plot_results.R /usr/local/bin/
ADD segment_ratio.R /usr/local/bin/
ADD zygosity.R /usr/local/bin/

# Heavily modified Adtex run script and post processing step
ADD ADTEx.py /usr/local/bin/
ADD basicDNAcopy.R /usr/local/bin/

ADD run_adtex /usr/local/bin/

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/usr/local/bin/run_adtex"]
CMD ["--help"]

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /home/bedtools-2.17.0


