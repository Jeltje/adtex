FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
	r-base 

# DNAcopy and wmtsa package (move to R library location)
# Note: cannot use ADD because this untars the file
COPY exec/DNAcopy_1.44.0.tar.gz /tmp/
RUN R CMD INSTALL /tmp/DNAcopy_1.44.0.tar.gz 
COPY exec/splus2R_1.2-0.tar.gz /tmp/
RUN R CMD INSTALL /tmp/splus2R_1.2-0.tar.gz
COPY exec/MASS_7.3-36.tar.gz /tmp/
RUN R CMD INSTALL /tmp/MASS_7.3-36.tar.gz
COPY exec/ifultools_2.0-1.tar.gz /tmp/
RUN R CMD INSTALL /tmp/ifultools_2.0-1.tar.gz
COPY exec/wmtsa_2.0-0.tar.gz /tmp/
RUN R CMD INSTALL /tmp/wmtsa_2.0-0.tar.gz

# bedtools
ADD exec/bedtools /usr/local/bin/
ADD exec/coverageBed /usr/local/bin/
ADD exec/intersectBed /usr/local/bin/

# adtex (add altered version instead of download)
ADD ADTEx.py /usr/local/bin/
ADD exec/base_cnv.R /usr/local/bin/
ADD exec/cnv_analyse.R /usr/local/bin/
ADD exec/extract_cnv.R /usr/local/bin/
ADD exec/plot_results.R /usr/local/bin/
ADD exec/RFunction.R /usr/local/bin/
ADD exec/segment_ratio.R /usr/local/bin/
ADD exec/zygosity.R /usr/local/bin/
ADD exec/getMeanCoverage.py /usr/local/bin/


# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* 

ENTRYPOINT ["python", "/usr/local/bin/ADTEx.py"]

