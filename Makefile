# Definitions
curpath = $(shell pwd)
runtime_fullpath = ${curpath}/runtime
build_output = runtime/exec/ADTEx.py
build_tool = runtime-container.DONE
nametag = jeltje/adtex

sources = test/normal.cov test/tumor.cov

# Steps
all: ${build_output} ${build_tool}

${build_output}: build/Dockerfile
	#cd build && docker build -t adtexbuild .
	docker run -v ${runtime_fullpath}:/data adtexbuild cp -rp exec /data

${build_tool}: ${build_output} runtime/Dockerfile
	cd runtime && docker build -t ${nametag} .
	docker rmi -f adtexbuild
	rm -rf runtime/exec
	touch ${build_tool}

${sources}: extract
extract: test/normal.cov.gz
	gunzip test/normal.cov.gz test/tumor.cov.gz

test: test/normal.cov
	docker run --rm -v ${curpath}/test:/data ${nametag} --normal normal.cov --tumor tumor.cov --targetbed targets.bed --centromeres centromeres.bed --baf input.baf --estimatePloidy --out test_out --sampleid testSample
	diff test/test_out/testSample.cnv expected_output/expected.cnv
	

