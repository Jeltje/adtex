# Definitions
build_path = build/
runtime_path = runtime/
build_output = ${runtime_path}/exec/ADTEx.py
runtime_fullpath = $(realpath ${runtime_path})
build_tool = runtime-container.DONE
nametag = jeltje/adtex

# Steps
all: ${build_output} ${build_tool}

${build_output}: ${build_path}/Dockerfile
	cd ${build_path} && docker build -t adtexbuild .
	docker run -v ${runtime_fullpath}:/data adtexbuild cp -rp exec /data

${build_tool}: ${build_output} ${runtime_path}/Dockerfile
	cd ${runtime_path} && docker build -t ${nametag} .
	docker rmi -f adtexbuild
	touch ${build_tool}

#push: all
#	# Requires ~/.dockercfg
#	docker push ${nametag}
