#!/bin/bash

if [ "$(uname -m)" = "aarch64" -o "$(uname -m)" = "arm64" ]; then
    docker build -t sa-eds . \
    --build-arg NUPACK_PK="nupack-4.0.1.8/package/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_aarch64.manylinux2014_aarch64.whl"
else
    docker build -t sa-eds . \
    --build-arg NUPACK_PK="nupack-4.0.1.8/package/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
fi

currentdir=$(pwd)
basename=`basename $currentdir`
docker run  --rm -it -v "${currentdir}:/home/user/${basename}" sa-eds 
