#!/bin/bash

# if [ "$(uname -m)" = "aarch64" -o "$(uname -m)" = "arm64" ]; then
#     docker build -t sa-eds . \
#     --build-arg NUPACK_PK="nupack-4.0.1.8/package/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_aarch64.manylinux2014_aarch64.whl"
# else
#     docker build -t sa-eds . \
#     --build-arg NUPACK_PK="nupack-4.0.1.8/package/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
# fi

currentdir=$(pwd)
basename=`basename $currentdir`
# docker run  --rm -it -v "${currentdir}:/home/user/${basename}" sa-eds /home/user/venv/bin/python3 /home/user/${basename}/${1}
docker run  --rm -v "${currentdir}:/home/user/${basename}" sa-eds /home/user/venv/bin/python3 /home/user/${basename}/${1} "${2}" "${3}" "${4}"
# docker run  --rm -v "${currentdir}:/home/user/${basename}" sa-eds /home/user/venv/bin/python3 /home/user/${basename}/${1} "${2}" "${3}"
# docker run  --rm -it -v "${currentdir}:/home/user/${basename}" sa-eds /bin/bash
# docker run --rm -it -v "${currentdir}:/home/user/SA-EDS" -w /home/user/SA-EDS/oxDNA/build ubuntu:latest /bin/bash -c "cmake .. && make -j4"
# docker run --rm -it -v "${currentdir}:/home/user/SA-EDS" -w /home/user/SA-EDS/oxDNA/ ubuntu:latest /bin/bash -c "apt-get update && apt-get install -y cmake make && mkdir build && cd build && echo pwd && cmake .. && make -j4"


