#!/bin/bash
docker build -t sa-eds . 
currentdir=$(pwd)
basename=`basename $currentdir`
docker run  --rm -it -v "${currentdir}:/home/user/${basename}" sa-eds 
