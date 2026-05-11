#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -J woltka

### If you want functions run the following lines (without the #) in your wol2 folder
# wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/function/
# wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/genomes/
# wget -r -np -nH --cut-dirs=2 http://ftp.microbio.me/pub/wol2/proteins/

#if not turn on --no-fun
#If you don't want to have the function please run the proteins line otherwise script will fail or hashtag out part of wol2sop.sh where .orf is created

#$1 = directory of Web of Life 2 database
#$2 = directory of QC'd data 
#./wol2sop.sh -d ./ -i ./align --no-tax --no-fun 
./wol2sop.sh -d ./wol2 -i align -o Output -f biom --no-tax --no-fun 
