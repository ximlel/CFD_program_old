#!/bin/bash  --login
shopt -s expand_aliases
source ~/.bash_aliases

for i in `ls -l | awk '/^d/{print $9}'`
do cd $i
preplot RHO.tec
preplot U.tec
preplot P.tec
cd .. 
done
