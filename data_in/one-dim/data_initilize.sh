#!/bin/bash

rm data_initilize.m 

for i in `ls -l | awk '/^d/{print $9}'`
do
echo "run ./$i/value_start.m;" >> data_initilize.m 
done

matlab -nojvm -nodisplay -nosplash -nodesktop <data_initilize.m
