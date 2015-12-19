#!/bin/bash

for i in `ls -l | awk '/^d/{print $9}'`
do
if [ -f "./$i/value_start.m" ]; then
echo "run ./$i/value_start.m;" >> data_initilize.m 
fi
done

matlab -nojvm -nodisplay -nosplash -nodesktop <data_initilize.m

rm data_initilize.m
