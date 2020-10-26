#!/bin/bash

for sr_cnt in 1 2 3 4 5
do
	echo "$sr_cnt"
	python3 main.py # --sr_cnt $sr_cnt
	/Applications/MATLAB_R2020b.app/bin/matlab -nodisplay -nosplash - nodesktop -r "run /Users/xiaofanyu/Github/robust-lora/relaxOpt/relaxedOpt.m; exit;"
done