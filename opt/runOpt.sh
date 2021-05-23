#!/bin/bash

for sr_cnt in 100 200 500 1000 2000
do
	echo "$sr_cnt"
	python3 ../alg/main.py --sr_cnt $sr_cnt
	#matlab -nodisplay -nosplash - nodesktop -r "run /home/xiaofan/Github/robust-lora/opt/optSolver.m; exit;"
done
