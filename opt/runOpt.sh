#!/bin/bash

for sr_cnt in 1 2 3 4 5 6 7 8 9 10
do
	echo "$sr_cnt"
	python3 ../alg/main.py #--sr_cnt $sr_cnt
	matlab -nodisplay -nosplash - nodesktop -r "run /home/xiaofan/Github/robust-lora/opt/optSolver.m; exit;"
done
