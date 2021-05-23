#!/bin/bash

for iter in {0..9}
do
	echo "$iter"
	python3 ../alg/main.py # --sr_cnt $sr_cnt
	matlab -nodisplay -nosplash - nodesktop -r "run /home/xiaofan/Github/robust-lora/opt/optSolver.m; exit;"
done
