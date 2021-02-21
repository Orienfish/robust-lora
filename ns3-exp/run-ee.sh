#!/bin/bash

for sr in 100 200 500 1000 2000
do
	for M in 1 2 3
	do
		cp "./bag-ee/sr_RGreedy_${M}_${sr}_0.txt" ./sr_loc.txt
		cp "./bag-ee/gw_RGreedy_${M}_${sr}_0.txt" ./gw_loc.txt
		cp "./bag-ee/pl_RGreedy_${M}_${sr}_0.txt" ./pl_mat.txt
		destlog=log
		start=`date +%s`
		./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
		end=`date +%s`
		runtime=$((end-start))
		echo "$runtime" >> "$destlog"
		#mv ./nodeData.txt ./res/nodeDataPL-GP.txt
		mv ./nodeEE.txt "./res/nodeEE_${M}_${sr}.txt"
		#mv ./globalPerformance.txt ./res/globalPerformancePL-GP.txt
		#mv ./phyPerformance.txt ./res/phyPerformancePL-GP.txt
	done
done

for sr in 100 200 500 1000 2000
do
	cp "./bag-ee/sr_ICIOT_7_${sr}_0.txt" ./sr_loc.txt
	cp "./bag-ee/gw_ICIOT_7_${sr}_0.txt" ./gw_loc.txt
	cp "./bag-ee/pl_ICIOT_7_${sr}_0.txt" ./pl_mat.txt
	destlog=log
	start=`date +%s`
	./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
	end=`date +%s`
	runtime=$((end-start))
	echo "$runtime" >> "$destlog"
	#mv ./nodeData.txt ./res/nodeDataPL-GP.txt
	mv ./nodeEE.txt "./res/nodeEE_ICIOT_${sr}.txt"
	#mv ./globalPerformance.txt ./res/globalPerformancePL-GP.txt
	#mv ./phyPerformance.txt ./res/phyPerformancePL-GP.txt
done