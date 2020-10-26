#!/bin/bash

cp ~/Github/robust-lora/bag/sr_RGreedy_1_100_0.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_1_100_0.txt ./gw_loc.txt
destlog=logG
start=`date +%s`
./waf --run "adr --MType=Confirmed" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataG.txt
mv ./nodeEE.txt ./res/nodeEEG.txt
mv ./globalPerformance.txt ./res/globalPerformanceG.txt
mv ./phyPerformance.txt ./res/phyPerformanceG.txt

cp ~/Github/robust-lora/res/sr_ICIOT_100_0.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_ICIOT_100_0.txt ./gw_loc.txt
destlog=logI
start=`date +%s`
./waf --run "adr --MType=Confirmed" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataI.txt
mv ./nodeEE.txt ./res/nodeEEI.txt
mv ./globalPerformance.txt ./res/globalPerformanceI.txt
mv ./phyPerformance.txt ./res/phyPerformanceI.txt

cp ~/Github/robust-lora/res/sr_relaxOpt.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_relaxOpt.txt ./gw_loc.txt
destlog=logO
start=`date +%s`
./waf --run "adr --MType=Confirmed" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataO.txt
mv ./nodeEE.txt ./res/nodeEEO.txt
mv ./globalPerformance.txt ./res/globalPerformanceO.txt
mv ./phyPerformance.txt ./res/phyPerformanceO.txt