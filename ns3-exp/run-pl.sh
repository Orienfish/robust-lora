#!/bin/bash

cp ./bag-pl/sr_RGreedy_1_264_0dp.txt ./sr_loc.txt
cp ./bag-pl/gw_RGreedy_1_264_0dp.txt ./gw_loc.txt
cp ./bag-pl/pl_RGreedy_1_264_0dp.txt ./pl_mat.txt
destlog=logPL-GP
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataPL-GP.txt
mv ./nodeEE.txt ./res/nodeEEPL-GP.txt
mv ./globalPerformance.txt ./res/globalPerformancePL-GP.txt
mv ./phyPerformance.txt ./res/phyPerformancePL-GP.txt

cp ./bag-pl/sr_RGreedy_1_264_0d.txt ./sr_loc.txt
cp ./bag-pl/gw_RGreedy_1_264_0d.txt ./gw_loc.txt
cp ./bag-pl/pl_RGreedy_1_264_0d.txt ./pl_mat.txt
destlog=logPL-G
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataPL-G.txt
mv ./nodeEE.txt ./res/nodeEEPL-G.txt
mv ./globalPerformance.txt ./res/globalPerformancePL-G.txt
mv ./phyPerformance.txt ./res/phyPerformancePL-G.txt

cp ./bag-pl/sr_ICIOT_7_264dp.txt ./sr_loc.txt
cp ./bag-pl/gw_ICIOT_7_264dp.txt ./gw_loc.txt
cp ./bag-pl/pl_ICIOT_7_264dp.txt ./pl_mat.txt
destlog=logPL-IP
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataPL-IP.txt
mv ./nodeEE.txt ./res/nodeEEPL-IP.txt
mv ./globalPerformance.txt ./res/globalPerformancePL-IP.txt
mv ./phyPerformance.txt ./res/phyPerformancePL-IP.txt

cp ./bag-pl/sr_ICIOT_8_264d.txt ./sr_loc.txt
cp ./bag-pl/gw_ICIOT_8_264d.txt ./gw_loc.txt
cp ./bag-pl/pl_ICIOT_8_264d.txt ./pl_mat.txt
destlog=logPL-I
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataPL-I.txt
mv ./nodeEE.txt ./res/nodeEEPL-I.txt
mv ./globalPerformance.txt ./res/globalPerformancePL-I.txt
mv ./phyPerformance.txt ./res/phyPerformancePL-I.txt