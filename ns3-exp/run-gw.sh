#!/bin/bash

cp ~/Github/robust-lora/res/sr_RGreedy_1_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_1_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_1_264_0dp.txt ./pl_mat.txt
destlog=log1
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData1.txt
mv ./nodeEE.txt ./res/nodeEE1.txt
mv ./globalPerformance.txt ./res/globalPerformance1.txt
mv ./phyPerformance.txt ./res/phyPerformance1.txt

cp ~/Github/robust-lora/res/sr_RGreedy_2_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_2_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_2_264_0dp.txt ./pl_mat.txt
destlog=log2
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData2.txt
mv ./nodeEE.txt ./res/nodeEE2.txt
mv ./globalPerformance.txt ./res/globalPerformance2.txt
mv ./phyPerformance.txt ./res/phyPerformance2.txt

cp ~/Github/robust-lora/res/sr_RGreedy_3_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_3_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_3_264_0dp.txt ./pl_mat.txt
destlog=log3
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData3.txt
mv ./nodeEE.txt ./res/nodeEE3.txt
mv ./globalPerformance.txt ./res/globalPerformance3.txt
mv ./phyPerformance.txt ./res/phyPerformance3.txt

cp ~/Github/robust-lora/res/sr_RGreedy_1_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_1_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_1_264_0dp.txt ./pl_mat.txt
destlog=log11
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=1 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData11.txt
mv ./nodeEE.txt ./res/nodeEE11.txt
mv ./globalPerformance.txt ./res/globalPerformance11.txt
mv ./phyPerformance.txt ./res/phyPerformance11.txt

cp ~/Github/robust-lora/res/sr_RGreedy_2_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_2_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_2_264_0dp.txt ./pl_mat.txt
destlog=log21
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=1 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData21.txt
mv ./nodeEE.txt ./res/nodeEE21.txt
mv ./globalPerformance.txt ./res/globalPerformance21.txt
mv ./phyPerformance.txt ./res/phyPerformance21.txt

cp ~/Github/robust-lora/res/sr_RGreedy_3_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_3_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_3_264_0dp.txt ./pl_mat.txt
destlog=log31
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=1 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData31.txt
mv ./nodeEE.txt ./res/nodeEE31.txt
mv ./globalPerformance.txt ./res/globalPerformance31.txt
mv ./phyPerformance.txt ./res/phyPerformance31.txt

cp ~/Github/robust-lora/res/sr_RGreedy_1_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_1_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_1_264_0dp.txt ./pl_mat.txt
destlog=log12
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=2 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData12.txt
mv ./nodeEE.txt ./res/nodeEE12.txt
mv ./globalPerformance.txt ./res/globalPerformance12.txt
mv ./phyPerformance.txt ./res/phyPerformance12.txt

cp ~/Github/robust-lora/res/sr_RGreedy_2_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_2_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_2_264_0dp.txt ./pl_mat.txt
destlog=log22
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=2 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData22.txt
mv ./nodeEE.txt ./res/nodeEE22.txt
mv ./globalPerformance.txt ./res/globalPerformance22.txt
mv ./phyPerformance.txt ./res/phyPerformance22.txt

cp ~/Github/robust-lora/res/sr_RGreedy_3_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_RGreedy_3_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_RGreedy_3_264_0dp.txt ./pl_mat.txt
destlog=log32
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=2 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData32.txt
mv ./nodeEE.txt ./res/nodeEE32.txt
mv ./globalPerformance.txt ./res/globalPerformance32.txt
mv ./phyPerformance.txt ./res/phyPerformance32.txt

cp ~/Github/robust-lora/res/sr_ICIOT_6_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_ICIOT_6_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_ICIOT_6_264_0dp.txt ./pl_mat.txt
destlog=logiciot
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot.txt
mv ./nodeEE.txt ./res/nodeEEiciot.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot.txt

cp ~/Github/robust-lora/res/sr_ICIOT_9_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_ICIOT_9_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_ICIOT_9_264_0dp.txt ./pl_mat.txt
destlog=logiciot1
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=1 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot1.txt
mv ./nodeEE.txt ./res/nodeEEiciot1.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot1.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot1.txt

cp ~/Github/robust-lora/res/sr_ICIOT_12_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_ICIOT_12_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_ICIOT_12_264_0dp.txt ./pl_mat.txt
destlog=logiciot2
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=2 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot2.txt
mv ./nodeEE.txt ./res/nodeEEiciot2.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot2.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot2.txt

