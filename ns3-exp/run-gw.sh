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
destlog=logiciot6
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot6.txt
mv ./nodeEE.txt ./res/nodeEEiciot6.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot6.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot6.txt

destlog=logiciot61
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=1 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot61.txt
mv ./nodeEE.txt ./res/nodeEEiciot61.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot61.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot61.txt

destlog=logiciot62
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=2 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot62.txt
mv ./nodeEE.txt ./res/nodeEEiciot62.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot62.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot62.txt

cp ~/Github/robust-lora/res/sr_ICIOT_9_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_ICIOT_9_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_ICIOT_9_264_0dp.txt ./pl_mat.txt
destlog=logiciot9
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot9.txt
mv ./nodeEE.txt ./res/nodeEEiciot9.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot9.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot9.txt

destlog=logiciot91
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=1 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot91.txt
mv ./nodeEE.txt ./res/nodeEEiciot91.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot91.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot91.txt

destlog=logiciot92
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=2 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot92.txt
mv ./nodeEE.txt ./res/nodeEEiciot92.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot92.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot92.txt

cp ~/Github/robust-lora/res/sr_ICIOT_12_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/res/gw_ICIOT_12_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/res/pl_ICIOT_12_264_0dp.txt ./pl_mat.txt
destlog=logiciot12
start=`date +%s`
./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot12.txt
mv ./nodeEE.txt ./res/nodeEEiciot12.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot12.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot12.txt

destlog=logiciot121
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=1 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot121.txt
mv ./nodeEE.txt ./res/nodeEEiciot121.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot121.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot121.txt

destlog=logiciot122
start=`date +%s`
./waf --run "adr --MType=Confirmed --nDownGateways=2 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot122.txt
mv ./nodeEE.txt ./res/nodeEEiciot122.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot122.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot122.txt