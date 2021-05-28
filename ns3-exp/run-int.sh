#!/bin/bash

cp ~/Github/robust-lora/alg/res/sr_RGreedy_1_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/alg/res/gw_RGreedy_1_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/alg/res/pl_RGreedy_1_264_0dp.txt ./pl_mat.txt
#destlog=log1
#start=`date +%s`
#./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
#end=`date +%s`
#runtime=$((end-start))
#echo "$runtime" >> "$destlog"
#mv ./nodeData.txt ./res/nodeData1.txt
#mv ./nodeEE.txt ./res/nodeEE1.txt
#mv ./globalPerformance.txt ./res/globalPerformance1.txt
#mv ./phyPerformance.txt ./res/phyPerformance1.txt

destlog=log1-130
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-130 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData11.txt
mv ./nodeEE.txt ./res/nodeEE11.txt
mv ./globalPerformance.txt ./res/globalPerformance11.txt
mv ./phyPerformance.txt ./res/phyPerformance11.txt

destlog=log1-124
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-124 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData12.txt
mv ./nodeEE.txt ./res/nodeEE12.txt
mv ./globalPerformance.txt ./res/globalPerformance12.txt
mv ./phyPerformance.txt ./res/phyPerformance12.txt

cp ~/Github/robust-lora/alg/res/sr_RGreedy_2_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/alg/res/gw_RGreedy_2_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/alg/res/pl_RGreedy_2_264_0dp.txt ./pl_mat.txt
#destlog=log2
#start=`date +%s`
#./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
#end=`date +%s`
#runtime=$((end-start))
#echo "$runtime" >> "$destlog"
#mv ./nodeData.txt ./res/nodeData2.txt
#mv ./nodeEE.txt ./res/nodeEE2.txt
#mv ./globalPerformance.txt ./res/globalPerformance2.txt
#mv ./phyPerformance.txt ./res/phyPerformance2.txt

destlog=log2-130
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-130 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData21.txt
mv ./nodeEE.txt ./res/nodeEE21.txt
mv ./globalPerformance.txt ./res/globalPerformance21.txt
mv ./phyPerformance.txt ./res/phyPerformance21.txt

destlog=log2-124
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-124 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData22.txt
mv ./nodeEE.txt ./res/nodeEE22.txt
mv ./globalPerformance.txt ./res/globalPerformance22.txt
mv ./phyPerformance.txt ./res/phyPerformance22.txt

cp ~/Github/robust-lora/alg/res/sr_RGreedy_3_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/alg/res/gw_RGreedy_3_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/alg/res/pl_RGreedy_3_264_0dp.txt ./pl_mat.txt
#destlog=log3
#start=`date +%s`
#./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
#end=`date +%s`
#runtime=$((end-start))
#echo "$runtime" >> "$destlog"
#mv ./nodeData.txt ./res/nodeData3.txt
#mv ./nodeEE.txt ./res/nodeEE3.txt
#mv ./globalPerformance.txt ./res/globalPerformance3.txt
#mv ./phyPerformance.txt ./res/phyPerformance3.txt

destlog=log3-130
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-130 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData31.txt
mv ./nodeEE.txt ./res/nodeEE31.txt
mv ./globalPerformance.txt ./res/globalPerformance31.txt
mv ./phyPerformance.txt ./res/phyPerformance31.txt

destlog=log3-124
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-124 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeData32.txt
mv ./nodeEE.txt ./res/nodeEE32.txt
mv ./globalPerformance.txt ./res/globalPerformance32.txt
mv ./phyPerformance.txt ./res/phyPerformance32.txt

cp ~/Github/robust-lora/alg/res/sr_ICIOT_1_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/alg/res/gw_ICIOT_1_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/alg/res/pl_ICIOT_1_264_0dp.txt ./pl_mat.txt
#destlog=logiciot
#start=`date +%s`
#./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
#end=`date +%s`
#runtime=$((end-start))
#echo "$runtime" >> "$destlog"
#mv ./nodeData.txt ./res/nodeDataiciot.txt
#mv ./nodeEE.txt ./res/nodeEEiciot.txt
#mv ./globalPerformance.txt ./res/globalPerformanceiciot.txt
#mv ./phyPerformance.txt ./res/phyPerformanceiciot.txt

destlog=logiciot6-130
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-130 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot6-130.txt
mv ./nodeEE.txt ./res/nodeEEiciot6-130.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot6-130.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot6-130.txt

destlog=logiciot6-124
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-124 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot6-124.txt
mv ./nodeEE.txt ./res/nodeEEiciot6-124.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot6-124.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot6-124.txt

cp ~/Github/robust-lora/alg/res/sr_ICIOT_2_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/alg/res/gw_ICIOT_2_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/alg/res/pl_ICIOT_2_264_0dp.txt ./pl_mat.txt
#destlog=logiciot9
#start=`date +%s`
#./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
#end=`date +%s`
#runtime=$((end-start))
#echo "$runtime" >> "$destlog"
#mv ./nodeData.txt ./res/nodeDataiciot9.txt
#mv ./nodeEE.txt ./res/nodeEEiciot9.txt
#mv ./globalPerformance.txt ./res/globalPerformanceiciot9.txt
#mv ./phyPerformance.txt ./res/phyPerformanceiciot9.txt

destlog=logiciot9-130
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-130 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot9-130.txt
mv ./nodeEE.txt ./res/nodeEEiciot9-130.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot9-130.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot9-130.txt

destlog=logiciot9-124
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-124 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot9-124.txt
mv ./nodeEE.txt ./res/nodeEEiciot9-124.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot9-124.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot9-124.txt

cp ~/Github/robust-lora/alg/res/sr_ICIOT_3_264_0dp.txt ./sr_loc.txt
cp ~/Github/robust-lora/alg/res/gw_ICIOT_3_264_0dp.txt ./gw_loc.txt
cp ~/Github/robust-lora/alg/res/pl_ICIOT_3_264_0dp.txt ./pl_mat.txt
#destlog=logiciot12
#start=`date +%s`
#./waf --run "adr --MType=Confirmed --MatrixPLEnabled" > "$destlog" 2>&1
#end=`date +%s`
#runtime=$((end-start))
#echo "$runtime" >> "$destlog"
#mv ./nodeData.txt ./res/nodeDataiciot12.txt
#mv ./nodeEE.txt ./res/nodeEEiciot12.txt
#mv ./globalPerformance.txt ./res/globalPerformanceiciot12.txt
#mv ./phyPerformance.txt ./res/phyPerformanceiciot12.txt

destlog=logiciot12-130
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-130 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot12-130.txt
mv ./nodeEE.txt ./res/nodeEEiciot12-130.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot12-130.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot12-130.txt

destlog=logiciot12-124
start=`date +%s`
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-124 --MatrixPLEnabled" > "$destlog" 2>&1
end=`date +%s`
runtime=$((end-start))
echo "$runtime" >> "$destlog"
mv ./nodeData.txt ./res/nodeDataiciot12-124.txt
mv ./nodeEE.txt ./res/nodeEEiciot12-124.txt
mv ./globalPerformance.txt ./res/globalPerformanceiciot12-124.txt
mv ./phyPerformance.txt ./res/phyPerformanceiciot12-124.txt