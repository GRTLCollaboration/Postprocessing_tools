#!/bin/bash

sed -n '/GRAMRLevel::/,$p' ./pout.0 > ./pout.0_temp
sed -e 's!GRAMRLevel::advance level !!' ./pout.0_temp > ./pout.0_temp2
sed -e 's/\(.\{1\}\)/\1,/' ./pout.0_temp2 > ./pout.0_temp
sed -e 's!at time !!' ./pout.0_temp > ./pout.0_temp2
sed -e 's/ (.*:/,/' ./pout.0_temp2 > ./pout.0_temp
sed -e 's/ \//,/' ./pout.0_temp > ./pout.0_parsed
rm ./pout.0_temp ./pout.0_temp2
