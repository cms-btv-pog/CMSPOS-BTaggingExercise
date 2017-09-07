#!/bin/bash -f

## run on Data samples

mkdir -p logs

for isHF in {1,0}; do
    for sample in {-100,-200,-300,2500,2300,2310,2514,2515,2600}  #Data & MC
    do
	root -b -q 'csvReweight.C+('$isHF','$sample')' > logs/log_$sample\_$isHF.log &
    done
    wait
done

wait 

echo "Finished producing all histograms."

