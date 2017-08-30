#!/bin/bash -f

## run on Data samples

for isHF in {1,0}
do
	for sample in {-100,-200,-300,2500,2300,2310,2514,2515,2600}  #Data & MC
	do
		root -b -q 'csvReweight.C+('$isHF','$sample')'
	done

done
echo "Finished producing all histograms."

