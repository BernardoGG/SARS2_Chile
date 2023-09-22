#!/bin/bash

# In parent folder execute: grep -v "last_sampling"  args_files/Gamma_TL_55_args.tsv | xargs -n5 -P 5 ./scripts/run_persistence_per_args.sh 

lineage=$1
comuna=$3
comuna="${comuna// /_}"
trees_folder=tree_files/${lineage}_minimal
evalTimes=$4
ancTimes=$5
out_folder=output/${lineage}/${comuna}
log_file=logs/${lineage}/

mkdir -p $out_folder
mkdir -p $log_file

log_file=${log_file}/${comuna}.txt

echo "Running for " $lineage >> $log_file
echo "evalTimes" $evalTimes >> $log_file
echo "ancTimes" $ancTimes >> $log_file
echo "Output to " $out_folder >> $log_file

total_trees=$(ls ${trees_folder}/*.tree | wc -l)
ctr=0
for tree in $(ls ${trees_folder}/*.tree)
do
	echo "Running tree ${ctr} out of ${total_trees} trees" >> $log_file
	outname=$(basename $tree)
	# Modify path to BEAST jar file accordingly
	java -Xms64m -Xmx2048m -cp "$HOME/Documents/code/beast-mcmc/build/dist/beast.jar" dr.app.tools.PersistenceSummarizer -evaluationTime $evalTimes  -ancestralTime $ancTimes -nodeStateAnnotation comuna $tree ${out_folder}/${outname/.tree/.tsv} > /dev/null 2>&1
	ctr=$((ctr+1))
done

