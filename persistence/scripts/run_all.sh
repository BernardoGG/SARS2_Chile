#!/bin/bash

cd ../

for i in $(ls args_files/); 
do 
	echo Running $i; 
	time grep -v "last_sampling" args_files/${i} | xargs -n5 -P 5 ./scripts/run_persistence_per_args.sh 
done
