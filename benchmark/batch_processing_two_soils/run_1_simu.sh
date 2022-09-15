#!/bin/sh
## 
for WATERDEPTH  in {0..30}
do

		ddn=depth_"$WATERDEPTH"
		mkdir $ddn
		cd    $ddn

		cp ../{sl_write_input.m,mesh.m,test.script} .
		sed -i s/FINDMENAME/$WATERDEPTH/g      sl_write_input.m
		sbatch test.script
		sleep 5
		cd ..
  
done


