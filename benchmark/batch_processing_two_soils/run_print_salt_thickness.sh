#!/bin/sh
## 
for WATERDEPTH  in {0..1}
do
		ddn=depth_"$WATERDEPTH"
		cd    $ddn
		let "EXCELROW=WATERDEPTH+1" #let zero water depth store in the first row
		
		cp ../{sl_salt_thickness.m,salt_thickness.script} .
		sed -i s/FINDMEROW/$EXCELROW/g      sl_salt_thickness.m

  		sbatch salt_thickness.script
		sleep 1
		cd ..
done
