#!/bin/sh
## 
for WATERDEPTH  in {0..30}
do
		ddn=depth_"$WATERDEPTH"
		cd    $ddn
		let "EXCELROW=WATERDEPTH+1" #let zero water depth store in the first row
		
		cp ../{sl_figure.m,sl_read.m,post_process.script} .
		sed -i s/FINDMEROW/$EXCELROW/g      sl_figure.m

  		sbatch post_process.script
		sleep 2
		cd ..
done
