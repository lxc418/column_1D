#!/bin/sh
## 
for WATERDEPTH  in {0..30}
do
		ddn=depth_"$WATERDEPTH"
		cd    $ddn		
		cp ../{sl_salt_thickness.m,salt_thickness.script} .

		cd ..
done


