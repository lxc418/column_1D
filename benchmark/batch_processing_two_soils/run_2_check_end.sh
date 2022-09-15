#!/bin/sh
## 

# used to check if the simulation terminated due to error
for WATERDEPTH  in {0..30}
do
		ddn=depth_"$WATERDEPTH"
		touch  check_simu_end.csv
		cd    $ddn
        tail -2 COLUMN.smy | head -n 1 >> ../check_simu_end.csv

		cd ..
done
