# This file is part of Private aggregate queries v0.1.
# 
# 
# Copyright (C) 2022, Warren Wnuck, Chitrabhanu Gupta, and Syed Mahbub Hafiz.
# 
# Private aggregate queries is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
# 
# Private aggregate queries is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Private aggregate queries. If not, see <http://www.gnu.org/licenses/>.

#!/bin/bash
DV=1
FILE="mimic_filtered"
echo "START" > "ResGPUThroughput_$FILE.txt"
#gpu
Q_MIN=$((7))
Q_MAX=$((7))
T_MAX=$((1))

cd ../../../../private_queries/
FILE_NAME="../../aggregate_queries/casestudies/mimic/bash_scripts/ResGPUThroughput_$FILE.txt"

for ((q=$Q_MIN; q<=$Q_MAX;q+=1));
do
    cd barrettCUDA

    ./genuint $((2**$q/32)) > uintX.h

    echo "making barrett"
    make barrett
    echo "after barrett"
    echo "q: 2^$q" >> $FILE_NAME
    for ((v =0; v<$T_MAX; v+=1));
    do				
        ./barrett ../../aggregate_queries/casestudies/mimic/ccs_files/val.22 ../../aggregate_queries/casestudies/mimic/ccs_files/out_$FILE.row ../../aggregate_queries/casestudies/mimic/ccs_files/out_$FILE.col >> $FILE_NAME	
    done
    cd ..
done
cd ../aggregate_queries/casestudies/mimic/bash_scripts/
echo "Finished Script"
