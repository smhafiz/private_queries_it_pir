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
string="mimic_filtered"
echo "START" > "BatchRes_$string.txt"
Q_MIN=$((7))
Q_MAX=$((7))
T_MAX=$((1))
files="mimic_new_query_5 mimic_new_query_6"
cp ../../../../private_queries/createindex ../ccs_files/
cp ../../../../private_queries/indexes.txt ../ccs_files/
cd ../ccs_files/
FILE_NAME="../bash_scripts/BatchRes_$string.txt"
for ((q=$Q_MIN;q<=$Q_MAX;q+=1));
do
    echo "start bash"
    echo "making batch"
    echo "q=2^$q in miliseconds:" >> $FILE_NAME
    for((t=0;t<$T_MAX;t+=1));
    do
        ./createindex $((2**$q/8)) indexes.txt $files >> $FILE_NAME
    done
    echo "trials done for 2^$q"
    mv -v out.col out_$string.col
    mv -v out.row out_$string.row
done
rm createindex
rm indexes.txt
cd ../bash_scripts
echo "Finished Script"
