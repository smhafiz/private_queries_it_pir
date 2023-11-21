# This file is part of Private aggregate queries v0.1.
# 
# 
# Copyright (C) 2022, Warren Wnuck and Syed Mahbub Hafiz.
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
echo "START" > "Fig5.txt"
#gpu
r=$((18))
p=$((10))
u=$((4))
q=$((8))
t=$((1))
files_min=$((1))
files=$((10))
prev=$((0))

cd ../../private_queries
FILE_NAME="../aggregate_queries/benchmarking/Fig5.txt"
for ((f=$files_min; f<=$files;f+=1));
do
	echo "start"
    
    for((i=0;i<2**$f;i+=1));
    do
        ./ccs $i $p $r $((2**$u))
        string+=$i" "
    done

    echo "files: $string" >> $FILE_NAME
    for ((v =0; v<$t; v+=1));
    do				
        ./createindex $((2**$q/8)) indexes.txt $string >> $FILE_NAME
        
    done
    mv -v val.22 "val_$i.22"

    for((k=0;k<2**$f;k+=1));
    do
        file1=$k".row"
        file2=$k".col"
        rm $file1
        rm $file2
    done
    unset string
    echo "removed $(($k+1)) files"
    mv -v out.col "out_$i.col"
    echo "done with set"
done
cd ../aggregate_queries/benchmarking/
echo "Finished Script"
