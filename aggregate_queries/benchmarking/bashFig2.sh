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
echo "START" > "Fig2.txt"
#gpu
R_MIN=$((16))
R_MAX=$((24))
p=$((16))
Q_MIN=$((7))
Q_MAX=$((9))
T_MAX=$((2))
u=$((2))
cd ../../private_queries
FILE_NAME="../../aggregate_queries/benchmarking/Fig2.txt"
for ((r=$R_MIN;r<=$R_MAX;r+=2));
do
    echo "start bash"
    ./ccs a $p $r $((2**$u))
    ./ccs b $p $r $((2**$u))		
    ./ccs c $p $r $((2**$u))
    ./ccs d $p $r $((2**$u))
    ./ccs e $p $r $((2**$u))
    ./ccs f $p $r $((2**$u))
    echo "created files for $r"
	
	for ((q=$Q_MIN; q<=$Q_MAX;q+=1));
	do
		cd barrettCUDA

		./genuint $((2**$q/32)) > uintX.h

		cd ..
		echo "making batch"
		./createindex $((2**$q/8)) indexes.txt a b c d e f
		
		cd barrettCUDA
		echo "making barrett"
		make barrett
		echo "after barrett"
		echo "q: 2^$q r: $r" >> $FILE_NAME
		for ((v =0; v<$T_MAX; v+=1));
		do				
			./barrett ../val.22 ../out.row ../out.col >> $FILE_NAME	
		done
		cd ..
	done
done
cd ../aggregate_queries/benchmarking/
echo "Finished Script"