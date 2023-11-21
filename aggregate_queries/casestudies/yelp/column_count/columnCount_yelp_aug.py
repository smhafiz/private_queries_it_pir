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

import random
from collections import Counter
files=['../ccs_files/out_yelp_aug.col']
for k in files:
    with open(k) as f:
        first_line = f.readline()
    first_line = list(map(int,first_line.split()))
    first_element = first_line[0]
    first_line = first_line[1:]
    #print(first_line)
    zcount =0
    nnzcount = 0
    zero_column_indexes = []
    for i in range(0,len(first_line)-1):
        sub=first_line[i+1]-first_line[i]
        #print("sub: ", sub)
        #if i == 50000:break
        if sub > 0:
            nnzcount+=1
        else:
            zcount+=1
            zero_column_indexes.append(zcount)
            #print(i)
    print("-----")
    print(k,": ")
    print("non zero: ", nnzcount)
    print("zero count: ", zcount)
    #print("zero_column_index: ", zero_column_indexes)
    file=k+".txt"
    my_str = ', '.join(str(item) for item in zero_column_indexes)
    with open(file, 'w') as g:
        g.write(my_str)
    print("number of columns: ", len(first_line[1:]))
    total = nnzcount+zcount

    if total == first_element:
        print("True")
        print("NonZero percent: ", ((nnzcount/total)*100))
    else:
        print("False")