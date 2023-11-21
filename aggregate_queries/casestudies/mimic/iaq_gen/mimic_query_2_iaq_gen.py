# This file is part of Private aggregate queries v0.1.
# 
# 
# Copyright (C) 2022, Chitrabhanu Gupta, Brijesh Vora, and Syed Mahbub Hafiz.
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

import numpy as np
import pandas as pd

import datetime
data_file = 'ADMISSIONS'
data_path = 'mimic_data/' + data_file + '.csv'
df = pd.read_csv(data_path)
database = df

filename = 'count_of_hispanic_admission'
filter_field = 'ADMISSION_TYPE'
aggregation_field = 'ETHNICITY'

#tracking time for index generation
times = []
for i in range(1):
    import timeit
    
  
    filter_field_list = database[filter_field].tolist() ## list of users
     
    users_map = {}
    index = 0
    start = timeit.default_timer()
    for i in filter_field_list:
        if i not in users_map:
            users_map[i] = index
            index+=1
    # print(users_map)
    hashmap_unique_users_tweet= {}
    hashmap_col = {}
    index = 0
    for i in range(len(filter_field_list)):
        # if (database[aggregation_field].loc[i][:8]=="HISPANIC"):
        # ### used for generating matrix
        #     if (users_map[filter_field_list[i]] not in hashmap_unique_users_tweet):
        #         hashmap_unique_users_tweet[users_map[filter_field_list[i]]] = [index] 
        #     elif((users_map[filter_field_list[i]] in hashmap_unique_users_tweet)):
        #         hashmap_unique_users_tweet[users_map[filter_field_list[i]]].append(index) 
        # else:
        #     if users_map[filter_field_list[i]] not in hashmap_unique_users_tweet:
        #         hashmap_unique_users_tweet[users_map[filter_field_list[i]]] = []
            # index-=1
        #### Not used for generating matrix  
        if (database[aggregation_field].iloc[i][:8]=="HISPANIC"):
            if index not in hashmap_col:
                hashmap_col[index] = [users_map[filter_field_list[i]]]
            else:
                hashmap_col[index].append(users_map[filter_field_list[i]])
        else:
            if index not in hashmap_col:
                hashmap_col[index] = []
            # index-=1
        index+=1

    unique_users = set(database[filter_field])

    ### calculate matrix
    # M = np.zeros((len(unique_users), len(database[filter_field])))
    # index = 0
    # for key, val in hashmap_unique_users_tweet.items():
    #     for i in val:
    #         M[key][i] = 1
    #     index+=1
    # stop = timeit.default_timer()
    # diff1 = (stop - start)

    # start = timeit.default_timer()
    ### Calculate col and row from hashmap directly
    # value_pairs = []
    index = 0
    my_row = []
    my_col = []
    for key, val in hashmap_col.items():
        for i in val:
            # value_pairs.append((i, key))
            my_row.append(i)
        my_col.append(index)
        index += len(val)
    my_col.append(index)


    # m = M.transpose()

    #generating compressed storage representation 
    # row = csr_matrix(m).indices.tolist()
    # col = csr_matrix(m).indptr.tolist()
   
    # print(row==my_row)
    # print(col==my_col)

    stop = timeit.default_timer()
    diff = (stop - start)

    a_row_file = [len(unique_users)] + my_row
    a_col_file = [len(database[filter_field])] + my_col
    
    a_row_file = [str(i) for i in a_row_file]
    a_col_file = [str(i) for i in a_col_file]
    
    # print("times 1:", diff1)
    # print("times 2:", diff2)
    times.append(diff)

print(sum(times)/len(times))
print(1.96*np.std(times)/np.sqrt(len(times)))

# Saving File
with open(filename+".row", 'w') as totxt_file:
   totxt_file.write(" ".join(a_row_file))#
with open(filename+".col", 'w') as totxt_file:
   totxt_file.write(" ".join(a_col_file))

print("P: " , len(unique_users))
print("R: ", len(database[filter_field]))
print("Row: ", len(my_row))
print(len(my_row)/(len(unique_users)* len(database[filter_field])))
