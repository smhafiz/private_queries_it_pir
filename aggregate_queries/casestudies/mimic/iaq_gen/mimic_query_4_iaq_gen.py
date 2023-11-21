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

import timeit
import numpy as np
import pandas as pd
start = timeit.default_timer()
filename = 'min_of_urgent_admissions'
filter_field = 'ADMISSION_TYPE'
aggregation_field = 'subject_id_by_at_URGENT'

data_file = 'ADMISSIONS'
data_path = 'mimic_data/' + data_file + '.csv'
df = pd.read_csv(data_path)
database = df

filter_field_list = database[filter_field].tolist()

times= []
for i in range(1):
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
    least_recent_time = {}
    unique_ad_type = database[filter_field].unique()
    for ad_type in unique_ad_type:
        least_recent_time[ad_type] = min(database[database[filter_field] == ad_type]["ADMITTIME"])
    for i in range(len(filter_field_list)):
        # if database["ADMITTIME"].loc[i] == least_recent_time[filter_field_list[i]]:
        #     if (users_map[filter_field_list[i]] not in hashmap_unique_users_tweet):
        #         hashmap_unique_users_tweet[users_map[filter_field_list[i]]] = [index] 
        #     elif((users_map[filter_field_list[i]] in hashmap_unique_users_tweet)):
        #         hashmap_unique_users_tweet[users_map[filter_field_list[i]]].append(index)
        # else:
        #     if users_map[filter_field_list[i]] not in hashmap_unique_users_tweet:
        #         hashmap_unique_users_tweet[users_map[filter_field_list[i]]] = []
        #### Not used for generating matrix  
        if database["ADMITTIME"].loc[i] == least_recent_time[filter_field_list[i]]:
            if index not in hashmap_col:
                hashmap_col[index] = [users_map[filter_field_list[i]]]
            else:
                hashmap_col[index].append(users_map[filter_field_list[i]])
        else:
            if index not in hashmap_col:
                hashmap_col[index] = []
        index+=1

    unique_users = set(database[filter_field])

    ### calculate matrix
    M = np.zeros((len(unique_users), len(database[filter_field])))
    index = 0
    for key, val in hashmap_unique_users_tweet.items():
        for i in val:
            M[key][i] = 1
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


    # import numpy as np
    # from scipy.sparse import csr_matrix


    # m = index_of_query_array
    # m = M.transpose()

    #print(sum([sum(i) for i in index_of_query]))
    #row = get_non_zero_indices(m)
    # row = csr_matrix(m).indices.tolist()
    # col = csr_matrix(m).indptr.tolist()
    # vals = csr_matrix(m).data.tolist()

    # print(row == my_row)
    # print(col == my_col)
    #print(vals)
    stop = timeit.default_timer()
    diff = (stop - start)
    # print("times 1:", diff1)
    # print("times 2:", diff2)
    times.append(diff)

    a_row_file = [len(unique_users)] + my_row
    a_col_file = [len(database[filter_field])] + my_col

    a_row_file = [str(i) for i in a_row_file]
    a_col_file = [str(i) for i in a_col_file]

with open(filename+".row", 'w') as totxt_file:
    totxt_file.write(" ".join(a_row_file))

with open(filename+".col", 'w') as totxt_file:
    totxt_file.write(" ".join(a_col_file))

print(sum(times)/len(times))
print(1.96*np.std(times)/np.sqrt(len(times)))
print("P: " , len(unique_users))
print("R: ", len(database[filter_field]))
print("Row: ", len(my_row))
print(len(my_row)/(len(unique_users)* len(database[filter_field])))

# stop = timeit.default_timer()