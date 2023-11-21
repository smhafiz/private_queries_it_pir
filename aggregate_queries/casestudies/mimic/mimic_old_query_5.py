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


data_file = 'ADMISSIONS'
# data_path = 'mimic_data/' + data_file + '.csv'
data_file2 = 'PRESCRIPTIONS'
# data_path1 = '/content/drive/MyDrive/Research/PIR/' + data_file + '.csv'
# data_path2 = '/content/drive/MyDrive/Research/PIR/' + data_file2 + '.csv'

# df = pd.read_csv(data_path1)
# df2 = pd.read_csv(data_path2)


from scipy.sparse import csr_matrix

import datetime

data_file = 'ADMISSIONS'
data_path = 'mimic_data/'+ data_file + '.csv'
df = pd.read_csv(data_path)


data_file = 'PRESCRIPTIONS'
data_path = 'mimic_data/' + data_file + '.csv'
df2 = pd.read_csv(data_path)

c = df.columns.tolist()+['DOSE_VAL_RX']

new_df = pd.merge(df2, df,  how='left', on=['SUBJECT_ID','HADM_ID'])[[i for i in c if i != 'ROW_ID']]


import timeit
start = timeit.default_timer()
filename = 'mimic_data/mimic_old_query_5'

filter_field = 'SUBJECT_ID'
aggregation_field = 'DOSE_VAL_RX'
database = new_df
filter_field_list = database[filter_field].tolist()


times= []
for i in range(100):

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
        # if (users_map[filter_field_list[i]] not in hashmap_unique_users_tweet):
        #     hashmap_unique_users_tweet[users_map[filter_field_list[i]]] = [index] 
        # elif((users_map[filter_field_list[i]] in hashmap_unique_users_tweet)):
        #     hashmap_unique_users_tweet[users_map[filter_field_list[i]]].append(index)
        #### Not used for generating matrix  
    
        if index not in hashmap_col:
            hashmap_col[index] = [users_map[filter_field_list[i]]]
        else:
            hashmap_col[index].append(users_map[filter_field_list[i]])
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


    a_row_file = [len(unique_users)] + my_row
    a_col_file = [len(database[filter_field])] + my_col

    a_row_file = [str(i) for i in a_row_file]
    a_col_file = [str(i) for i in a_col_file]

    stop = timeit.default_timer()
    diff  = (stop - start)
    # print("times 1:", diff1)
    # print("times 2:", diff2)
    times.append(diff)


print(sum(times)/len(times))
print(1.96*np.std(times)/np.sqrt(len(times)))
print("P: " , len(unique_users))
print("R: ", len(database[filter_field]))
print("Row: ", len(my_row))
print(len(my_row)/(len(unique_users)* len(database[filter_field])))
# with open(filename+".row", 'w') as totxt_file:
#     totxt_file.write(" ".join(a_row_file))

# with open(filename+".col", 'w') as totxt_file:
#     totxt_file.write(" ".join(a_col_file))


# stop = timeit.default_timer()