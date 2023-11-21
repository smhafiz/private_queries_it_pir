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

data_file = 'yelp_academic_dataset_business'
data_path = 'yelp_data/' + data_file + '.json'
df = pd.read_json(data_path, lines=True)
frames = [df, df, df, df]
df = pd.concat(frames, axis =0, ignore_index=True)
df['categories'] = df['categories'].replace(np.nan, 'None')

split_categories = [i.replace(" ", "").split(',') for i in df['categories'].tolist()]

all_categories = []
unique_categories = {}
most_stars= {}

for i in range(len(split_categories)):
    for category in split_categories[i]:
        if category not in unique_categories:
            unique_categories[category] = 1
            all_categories.append(category)
        if category not in most_stars:
            most_stars[category] = df["stars"].loc[i]
        else:
            most_stars[category] = max(most_stars[category], df["stars"].loc[i])

filename = 'max_yelp_categories_exhaustive_query_2'
filter_field = 'categories'

database = df
all_categories = []
unique_categories = {}
for sub_list in split_categories:
    for category in sub_list:
        if category not in unique_categories:
            unique_categories[category] = 1
            all_categories.append(category)

filter_field_list = all_categories

times= []
for i in range(1):

    users_map = {}
    index = 0
    start = timeit.default_timer()
    for i in split_categories:
        for j in i:
            if j not in users_map:
                users_map[j] = index
                index+=1
    # print(users_map)

    hashmap_unique_users_tweet= {}
    hashmap_col = {}
    index = 0
    first_occurence = set()
    for i in range(len(split_categories)):
        for category in split_categories[i]:
            if most_stars[category] == df["stars"].loc[i] and category not in first_occurence:
                # if (users_map[category] not in hashmap_unique_users_tweet):
                #     hashmap_unique_users_tweet[users_map[category]] = set([index])
                # elif((users_map[category] in hashmap_unique_users_tweet)):
                #     hashmap_unique_users_tweet[users_map[category]].add(index)
                #### Not used for generating matrix  
            
                if index not in hashmap_col:
                    hashmap_col[index] = set([users_map[category]])
                else:
                    hashmap_col[index].add(users_map[category])
                first_occurence.add(category)
        index+=1

    for i in range(database.shape[0]):
        if i not in hashmap_col:
            hashmap_col[i] = set()

    hashmap_unique_users_tweet = dict(sorted(hashmap_unique_users_tweet.items()))
    hashmap_col = dict(sorted(hashmap_col.items()))
    unique_users = all_categories

    ### calculate matrix
    # M = np.zeros((len(unique_users), database.shape[0]))
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
        val = sorted(val)
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