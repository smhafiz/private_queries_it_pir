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

data_files = ['twitter_1050_11','twitter_1051_11','twitter_1053_11','twitter_1054_11',
                'twitter_1055_11','twitter_1056_11','twitter_1057_11','twitter_1058_11',
                'twitter_1059_11', 'twitter_1060_11', 'twitter_1061_11','twitter_1062_11']

df1 = pd.read_csv("twitter_data/"+data_files[0]+'.csv',  lineterminator='\n')

df2 = pd.read_csv("twitter_data/"+ data_files[1]+'.csv',  lineterminator='\n')

df3 = pd.read_csv("twitter_data/"+ data_files[2]+'.csv',  lineterminator='\n')
df4 = pd.read_csv("twitter_data/"+ data_files[3]+'.csv',  lineterminator='\n')
df5 = pd.read_csv("twitter_data/"+ data_files[4]+'.csv',  lineterminator='\n')
df6 = pd.read_csv("twitter_data/"+ data_files[5]+'.csv',  lineterminator='\n')
df7 = pd.read_csv("twitter_data/"+ data_files[6]+'.csv',  lineterminator='\n')
df8 = pd.read_csv("twitter_data/"+ data_files[7]+'.csv',  lineterminator='\n')
df9 = pd.read_csv("twitter_data/"+ data_files[8]+'.csv',  lineterminator='\n')
df10= pd.read_csv("twitter_data/"+ data_files[9]+'.csv',  lineterminator='\n')
df11 = pd.read_csv("twitter_data/"+ data_files[10]+'.csv',  lineterminator='\n')
df12 = pd.read_csv("twitter_data/"+ data_files[11]+'.csv',  lineterminator='\n')

main_df = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12], axis =0, ignore_index=True)
main_df = main_df.dropna(subset=['user'])



# database = df  #trimming dataset
database = main_df
filtered_df = database[database['likeCount']>100]
#tracking time for index generation
times = []
for i in range(100):
    import timeit
    
    filter_field = 'user'   #each row of index of queries matrix will be unique value of this field
    # aggregation_field = 'likeCount'  
    aggregation_field = 'retweetCount'
    
    # agg_list = database[aggregation_field].tolist()
    filter_field_list = filtered_df[filter_field].tolist() ## list of users
    # possible_filter_vals = list(set(filter_field_list)) ## list of unique users
    

    #### Brijesh's 
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
    for i in range(filtered_df.shape[0]):
        # if (filtered_df[aggregation_field].iloc[i]==0.0):
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
        if (filtered_df[aggregation_field].iloc[i]==0.0):
            if index not in hashmap_col:
                hashmap_col[index] = [users_map[filter_field_list[i]]]
            else:
                hashmap_col[index].append(users_map[filter_field_list[i]])
        else:
            if index not in hashmap_col:
                hashmap_col[index] = []
            # index-=1
        index+=1

    unique_users = set(filtered_df['user'])

    ### calculate matrix
    # M = np.zeros((len(unique_users), len(database['user'])))
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
    length = len(my_col)
    my_col = my_col + [my_col[-1] for i in range(len(database[filter_field])-len(my_col)+1)]


    
    # m = index_of_query_array
    # m = M.transpose()
    

    #generating compressed storage representation 
    # row = csr_matrix(m).indices.tolist()
    # col = csr_matrix(m).indptr.tolist()
    # vals = csr_matrix(m).data.tolist()
    # print(row==my_row)
    # print(col==my_col)
    # a_row_file = [len(unique_users)] + my_row
    # a_col_file = [len(database['user'])] + my_col
    

    # a_row_file = [str(i) for i in a_row_file]
    # a_col_file = [str(i) for i in a_col_file]
    
    stop = timeit.default_timer()
    diff  = (stop - start)
    # print("times 1:", diff1)
    # print("times 2:", diff2)
    times.append(diff)

print(sum(times)/len(times))
print(1.96*np.std(times)/np.sqrt(len(times)))

filename = "twitter_filtered_retweet_"+str(len(unique_users))+"_"+str((database.shape[0]))
# Saving File
# with open(filename+".row", 'w') as totxt_file:
#    totxt_file.write(" ".join(a_row_file))#
# with open(filename+".col", 'w') as totxt_file:
#    totxt_file.write(" ".join(a_col_file))

print("P: " , len(unique_users))
print("R: ", len(database[filter_field]))
print("Row: ", len(my_row))
print(len(my_row)/(len(unique_users)* len(database[filter_field])))