#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import random
from collections import Counter
import math


# In[2]:


dataset = pd.read_csv('cold database.csv')
Yama_dataset = dataset.copy()


# In[3]:


Yama_dataset = Yama_dataset.drop(['Unnamed: 0', 'Unnamed: 0.1'], axis = 1)
Yama_1hr = Yama_dataset.iloc[:, [0, 1, 2]]

#1. 把感性的drop掉
drop_list = []
for i in range(len(Yama_1hr['gene'])):
    if Yama_1hr.iloc[i][1] < 0:
        drop_list.append(i)
Yama_1hr_drop = Yama_1hr.drop(drop_list)


# In[6]:


#轉換有drop的
#2. p轉換成score經過-log(p)
p_score_1hr_drop = []
p_1hr = list(Yama_1hr_drop['logFC 1hr/0hr p-value'])
for i in range(len(p_1hr)):
    p_score_1hr_drop.append(-math.log(p_1hr[i], 10))
Yama_1hr_drop['score from p'] = p_score_1hr_drop


# In[8]:


Yama_1hr_drop = Yama_1hr_drop.reset_index().drop(['index'], axis = 1)


# In[10]:


Yama_1hr_drop


# In[12]:


def prioritized_external(data, boots_round, permut_round, ranks, file_num):  #有file_num
    bootstrap_genes_record = []  #紀錄bootstrap取出的基因頻率用

    #創立一個變數存入原始資料庫
    ct_genes_score = data.iloc[:, [0, -1]]

    #建立所有基因的字典以及bootstrap超越permuted的次數，一開始的次數接設定為0
    all_genes = list(ct_genes_score['gene'])
    original_time = [0 for i in range(len(ct_genes_score['gene']))]
    boots_over_permute_times = dict(zip(all_genes, original_time))

    for bootstrap_times in range(boots_round):  
        #bootstrap
        bootstrap_genes = ct_genes_score.sample(n = len(ct_genes_score.index), replace = True)
        #排序bootstrap結果分數
        bootstrap_score = bootstrap_genes.sort_values(['score from p'], ascending = False)[0:ranks]

        #紀錄這一次bootstrap拿到的基因
        for bootstrap_gene in bootstrap_score['gene']:
            bootstrap_genes_record.append(bootstrap_gene)

        for permutation_times in range(permut_round):  #每一次的bootstrap中的permutation過程
            #先複製原先的Dataframe
            new_gene_score = ct_genes_score.copy()
            #把gene這一欄取出來並轉換成list，重新排列
            permuted_genes = random.sample(list(ct_genes_score['gene']), len(ct_genes_score['gene']))
            #再存回去取代原本的gene欄位
            new_gene_score['gene'] = permuted_genes
            #改變欄位成permuted後權重分數
            new_gene_score.columns = ['gene', 'permuted score']

            #找出bootstrap樣本出現的基因在permuted之後的對應分數
            score = []
            permu_score_dict = dict(zip(list(new_gene_score['gene']), list(new_gene_score['permuted score'])))
            for i in range(len(bootstrap_score['gene'])):
                score.append(permu_score_dict[list(bootstrap_score['gene'])[i]])

            #合併儲存格，可看到permuted 和bootstrap的分數
            final_genes_2score = bootstrap_score
            final_genes_2score['permuted'] = score

            #開始篩選，假設加權基因大於permuted gene 就把所在字典的那個基因對應到的超過次數加一
            for i in range(len(final_genes_2score['gene'])):
                if final_genes_2score.iloc[i, :]['score from p'] > final_genes_2score.iloc[i, :]['permuted']:
                    boots_over_permute_times[final_genes_2score.iloc[i, :]['gene']] = boots_over_permute_times[final_genes_2score.iloc[i, :]['gene']] + 1 
        print('finished round {} bootstrap'.format(bootstrap_times + 1))
    
    #存檔
    all_result = pd.DataFrame()
    all_result['gene'] = list(boots_over_permute_times.keys())
    all_result['boots over times'] = list(boots_over_permute_times.values())
                       
    gene_frequency = sorted(Counter(bootstrap_genes_record).items(), key = lambda x : x[1], reverse = True)
        
    gene2 = []
    frequency = []
    for i in range(len(gene_frequency)):
        gene2.append(gene_frequency[i][0])
        frequency.append(gene_frequency[i][1])
        
    file2 = pd.DataFrame()
    file2['gene'] = gene2
    file2['frequency'] = frequency
    file2.to_csv('1hr_44 freq drop' + file_num + '.csv')
    all_result.to_csv('1hr_44 all drop' + file_num + '.csv')


# In[ ]:


prioritized_external(Yama_1hr_drop, 50, 10000, 44, '1')


# bootstrap_genes_record = []  #紀錄bootstrap取出的基因頻率用
# 
# #創立一個變數存入原始資料庫
# ct_genes_score = Yama_1hr_drop.iloc[:, [0, -1]]
# 
# 
# #建立所有基因的字典以及bootstrap超越permuted的次數，一開始的次數接設定為0
# all_genes = list(ct_genes_score['gene'])
# original_time = [0 for i in range(len(ct_genes_score['gene']))]
# boots_over_permute_times = dict(zip(all_genes, original_time))
# 
# for bootstrap_times in range(5):  
#     #bootstrap
#     bootstrap_genes = ct_genes_score.sample(n = len(ct_genes_score.index), replace = True)
#     #排序bootstrap結果分數
#     bootstrap_score = bootstrap_genes.sort_values(['score from p'], ascending = False)[0:44]
#     
#     #紀錄這一次bootstrap拿到的基因
#     for bootstrap_gene in bootstrap_score['gene']:
#         bootstrap_genes_record.append(bootstrap_gene)
#     
#     for permutation_times in range(5):  #每一次的bootstrap中的permutation過程
#         #先複製原先的Dataframe
#         new_gene_score = ct_genes_score.copy()
#         #把gene這一欄取出來並轉換成list，重新排列
#         permuted_genes = random.sample(list(ct_genes_score['gene']), len(ct_genes_score['gene']))
#         #再存回去取代原本的gene欄位
#         new_gene_score['gene'] = permuted_genes
#         #改變欄位成permuted後權重分數
#         new_gene_score.columns = ['gene', 'permuted score']
#         
#         #找出bootstrap樣本出現的基因在permuted之後的對應分數
#         score = []
#         permu_score_dict = dict(zip(list(new_gene_score['gene']), list(new_gene_score['permuted score'])))
#         for i in range(len(bootstrap_score['gene'])):
#             score.append(permu_score_dict[list(bootstrap_score['gene'])[i]])
#         
#         #合併儲存格，可看到permuted 和bootstrap的分數
#         final_genes_2score = bootstrap_score
#         final_genes_2score['permuted'] = score
#         
#         #開始篩選，假設加權基因大於permuted gene 就把所在字典的那個基因對應到的超過次數加一
#         for i in range(len(final_genes_2score['gene'])):
#             if final_genes_2score.iloc[i, :]['score from p'] > final_genes_2score.iloc[i, :]['permuted']:
#                 boots_over_permute_times[final_genes_2score.iloc[i, :]['gene']] = boots_over_permute_times[final_genes_2score.iloc[i, :]['gene']] + 1 
#     print('finished round {} bootstrap'.format(bootstrap_times + 1))

# all_result_drop = pd.DataFrame()
# all_result_drop['gene'] = list(boots_over_permute_times.keys())
# all_result_drop['boots over times'] = list(boots_over_permute_times.values())

# boots_over = sorted(boots_over_permute_times.items(), key = lambda x : x[1], reverse = True)[0:44]
# gene = []
# times_over_permute = []
# for i in range(len(boots_over)):
#     gene.append(boots_over[i][0])
#     times_over_permute.append(boots_over[i][1])

# file_drop = pd.DataFrame()
# file_drop['gene'] = gene
# file_drop['times over permute'] = times_over_permute

# gene_frequency = sorted(Counter(bootstrap_genes_record).items(), key = lambda x : x[1], reverse = True)
# gene2 = []
# frequency = []
# for i in range(len(gene_frequency)):
#     gene2.append(gene_frequency[i][0])
#     frequency.append(gene_frequency[i][1])

# file2_drop = pd.DataFrame()
# file2_drop['gene'] = gene2
# file2_drop['frequency'] = frequency

# file_drop.to_csv('short-term cold external (drop) boots.csv')
# file2_drop.to_csv('short-term cold external (drop) freq.csv')
# all_result_drop.to_csv('short-term cold external (drop) all.csv')

# In[ ]:


# with pd.ExcelWriter('short-term cold external (drop).xlsx') as writer:
#     file_drop.to_excel(writer, sheet_name = 'boots樣本超越次數')
#     file2_drop.to_excel(writer, sheet_name = 'appeared frequency of genes')
#     all_result_drop.to_excel(writer, sheet_name = 'whole results')

