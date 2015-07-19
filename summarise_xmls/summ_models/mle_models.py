
# coding: utf-8

# ## Bacterial rates analyses summary
# 
# - Pull marginal likelihoods from out_data.txt, which is obtained by running get_mle.sh in shell

# In[1]:

import pandas as pd
import numpy as np
import re, os, sys, subprocess
import matplotlib.pyplot as plt

get_ipython().magic(u'pylab inline')


# ### Pull mles: in the directory, run ./get_mles.sh The output is in out_data.txt

# In[2]:

mles_raw_file = open('out_data.txt', 'r').readlines()
mles_raw_file[:10]


# In[3]:

mle_data = pd.DataFrame(np.empty(shape = (1, 3)))

for i in range(len(mles_raw_file)):
    print mles_raw_file[i]
    if re.match('^log.+', mles_raw_file[i]) is None:
        mles_raw_file[i]
        temp_name = re.sub('[.]|mle[.]|log|\n', '', mles_raw_file[i])
#        print re.findall('=.+', mles_raw_file[i+1])
        temp_ps = float(re.sub('=| ', '', re.findall('=.+', mles_raw_file[i+1])[0]))
        temp_ss = float(re.sub('=| ', '', re.findall('=.+', mles_raw_file[i+2])[0]))
        temp_frame = pd.DataFrame([temp_name, temp_ps, temp_ss]).transpose()
        mle_data = pd.concat([mle_data, temp_frame], ignore_index=True, axis=0)

mle_data = mle_data.ix[1:, ]
mle_data.index = [i for i in range(mle_data.shape[0])]
mle_data.columns = ['file_name', 'ps', 'ss']


# In[4]:

mle_data.head()


# In[5]:

rates_data = np.empty(shape = (1, 3))
for i in mle_data['file_name']:
    temp_file = pd.read_csv(i+'.log', comment='#', sep = '\t')
    temp_rate_name = [c_ for c_ in temp_file.keys() if not re.match('ucld[.]mean|clock[.]rate', c_) is None]
    temp_data = np.concatenate([np.mean(temp_file[temp_rate_name]), np.percentile(temp_file[temp_rate_name], [2.5, 97.5])])
    rates_data = np.vstack([rates_data, temp_data])
rates_data = pd.DataFrame(rates_data[1:])


# In[6]:

print rates_data.shape
print mle_data.shape


# In[7]:

mle_data.head()


# In[8]:

summary_data = pd.concat([mle_data, rates_data], ignore_index=True, axis = 1)


# In[9]:

summary_data.columns = ['file_name', 'ps', 'ss', 'mean_rate', 'lowerHPD', 'higherHPD']
summary_data.head()


# In[10]:

#Next: select each bacteria type. For each type select best, and set up randomisations to validate the estimates.
summary_data.shape


# In[11]:

summary_data = pd.concat([summary_data, pd.DataFrame(np.empty(shape = (summary_data.shape[0], 2)))], axis = 1, ignore_index=True)
summary_data.columns = ['file_name', 'ps', 'ss', 'mean_rate', 'lowerHPD', 'higherHPD', 'genus', 'best_ML']


# In[12]:

summary_data.head()


# In[13]:

for i in range(summary_data.shape[0]):
    summary_data.ix[i, 'genus'] = '_'.join(re.split('_', summary_data.ix[i, 'file_name'])[:2])
summary_data


# ### Find the model with the highest ML 
# 
# - These should be run for the date randomisation test

# In[14]:

genus_set = set(summary_data.ix[:, 'genus'])
print genus_set
for g in genus_set:
    ps_temp = np.max(summary_data.ix[summary_data.ix[:, 'genus'] == g, 'ps'])
    print ps_temp
    
    #Find location of best mle for genus
    genus_best = []
    for r in range(summary_data.shape[0]):    
        if (summary_data.ix[r, 'ps'] ==  ps_temp) and (summary_data.ix[r, 'genus'] == g):
            genus_best.append(r)

#    max_loc = np.array([(summary_data.ix[:, 'ps'] ==  ps_temp), (summary_data.ix[:, 'genus'] == g)])

#    max_loc = max_loc.transpose()
#    max_loc = np.where([all(max_loc[i, :]) for i in range(max_loc.shape[0])])[0][0]
    
    summary_data.ix[genus_best, 'best_ML'] = 'BEST'


# In[14]:




# In[15]:

summary_data.ix[summary_data.ix[:, 'best_ML'] == 'BEST', :]


# In[16]:

summary_data.ix[summary_data.ix[:, 'best_ML'] == 'BEST', :].to_csv('res_10_datasets.csv', index = False)


# In[18]:

summary_data.ix[summary_data.ix[:, 'best_ML'] == 'BEST', :].to_csv('best_models.csv')


# In[ ]:



