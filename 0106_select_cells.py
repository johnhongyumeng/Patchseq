#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


count_data = pd.read_csv("20200513_Mouse_PatchSeq_Release_count.v2.csv",
                        index_col = 0,)
count_data # genes in row, samples/cells in column, 4435 cells, 45768 genes


# In[3]:


metadata = pd.read_csv("20200711_patchseq_metadata_mouse.csv")
metadata.head()


# In[4]:


metadata.rename(columns={'Unnamed: 21':'T-type_Label'}, inplace=True)
metadata.head()


# In[5]:


# select cells whose T-type_Label is "Sst" and is of 'VISp2/3'
sst_L23_metadata = metadata.loc[(metadata["T-type_Label"].str.startswith("Sst")) & (metadata["structure"] == "VISp2/3"), :]

sst_L23_metadata # 172 cells


# In[6]:


# select cells whose T-type_Label is "Sst" and is of 'VISp2/3'
pv_L23_metadata = metadata.loc[(metadata["T-type_Label"].str.startswith("Pvalb")) & (metadata["structure"] == "VISp2/3"), :]

pv_L23_metadata # 145 cells


# In[ ]:





# In[7]:


sst_L23_counts = count_data.loc[:, count_data.columns.isin(sst_L23_metadata["transcriptomics_sample_id"].tolist())]
sst_L23_counts


# In[8]:


pv_L23_counts = count_data.loc[:, count_data.columns.isin(pv_L23_metadata["transcriptomics_sample_id"].tolist())]
pv_L23_counts


# In[ ]:





# In[9]:


sst_L23_counts.to_csv('sst_L23_counts.csv', index=True)
pv_L23_counts.to_csv('pv_L23_counts.csv', index=True)
# then use R to do basic quality control of their scRNA-seq data


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




