#!/usr/bin/env python
# coding: utf-8

# In[15]:


import numpy as np
import pandas as pd
import dandi


# In[ ]:





# In[16]:


metadata = pd.read_csv("20200711_patchseq_metadata_mouse.csv")
metadata.rename(columns={'Unnamed: 21':'T-type_Label'}, inplace=True)


# In[17]:


metadata.insert(loc=0, column="index", value = range(len(metadata))) # add an column of original index


# In[18]:


# metadata


# In[ ]:





# In[ ]:





# In[ ]:





# In[19]:


sst_L23_metadata = metadata.loc[(metadata["T-type_Label"].str.startswith("Sst")) & (metadata["structure"] == "VISp2/3"), :]
pv_L23_metadata = metadata.loc[(metadata["T-type_Label"].str.startswith("Pvalb")) & (metadata["structure"] == "VISp2/3"), :]

# no cells are filtered out by Seurat quality control,tid above can be used directly.
# some genes are filtered


# In[ ]:





# In[20]:


# if some cells are filtered out and lists of cells with good quality are generated, then:
# filtered_pv_L23_tid = pd.read_csv("filtered_pv_L23_cells.csv",
#                                    index_col = 0,)

# filtered_sst_L23_tid = pd.read_csv("filtered_sst_L23_cells.csv",
#                                     index_col = 0,)


# In[ ]:





# In[23]:


# L23_pv_metadata = pv_L23_metadata.loc[pv_L23_metadata["transcriptomics_sample_id"].isin(filtered_pv_L23_tid["x"]), :]
# L23_sst_metadata = sst_L23_metadata.loc[sst_L23_metadata["transcriptomics_sample_id"].isin(filtered_sst_L23_tid["x"]), :]
L23_pv_metadata = pv_L23_metadata
L23_sst_metadata = sst_L23_metadata


# In[24]:


L23_pv_specimen_id = L23_pv_metadata.loc[:, "cell_specimen_id"]
L23_sst_specimen_id = L23_sst_metadata.loc[:, "cell_specimen_id"]


# In[ ]:





# In[25]:


# this csv file is rewrited from xlsx file by R. Code is in 0929_demo.Rmd
file_manifest = pd.read_csv("2021-09-13_mouse_file_manifest.csv",index_col = 0,)
# file_manifest.head()


# In[26]:


# get download links from file_manifest mapped by filtered metadata's cell_specimen_id
nwb_pv_L23_urls = file_manifest.loc[
    (file_manifest["cell_specimen_id"].isin(L23_pv_specimen_id))  &
    (file_manifest["file_type"] == "nwb"),
    :]


# In[27]:


# get download links from file_manifest mapped by filtered metadata's cell_specimen_id
nwb_sst_L23_urls = file_manifest.loc[
    (file_manifest["cell_specimen_id"].isin(L23_sst_specimen_id))  &
    (file_manifest["file_type"] == "nwb"),
    :]


# In[ ]:





# In[28]:


# nwb_pv_L23_urls["archive_uri"].values[0]


# In[29]:


nwb_pv_L23_downloads = pd.DataFrame(nwb_pv_L23_urls["archive_uri"] + nwb_pv_L23_urls["file_name"])
nwb_pv_L23_downloads.to_csv('nwb_pv_L23_downloads.csv', index=False, header = 0, sep = '\n')

nwb_sst_L23_downloads = pd.DataFrame(nwb_sst_L23_urls["archive_uri"] + nwb_sst_L23_urls["file_name"])
nwb_sst_L23_downloads.to_csv('nwb_sst_L23_downloads.csv', index=False, header = 0, sep = '\n')

# use these files to download corresponding nwb files for each cell type


# In[ ]:





# In[ ]:





# In[32]:


# prepare nwb file list with the %0D on the HPC downloaded files
nwb_pv_L23_filelist = pd.DataFrame(nwb_pv_L23_urls["file_name"] + '%0D')
nwb_pv_L23_filelist.to_csv('nwb_pv_L23_filelist.csv', index=True)

nwb_sst_L23_filelist = pd.DataFrame(nwb_sst_L23_urls["file_name"] + '%0D')
nwb_sst_L23_filelist.to_csv('nwb_sst_L23_filelist.csv', index=True)


# In[14]:





# In[ ]:





# In[35]:





# In[40]:


# merge the nwb file related information to the metadata 
nwb_related = file_manifest[['cell_specimen_id', 'file_name', 'file_type', 'archive', 'archive_uri']]
nwb_related = nwb_related.loc[(file_manifest["file_type"] == "nwb"),:]

L23_pv_metadata_nwbs = pd.merge(L23_pv_metadata, nwb_related, on = "cell_specimen_id", sort = False, how = 'left')
L23_sst_metadata_nwbs = pd.merge(L23_sst_metadata, nwb_related, on = "cell_specimen_id", sort = False, how = 'left')


# In[43]:


# save merged metadata as csv files
L23_pv_metadata_nwbs.to_csv('L23_pv_metadata_nwbs.csv', index=True)
L23_sst_metadata_nwbs.to_csv('L23_sst_metadata_nwbs.csv', index=True)


# In[42]:


# L23_pv_metadata_nwbs


# In[ ]:


# L23_sst_metadata_nwbs


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




