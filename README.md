# Patchseq
This folder is to save and update the codes that used to analyze Patchseq data from Allen Institue.  
Developed by John Hongyu Meng and Yijie Kang.  
01/05/2022  


## Data Sources
**Data instruction, metadata and file manifest:**  
[ALLEN BRAIN MAP: Multimodal Characterization in Mouse Visual Cortex](https://portal.brain-map.org/explore/classes/multimodal-characterization/multimodal-characterization-mouse-visual-cortex) [(Gouwens, et al, 2020)](https://pubmed.ncbi.nlm.nih.gov/33186530/).  
**Transcriptomic data download index:**  
http://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/ (20200513_Mouse_PatchSeq_Release_count.v2.csv.tar is used here)


## Workflow
1. select cells according to specific layer and cell types
2. basic quality control based on selected cells' normalized (CPM) count matrix - **cell type annotation contamination issue unsolved, 01/06/2022**
3. for cells that pass the quality control, find their electrophysiological data's url (of nwb files) from file manifest
4. download nwb files
5. Requires ipfx package from Allen. Needs to be installed with Python 3.6 or 3.8. Seems 3.6 is a bit more robust. (No error on John's Windows Laptop, vs. 1 error with Python 3.8)
