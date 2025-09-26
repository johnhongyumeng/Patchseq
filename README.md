# Patchseq
This folder is to save and update the codes that used to analyze Patchseq data from Allen Institue.  The code is used in generating the figures for paper:

In Search of Transcriptomic Correlates of Neuronal Firing-Rate Adaptation across Subtypes, Regions and Species: A Patch-seq Analysis

https://www.biorxiv.org/content/10.1101/2024.12.05.627057v2

For specific figure reproduction, see the correponding R code. 

For HH model, see code under ./HHmodel/

Developed by John Hongyu Meng and Yijie Kang.  
09/26/2025  


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
5. Requires ipfx package from Allen. Needs to be installed with Python 3.6 or 3.8. It seems that build a new environment with Python 3.6 is most robust way using the package.

## To analyze the electronical features, using AllenAI_Loop_* to generate a combined analyze. 
1. The code AllenAI_cell_*, AllenAI_ind_* are used to check individual cell or individual sweep, respectively.
2. To analyze batch data, run code ./EphysCode/Allen_V1_general_nowaring.py to get most of E-feature, including half-width (HW). Please note that the AI generated from the code is not used anymore. We keep it here for historical reason.
3. To generate measurement about firing rate adaptation and onset bursting (i.e., dAdap and p-value for onser bursting.), run code ./EphysCode/Allen_V1_getrheo.py. This code using the previous generated AP info directly, so can be used independent of ipfx package from Allen
