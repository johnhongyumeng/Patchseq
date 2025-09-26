# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 11:01:32 2025
Calculate the abundancy

@author: yawnt
import os
print(os.getcwd())
"""

import pandas as pd
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

dataset_str='AllenV1PV_ephys_mega0701'
# dataset_str='AllenV1SST_ephys_mega0701'
# 
# df = pd.read_excel('AllenV1SST_ehys_mega0623.xlsx', sheet_name='Relevant')

df = pd.read_excel(dataset_str+'.xlsx', sheet_name='Relevant')

df.head(10)      # Shows the first 10 rows
df.columns



from matplotlib import rcParams, font_manager

#Figure Parameters 
rcParams['figure.figsize'] = (11,7.33) #figure size in inches
rcParams['font.family'] = "sans-serif" 
rcParams['font.sans-serif'] = ['Calibri']
rcParams['font.weight'] = "roman" 
rcParams['font.style'] = "normal" 
rcParams['font.size'] = 32 
rcParams['pdf.fonttype'] = 42
rcParams['axes.linewidth'] = 2 #thickness of the border
rcParams['xtick.major.width'] = rcParams['axes.linewidth'] 
rcParams['ytick.major.width'] = rcParams['axes.linewidth']
rcParams['axes.spines.right'] = False #hides right border
rcParams['axes.spines.top'] = False #hides top
rcParams['legend.frameon'] = False #hides box around the legend
rcParams['legend.fontsize'] = 18 #font size in pt
rcParams['axes.labelsize'] = 32 
rcParams['xtick.labelsize'] = 24
rcParams['ytick.labelsize'] = rcParams['xtick.labelsize']
rcParams['lines.linewidth'] = 3
rcParams['xtick.major.size'] = 5
rcParams['lines.markersize'] = 16
rcParams['ytick.major.size'] = rcParams['xtick.major.size'] 
rcParams['lines.color'] = 'black'
rcParams['axes.prop_cycle'] = plt.cycler(color=['black', 'red', 'blue', 'green', 'purple', 'brown'])





# --- 1. keep only “hero” rows ----------------------------------------
df_hero = df[df['hero or Max'].str.lower() == 'hero'].copy()

# --- 2A. BH correction (global) --------------------------------------
group_temp = df_hero.copy()  # avoid SettingWithCopyWarning**
group_temp['pvalue_ob'] = group_temp['pvalue_ob'].fillna(1)

_, p_adj_global, _, _ = multipletests(group_temp['pvalue_ob'],
                                      method='fdr_bh')
df_hero['p_adj_global'] = p_adj_global

# --- 2B. BH correction (per T-type) ----------------------------------
def bh_per_group(group):
    group = group.copy()  # avoid SettingWithCopyWarning**
    group['pvalue_ob'] = group['pvalue_ob'].fillna(1)
    _, p_adj, _, _ = multipletests(group['pvalue_ob'], method='fdr_bh')
    group['p_adj_Ttype'] = p_adj
    return group

df_hero = (
    df_hero
    .groupby('T-type_Label', group_keys=False, sort=False)
    .apply(bh_per_group)
)

# --- 3. Save ----------------------------------------------------------
out_file = dataset_str+'_BH_corrected.xlsx'
df_hero.to_excel(out_file, index=False)


# --- Compute OB sensitivity ------------------------------------------
df_hero['OB_global_sensitive'] = (df_hero['OBrate_rec'] > 10) & (df_hero['p_adj_global'] < 0.05)

# --- Count sensitive vs. not -----------------------------------------
counts = df_hero['OB_global_sensitive'].value_counts()
sizes = [counts.get(True, 0), counts.get(False, 0)]
labels = [f'n={sizes[0]}', f'n={sizes[1]}']
total = sum(sizes)

# --- Plot pie chart --------------------------------------------------
fig, ax = plt.subplots(figsize=(4, 4))
ax.pie(sizes, labels=labels, autopct='%1.1f%%',
       startangle=90, counterclock=False, colors=['red', 'lightgray'])
ax.set_title(f'n={total}')
ax.axis('equal')
plt.show()
fig.savefig('./fig/'+dataset_str+'CellType_abundancy'+'.jpg',bbox_inches='tight')
fig.savefig('./fig/'+dataset_str+'CellType_abundancy'+'.svg',bbox_inches='tight')



# Breakdown to sub-type
df_hero['OB_Ttype_sensitive'] = (df_hero['OBrate_rec'] > 10) & (df_hero['p_adj_Ttype'] < 0.05)

# Group by T-type and compute counts
summary_df = (
    df_hero.groupby('T-type_Label')
    .agg(
        total_cells=('OB_Ttype_sensitive', 'count'),
        num_OB_sensitive=('OB_Ttype_sensitive', 'sum')
    )
    .reset_index()
)

# Add fraction
summary_df['fraction_OB_sensitive'] = summary_df['num_OB_sensitive'] / summary_df['total_cells']

# Optional: sort by most OB-sensitive T-types
summary_df = summary_df.sort_values('fraction_OB_sensitive', ascending=False)
out_file = dataset_str+'_TtypeAbundancy.xlsx'
summary_df.to_excel(out_file, index=False)



# # --- Select top 6 T-types by OB-sensitive cell count --------------------
# top6_df = summary_df.sort_values('fraction_OB_sensitive', ascending=False).head(6)



# for _, row in top6_df.iterrows():
#     ttype   = row['T-type_Label']
#     total   = int(row['total_cells'])
#     n_sig   = int(row['num_OB_sensitive'])
#     n_nonsig = total - n_sig

#     # -------- Pie chart ------------------------------------------------------
#     fig, ax = plt.subplots(figsize=(6, 6))
#     ax.pie([n_sig, n_nonsig],
#            labels=['', ''],
#            autopct='%1.1f%%', startangle=90, counterclock=False,
#            colors=['red', 'lightgray'])
#     ax.set_title(f'{ttype} (n={total})')
#     ax.axis('equal')

#     # -------- Safe filename  -------------------------------------------------
#     safe_ttype = ttype.replace(' ', '_').replace('/', '_')
#     fname_base = './fig/'+dataset_str+f'Ttype_{safe_ttype}_abundancy'

#     fig.savefig(fname_base + '.jpg',  bbox_inches='tight')
#     fig.savefig(fname_base + '.svg',  bbox_inches='tight')
#     plt.close(fig)          # free memory if you’re looping over many types