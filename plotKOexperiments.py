import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse, ast
from argparse import RawTextHelpFormatter
from os.path import isfile, join
import re
import pandas as pd
import numpy as np
import itertools
import pdb
import matplotlib.pyplot as plt
import seaborn as sns
import math

swap_names = {
	'TUMOR-11-28-N-GEM-00-C-L-09-Mean-010203':'C11-28-L-09-wGF-w/oDrug',
	'TUMOR-11-28-N-GEM-00-C-L-09-Mean-010203-MK-2206-Conc1000.0nM': 'C11-28-L-09-wGF-wMK-2206',
	'TUMOR-11-28-N-GEM-00-C-L-09-Mean-010203-Wortmannin-Conc5000.0nM': 'C11-28-L-09-wGF-wWortmannin',
	'TUMOR-11-28-N-GEM-00-C-L-13-Mean-010203': 'C11-28-L-13-w/oGF-w/oDrug',
	'TUMOR-11-28-N-GEM-00-C-L-13-Mean-010203-MK-2206-Conc1000.0nM': 'C11-28-L-13-w/oGF-wMK-2206',
	'TUMOR-11-28-N-GEM-00-C-L-13-Mean-010203-Wortmannin-Conc5000.0nM':'C11-28-L-13-w/oGF-wWortmannin',
	'TUMOR-13-29-N-GEM-00-C-L-01-Mean-010203': 'C13-29-L-01-wGF-w/oDrug',
	'TUMOR-13-29-N-GEM-00-C-L-01-Mean-010203-MK-2206-Conc1000.0nM': 'C13-29-L-01-wGF-wMK-2206',
	'TUMOR-13-29-N-GEM-00-C-L-01-Mean-010203-Wortmannin-Conc5000.0nM': 'C13-29-L-01-wGF-wWortmannin',
	'TUMOR-13-29-N-GEM-00-C-L-07-Mean-010203': 'C13-29-L-07-w/oGF-w/oDrug',
	'TUMOR-13-29-N-GEM-00-C-L-07-Mean-010203-MK-2206-Conc1000.0nM':	'C13-29-L-07-w/oGF-wMK-2206',
	'TUMOR-13-29-N-GEM-00-C-L-07-Mean-010203-Wortmannin-Conc5000.0nM': 'C13-29-L-07-w/oGF-wWortmannin',
}

file = sys.argv[1]

df = pd.read_table(file, sep = '\t', low_memory=False)
temp = df[(df['fold_change']>1.3) | (df['fold_change']<0.8)][['Knockout of', 'Sample', 'fold_change']].sort_values(by=['fold_change'])
zeros = {0:10e-8}
temp = temp.replace(zeros)
# temp = df[['Knockout of', 'Sample', 'fold_change']].sort_values(by=['fold_change'])
# pdb.set_trace()
dic = temp.to_dict()
newi = set(dic['Sample'].values())
newc = set(dic['Knockout of'].values())
z = {k:{v: np.log2(temp[(temp['Sample']==v) & (temp['Knockout of']==k)]['fold_change'].values[0]) if temp[(temp['Sample']==v) & (temp['Knockout of']==k)]['fold_change'].values else 1.0 for v in newi} for k in newc}

df_z = pd.DataFrame(z)
df_z=df_z.reindex_axis(df_z.mean().sort_values().index, axis=1)

ax = sns.heatmap(df_z, linewidths=.5, cmap='RdBu_r', center=0, vmin = -10, vmax=10)
ax.set_yticklabels(ax.get_yticklabels(),rotation=0) 
ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
ax.collections[0].colorbar.set_label("log2(fold change)")



plt.savefig(file[:file.rindex('/')+1] +'figure_'+ file[file.rindex('/')+1:file.rindex('.')] +'.png', bbox_inches='tight', dpi=200)
