import pandas as pd
import numpy as np
from collections import OrderedDict as Odict
import matplotlib.pyplot as plt
import collections
import sys 
import pdb
import random
import os
'''
Add noise to input data - experiment table
- shuffle columns internaly
- shuffle rows RPKM values
- output 4 combinations, row, cols and mixed of both
	and one repeated shuffle.

'''
def plot_box(df, path, name, axis):

	print 'Generating boxplot....'
	med = df.median()
	med.sort_values(inplace=True)
	newdf = df[med.index]
	fig=plt.figure(figsize=(1800/100, 1200/100), dpi=100)
	newdf.boxplot(rot=90,  patch_artist=True)
	plt.draw()
	# plt.show()
	
	file_name ='boxplot_'+ axis + '_' + name +'.png'
	figure_file_path = os.path.join(path, file_name)
	fig.savefig(figure_file_path, dpi=100, bbox_inches="tight")
	plt.cla()	
	plt.close(fig)

	return


def save_tables(column_shuff, row_shuff, all_shuff, all_shuff2, path, name):

	print 'Saving tables....'
	column_shuff.to_csv(os.path.join(path, name + '_column_shuff.csv'), sep='\t')
	row_shuff.to_csv(os.path.join(path, name + '_row_shuff.csv'), sep='\t')
	all_shuff.to_csv(os.path.join(path, name + '_all_shuff.csv'), sep='\t')
	all_shuff2.to_csv(os.path.join(path, name + '_all2_shuff.csv'), sep='\t')

	return

def shuffle_row(df):
	df = df.copy()
	for index, row in df.iterrows():
		df.loc[index] = np.random.permutation(df.loc[index])
	return df

def main(exp_table_file, task):

	# exp_table = pd.DataFrame({'A':np.random.rand(10), 'B':range(10), 'C':np.random.rand(10) })
	exp_table = pd.read_table(exp_table_file, index_col=0, low_memory=False)
	index_org = list(exp_table.index)

	index = [i for i in index_org if 'RPKM' in i]
	exp_table_rpkm = exp_table.loc[index]

	# pdb.set_trace()

	path = os.path.dirname(os.path.abspath(exp_table_file))
	name = os.path.basename(exp_table_file)[:-4]
	exp_table_rpkm.to_csv(os.path.join(path, name + '_rpkm.csv'), sep='\t')

	if 'generate' in task:

		column_shuff = exp_table_rpkm.sample(frac=1)
		column_shuff.index = index
		column_shuff.index.name = 'ID'

		row_shuff = shuffle_row(exp_table_rpkm)

		all_shuff = row_shuff.sample(frac=1)
		all_shuff.index = index
		all_shuff.index.name = 'ID'

		all_shuff2 = shuffle_row(all_shuff)
		all_shuff2 = all_shuff2.sample(frac=1)
		all_shuff2.index = index
		all_shuff2.index.name = 'ID'

		save_tables(column_shuff, row_shuff, all_shuff, all_shuff2, path, name)

		if 'plot' in task:
			plot_box(exp_table_rpkm, path, name, 'sample')
			tmp = exp_table_rpkm.transpose()
			plot_box(tmp, path, name, 'rpkm')
			plot_box(column_shuff, path, 'exp_table_column_shuff', 'sample')
			tmp = column_shuff.transpose()
			plot_box(tmp, path, 'exp_table_column_shuff', 'rpkm')
			plot_box(row_shuff, path, 'exp_table_row_shuff', 'sample')
			tmp = row_shuff.transpose()
			plot_box(tmp, path, 'exp_table_row_shuff', 'rpkm')
			plot_box(all_shuff, path, 'exp_table_all_shuff', 'sample')
			tmp = all_shuff.transpose()
			plot_box(tmp, path, 'exp_table_all_shuff', 'rpkm')
			plot_box(all_shuff2, path, 'exp_table_all2_shuff', 'sample')
			tmp = all_shuff2.transpose()
			plot_box(tmp, path, 'exp_table_all2_shuff', 'rpkm')
	
	elif task == 'plot':
		plot_box(exp_table_rpkm, path, name, 'sample')
		tmp = exp_table_rpkm.transpose()
		plot_box(tmp, path, name, 'rpkm')

	else:
		print 'Please select as input one of the different tasks: generate, plot or generate_plot!'
	# pdb.set_trace()
	return

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])