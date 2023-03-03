import subprocess
import pandas as pd


filtered_df = pd.read_csv('filtered_df_grouped.csv')


executable_list = []


for index, row_data in filtered_df.iterrows():

	row_data = row_data.to_frame()  # series to frame
	row_data = row_data.T

	pdb = row_data['pdb_id'].item().lower()
	target = row_data['target'].item()

	import re
	s = re.sub(r'[()\[\],\'"]', ' ', target)
	target_list = s.split(' ')
	lst = [s for s in target_list if s]

	# print(pdb, lst)
	#
	executable_path = "../foldx5MacC11.tar_/./foldx_20231231"
	argument1 = "--command=PositionScan"
	argument2 = "--pdb="+pdb+".pdb"
	argument3 = "--pdb-dir=/Users/natalieso/Downloads/20230217_0084705/"
	argument4 = "--positions="+','.join(lst)

	executable_list.append([executable_path, argument1, argument2, argument3, argument4])
# print(executable_list)
for executable in executable_list:
	subprocess.run(executable)
