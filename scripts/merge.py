import pandas as pd
import argparse
import os
import glob

parser = argparse.ArgumentParser(description="Merge files on first column.")
parser.add_argument('-i', '--input', dest='input', help='Folder containing files to merge.')
parser.add_argument('-s', '--suf', dest='suf', help='Suffix')
parser.add_argument('-o', '--output', dest='output', help='Output.')
#parser.add_argument('-m', '--merge', dest='merge', help='Inner or outer.')
args = parser.parse_args()

files=glob.glob(os.path.join(args.input, "*"+args.suf))
print(files)
data=pd.DataFrame()
for f in files:
    ff=f.split('/')[-1]
    data2 = pd.read_table(f,sep='\t',header=0,index_col=None)
    data3 = data2.transpose()
    data3.columns=[ff]
    data=data.merge(data3, how='outer', left_index=True, right_index=True).fillna('NA')

data.index.name='Rows'
data.to_csv(args.output+'metrics.txt',sep='\t',header=True,index=True)
