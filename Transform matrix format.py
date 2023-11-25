import os
import sys
import numpy as np
import pandas as pd

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')


dirs = os.listdir("./input/")

print(dirs)

for s in dirs:
    path = './input/'+s
    df = pd.read_csv(path,header=0)


    df["seqnames_start_end"] =df["seqnames"].map(str)+"_"+ df["start"].map(str)+"_"+ df["end"].map(str)

    del df[df.keys()[0]]
#df.drop(axis=0)
    del df['width']
    del df['strand']
    del df['idx']
    del df['seqnames']
    del df['start']
    del df['end']

    d = df.pop('seqnames_start_end')
    df.insert(0,'seqnames_start_end',d)

    print(df)

    df.to_csv('./output/'+s,index=False)
