import pandas as pd
import glob

from os import path

def load(pattern='*LT*.csv', directory='results'):
    '''load df's from directory, average the columns and return a new df'''
    pattern = path.join(directory, pattern)
    print(glob.glob(pattern))
    return
    results = list()
    for df in [pd.read_csv(f) for f in glob.glob(pattern)]:
        result = dict()
        for column in df:
            result[column] = df[column].mean()
        results.append(result)
    return pd.DataFrame(results)

load()
