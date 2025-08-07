import pandas as pd, numpy as np
import os
import glob
mdpath = "data/gamb/metadata"

sets = [i for i in os.listdir(f'{mdpath}/general') if os.path.isdir(f'{mdpath}/general/{i}')]
metadata = pd.DataFrame()
for s in sets:
    df = pd.read_csv(f'{mdpath}/general/{s}/samples.meta.csv')
    df['set'] = s
    metadata = pd.concat([metadata, df])
from functools import reduce
cohorts = ['data/gamb/metadata/cohorts_20250502']
cohortdata = pd.DataFrame()
for c in cohorts:
    for s in os.listdir(c):
        df = reduce(lambda left, right: pd.merge(left, right, on='sample_id', how='outer'), 
                    [pd.read_csv(f'{c}/{s}/{i}') for i in os.listdir(f'{c}/{s}') if i.endswith('.csv')])
        cohortdata = pd.concat([cohortdata, df])
pd.merge(metadata, cohortdata, on='sample_id', how='outer').to_csv('data/gamb/metadata/merged_metadata.csv', index=False
