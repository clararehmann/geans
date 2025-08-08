import pandas as pd, numpy as np
import os, shutil, glob

# anopheles gambiae complex
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
pd.merge(metadata, cohortdata, on='sample_id', how='outer').to_csv('data/gamb/metadata/merged_metadata.csv', index=False)

# now separate out cohorts and admin groups
# gambiae
cohortdata = pd.read_csv('data/gamb/metadata/merged_metadata.csv')
gambiae_data = cohortdata[cohortdata['taxon'] == 'gambiae']
for s in set(gambiae_data['set']):
    df = gambiae_data[gambiae_data['set'] == s]
    df.to_csv(f'data/gamb/metadata/set_{s}_gamb.txt', index=False, sep='\t')
for c in set(gambiae_data['cohort_admin1_year']):
    df = gambiae_data[gambiae_data['cohort_admin1_year'] == c]
    df.to_csv(f'data/gamb/metadata/admin1_year_{c}.txt', index=False, sep='\t')
# coluzzii
os.makedirs('data/colu/metadata', exist_ok=True)
colu_data = cohortdata[cohortdata['taxon']=='coluzzii']
colu_data.to_csv('data/colu/metadata/colu_metadata_all.txt', index=False, sep='\t')
for s in set(colu_data['set']):
    df = colu_data[colu_data['set'] == s]
    df.to_csv(f'data/colu/metadata/set_{s}_colu.txt', index=False, sep='\t')
for c in set(colu_data['cohort_admin1_year']):
    df = colu_data[colu_data['cohort_admin1_year'] == c]
    df.to_csv(f'data/colu/metadata/admin1_year_{c}.txt', index=False, sep='\t')
# arabiensis
os.makedirs('data/arab/metadata', exist_ok=True)
arab_data = cohortdata[cohortdata['taxon']=='arabiensis']
arab_data.to_csv('data/arab/metadata/arab_metadata_all.txt', index=False, sep='\t')
for s in set(arab_data['set']):
    df = arab_data[arab_data['set'] == s]
    df.to_csv(f'data/arab/metadata/set_{s}_arab.txt', index=False, sep='\t')
for c in set(arab_data['cohort_admin1_year']):
    df = arab_data[arab_data['cohort_admin1_year'] == c]
    df.to_csv(f'data/arab/metadata/admin1_year_{c}.txt', index=False, sep='\t')

# anopheles funestus
os.makedirs('data/afun/metadata', exist_ok=True)
mdpath = "/sietch_colab/crehmann/vo_afun_release/v1.0/metadata"
sets = [i for i in os.listdir(f'{mdpath}') if os.path.isdir(f'{mdpath}/{i}')]
metadata = pd.DataFrame()
for s in sets:
    df = pd.read_csv(f'{mdpath}/{s}/samples.meta.csv')
    df['set'] = s
    metadata = pd.concat([metadata, df])
metadata.to_csv('data/afun/metadata/afun_metadata_all.txt', index=False, sep='\t')
# now separate out cohorts and admin groups
for s in set(metadata['set']):
    df = metadata[metadata['set'] == s]
    df.to_csv(f'data/afun/metadata/set_{s}_afun.txt', index=False, sep='\t')

# anopheles minimus
os.makedirs('data/amin/metadata', exist_ok=True)
mdpath = "/sietch_colab/crehmann/vo_amin_release/v1/metadata"
shutil.copyfile(f'{mdpath}/samples.meta.csv', 'data/amin/metadata/amin_metadata_all.txt')