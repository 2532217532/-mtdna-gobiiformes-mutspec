## Plot output stacked barplots for mSigAct

import pandas as pd
from pandas.api.types import CategoricalDtype
from plotActivity import plotActivity
import re

import warnings
warnings.filterwarnings("ignore")

set_order = ['high-Actinopteri', 'low-Actinopteri', 'diff-Actinopteri', 
             'high-Amphibia', 'low-Amphibia', 'diff-Amphibia',
             'high-Lepidosauria', 'low-Lepidosauria', 'diff-Lepidosauria', 
             'high-Mammalia', 'low-Mammalia', 'diff-Mammalia',
             'high-Aves', 'low-Aves', 'diff-Aves'] ### order in plot

cat_type = CategoricalDtype(categories=set_order, ordered=True) # sorting help


pt_del = r'__Ts\.\.\.Tv|__Ts\.only' ### set pattern do delete from samples

# without priors
# total Ts and Tv
inpath  = "./data/mSigAct/output/all_relatable_sbs_prop1_Activities.txt"
outpath = "./data/mSigAct/output/figures/all_relatable_sbs_prop1_Activities.pdf"

df = pd.read_csv(inpath, sep='\t')
print(df.shape)

df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])
df['WithTv'] = df.Samples.str.contains('Tv')

df['Samples_Sort'] = df['Samples'].str.replace(pt_del, '', regex=True).str.replace('_', '-')
df['Samples'] = df['Samples'].str.replace('.only', '').str.replace(r'\.\.\.', ' & ')

df['Samples_Sort'] =df['Samples_Sort'].astype(cat_type)
df = df.sort_values(['Samples_Sort', 'Set']).drop(['Set', 'WithTv', 'Samples_Sort'], 1)
df = df.set_index('Samples').astype(int)

# drop signatures with low fraction (< 0.5%)
df = df.loc[:, (df.T / df.sum(1)).T.mean() >= 0.005]
print(df.shape)
df.to_csv(inpath.replace('.txt', '') + '.sorted.txt', sep='\t')

plotActivity(
    inpath.replace('.txt', '') + '.sorted.txt', 
    outpath, 
    bin_size=30, 
    delimiter_step=6, delimiter_size=1,
)



# without priors
# only Ts
inpath  = "./data/mSigAct/output/all_relatable_sbs_prop1_Activities.txt"
outpath = "./data/mSigAct/output/figures/all_relatable_sbs_prop1_Activities_Ts.pdf"

df = pd.read_csv(inpath, sep='\t')
print(df.shape)
df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])
df['WithTv'] = df.Samples.str.contains('Tv')
df['Samples'] = df['Samples'].str.replace('__Ts.only', '').str.replace('_', '-')
df = df[~df['Samples'].str.contains(r'\.\.\.')]

df['Samples'] =df['Samples'].astype(cat_type)
df = df.sort_values(['Samples']).drop(['Set', 'WithTv'], 1)

df = df.set_index('Samples').astype(int)

# drop signatures with low fraction (< 0.5%)
df = df.loc[:, (df.T / df.sum(1)).T.mean() >= 0.005]
print(df.shape)
df.to_csv(inpath.replace('.txt', '') + '.sorted.ts.txt', sep='\t')

plotActivity(
    inpath.replace('.txt', '') + '.sorted.ts.txt', 
    outpath, 
    bin_size=50, 
    delimiter_step=3, delimiter_size=1,
    rename=True, figure_width=4,
)



### When using priors from SigProfiller output

inpath  = "./data/mSigAct/output/custom_prop_Activities.txt"
outpath = "./data/mSigAct/output/figures/custom_prop_Activities.pdf"

df = pd.read_csv(inpath, sep='\t')
df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])

df['Samples_Sort'] = df['Samples'].str.replace(pt_del, '', regex=True).str.replace('_', '-')
df['Samples'] = df['Samples'].str.replace('.only', '').str.replace(r'\.\.\.', ' & ')
df['Samples_Sort'] =df['Samples_Sort'].astype(cat_type)

df = df.sort_values(['Samples_Sort', 'Set']).drop(['Set','Samples_Sort'], 1)
df = df.set_index('Samples', append=True).astype(int).reset_index(1)

df.to_csv(inpath.replace('.txt', '') + '.sorted.txt', index=False, sep='\t')

plotActivity(
    inpath.replace('.txt', '') + '.sorted.txt', 
    outpath, 
    bin_size=30, 
    delimiter_step=6, delimiter_size=1,
)



### When using priors from first mSigAct output without priors

inpath  = "./data/mSigAct/output/custom_from1_Activities.txt"
outpath = "./data/mSigAct/output/figures/custom_from1_Activities.pdf"

df = pd.read_csv(inpath, sep='\t')
df['Set'] = df.Samples.apply(lambda x: x.split('_')[0])

df['Samples_Sort'] = df['Samples'].str.replace(pt_del, '', regex=True).str.replace('_', '-')
df['Samples'] = df['Samples'].str.replace('.only', '').str.replace(r'\.\.\.', ' & ')
df['Samples_Sort'] =df['Samples_Sort'].astype(cat_type)

df = df.sort_values(['Samples_Sort', 'Set']).drop(['Set','Samples_Sort'], 1)
df = df.set_index('Samples', append=True).astype(int).reset_index(1)

df.to_csv(inpath.replace('.txt', '') + '.sorted.txt', index=False, sep='\t')

plotActivity(
    inpath.replace('.txt', '') + '.sorted.txt', 
    outpath, 
    bin_size=30, 
    delimiter_step=6, delimiter_size=3,
)