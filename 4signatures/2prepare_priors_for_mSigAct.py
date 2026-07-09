import pandas as pd

df = pd.read_csv('./data/SigProfilerAssignment/output/Assignment_Solution_Activities.txt', sep='\t', index_col=0)

# calculate freqs of all SBS in SigProfiller output averaged for all runs
sbs_priors = df.sum(axis=0)
sbs_priors = (sbs_priors[sbs_priors > 0] / sbs_priors.sum()).round(2)
sbs_priors = sbs_priors[sbs_priors >= 0.01] # frequent then 1%

print('Priors based on SigProfiller output')
print(sbs_priors.sort_values())
