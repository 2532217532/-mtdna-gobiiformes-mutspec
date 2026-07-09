#!/usr/bin/env
rm(list=ls(all=TRUE))

library(ICAMS)
library(mSigAct)
library(cosmicsig)
library(dplyr)

path = './data/SigProfilerAssignment/input'
samples_names = list.files(path, pattern='_samples.txt', full.names=T)
samples_names

diff_df = read.table(samples_names[1], sep='\t', header =T)
high_df = read.table(samples_names[2], sep='\t', header =T)
low_df = read.table(samples_names[3], sep='\t', header =T)

diff_df = diff_df %>% 
          rename_at(2:ncol(diff_df),~ paste0('diff', '_', .))

high_df = high_df %>% 
  rename_at(2:ncol(high_df),~ paste0('high', '_', .)) 

low_df = low_df %>% 
  rename_at(2:ncol(low_df),~ paste0('low', '_', .)) 

low_df = low_df[,-1]
high_df = high_df[,-1]

final_catalog = cbind(diff_df, high_df, low_df)

write.table(final_catalog, './data/mSigAct/input/samples_mSigAct.txt', sep='\t', row.names = F)

path_full_catalog = './data/mSigAct/input/samples_mSigAct.txt'
input_catalog <- ICAMS::ReadCatalog(file = path_full_catalog)


sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96

sig_use= c('SBS2', 'SBS26', 'SBS19', 'SBS23', 'SBS44', 'SBS5', 'SBS12', 'SBS30')

sigs_to_use <- sigs[, colnames(sigs) %in% sig_use]


### Proportions are taken from SigProfiler

output_home <- "./data/mSigAct/output/raw_output/custom_prop"

sig_prop = c(0.01, 0.02, 0.03, 0.03, 0.03, 0.23, 0.31, 0.35)
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 1, 
                             mc.cores.per.sample = 4)

### Take all relatable SBS with 1 prop 
output_home <- "./data/mSigAct/output/raw_output/all_relatable_sbs_prop1"


sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96

sig_to_delete= c('SBS32', 'SBS11', 'SBS25', 'SBS31', 'SBS32', 'SBS35', 'SBS86',
                 'SBS87', 'SBS90', 'SBS22', 'SBS88', 'SBS27', 'SBS43', 'SBS45',
                 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52',
                 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59',
                 'SBS60', 'SBS95', 'SBS9', 'SBS84', 'SBS85', 'SBS4', 'SBS29',
                 'SBS92', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS38')

sigs_to_use = sigs[, !(colnames(sigs) %in% sig_to_delete)] 

sig_prop = rep(1, ncol(sigs_to_use))
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             max.subsets = 100000,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 4, 
                             mc.cores.per.sample = 8)

### Before running this part go to ./4aggregate_mSigAct_outputs.ipynb script and calculate proportions from mSigAct 1 prop
### Proportions from mSigAct with 1 prop

output_home <- "./data/mSigAct/output/raw_output/custom_from1"

sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96

sig_use= c('SBS12', 'SBS30', 'SBS23', 'SBS26', 'SBS2', 'SBS42', 'SBS21')

sigs_to_use <- sigs[, colnames(sigs) %in% sig_use]

sig_prop = c(0.319994, 0.292430, 0.197766, 0.052957, 0.033872, 0.019849, 0.017117)
names(sig_prop) = colnames(sigs_to_use)

retval <-
  mSigAct::MAPAssignActivity(spectra = input_catalog, 
                             sigs = sigs_to_use,
                             sigs.presence.prop = sig_prop,
                             output.dir = output_home,
                             max.level = ncol(sigs_to_use) - 1,
                             max.subsets = 100000,
                             p.thresh = 0.05 / ncol(sigs_to_use), 
                             num.parallel.samples = 1, 
                             mc.cores.per.sample = 24)

