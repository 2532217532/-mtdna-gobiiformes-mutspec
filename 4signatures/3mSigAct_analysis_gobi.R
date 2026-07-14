#!/usr/bin/env
# mSigAct analysis for Gobiiformes mtDNA mutational spectra
# Adapted from the chordata project template.
#
# PREREQUISITE: Step 3 notebook (1signatures_analysis_gobi.ipynb) must be run first
# to generate the SigProfilerAssignment input files and output.
#
# This script runs mSigAct in 3 modes:
#   1. custom_prop: priors from SigProfilerAssignment output
#   2. all_relatable_sbs_prop1: uniform priors on all relatable (non-artifact) SBS
#   3. custom_from1: priors from mSigAct prop1 output (run after 4aggregate notebook)

rm(list=ls(all=TRUE))

library(ICAMS)
library(mSigAct)
library(cosmicsig)
library(dplyr)

# ============================================================
# Load input spectra (96-component, renormalized to human genome)
# ============================================================
path = './data/SigProfilerAssignment/input'
samples_names = list.files(path, pattern='_samples.txt', full.names=T)
samples_names

diff_df = read.table(samples_names[1], sep='\t', header=T)
high_df = read.table(samples_names[2], sep='\t', header=T)
low_df  = read.table(samples_names[3], sep='\t', header=T)

diff_df = diff_df %>%
          rename_at(2:ncol(diff_df),~ paste0('diff', '_', .))

high_df = high_df %>%
  rename_at(2:ncol(high_df),~ paste0('high', '_', .))

low_df = low_df %>%
  rename_at(2:ncol(low_df),~ paste0('low', '_', .))

low_df  = low_df[,-1]
high_df = high_df[,-1]

final_catalog = cbind(diff_df, high_df, low_df)

write.table(final_catalog, './data/mSigAct/input/samples_mSigAct.txt', sep='\t', row.names=F)

path_full_catalog = './data/mSigAct/input/samples_mSigAct.txt'
input_catalog <- ICAMS::ReadCatalog(file=path_full_catalog)


# ============================================================
# Run 1: custom_prop — priors from SigProfilerAssignment
# ============================================================
# *** UPDATE AFTER running Step 3 notebook + Step 4 priors script ***
# These signature IDs and proportions come from 2prepare_priors_for_mSigAct.py output.
# Replace with the actual output after running that script.
#
# Template (from chordata project):
#   sig_use  = c('SBS2', 'SBS26', 'SBS19', 'SBS23', 'SBS44', 'SBS5', 'SBS12', 'SBS30')
#   sig_prop = c(0.01, 0.02, 0.03, 0.03, 0.03, 0.23, 0.31, 0.35)
# *** UPDATE ABOVE after Step 4 ***

sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96

# TODO: Replace with actual Gobiiformes SigProfilerAssignment output
sig_use  = c('SBS2', 'SBS26', 'SBS19', 'SBS23', 'SBS44', 'SBS5', 'SBS12', 'SBS30')
sigs_to_use <- sigs[, colnames(sigs) %in% sig_use]

output_home <- "./data/mSigAct/output/raw_output/custom_prop"

# TODO: Replace with actual proportions from Step 4
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


# ============================================================
# Run 2: all_relatable_sbs_prop1 — uniform priors
# ============================================================
# Excludes artifact, treatment, immunosuppressant, lymphoid, colibactin,
# tobacco, and UV signatures — these are not expected in mtDNA.
# Same exclusion list as the chordata template.
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


# ============================================================
# Run 3: custom_from1 — priors from mSigAct prop1 output
# ============================================================
# *** RUN AFTER Step 6 aggregation notebook which computes top SBS proportions ***
# *** UPDATE these values after running 4aggregate_mSigAct_outputs_gobi.ipynb ***
#
# Template (from chordata project):
#   sig_use  = c('SBS12', 'SBS30', 'SBS23', 'SBS26', 'SBS2', 'SBS42', 'SBS21')
#   sig_prop = c(0.319994, 0.292430, 0.197766, 0.052957, 0.033872, 0.019849, 0.017117)
# *** UPDATE ABOVE after 4aggregate notebook ***

output_home <- "./data/mSigAct/output/raw_output/custom_from1"

sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96

# TODO: Replace with actual top SBS from prop1 output
sig_use = c('SBS12', 'SBS30', 'SBS23', 'SBS26', 'SBS2', 'SBS42', 'SBS21')
sigs_to_use <- sigs[, colnames(sigs) %in% sig_use]

# TODO: Replace with actual proportions
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
