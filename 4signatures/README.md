# Decomposition of mtDNA spectra to COSMIC signatures

## Workflow 

1. Need to renormalize chordates classes spectra to human genome (see details in article methods) - [py script](./0count_human_triplets_freqs.py)
2. Run SigProfilerAssignment cosmic_fit function to run signatures assignment, aggregate results and plot figures - [ipy notebook](./1signatures_analysis_sigpro.ipynb)

We splitted 192-component spectra to parts based on transitions freqs: Low High Diff (see details in article methods). Besides we renormalized spectra to emulate human genome mutation counts

3. To run another tool mSigAct we prepared priors for signatures distribution in the data based on SigProfilerAssignment results - [py script](./2prepare_priors_for_mSigAct.py)
4. Run mSigAct MAPAssignActivity func to decompose spectra by another approach - [R script](./3mSigAct_analysis.R). Before running custom_from1 way with mSigAct use [ipy notebook](./4aggregate_mSigAct_outputs.ipynb) to calculate proportions of SBS and take 6 most coomon of them. 

We executed mSigAct in several ways to get stable results:
- custom_prop - using priors for signatures calculated using SigProfilerAssignment output (see previous step)
- prop1 - priors for signatures are uniform and equal to 1
- custom_from1 - using priors for valuable signatures calculated using *prop1* output

5. Aggregation and little EDA of mSigAct results - [ipy notebook](./4aggregate_mSigAct_outputs.ipynb)
6. Plot figures for mSigAct output - [py script](./5plot_output_for_mSigAct.py)


## Scripts and notebooks

- [0count_human_triplets_freqs.py](./0count_human_triplets_freqs.py) - renormalize chordates classes spectra to human genome
- [./1signatures_analysis_sigpro.ipynb](./1signatures_analysis_sigpro.ipynb) - run SigProfiletAssignment analysis
- [2prepare_priors_for_mSigAct.py](./2prepare_priors_for_mSigAct.py) - prepare SigProfiler outputs for mSigAct running 
- [3mSigAct_analysis.R](./3mSigAct_analysis.R) - run mSigAct on renormolized spectra in three ways
- [4aggregate_mSigAct_outputs.ipynb](./4aggregate_mSigAct_outputs.ipynb) - calculate prop of the SBS from prop1 in mSigAct and then make table of all test from mSigAct
- [5plot_output_for_mSigAct.py](./5plot_output_for_mSigAct.py) - draw mSigAct output in COSMIC format
- [plotActivity.py](./plotActivity.py) - visualization of barplots in the SigProfiler way
- [utils.py](./utils.py) - some functions for data processing

## Materials in the article

**Figures:**

- Figure 3c [1](./data/SigProfilerAssignment/output/only_Ts.pdf) and [2](./data/mSigAct/output/figures/all_relatable_sbs_prop1_Activities_Ts.pdf)
- Supplementary Figure 12 [a](./data/SigProfilerAssignment/output/total.pdf) and [b](./data/mSigAct/output/figures/all_relatable_sbs_prop1_Activities.pdf)

**Tables:**

- [Supplementary Table 6](./data/SigProfilerAssignment/output/Solution_Stats.txt)
- [Supplementary Table 7](./data/mSigAct/output/Distances_mSigAct.csv)

## Used Data

- [192-comp mutational spectra data](../1data_derivation/dataset/MutSpecVertebrates192.csv.gz)

### External Data

- [Human Genome](./data/human_genome/download_human_genome.sh) - bash script to download human genome