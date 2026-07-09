from collections import defaultdict
import json

import numpy as np
import pandas as pd

from pymutspec.annotation import rev_comp
from pymutspec.constants import possible_sbs192, possible_sbs12, possible_sbs96, possible_codons


def _calc_mutspec_class_legacy(df: pd.DataFrame, group_col="Class"):
    ms_cls = df.groupby([group_col, 'Mut'])['RawMutSpec'].sum().reset_index()
    ms_cls["RawMutSpecSum"] = ms_cls.Class.map(ms_cls.groupby("Class").RawMutSpec.sum().to_dict())
    ms_cls['MutSpec'] = ms_cls.RawMutSpec / ms_cls.RawMutSpecSum
    ms_cls = ms_cls.drop(['RawMutSpec', 'RawMutSpecSum'], axis=1)
    return ms_cls


def calc_mutspec_class(df: pd.DataFrame, group_col="Class"):
    assert not df.MutSpec.isna().any()
    ms_cls = df.groupby([group_col, 'Mut'])['MutSpec'].mean().reset_index()
    return ms_cls


def collapse_sbs192(df: pd.DataFrame, to=12):
    assert (df.columns == possible_sbs192).all()
    df = df.copy()
    if to == 12:
        for sbs192 in possible_sbs192:
            sbs12 = sbs192[2:5]
            if sbs12 in df.columns.values:
                df[sbs12] += df[sbs192]
            else:
                df[sbs12] = df[sbs192]

        return df[possible_sbs12]
    elif to == 96:
        for sbs96 in possible_sbs96:
            sbs96_rev = rev_comp(sbs96)
            df[sbs96] = df[sbs96] + df[sbs96_rev]
        return df[possible_sbs96]
    else:
        raise NotImplementedError()


def complete_sbs_columns(df: pd.DataFrame, ncomp):
    if ncomp == 96:
        assert df.shape[1] < 96
        possible_sbs = possible_sbs96
    elif ncomp == 192:
        possible_sbs = possible_sbs192
    else:
        raise NotImplementedError
    
    df = df.copy()
    if len(df.columns) != 192:
        for sbs in set(possible_sbs).difference(df.columns.values):
            df[sbs] = 0.
    df = df[possible_sbs]
    return df


effect2sbs = {
    'MMR_deficiency': ['SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44', 'DBS7', 'DBS10', 'ID7'], 
    'POL_deficiency': ['SBS10a', 'SBS10b', 'SBS10c', 'SBS10d', 'SBS28', 'DBS3'], 
    'HR_deficiency': ['SBS3', 'ID6'], 
    'BER_deficiency': ['SBS30', 'SBS36'], 
    'Chemotherapy': ['SBS11', 'SBS25', 'SBS31', 'SBS35', 'SBS86', 'SBS87', 'SBS90', 'DBS5'], 
    'Immunosuppressants': ['SBS32'], 
    'Treatment': ['SBS11', 'SBS25', 'SBS31', 'SBS32', 'SBS35', 'SBS86', 'SBS87', 'SBS90', 'DBS5'], 
    'APOBEC': ['SBS2', 'SBS13'], 
    'Tobacco': ['SBS4', 'SBS29', 'SBS92', 'DBS2', 'ID3'], 
    'UV': ['SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS38', 'DBS1', 'ID13'], 
    'AA': ['SBS22'], 
    'Colibactin': ['SBS88', 'ID18'], 
    'Artifact': ['SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'SBS95'], 
    'Lymphoid': ['SBS9', 'SBS84', 'SBS85']
}
sbs2effect = {sbs: ef for ef,sbs_lst in effect2sbs.items() for sbs in sbs_lst}
sbs2effect["SBS18"] = "ROS"
sbs2effect["SBS1"] = "C_deamination"
sbs2effect["SBS5"] = "UNK_clock_like"
