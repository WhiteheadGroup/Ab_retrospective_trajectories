#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 09:46:14 2024

@author: siobhan
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import os

file1 = '/Users/siobhan/PycharmProjects/MAGMA-seq_proccessing/Outputs_Raw MLE/YL009_MLEoutput_variants_bright_adjust.csv'
file2 = '/Users/siobhan/PycharmProjects/MAGMA-seq_proccessing/Outputs_Raw MLE/YL011_MLEoutput_variants_brightAdj.csv'

file1_df = pd.read_csv(file1)
file2_df = pd.read_csv(file2)

successFilter = True
countFilter = True
outlierFilter = False
commonVars = False
CIfilter = False
Kdfilter = False
TrueVarFilter = True

antibodies_with_allowed_mutations = {
    "Ab_1-20_UCA": {
        "VH": ["F27L", "Y58F", "K109Q", "WT"],
        "VL": ["WT", "A25A"]
    },
    "002-S21F2_UCA": {
        "VH": ["G25A", "S29N", "T31A", "S32Y", "T104A", "WT", "Q112Q", "S119S", "L19L", "P15P", "Q1G", "E2Q", "G9G",
               "A10A", "Q40Q", "S62S", "V6V"],
        "VL": ["F73L", "D92K", "N93I", "WT", "L94L", "T5T", "G102G", "K104K", "Y87Y", "Q37Q", "F99F", "C23C", "E106E",
               "T96T"]
    },
    "Ab_2-15_UCA": {
        "VH": ["K19R", "N54I", "G56D", "A79V", "S107N", "V118I", "V125I", "WT", "A9A", "G15G", "K13K", "G49G", "A97A",
               "K23K", "R85R", "S84S", "Y109Y", "K63K", "E82E", "V20V", "A16A", "A24A", "F116F", "A40A", "G122G",
               "L83L", "L86L", "R72R", "S21S", "T126T", "G100G", "V18V"],
        "VL": ["Y34F", "E52D", "N55K", "Y88C", "L99F", "WT", "D28D", "S13S", "R56R", "G79G", "T108T", "N33N", "A82A",
               "P57P", "P14P", "V53V", "L4L", "Q39Q", "S36S", "V35V", "L80L", "S21S", "A8A", "T25T", "A73A", "L99L",
               "G24G", "S9S", "T103T", "S69S", "S91S", "I18I", "S2S", "A45A", "V100V", "T19T", "S78S", "Q40Q", "S67S",
               "S11S", "G70G", "D84D", "P46P", "L75L", "F64F"]
    },
    "Ab_2-7_UCA": {
        "VH": ["WT"],
        "VL": ["A15T", "G31A", "E52D", "N55K", "S95T", "WT", "S13S"]
    },
    "Ab_1-57_UCA": {
        "VH": ["Q3H", "S77A", "T108I", "G109N", "F117S", "WT"],
        "VL": ["WT"]
    },
    "C118_UCA": {
        "VH": ["S31N", "V93I", "S101T", "M114L", "WT", "G120G", "G109G", "V107V", "A97A", "S126S"],
        "VL": ["L2P", "S52T", "WT"]
    },
    "CC631_UCA": {
        "VH": ["K23M", "N52S", "S57G", "M70L", "WT", "G49G"],
        "VL": ["V105L", "WT", "A25A", "N31N", "F99F", "T103T", "K104K", "S63S"]
    },
    "Ab_5-7_UCA": {
        "VH": ["I50V", "Q62E", "Q65R", "G120A", "M121L", "K126Q", "WT"],
        "VL": ["K45E", "S93T", "WT"]
    },
    "CC121_UCA": {
        "VH": ["F27L", "Y58F", "M104L", "WT"],
        "VL": ["Q3V", "L4M", "D107E", "WT", "A25A", "G101G", "Q38Q"]
    },
    # Add more antibodies and their allowed mutations here
}

silent_mutations = {
    "Ab_1-20_UCA": {
        "VH": [],
        "VL": ["A25A"]
    },
    "002-S21F2_UCA": {
        "VH": ["Q112Q", "S119S", "L19L", "P15P", "Q1G", "E2Q", "G9G",
               "A10A", "Q40Q", "S62S", "V6V"],
        "VL": ["L94L", "T5T", "G102G", "K104K", "Y87Y", "Q37Q", "F99F", "C23C", "E106E",
               "T96T"]
    },
    "Ab_2-15_UCA": {
        "VH": ["A9A", "G15G", "K13K", "G49G", "A97A",
               "K23K", "R85R", "S84S", "Y109Y", "K63K", "E82E", "V20V", "A16A", "A24A", "F116F", "A40A", "G122G",
               "L83L", "L86L", "R72R", "S21S", "T126T", "G100G", "V18V"],
        "VL": ["D28D", "S13S", "R56R", "G79G", "T108T", "N33N", "A82A",
               "P57P", "P14P", "V53V", "L4L", "Q39Q", "S36S", "V35V", "L80L", "S21S", "A8A", "T25T", "A73A", "L99L",
               "G24G", "S9S", "T103T", "S69S", "S91S", "I18I", "S2S", "A45A", "V100V", "T19T", "S78S", "Q40Q", "S67S",
               "S11S", "G70G", "D84D", "P46P", "L75L", "F64F"]
    },
    "Ab_2-7_UCA": {
        "VH": [],
        "VL": ["S13S"]
    },
    "Ab_1-57_UCA": {
        "VH": [],
        "VL": []
    },
    "C118_UCA": {
        "VH": ["G120G", "G109G", "V107V", "A97A", "S126S"],
        "VL": []
    },
    "CC631_UCA": {
        "VH": ["G49G"],
        "VL": ["A25A", "N31N", "F99F", "T103T", "K104K", "S63S"]
    },
    "Ab_5-7_UCA": {
        "VH": [],
        "VL": []
    },
    "CC121_UCA": {
        "VH": [],
        "VL": ["A25A", "G101G", "Q38Q"]
    },
    # Add more antibodies and their allowed mutations here
}


def parse_mutations(variant):
    """Extract VH and VL mutations from the variant string."""
    vh_mutations = []
    vl_mutations = []
    sections = variant.split('|')

    for section in sections:
        if '>VH:' in section:
            vh_part = section.split('>VH:')[1]
            for mut in vh_part.split(';'):
                vh_mutations.append(mut.split('-')[0])
        elif '>VL:' in section:
            vl_part = section.split('>VL:')[1]
            for mut in vl_part.split(';'):
                vl_mutations.append(mut.split('-')[0])

    return vh_mutations, vl_mutations


def remove_silent_mutations(vh_mutations, vl_mutations, antibody):
    """Remove silent mutations from the VH and VL mutation lists."""
    if antibody in silent_mutations:
        vh_mutations = [mut for mut in vh_mutations if mut not in silent_mutations[antibody]['VH']]
        vl_mutations = [mut for mut in vl_mutations if mut not in silent_mutations[antibody]['VL']]
    return vh_mutations, vl_mutations

def is_variant_allowed(variant, antibody, allowed_mutations):
    """Check if the variant contains only allowed VH and VL mutations for the given antibody."""
    vh_mutations, vl_mutations = parse_mutations(variant)
    vh_mutations, vl_mutations = remove_silent_mutations(vh_mutations, vl_mutations, antibody)
    vh_allowed = all(mutation in allowed_mutations[antibody]['VH'] for mutation in vh_mutations)
    vl_allowed = all(mutation in allowed_mutations[antibody]['VL'] for mutation in vl_mutations)
    return vh_allowed and vl_allowed

def filter_csv(dataFrame, successFilter, countFilter, CIfilter, Kdfilter, TrueVarFilter, filt):
    if successFilter:
        dataFrame = dataFrame[dataFrame['Success'] != False]

    if countFilter:
        dataFrame = dataFrame[dataFrame['Avg_counts'] >= 100]

    if CIfilter:
        dataFrame = dataFrame[dataFrame['95% CI high'] != 0]

    if Kdfilter:
        dataFrame = dataFrame[dataFrame['Kd'] < 3000]

    if TrueVarFilter:
        def variant_allowed(row):
            antibody = row['Ab']
            variant = row['Variant']
            if antibody in antibodies_with_allowed_mutations:
                return is_variant_allowed(variant, antibody, antibodies_with_allowed_mutations)
            return False

        dataFrame = dataFrame[dataFrame.apply(variant_allowed, axis=1)]

    return dataFrame


file1_filtered = filter_csv(file1_df, successFilter, countFilter, CIfilter, Kdfilter, TrueVarFilter)

file2_filtered = filter_csv(file2_df, successFilter, countFilter, CIfilter, Kdfilter, TrueVarFilter)

# Merge filtered DataFrames
common_variants = pd.merge(file1_filtered, file2_filtered, on='Variant', suffixes=('_file1', '_file2'))

# Calculate the Spearman rank correlation
spearman_corr, _ = spearmanr(common_variants['Kd_file1'], common_variants['Kd_file2'])

# Plot Kd values for common variants on log axes
plt.figure(figsize=(10, 6))
plt.scatter(common_variants['Kd_file1'], common_variants['Kd_file2'], alpha=0.7, label='Common Variants')

# Find rows where Kd values differ by more than 10x
diff_variants = common_variants[(common_variants['Kd_file1'] > 3 * common_variants['Kd_file2']) |
                                (common_variants['Kd_file2'] > 3 * common_variants['Kd_file1'])]

# Highlight the differing variants with red dots
plt.scatter(diff_variants['Kd_file1'], diff_variants['Kd_file2'], color='red', alpha=0.7,
            label='Variants > 10x Kd Difference')

# Add a line y = x
max_kd = max(common_variants['Kd_file1'].max(), common_variants['Kd_file2'].max())
min_kd = min(common_variants['Kd_file1'].min(), common_variants['Kd_file2'].min())
plt.plot([min_kd, max_kd], [min_kd, max_kd], color='black', linestyle='--')

plt.xscale('log')  # Set x-axis to logarithmic scale
plt.yscale('log')  # Set y-axis to logarithmic scale

#plt.xlabel('Kd from ' + f1_name + '')
#plt.ylabel('Kd from ' + f2_name + '')
plt.title(f'Comparison of Kd values for common variants\nSpearman correlation = {spearman_corr:.2f}')
plt.legend()
plt.grid(True)

plt.show()

def rem_outliers(dataFrame, outliers, filt):
    merged_df = dataFrame.merge(outliers, on=['Variant'], how='outer', indicator=True, suffixes=('', '_outliers'))

    no_outliers = merged_df[merged_df['_merge'] == 'left_only']

    filtered_columns = [col for col in no_outliers.columns if
                        not col.endswith('_file1') and not col.endswith('_file2') and not col.endswith('_merge')]

    # Filter the DataFrame to keep only the selected columns
    df_filtered = no_outliers[filtered_columns]

    filt = [filt, '_remOut']

    return df_filtered

if outlierFilter:
    file1_filtered = rem_outliers(file1_filtered, diff_variants, filt)

    file2_filtered = rem_outliers(file2_filtered, diff_variants, filt)

def common_only(file1, file2):
    Vars = file2[['Variant']]

    merged_df = file1[file1['Variant'].isin(Vars['Variant'])]

    #filt = [filt, '_remOut']

    return merged_df


if commonVars:
    file1_filtered = common_only(file1_filtered, file2_filtered)

    file2_filtered = common_only(file2_filtered, file1_filtered)

file1_filtered.to_csv('YL009_MLE_BA_trueVar-SM_succ_count100.csv')

file2_filtered.to_csv('YL011_MLE_BA_trueVar-SM_succ_count100.csv')
