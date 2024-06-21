import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from matplotlib.colors import ListedColormap
import seaborn as sns
import numpy as np
import os

plt.clf()

file1 = '/Users/siobhan/PycharmProjects/MAGMA-seq_proccessing/Outputs_Raw MLE/YL009_MLEoutput_variants_bright_adjust.csv'
file2 = '/Users/siobhan/PycharmProjects/MAGMA-seq_proccessing/Outputs_Raw MLE/YL011_MLEoutput_variants_brightAdj.csv'

filt = ''

#f1_name = os.path.basename(file1)
#f2_name = os.path.basename(file2)

file1_df = pd.read_csv(file1)
file2_df = pd.read_csv(file2)

# #Filter settings for Kd rep1 vs Kd rep2
# successFilter = True
# countFilter = True
#count_limit=500
# outlierFilter = False
# commonVars = False
# CIfilter = True
# Kdfilter = True
#Kd_upper_lim=3000 #need to adjust > vs <
# TrueVarFilter = True

#Filter settings for box plots
successFilter = False
countFilter = True
count_limit = 250
outlierFilter = False
commonVars = False
CIfilter = False
Kdfilter = True
Kd_lower_lim = 0  #need to adjust > vs <
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


def is_variant_allowed(variant, antibody, allowed_mutations):
    """Check if the variant contains only allowed VH and VL mutations for the given antibody."""
    vh_mutations, vl_mutations = parse_mutations(variant)
    vh_allowed = all(mutation in allowed_mutations[antibody]['VH'] for mutation in vh_mutations)
    vl_allowed = all(mutation in allowed_mutations[antibody]['VL'] for mutation in vl_mutations)
    return vh_allowed and vl_allowed


def filter_csv(dataFrame, successFilter, countFilter, CIfilter, Kdfilter, TrueVarFilter, filt):
    if successFilter:
        dataFrame = dataFrame[dataFrame['Success'] != False]
        filt = [filt, '_succ']

    if countFilter:
        dataFrame = dataFrame[dataFrame['Avg_counts'] >= count_limit]
        filt = [filt, '_ct']

    if CIfilter:
        dataFrame = dataFrame[dataFrame['95% CI high'] > 0]

    if Kdfilter:
        dataFrame = dataFrame[dataFrame['Kd'] > Kd_lower_lim]

    if TrueVarFilter:
        def variant_allowed(row):
            antibody = row['Ab']
            variant = row['Variant']
            if antibody in antibodies_with_allowed_mutations:
                return is_variant_allowed(variant, antibody, antibodies_with_allowed_mutations)
            return False

        dataFrame = dataFrame[dataFrame.apply(variant_allowed, axis=1)]

    filt = [filt, '_trueVar']

    return dataFrame


file1_filtered = filter_csv(file1_df, successFilter, countFilter, CIfilter, Kdfilter, TrueVarFilter, filt)

file2_filtered = filter_csv(file2_df, successFilter, countFilter, CIfilter, Kdfilter, TrueVarFilter, filt)

# Merge filtered DataFrames
common_variants = pd.merge(file1_filtered, file2_filtered, on='Variant', suffixes=('_YL009', '_YL011'))

# Calculate the Spearman rank correlation
spearman_corr, _ = spearmanr(common_variants['Kd_YL009'], common_variants['Kd_YL011'])

# Define unique categories
categories = common_variants['Ab_YL009'].unique()

# Define a colormap
# Use a seaborn palette
# Define a list of 7 unique hexadecimal color codes
hex_colors = ['#FFC305', '#838383', '#45ACB3', '#F27F34', '#3954B2', '#E24249', '#A451C2']

# Create a custom colormap
custom_cmap = ListedColormap(hex_colors)

# Plot each category with its assigned color from the custom colormap
plt.figure(figsize=(10, 6))
for i, category in enumerate(categories):
    subset = common_variants[common_variants['Ab_YL009'] == category]
    plt.scatter(subset['Kd_YL009'], subset['Kd_YL011'],
                label=category,
                color=custom_cmap(i))

# Add a line y = x
max_kd = max(common_variants['Kd_YL009'].max(), common_variants['Kd_YL011'].max())
min_kd = min(common_variants['Kd_YL009'].min(), common_variants['Kd_YL011'].min())
plt.plot([min_kd, max_kd], [min_kd, max_kd], color='black', linestyle='--')

plt.xscale('log')  # Set x-axis to logarithmic scale
plt.yscale('log')  # Set y-axis to logarithmic scale

# Define font properties
font_properties = {'family': 'sans-serif', 'weight': 'normal', 'size': 14}

# Add labels, title, and legend with custom font
plt.xlabel('KD - Rep 1 [nM]', fontdict=font_properties)
plt.xticks(fontproperties=font_properties)
plt.yticks(fontproperties=font_properties)
plt.ylabel('Kd - Rep 2 [nM]', fontdict=font_properties)
plt.title(f'Comparison of Kd values for common variants\nSpearman correlation = {spearman_corr:.2f}')
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.legend(framealpha=1, prop=font_properties, frameon=False)

# Save the plot as a PNG file
plt.savefig('Fig2B_KdComp.png', dpi=600, bbox_inches='tight')

plt.show()

final_commonVariants = common_variants[['Variant', 'Kd_YL009', 'Kd_YL011']]

final_commonVariants.to_csv('Kd_repComparison_data.csv', index=False)


def parse_and_count_valid_mutations(variant_str, exclusion_list):
    """Extract mutations from the variant string and count valid mutations."""
    vh_mutations = []
    vl_mutations = []
    sections = variant_str.split('|')

    for section in sections:
        if '>VH:' in section:
            vh_part = section.split('>VH:')[1]
            for mut in vh_part.split(';'):
                mutation = mut.split('-')[0]
                # Check if mutation is not 'WT', not in exclusion list, and does not have same first and last letter
                if mutation != 'WT' and mutation not in exclusion_list and len(mutation) >= 2 and mutation[0] != \
                        mutation[-1]:
                    vh_mutations.append(mutation)
        elif '>VL:' in section:
            vl_part = section.split('>VL:')[1]
            for mut in vl_part.split(';'):
                mutation = mut.split('-')[0]
                # Check if mutation is not 'WT', not in exclusion list, and does not have same first and last letter
                if mutation != 'WT' and mutation not in exclusion_list and len(mutation) >= 2 and mutation[0] != \
                        mutation[-1]:
                    vl_mutations.append(mutation)

    return len(vh_mutations) + len(vl_mutations)


# Example usage:
mutations_to_remove = {"Q1G", "E2Q"}

# Apply the function to each row in file1_filtered
file1_filtered['mutation_count'] = file1_filtered['Variant'].apply(
    lambda x: parse_and_count_valid_mutations(x, mutations_to_remove))

file2_filtered['mutation_count'] = file2_filtered['Variant'].apply(lambda x: parse_and_count_valid_mutations(x, mutations_to_remove))

#Violin Plots

file1_filtered = file1_filtered.assign(Source='YL009')
file2_filtered = file2_filtered.assign(Source='YL011')

df_combined = pd.concat([file1_filtered, file2_filtered])

df_combined.loc[df_combined['Kd'] > 2000, 'Kd'] = 2500

#Plot 002-S21F2
Ab002_df = df_combined[df_combined['Ab'] == '002-S21F2_UCA']

# Step 1: Get unique mutation_count categories from Ab002_df
mutation_counts = sorted(Ab002_df['mutation_count'].unique())

# Step 2: Define the height and width ratio
subplot_height = 1  # Adjust the height of each subplot
subplot_width = 12  # Calculate the width based on the aspect ratio 1:3

# Step 3: Create subplots with the specified size and layout
fig, axes = plt.subplots(nrows=len(mutation_counts), ncols=1,
                         figsize=(subplot_width, subplot_height * len(mutation_counts)), sharex=True)

# Define logarithmically spaced bin edges

# Define logarithmically spaced bin edges

# Define subplot layout
subplot_height = 2  # Adjust as needed
subplot_width = 12  # Adjust as needed
fig, axes = plt.subplots(nrows=len(mutation_counts), ncols=1,
                         figsize=(subplot_width, subplot_height * len(mutation_counts)), sharex=True)


# Iterate over each mutation_count category and plot the corresponding histogram
for i, count in enumerate(mutation_counts):
    ax = axes[i]
    test=Ab002_df[Ab002_df['mutation_count'] == count]

    ax.set_xscale('log')

    sns.histplot(data=Ab002_df[Ab002_df['mutation_count'] == count], color='#FFC305', edgecolor='#A88103', x='Kd', ax=ax, bins=[1,1.07931034482759,1.15862068965517,1.23793103448276,1.31724137931035,1.39655172413793,1.47586206896552,1.5551724137931,1.63448275862069,1.71379310344828,1.79310344827586,1.87241379310345,1.95172413793104,2.03103448275862,2.11034482758621,2.18965517241379,2.26896551724138,2.34827586206897,2.42758620689655,2.50689655172414,2.58620689655172,2.66551724137931,2.7448275862069,2.82413793103448,2.90344827586207,2.98275862068965,3.06206896551724,3.14137931034483,3.22068965517241,3.3, 3.379310345,3.45862069])

    # Set x-axis to log scale


    # Remove right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_ylabel(count, rotation=0, ha='right', fontsize=30, labelpad=20)

    # Set font size for x-axis label
    if i == len(mutation_counts) - 1:
        ax.set_xlabel('Kd', fontsize=30)
    else:
        ax.set_xlabel('')
        ax.xaxis.set_tick_params(labelbottom=False)

    # Set font size and padding for tick labels
    ax.tick_params(axis='both', which='major', labelsize=25, pad=10)

    # Ensure tight layout for each subplot
    ax.figure.tight_layout()

# Adjust overall layout of the figure
plt.subplots_adjust(hspace=0.2)  # Adjust vertical spacing between subplots if needed

#plt.savefig('002_mutationalSteps.png', dpi=1200, bbox_inches='tight')
# Show the plot
plt.show()

#Plot 2-15
Ab215_df = df_combined[df_combined['Ab'] == 'Ab_2-15_UCA']

# Step 1: Get unique mutation_count categories from 2-15_df
mutation_counts = sorted(Ab215_df['mutation_count'].unique())

# Define logarithmically spaced bin edges

# Define logarithmically spaced bin edges

# Define subplot layout
subplot_height = 2  # Adjust as needed
subplot_width = 12  # Adjust as needed
fig, axes = plt.subplots(nrows=len(mutation_counts), ncols=1,
                         figsize=(subplot_width, subplot_height * len(mutation_counts)), sharex=True)


# Iterate over each mutation_count category and plot the corresponding histogram
for i, count in enumerate(mutation_counts):
    ax = axes[i]
    test=Ab215_df[Ab215_df['mutation_count'] == count]

    ax.set_xscale('log')

    sns.histplot(data=Ab215_df[Ab215_df['mutation_count'] == count], color='#F27F34', edgecolor='#9A5120', x='Kd', ax=ax, bins=[1,1.07931034482759,1.15862068965517,1.23793103448276,1.31724137931035,1.39655172413793,1.47586206896552,1.5551724137931,1.63448275862069,1.71379310344828,1.79310344827586,1.87241379310345,1.95172413793104,2.03103448275862,2.11034482758621,2.18965517241379,2.26896551724138,2.34827586206897,2.42758620689655,2.50689655172414,2.58620689655172,2.66551724137931,2.7448275862069,2.82413793103448,2.90344827586207,2.98275862068965,3.06206896551724,3.14137931034483,3.22068965517241,3.3, 3.379310345,3.45862069])

    # Set x-axis to log scale


    # Remove right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_ylabel(count, rotation=0, ha='right', fontsize=30, labelpad=20)

    # Set font size for x-axis label
    if i == len(mutation_counts) - 1:
        ax.set_xlabel('Kd', fontsize=30)
    else:
        ax.set_xlabel('')
        ax.xaxis.set_tick_params(labelbottom=False)

    # Set font size and padding for tick labels
    ax.tick_params(axis='both', which='major', labelsize=25, pad=10)

    # Ensure tight layout for each subplot
    ax.figure.tight_layout()

# Adjust overall layout of the figure
plt.subplots_adjust(hspace=0.2)  # Adjust vertical spacing between subplots if needed

#plt.savefig('215_mutationalSteps.png', dpi=1200, bbox_inches='tight')
# Show the plot
plt.show()



