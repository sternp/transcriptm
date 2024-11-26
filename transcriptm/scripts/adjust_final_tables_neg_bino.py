# Importing necessary libraries
import pandas as pd
import os
import numpy as np
from scipy.stats import trim_mean, skew
from statsmodels.formula.api import glm
import statsmodels.api as sm

############################################################################

def skewness_based_trim_mean_top(values):
    data_skewness = skew(values)
    if data_skewness < 1:
        trim_pct = 0.10  # Mild skewness, trim 10%
    elif 1 <= data_skewness < 10:
        trim_pct = 0.20  # Moderate skewness, trim 20%
    else:
        trim_pct = 0.40  # High skewness, trim 40%

    sorted_values = np.sort(values)
    num_to_keep = int(len(sorted_values) * (1 - trim_pct))
    trimmed_values = sorted_values[:num_to_keep]
    return np.mean(trimmed_values)

###########################################################################

# Function to perform negative binomial regression and apply skewness-based trimmed mean per contig
def calculate_nb_trim_mean(raw_counts, contig):
    df = pd.DataFrame({'raw_count': raw_counts, 'contig': contig})
    df['raw_count'] = pd.to_numeric(df['raw_count'], errors='coerce')
    df = df.dropna(subset=['raw_count'])

    # Set alpha to 1.0 for the negative binomial regression
    alpha = 1.0
    model = glm('raw_count ~ 1', data=df, family=sm.families.NegativeBinomial(alpha=alpha)).fit()

    df['predicted'] = model.predict()

    # Group by 'contig' and calculate the skewness-based trimmed mean for predicted values
    trim_mean_nb_per_contig = df.groupby('contig')['predicted'].apply(skewness_based_trim_mean_top).reset_index()
    trim_mean_nb_per_contig.rename(columns={'predicted': 'nb_trimmed_mean'}, inplace=True)

    return trim_mean_nb_per_contig

###########################################################################

# Load intergenic_rawcount_tpm_table.txt
intergenic_rawcount_tpm_table_path = 'final_table/intergenic_rawcount_fwd_tpm_table.txt'
with open(intergenic_rawcount_tpm_table_path, 'r') as file:
    intergenic_rawcount_tpm_table_content = file.readlines()

# Load rawcount_tpm_table.txt
rawcount_tpm_table_path = 'final_table/rawcount_fwd_tpm_table.txt'
with open(rawcount_tpm_table_path, 'r') as file:
    rawcount_tpm_table_content = file.readlines()

# Define the columns for both tables
intergenic_columns_new = ["gene_id", "contig", "start", "end", "strand", "length", "raw_count", "rpk", "tpm"]
feature_columns_new = ["gene_id", "contig", "start", "end", "strand", "length", "gene", "product", "db_xref", "Dbxref", "inference", "UniProtKB", "raw_count", "rpk", "tpm"]

# Parsing intergenic_rawcount_tpm_table.txt
intergenic_data_new = [line.strip().split('\t') for line in intergenic_rawcount_tpm_table_content[1:]]
intergenic_df_new = pd.DataFrame(intergenic_data_new, columns=intergenic_columns_new)

# Parsing rawcount_tpm_table.txt
feature_data_new = [line.strip().split('\t') for line in rawcount_tpm_table_content[1:]]
feature_df_new = pd.DataFrame(feature_data_new, columns=feature_columns_new)

# Convert relevant columns to numeric for calculations
intergenic_df_new['rpk'] = pd.to_numeric(intergenic_df_new['rpk'])
feature_df_new['rpk'] = pd.to_numeric(feature_df_new['rpk'])
feature_df_new['length'] = pd.to_numeric(feature_df_new['length'])

# Group by 'contig' and calculate the trimmed mean for intergenic counts
trimmed_mean_rpk_per_contig_new = intergenic_df_new.groupby('contig')['rpk'].apply(skewness_based_trim_mean_top).reset_index()

# Merge the trimmed mean RPK values with the feature count table on the 'contig' column
merged_df_new = pd.merge(feature_df_new, trimmed_mean_rpk_per_contig_new, on='contig', how='left', suffixes=('', '_intergenic_trimmed_mean'))

# Subtract the trimmed mean RPK value from each RPK value in the feature count table to get rpk_corrected
merged_df_new['rpk_corrected'] = merged_df_new['rpk'] - merged_df_new['rpk_intergenic_trimmed_mean']
merged_df_new['rpk_corrected'] = merged_df_new['rpk_corrected'].apply(lambda x: max(x, 0))

# Calculate the raw_count_corrected based on rpk_corrected and length
merged_df_new['raw_count_corrected'] = merged_df_new['rpk_corrected'] * merged_df_new['length'] / 1000

# Correct TPM calculation approach
merged_df_new['length_kb'] = merged_df_new['length'] / 1000
merged_df_new['raw_count_corrected_per_kb'] = merged_df_new['raw_count_corrected'] / merged_df_new['length_kb']
sum_normalized_counts = merged_df_new['raw_count_corrected_per_kb'].sum()
merged_df_new['tpm_corrected'] = (merged_df_new['raw_count_corrected_per_kb'] * 1e6) / sum_normalized_counts
merged_df_new.drop(columns=['length_kb', 'raw_count_corrected_per_kb'], inplace=True)

# Apply negative binomial regression and add skewness-based trimmed mean NB values per contig
nb_trim_mean_per_contig = calculate_nb_trim_mean(feature_df_new['raw_count'], feature_df_new['contig'])
merged_df_new = pd.merge(merged_df_new, nb_trim_mean_per_contig, on='contig', how='left')

# Ensure that raw_count and nb_trimmed_mean in merged_df_new are numeric before subtraction
merged_df_new['raw_count'] = pd.to_numeric(merged_df_new['raw_count'], errors='coerce')
merged_df_new['nb_trimmed_mean'] = pd.to_numeric(merged_df_new['nb_trimmed_mean'], errors='coerce')

# Create raw_count_adjusted by subtracting nb_trimmed_mean from raw_count
merged_df_new['raw_count_adjusted'] = merged_df_new['raw_count'] - merged_df_new['nb_trimmed_mean']

# Replace any NaN values in raw_count_adjusted with 0
merged_df_new['raw_count_adjusted'] = merged_df_new['raw_count_adjusted'].fillna(0)

# Save the final table with the new columns
output_directory = 'final_table_corrected'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

output_file_path = os.path.join(output_directory, 'corrected_feature_count_table_with_nb_trim_mean.tsv')
merged_df_new.to_csv(output_file_path, index=False, sep='\t')

print(f"\nThe updated table, including NB regression trimmed mean and adjusted raw count values, is saved at: {output_file_path}")
