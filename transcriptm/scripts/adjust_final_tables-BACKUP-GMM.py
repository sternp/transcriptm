# Import necessary libraries
import pandas as pd
import os
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import mode

# Function to calculate the mode of the lowest Gaussian component across the entire contig
def gmm_mode_lowest_component(contig, df, n_components=2):
    # Filter for rows within the specified contig
    contig_df = df[df['contig'] == contig]
    
    # Extract 'rpk' values after ensuring there are values within the contig
    values = contig_df['rpk'].dropna().values.reshape(-1, 1)
    if len(values) < n_components:
        return np.nan  # Return NaN if there are not enough values to fit the GMM
    
    # Fit a Gaussian Mixture Model with multiple components
    gmm = GaussianMixture(n_components=n_components, random_state=0)
    gmm.fit(values)
    
    # Identify the component with the lowest mean
    component_means = gmm.means_.flatten()
    lowest_component_index = np.argmin(component_means)
    
    # Get the data points most likely to belong to the lowest component
    labels = gmm.predict(values)
    lowest_component_values = values[labels == lowest_component_index].flatten()
    
    # Calculate the mode of the values in the lowest component, handling scalar output
    if len(lowest_component_values) > 0:
        mode_result = mode(lowest_component_values)
        lowest_mode = mode_result.mode[0] if isinstance(mode_result.mode, np.ndarray) else mode_result.mode
    else:
        lowest_mode = np.nan
    
    return lowest_mode

########################################################################################################################
#REV


# Load intergenic_rawcount_tpm_table.txt with genomic coordinates
intergenic_rawcount_tpm_table_path = 'final_table/intergenic_rawcount_tpm_table.txt'
intergenic_df_new = pd.read_csv(intergenic_rawcount_tpm_table_path, sep='\t')

# Load rawcount_tpm_table.txt
rawcount_tpm_table_path = 'final_table/rawcount_tpm_table.txt'
feature_df_new = pd.read_csv(rawcount_tpm_table_path, sep='\t')

# Convert relevant columns to numeric, forcing errors to NaN
for col in ['rpk', 'raw_count', 'length', 'start', 'end']:
    intergenic_df_new[col] = pd.to_numeric(intergenic_df_new[col], errors='coerce')
    feature_df_new[col] = pd.to_numeric(feature_df_new[col], errors='coerce')

# Calculate the GMM-based mode for the lowest component across each contig
gmm_mode_per_contig_new = []
for contig in feature_df_new['contig'].unique():
    gmm_mode = gmm_mode_lowest_component(contig, intergenic_df_new)
    gmm_mode_per_contig_new.append([contig, gmm_mode])

# Create a DataFrame for the GMM-based mode values
gmm_mode_df = pd.DataFrame(gmm_mode_per_contig_new, columns=['contig', 'rpk_intergenic_gmm_mode'])

# Merge the GMM-based mode values with the feature count table on 'contig'
merged_df_new = pd.merge(feature_df_new, gmm_mode_df, on='contig', how='left')

# Interpolate missing values in 'rpk_intergenic_gmm_mode' within each contig group
merged_df_new['rpk_intergenic_gmm_mode'] = (
    merged_df_new.groupby('contig')['rpk_intergenic_gmm_mode']
    .apply(lambda group: group.interpolate(method='linear'))
    .reset_index(level=0, drop=True)  # Ensure the index matches the original DataFrame
)

# Assign 0 to any remaining NaN values in 'rpk_intergenic_gmm_mode' after interpolation
merged_df_new['rpk_intergenic_gmm_mode'] = merged_df_new['rpk_intergenic_gmm_mode'].fillna(0)

# Add column to indicate whether each value was interpolated
merged_df_new['was_interpolated'] = merged_df_new['rpk_intergenic_gmm_mode'].isna() | (
    merged_df_new['rpk_intergenic_gmm_mode'] != merged_df_new['rpk_intergenic_gmm_mode'].ffill()
)
merged_df_new['was_interpolated'] = merged_df_new['was_interpolated'].astype(int)  # 1 for interpolated, 0 otherwise

# Subtract the GMM-based mode from each RPK value to get rpk_corrected, handling NaN values
merged_df_new['rpk_corrected'] = merged_df_new['rpk'] - merged_df_new['rpk_intergenic_gmm_mode']
merged_df_new['rpk_corrected'] = merged_df_new['rpk_corrected'].apply(lambda x: max(x, 0) if pd.notnull(x) else 0)

# Calculate raw_count_corrected based on rpk_corrected and length
merged_df_new['raw_count_corrected'] = merged_df_new['rpk_corrected'] * merged_df_new['length'] / 1000

# Calculate corrected TPM values
merged_df_new['length_kb'] = merged_df_new['length'] / 1000
merged_df_new['raw_count_corrected_per_kb'] = merged_df_new['raw_count_corrected'] / merged_df_new['length_kb']
sum_normalized_counts = merged_df_new['raw_count_corrected_per_kb'].sum()
merged_df_new['tpm_corrected'] = (merged_df_new['raw_count_corrected_per_kb'] * 1e6) / sum_normalized_counts
merged_df_new.drop(columns=['length_kb', 'raw_count_corrected_per_kb'], inplace=True)

# Save the final table without 'raw_count_adjusted' column
output_directory = 'final_table_corrected'
os.makedirs(output_directory, exist_ok=True)
output_file_path = os.path.join(output_directory, 'corrected_feature_count_table.tsv')
merged_df_new.to_csv(output_file_path, index=False, sep='\t')




########################################################################################################################
#FWD

# Load intergenic_rawcount_tpm_table.txt with genomic coordinates
intergenic_rawcount_tpm_table_path = 'final_table/intergenic_rawcount_fwd_tpm_table.txt'
intergenic_df_new = pd.read_csv(intergenic_rawcount_tpm_table_path, sep='\t')

# Load rawcount_tpm_table.txt
rawcount_tpm_table_path = 'final_table/rawcount_tpm_table.txt'
feature_df_new = pd.read_csv(rawcount_tpm_table_path, sep='\t')

# Convert relevant columns to numeric, forcing errors to NaN
for col in ['rpk', 'raw_count', 'length', 'start', 'end']:
    intergenic_df_new[col] = pd.to_numeric(intergenic_df_new[col], errors='coerce')
    feature_df_new[col] = pd.to_numeric(feature_df_new[col], errors='coerce')

# Calculate the GMM-based mode for the lowest component across each contig
gmm_mode_per_contig_new = []
for contig in feature_df_new['contig'].unique():
    gmm_mode = gmm_mode_lowest_component(contig, intergenic_df_new)
    gmm_mode_per_contig_new.append([contig, gmm_mode])

# Create a DataFrame for the GMM-based mode values
gmm_mode_df = pd.DataFrame(gmm_mode_per_contig_new, columns=['contig', 'rpk_intergenic_gmm_mode'])

# Merge the GMM-based mode values with the feature count table on 'contig'
merged_df_new = pd.merge(feature_df_new, gmm_mode_df, on='contig', how='left')

# Interpolate missing values in 'rpk_intergenic_gmm_mode' within each contig group
merged_df_new['rpk_intergenic_gmm_mode'] = (
    merged_df_new.groupby('contig')['rpk_intergenic_gmm_mode']
    .apply(lambda group: group.interpolate(method='linear'))
    .reset_index(level=0, drop=True)  # Ensure the index matches the original DataFrame
)

# Assign 0 to any remaining NaN values in 'rpk_intergenic_gmm_mode' after interpolation
merged_df_new['rpk_intergenic_gmm_mode'] = merged_df_new['rpk_intergenic_gmm_mode'].fillna(0)

# Add column to indicate whether each value was interpolated
merged_df_new['was_interpolated'] = merged_df_new['rpk_intergenic_gmm_mode'].isna() | (
    merged_df_new['rpk_intergenic_gmm_mode'] != merged_df_new['rpk_intergenic_gmm_mode'].ffill()
)
merged_df_new['was_interpolated'] = merged_df_new['was_interpolated'].astype(int)  # 1 for interpolated, 0 otherwise

# Subtract the GMM-based mode from each RPK value to get rpk_corrected, handling NaN values
merged_df_new['rpk_corrected'] = merged_df_new['rpk'] - merged_df_new['rpk_intergenic_gmm_mode']
merged_df_new['rpk_corrected'] = merged_df_new['rpk_corrected'].apply(lambda x: max(x, 0) if pd.notnull(x) else 0)

# Calculate raw_count_corrected based on rpk_corrected and length
merged_df_new['raw_count_corrected'] = merged_df_new['rpk_corrected'] * merged_df_new['length'] / 1000

# Calculate corrected TPM values
merged_df_new['length_kb'] = merged_df_new['length'] / 1000
merged_df_new['raw_count_corrected_per_kb'] = merged_df_new['raw_count_corrected'] / merged_df_new['length_kb']
sum_normalized_counts = merged_df_new['raw_count_corrected_per_kb'].sum()
merged_df_new['tpm_corrected'] = (merged_df_new['raw_count_corrected_per_kb'] * 1e6) / sum_normalized_counts
merged_df_new.drop(columns=['length_kb', 'raw_count_corrected_per_kb'], inplace=True)

# Save the final table without 'raw_count_adjusted' column
output_directory = 'final_table_corrected'
os.makedirs(output_directory, exist_ok=True)
output_file_path = os.path.join(output_directory, 'FWD-corrected_feature_count_table.tsv')
merged_df_new.to_csv(output_file_path, index=False, sep='\t')






########################################################################################################################







# Finish
os.system("touch final_table_corrected/done")

