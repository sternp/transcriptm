# Import necessary libraries
import pandas as pd
import os
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import mode

# Function to find the optimal number of components using BIC and calculate the lowest mode across all components
def gmm_lowest_mode_across_all_components(contig, df, max_components=5):
    # Filter for rows within the specified contig
    contig_df = df[df['contig'] == contig]
    
    # Extract 'rpk' values after ensuring there are values within the contig
    values = contig_df['rpk'].dropna().values.reshape(-1, 1)
    if len(values) < 2:
        return np.nan, 1  # Return NaN if there are not enough values to fit the GMM
    
    # Determine the optimal number of components based on BIC
    lowest_bic = np.inf
    best_gmm = None
    optimal_components = 1
    
    for n_components in range(1, max_components + 1):
        gmm = GaussianMixture(n_components=n_components, random_state=0)
        gmm.fit(values)
        bic = gmm.bic(values)
        if bic < lowest_bic:
            lowest_bic = bic
            best_gmm = gmm
            optimal_components = n_components
    
    # Calculate the mode of the values in each component and identify the lowest mode
    labels = best_gmm.predict(values)
    component_modes = []

    for i in range(optimal_components):
        component_values = values[labels == i].flatten()
        if len(component_values) > 0:
            component_mode_result = mode(component_values)
            component_mode = component_mode_result.mode[0] if isinstance(component_mode_result.mode, np.ndarray) else component_mode_result.mode
            component_modes.append(component_mode)
    
    # Find the lowest mode across all components
    lowest_mode_across_all_components = min(component_modes) if component_modes else np.nan
    
    return lowest_mode_across_all_components, optimal_components


#################################################################### REVERSE


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

# Calculate the GMM-based lowest mode across all components for each contig
gmm_modes_per_contig_new = []
for contig in feature_df_new['contig'].unique():
    lowest_mode, optimal_components = gmm_lowest_mode_across_all_components(contig, intergenic_df_new)
    gmm_modes_per_contig_new.append([contig, lowest_mode, optimal_components])

# Create a DataFrame for the GMM-based lowest mode and optimal component count
gmm_modes_df = pd.DataFrame(gmm_modes_per_contig_new, columns=['contig', 'rpk_intergenic_gmm_mode_lowest', 'optimal_components'])

# Merge the GMM-based lowest mode and component count with the feature count table on 'contig'
merged_df_new = pd.merge(feature_df_new, gmm_modes_df, on='contig', how='left')

# Interpolate missing values in 'rpk_intergenic_gmm_mode_lowest' within each contig group
merged_df_new['rpk_intergenic_gmm_mode_lowest'] = (
    merged_df_new.groupby('contig')['rpk_intergenic_gmm_mode_lowest']
    .apply(lambda group: group.interpolate(method='linear'))
    .reset_index(level=0, drop=True)  # Ensure the index matches the original DataFrame
)

# Assign 0 to any remaining NaN values after interpolation
merged_df_new['rpk_intergenic_gmm_mode_lowest'] = merged_df_new['rpk_intergenic_gmm_mode_lowest'].fillna(0)

# Calculate corrected values using the lowest mode across all components
# Calculate rpk_corrected based on the lowest mode
merged_df_new['rpk_corrected'] = merged_df_new['rpk'] - merged_df_new['rpk_intergenic_gmm_mode_lowest']
merged_df_new['rpk_corrected'] = merged_df_new['rpk_corrected'].apply(lambda x: max(x, 0) if pd.notnull(x) else 0)

# Calculate raw_count_corrected and tpm_corrected based on the lowest mode
merged_df_new['raw_count_corrected'] = merged_df_new['rpk_corrected'] * merged_df_new['length'] / 1000
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



#################################################################### FORWARD


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

# Calculate the GMM-based lowest mode across all components for each contig
gmm_modes_per_contig_new = []
for contig in feature_df_new['contig'].unique():
    lowest_mode, optimal_components = gmm_lowest_mode_across_all_components(contig, intergenic_df_new)
    gmm_modes_per_contig_new.append([contig, lowest_mode, optimal_components])

# Create a DataFrame for the GMM-based lowest mode and optimal component count
gmm_modes_df = pd.DataFrame(gmm_modes_per_contig_new, columns=['contig', 'rpk_intergenic_gmm_mode_lowest', 'optimal_components'])

# Merge the GMM-based lowest mode and component count with the feature count table on 'contig'
merged_df_new = pd.merge(feature_df_new, gmm_modes_df, on='contig', how='left')

# Interpolate missing values in 'rpk_intergenic_gmm_mode_lowest' within each contig group
merged_df_new['rpk_intergenic_gmm_mode_lowest'] = (
    merged_df_new.groupby('contig')['rpk_intergenic_gmm_mode_lowest']
    .apply(lambda group: group.interpolate(method='linear'))
    .reset_index(level=0, drop=True)  # Ensure the index matches the original DataFrame
)

# Assign 0 to any remaining NaN values after interpolation
merged_df_new['rpk_intergenic_gmm_mode_lowest'] = merged_df_new['rpk_intergenic_gmm_mode_lowest'].fillna(0)

# Calculate corrected values using the lowest mode across all components
# Calculate rpk_corrected based on the lowest mode
merged_df_new['rpk_corrected'] = merged_df_new['rpk'] - merged_df_new['rpk_intergenic_gmm_mode_lowest']
merged_df_new['rpk_corrected'] = merged_df_new['rpk_corrected'].apply(lambda x: max(x, 0) if pd.notnull(x) else 0)

# Calculate raw_count_corrected and tpm_corrected based on the lowest mode
merged_df_new['raw_count_corrected'] = merged_df_new['rpk_corrected'] * merged_df_new['length'] / 1000
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





# Finish
os.system("touch final_table_corrected/done")

