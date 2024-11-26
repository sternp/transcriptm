# Import necessary libraries
import pandas as pd
import os
import numpy as np

# Function to calculate the median within 50,000 bases
def median_within_range(query_start, query_end, df, contig, base_range=50000):
    # Filter for rows within the specified range around the query
    nearby_df = df[(df['contig'] == contig) & 
                   (df['start'] <= query_end + base_range) & 
                   (df['end'] >= query_start - base_range)]
    
    # Extract 'rpk' values after ensuring there are values within range
    values = nearby_df['rpk'].dropna().values
    if len(values) == 0:
        return np.nan  # Return NaN if there are no values within range
    
    # Calculate the median
    return np.median(values)

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

# Calculate the median RPK for intergenic counts within 50,000 bases by contig
median_rpk_per_contig_new = []
for _, row in feature_df_new.iterrows():
    contig = row['contig']
    query_start = row['start']
    query_end = row['end']
    median_rpk = median_within_range(query_start, query_end, intergenic_df_new, contig)
    median_rpk_per_contig_new.append([contig, query_start, query_end, median_rpk])

# Create a DataFrame for the median values
median_rpk_df = pd.DataFrame(median_rpk_per_contig_new, columns=['contig', 'start', 'end', 'rpk_intergenic_median'])

# Merge the median RPK values with the feature count table on 'contig', 'start', and 'end'
merged_df_new = pd.merge(feature_df_new, median_rpk_df, on=['contig', 'start', 'end'], how='left')

# Interpolate missing values in 'rpk_intergenic_median' within each contig group
merged_df_new['rpk_intergenic_median'] = (
    merged_df_new.groupby('contig')['rpk_intergenic_median']
    .apply(lambda group: group.interpolate(method='linear'))
    .reset_index(level=0, drop=True)  # Ensure the index matches the original DataFrame
)

# Assign 0 to any remaining NaN values in 'rpk_intergenic_median' after interpolation
merged_df_new['rpk_intergenic_median'] = merged_df_new['rpk_intergenic_median'].fillna(0)

# Add column to indicate whether each value was interpolated
merged_df_new['was_interpolated'] = merged_df_new['rpk_intergenic_median'].isna() | (
    merged_df_new['rpk_intergenic_median'] != merged_df_new['rpk_intergenic_median'].ffill()
)
merged_df_new['was_interpolated'] = merged_df_new['was_interpolated'].astype(int)  # 1 for interpolated, 0 otherwise

# Subtract median RPK from each RPK value to get rpk_corrected, handling NaN values
merged_df_new['rpk_corrected'] = merged_df_new['rpk'] - merged_df_new['rpk_intergenic_median']
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

# Finish
os.system("touch final_table_corrected/done")

