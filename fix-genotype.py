# Import important libraries
import pandas as pd

# Define the file path
fileInput = 'input.vcf'

# Initialize empty lists to store header lines
header_lines = []

# Open the VCF file and read it line by line
with open(fileInput, 'r') as file:
    for line in file:
        line = line.strip()
        # Check if the line starts with #CHROM, indicating the end of the header
        if line.startswith("#CHROM"):
            break
        header_lines.append(line)

# Removing additional information which is not important
with open(fileInput) as f:
  total = f.readlines()
  for elem in total:
    if elem.find('#CHROM') > -1:
      skip_value = total.index(elem)

# Import vcf file
body_df = pd.read_csv(fileInput, skiprows=skip_value, sep='\t', low_memory=False)

# Assign a name to the column at position 10 (0-based index)
column_position = 9  # This is position 10 (0-based index)
column_variable = body_df.columns[column_position]

# Remove rows where 'column_variable' starts with '.' or NaN or '0/0'
body_df = body_df[~(body_df[column_variable].astype(str).str.startswith('.') | body_df[column_variable].isna())]

# Filter rows where 'column_variable' starts with '0/0' or '0|0' and set 'ALT' equal to 'REF'
body_df.loc[body_df[column_variable].str.startswith(('0/0', '0|0')), 'ALT'] = body_df['REF']

# Separated sample value
body_df['GenTypCh'] = body_df[column_variable].str.split(':').str[0]

# Separated sample value
body_df['GenCvrgCh'] = body_df[column_variable].str.split(':').str[1]

# Split the 'GenTypCh' column into three separate columns
body_df[['GenTypChVal1', 'GenTypChValSep', 'GenTypChVal2']] = body_df['GenTypCh'].str.extract(r'(\d)([/|])(\d)')

# Define a custom function to calculate the GenCvrgChVal1 values
def process_GenCvrgChVal1(row):
    GenCvrgCh_values = row['GenCvrgCh'].split(',')
    GenTypChVal1_value = int(row['GenTypChVal1'])

    if len(GenCvrgCh_values) == 1:
        return row['GenCvrgCh']
    elif GenTypChVal1_value <= len(GenCvrgCh_values):
        return GenCvrgCh_values[GenTypChVal1_value]
    else:
        return ''

body_df['GenCvrgChVal1'] = body_df.apply(process_GenCvrgChVal1, axis=1)

# Define a custom function to calculate the GenCvrgChVal2 values
def process_GenCvrgChVal2(row):
    GenCvrgCh_values = row['GenCvrgCh'].split(',')
    GenTypChVal2_value = int(row['GenTypChVal2'])

    if len(GenCvrgCh_values) == 1:
        return row['GenCvrgCh']
    elif GenTypChVal2_value <= len(GenCvrgCh_values):
        return GenCvrgCh_values[GenTypChVal2_value]
    else:
        return ''

body_df['GenCvrgChVal2'] = body_df.apply(process_GenCvrgChVal2, axis=1)

# Assuming you have a DataFrame called body_df with 'GenCvrgChVal1' and 'GenCvrgChVal2'
body_df['GenCvrgChX'] = body_df['GenCvrgChVal1'] + ',' + body_df['GenCvrgChVal2']

# Replace the second position value in column_variable with GenCvrgChX value separated by ':'
body_df[column_variable] = body_df.apply(lambda row: row[column_variable][:row[column_variable].find(':')+1] + row['GenCvrgChX'] + row[column_variable][row[column_variable].find(':', row[column_variable].find(':')+1):], axis=1)

# Define a custom function to calculate the altVal2 values
def process_altVal2(row):
    alt_values = row['ALT'].split(',')
    GenTypChVal2_value = int(row['GenTypChVal2'])

    if len(alt_values) == 1:
        return row['ALT']
    elif GenTypChVal2_value <= len(alt_values):
        return alt_values[GenTypChVal2_value - 1]
    else:
        return ''

body_df['altVal2'] = body_df.apply(process_altVal2, axis=1)

# Define a custom function to calculate the altVal1 values
def process_altVal1(row):
    if row['GenTypChVal1'] == '0':
        return None
    else:
        alt_values = row['ALT'].split(',')
        GenTypChVal1_value = int(row['GenTypChVal1'])
        if GenTypChVal1_value <= len(alt_values):
            return alt_values[GenTypChVal1_value - 1]
        else:
            return None

# Apply the custom function to create the altVal1 column
body_df['altVal1'] = body_df.apply(process_altVal1, axis=1)

# Use a conditional statement to create the new column
body_df['altValX'] = body_df.apply(lambda row: row['altVal1'] + ',' + row['altVal2'] if row['altVal1'] is not None else row['altVal2'], axis=1)

# Function to apply to each row
def process_altValX(row):
    values = row['altValX'].split(',')

    if len(values) == 1:
        return '0,1'
    elif len(values) == 2:
        if values[0] == values[1]:
            return '1,1'
        else:
            return '1,2'
    else:
        return None  # You can decide what to do with other cases

# Apply the function to create a new column
body_df['GenTypChValX'] = body_df.apply(process_altValX, axis=1)

# Function to remove adjacent duplicate values within a cell
def remove_duplicates(cell_value):
    values = cell_value.split(',')
    result = [values[0]]  # Initialize the result with the first value
    for i in range(1, len(values)):
        if values[i] != values[i - 1]:
            result.append(values[i])
    return ','.join(result)

# Apply the function to the specified column
body_df['altValNew'] = body_df['altValX'].apply(remove_duplicates)

# Function to combine GenTypChValSep and GenTypChValX with the corresponding punctuation
def combine_columns(row):
    GenTypChValSep_value = row['GenTypChValSep']
    GenTypChValX_value = row['GenTypChValX']

    # Split GenTypChValX by comma
    GenTypChValX_parts = GenTypChValX_value.split(',')

    # Join GenTypChValSep and GenTypChValX parts with the same punctuation from GenTypChValSep
    GenTypChValNew_value = GenTypChValSep_value.join(GenTypChValX_parts)

    return GenTypChValNew_value

# Apply the function to create GenTypChValNew
body_df['GenTypChValNew'] = body_df.apply(combine_columns, axis=1)

# Filter rows where 'GenTypCh' is either '0/0' or '0|0' and set 'GenTypChValNew' equal to 'GenTypCh'
body_df.loc[body_df['GenTypCh'].isin(['0/0', '0|0']), 'GenTypChValNew'] = body_df['GenTypCh']

body_df = body_df[body_df['altValNew'] != '*']

# Replace ALT values with altValNew values
body_df['ALT'] = body_df['altValNew']

# Replace the first position value in column_variable with GenTypChValNew value separated by ':'
body_df[column_variable] = body_df.apply(lambda row: row['GenTypChValNew'] + row[column_variable][row[column_variable].find(':'):], axis=1)

# Determine the number of columns to keep (all except the last 10)
columns_to_keep = body_df.shape[1] - 14

# Drop the last 10 columns
body_df = body_df.iloc[:, :columns_to_keep]

# Create an empty dictionary to store the output data
output_data = {}

# Iterate through columns and assign numbers as keys
for i, col_name in enumerate(body_df.columns):
    output_data[str(i)] = [col_name] + body_df[col_name].tolist()

# Create the update_body_df DataFrame
update_body_df = pd.DataFrame(output_data)

# Combine the header lines into a single string
header_string = '\n'.join(header_lines)

# Define the output file path
output_file = 'output.vcf'

# Write the combined header and data to the VCF file
with open(output_file, 'w') as vcf_file:
    vcf_file.write(header_string + '\n')
    update_body_df.to_csv(vcf_file, sep='\t', index=False, header=False, mode='a')
