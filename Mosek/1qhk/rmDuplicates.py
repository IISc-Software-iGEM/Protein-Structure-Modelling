import pandas as pd

# Load the data from a text file
input_file = 'testfile.txt'
output_file = 'output.txt'

# Read the data into a DataFrame
df = pd.read_csv(input_file, delim_whitespace=True, header=None)

# Remove duplicate rows
df = df.drop_duplicates()

# Save the cleaned data back to a text file
df.to_csv(output_file, sep=' ', header=False, index=False)

# Print the cleaned data
print(df)
