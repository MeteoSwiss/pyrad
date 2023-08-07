
# %%
from pathlib import Path
import pandas as pd
mainpath = Path(__file__).resolve().parent.parent
OUT_DIRECTORY = str(Path(mainpath, 'doc', 'mappings'))


def dataframe_to_restructured_text_table(df):
    # Convert the DataFrame to a tabular string representation
 # Determine the maximum width for each column
    max_widths = [max(len(str(value)) for value in df[col])
                  for col in df.columns]

    # Format the DataFrame columns to have the same width
    df_formatted = df.copy()
    for col, max_width in zip(df.columns, max_widths):
        df_formatted[col] = df[col].apply(lambda x: str(x).ljust(max_width))

    # Convert the formatted DataFrame to a tabular string representation
    table_str = df_formatted.to_string(index=False, col_space=max_widths)

    # Split the table string into lines
    lines = table_str.split('\n')

    # Extract column headers and the horizontal separator
    headers = lines[0].split()

    # Create reST table headers
    rst_table = ['| ' + ' | '.join(headers) + ' |']
    # Add rows to the reST table
    for line in lines[1:]:
        rst_table.append('| ' + ' | '.join(line.split()) + ' |')

    # Join the reST table lines
    rst_table_str = '\n'.join(rst_table)

    return rst_table_str


df = pd.read_csv(Path(OUT_DIRECTORY, 'pyrad_to_pyart.txt'))

# Convert DataFrame to reStructuredText table
rst_table = dataframe_to_restructured_text_table(df)

print(rst_table)


# %%
