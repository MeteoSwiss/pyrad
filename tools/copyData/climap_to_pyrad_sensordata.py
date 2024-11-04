import os
import pandas as pd
import argparse

"""
This script processes a single CSV file containing meteorological data extracted from Climap
and organizes the data into a specified directory structure based on date and station identifier
as is required by prad

For each unique date in the data,
it creates a file with all records for that date and station in the format "YYYYMMDD_abbr.csv", where 
"YYYYMMDD" represents the year, month, and day, and "abbr" is the station identifier.

The script takes two command-line arguments: the path to the input CSV file and the output directory. 
The output structure is created under the specified output directory in a nested format by year and 
month (i.e., "/output_dir/YYYYMM/"). Each daily CSV file includes columns "StationID", "DateTime", 
and "Value", with "StationID" representing the station abbreviation.

## CSV Input Format
The input CSV file should have at least the following columns:
    - "abbr": a string representing the station abbreviation.
    - "time": an int64 date and time field in the format "YYYYMMDDHHMM".
    - "rre150z0": precipitation
    as well as any additional columns that correspond to climap variables
    (rre150z0, tre200s0, etc)

## Directory and Filename Structure
The output directory is organized as follows:
    /output_dir/
        └── YYYYMM/
            └── YYYYMMDD_abbr.csv

## Example
Given an input CSV with the following data:
    abbr,time,rre150z0
    BRZ,202408120000,0.0
    BRZ,202408120010,0.1
    ABC,202408120000,0.2

The output structure would be:
    /output_dir/202408/
        ├── 20240812_BRZ.csv
        └── 20240812_ABC.csv

Each CSV file includes the columns:
    - "StationID": the station abbreviation (e.g., "BRZ").
    - "DateTime": the timestamp (e.g., "202408120000").
    - "Value": the measurement (e.g., 0.0).

## Functions
- `process_file(input_csv, output_dir)`: Reads the input CSV file, organizes records by unique date
  and station identifier, and writes each set of records to a new CSV file under the specified output
  directory.
- `main()`: Parses command-line arguments and calls `process_file` to initiate the file processing.

## Usage
Run the script from the command line:
```bash
python process_csv.py input.csv /path/to/output_directory

"""

def process_file(input_csv, output_dir):
    # Read the input CSV file
    df = pd.read_csv(input_csv, dtype={'time': str})
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True,)

    # Iterate over unique dates in the 'time' column
    for date in df['time'].str[:8].unique():
        # Filter rows for the specific date
        daily_df = df[df['time'].str.startswith(date)].copy()
        
        # Get unique abbrv
        station_ids = set(daily_df['abbr'])
        
        # Create the directory path YYYYMM under the output directory
        year_month = date[:6]
        year_month_dir = os.path.join(output_dir, year_month)
        os.makedirs(year_month_dir, exist_ok=True)
        
        for sid in station_ids:
            station_df = daily_df[daily_df['abbr'] == sid]
            
            # Update the column names to match the required output
            station_df.columns = ['StationID', 'DateTime', 'Value']
            
            # Define the output filename YYYYMMDD_abbr.csv
            output_filename = f"{date}_{sid}.csv"
            output_filepath = os.path.join(year_month_dir, output_filename)
            
            # Write the filtered data to the output file
            station_df.to_csv(output_filepath, index=False)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process a CSV file and save rows by date and station.")
    parser.add_argument('input_csv', type=str, help='Path to the input CSV file')
    parser.add_argument('output_dir', type=str, help='Path to the output directory')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process the file
    process_file(args.input_csv, args.output_dir)

if __name__ == '__main__':
    main()
