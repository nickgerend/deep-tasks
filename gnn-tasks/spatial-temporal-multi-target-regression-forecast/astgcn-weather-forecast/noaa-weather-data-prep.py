#%% NOAA weather data: https://www.ncei.noaa.gov/pub/data/ghcn/daily/

# Downoload Files:
# - readme.txt: contains the descriptions of the data and instructions
# - ghcnd-stations.txt: contains the metadata for the stations
# - ghcnd_all.tar.gz: contains the daily data for all stations

# Extract the daily data:

import tarfile
from tqdm import tqdm
import os
file_path = '.\\noaa\\ghcnd_all.tar.gz'
try:
    with tarfile.open(file_path, 'r:gz') as tar:
        members = tar.getmembers()  # Get a list of all members in the archive
        total_files = len(members)  # Total number of files to extract

        # Use tqdm for a progress bar
        with tqdm(total=total_files, desc="Extracting files", unit="file") as progress_bar:
            for member in members:
                try:
                    tar.extract(member, path='extracted_files')  # Extract the file
                except Exception as e:
                    print(f"Error extracting {member.name}: {e}")
                progress_bar.update(1)  # Update the progress bar
    print('Extraction completed successfully.')
except Exception as e:
    print(f"Failed to open the tar file: {e}")

#%% Format the data

# per the instructions in the readme.txt file, the data is in fixed-width format
def read_dly_file(file_path):

    # Define column specifications based on the NOAA .dly format
    colspecs = [(0, 11), (11, 15), (15, 17), (17, 21)]  # ID, YEAR, MONTH, ELEMENT
    colspecs += [(21 + i * 8, 26 + i * 8) for i in range(31)]  # VALUE1 to VALUE31
    colspecs += [(27 + i * 8, 27 + i * 8) for i in range(31)]  # MFLAG1 to MFLAG31
    colspecs += [(28 + i * 8, 28 + i * 8) for i in range(31)]  # QFLAG1 to QFLAG31
    colspecs += [(29 + i * 8, 29 + i * 8) for i in range(31)]  # SFLAG1 to SFLAG31

    # Define column names
    columns = ["ID", "YEAR", "MONTH", "ELEMENT"]
    columns += [f"VALUE{i + 1}" for i in range(31)]
    columns += [f"MFLAG{i + 1}" for i in range(31)]
    columns += [f"QFLAG{i + 1}" for i in range(31)]
    columns += [f"SFLAG{i + 1}" for i in range(31)]

    # Read the .dly file into a pandas DataFrame
    df = pd.read_fwf(file_path, colspecs=colspecs, names=columns)

    # Replace missing values (-9999) with NaN
    for i in range(31):
        df[f"VALUE{i + 1}"] = df[f"VALUE{i + 1}"].replace(-9999, pd.NA)

    # Reshape the DataFrame to long format
    long_format = df.melt(
        id_vars=["ID", "YEAR", "MONTH", "ELEMENT"],
        value_vars=[f"VALUE{i + 1}" for i in range(31)],
        var_name="DAY",
        value_name="VALUE",
    )

    # Extract the day number from the column name (e.g., VALUE1 → 1)
    long_format["DAY"] = long_format["DAY"].str.extract(r'(\d+)').astype(int)

    # Convert to a datetime format
    long_format["DATE"] = pd.to_datetime(
        long_format[["YEAR", "MONTH", "DAY"]], errors="coerce"
    )

    # Drop rows with invalid dates
    long_format = long_format.dropna(subset=["DATE"])

    # Merge the flags (MFLAG, QFLAG, SFLAG) for each day
    long_format["MFLAG"] = df.melt(
        id_vars=["ID", "YEAR", "MONTH", "ELEMENT"],
        value_vars=[f"MFLAG{i + 1}" for i in range(31)],
        var_name="DAY",
        value_name="MFLAG",
    )["MFLAG"]

    long_format["QFLAG"] = df.melt(
        id_vars=["ID", "YEAR", "MONTH", "ELEMENT"],
        value_vars=[f"QFLAG{i + 1}" for i in range(31)],
        var_name="DAY",
        value_name="QFLAG",
    )["QFLAG"]

    long_format["SFLAG"] = df.melt(
        id_vars=["ID", "YEAR", "MONTH", "ELEMENT"],
        value_vars=[f"SFLAG{i + 1}" for i in range(31)],
        var_name="DAY",
        value_name="SFLAG",
    )["SFLAG"]

    # Drop unnecessary columns and reorder
    long_format = long_format[["ID", "DATE", "ELEMENT", "VALUE", "MFLAG", "QFLAG", "SFLAG"]]

    return long_format

# process the dly files and add some filters, save to csv
def process_dly_files(directory_path, output_csv, file_list):
   
   # process specific .dly files in a directory and save a filtered CSV file

    consolidated_data = []

    # Iterate through the specified files in the file_list
    for file_name in file_list:
        file_path = os.path.join(directory_path, file_name)
        if os.path.exists(file_path) and file_name.endswith(".dly"):
            print(f"Processing: {file_name}")

            # Read and process the .dly file
            df = read_dly_file(file_path)

            # Filter for files containing the core elements
            core_elements = {"PRCP", "TMAX", "TMIN"}
            available_elements = set(df["ELEMENT"].unique())
            if not core_elements.issubset(available_elements):
                print(f"Skipping {file_name}: Missing core elements.")
                continue

            # Filter for dates in 2025 or later
            df_filtered = df[df["DATE"].dt.year >= 2020]

            # Add to the consolidated data if there are valid rows
            if not df_filtered.empty:
                consolidated_data.append(df_filtered)
        else:
            print(f"Skipping {file_name}: File not found or invalid format.")

    # Combine all filtered data
    if consolidated_data:
        final_df = pd.concat(consolidated_data, ignore_index=True)
        final_df.to_csv(output_csv, index=False)
        print(f"Consolidated CSV saved to: {output_csv}")
    else:
        print("No valid data found for the filter years with all filter elements.")

# load the stations file
def load_stations(file_path):
    
    # Loads the stations file (ghcnd-stations.txt) and returns a DataFrame.

    # Define the column specifications and names based on the file format
    colspecs = [
        (0, 11),   # ID
        (12, 20),  # LATITUDE
        (21, 30),  # LONGITUDE
        (31, 37),  # ELEVATION
        (38, 40),  # STATE
        (41, 71),  # NAME
        (72, 75),  # GSN FLAG
        (76, 79),  # HCN/CRN FLAG
        (80, 85),  # WMO ID
    ]
    column_names = [
        "ID", "LATITUDE", "LONGITUDE", "ELEVATION", "STATE",
        "NAME", "GSN_FLAG", "HCN_CRN_FLAG", "WMO_ID"
    ]

    # Read the file into a DataFrame
    df = pd.read_fwf(file_path, colspecs=colspecs, names=column_names)

    # Replace missing values for ELEVATION and WMO ID
    df["ELEVATION"] = df["ELEVATION"].replace(-999.9, pd.NA)
    df["WMO_ID"] = df["WMO_ID"].replace("     ", pd.NA)  # Blank WMO IDs to NaN

    return df

# df date and element filters
def filter_data(df, date_range, required_elements):

    # Filter data within the specified date range
    start_date, end_date = pd.to_datetime(date_range[0]), pd.to_datetime(date_range[1])
    df_filtered = df[(df["DATE"] >= start_date) & (df["DATE"] <= end_date)]

    # Check IDs that are present for all days in the date range
    date_range_days = pd.date_range(start=start_date, end=end_date)
    ids_with_all_days = (
        df_filtered.groupby("ID")["DATE"]
        .apply(lambda dates: set(date_range_days).issubset(set(dates)))
        .reset_index()
    )
    ids_with_all_days = ids_with_all_days[ids_with_all_days["DATE"] == True]["ID"]

    # Filter to only include IDs that are present for all days
    df_filtered = df_filtered[df_filtered["ID"].isin(ids_with_all_days)]

    # Check IDs that contain all required elements
    ids_with_all_elements = (
        df_filtered.groupby("ID")["ELEMENT"]
        .apply(lambda elements: set(required_elements).issubset(set(elements)))
        .reset_index()
    )
    ids_with_all_elements = ids_with_all_elements[ids_with_all_elements["ELEMENT"] == True]["ID"]

    # Filter to only include IDs that have all required elements
    df_filtered = df_filtered[df_filtered["ID"].isin(ids_with_all_elements)]

    return df_filtered

# filter the stations
df_nodes = load_stations('.\\ghcnd-stations.txt')
lower_48_states = [
    'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'IA', 'ID',
    'IL', 'IN', 'KS', 'KY', 'LA', 'MA', 'MD', 'ME', 'MI', 'MN', 'MO',
    'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ', 'NM', 'NV', 'NY', 'OH',
    'OK', 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VA', 'VT',
    'WA', 'WI', 'WV', 'WY'
]
df_L48_nodes = df_nodes[df_nodes["STATE"].isin(lower_48_states)]
stations_L48 = [item + '.dly' for item in df_L48_nodes['ID'].unique()]

# dly files from selected stations to csv (can exclude US1 files, etc. from the dir for US only)
directory = '.\\extracted_files\ghcnd_all\\'
output_file = "dly_files_filtered.csv"
process_dly_files(directory, output_file, stations_L48)

#%% Setup data for modeling

# apply refined filters
df = pd.read_csv('dly_files_filtered.csv', parse_dates=['DATE'])
df = df.loc[df['ELEMENT'].isin(["PRCP", "SNOW", "SNWD", "TMAX", "TMIN"]), 
    ['ID', 'DATE', 'ELEMENT', 'VALUE']]
date_range = ("2020-01-01", "2025-01-22")
required_elements = ["PRCP", "SNOW", "SNWD", "TMAX", "TMIN"] 
filtered_data = filter_data(df, date_range, required_elements)

# pivot and save
df_pivoted = filtered_data.pivot_table(
    index=["ID", "DATE"],  # Keep ID and DATE as the index
    columns="ELEMENT",    # Pivot on the ELEMENT column
    values="VALUE",       # Use VALUE as the data to fill the columns
    aggfunc="first"       # Use the first value if duplicates exist (shouldn't happen for clean data)
).reset_index()

# Flatten the column index after pivoting
df_pivoted.columns.name = None  # Remove the column name
df_pivoted.columns = [col if isinstance(col, str) else col for col in df_pivoted.columns]
df_pivoted = df_pivoted.groupby("ID").filter(lambda group: not group.isna().any().any())
df_pivoted.to_csv('df_events.csv', index=False)

# align the stations to the filtered dly data
df_events = pd.read_csv('df_events.csv', parse_dates=['DATE'])
stations = df_events['ID'].unique()
df_nodes = df_nodes[df_nodes["ID"].isin(stations)]
df_nodes.to_csv('df_nodes.csv', index=False)