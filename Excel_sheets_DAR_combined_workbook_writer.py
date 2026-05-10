import pandas as pd
import glob
import os

# Folder containing your Excel files
folder_path = "pathtoDARexcelsheets"

# Get all Excel files in the folder
excel_files = glob.glob(os.path.join(folder_path, "*.xlsx"))
# Mapping of old column names to new ones
header_mapping = {"Sequence Name": "Protein Name"}

# Output file
output_file = "combined_workbook_ADC.xlsx"

# Create a writer object
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    for file in excel_files:
        # Read the Excel file
        df = pd.read_excel(file)
        # Rename columns
        df.rename(columns=header_mapping, inplace=True)
        # Insert new columns at the beginning
        df.insert(0, "Chain", "")  # Column A
        df.insert(1, "DAR", "")  # Column B

        # Use the file name (without extension) as sheet name
        file_name = os.path.splitext(os.path.basename(file))[0]

        # Get the last 10 characters of the file name
        sheet_name = file_name[-10:]

        # Write to the workbook
        df.to_excel(writer, sheet_name=sheet_name, index=False)

print(f"✅ Combined {len(excel_files)} files into '{output_file}' successfully!")
