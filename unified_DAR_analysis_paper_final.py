########################
#python code to automatize mean DAR/Structural Integrity calculation for different DMPK studies and method applicability using BiopharmaFinder software by Andrea Di Ianni and LC-HRMS lab Merck Ivrea
########################

################################ PACKAGES #######################################################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import openpyxl

#### GUI ####

import tkinter
from tkinter import messagebox
from tkinter import filedialog

main_win = tkinter.Tk()
main_win.geometry("600x700")
main_win.title('DAR Analysis')
main_win.sourceFile = ''
main_win.Mol_ID = ''
main_win.x_axis = ''
main_win.studyID = ''
main_win.plasma = ''
main_win.studytype = ''
main_win.output_value =''

#### Create a LABEL for the entry widget Mol_ID
label_mol = tkinter.Label(main_win, text="Enter the molecule ID:")
label_mol.pack(padx=10, pady=10)

# Create an entry widget for Mol_ID
entry_mol = tkinter.Entry(main_win)
entry_mol.pack()
entry_mol.focus_set()

#### Create a LABEL for the entry widget x_axis
label_x_axis = tkinter.Label(main_win, text="Define your x axis (conc./time):")
label_x_axis.pack(padx=10, pady=10)

# Create an entry widget x_axis
entry_x_axis = tkinter.Entry(main_win)
entry_x_axis.pack(pady=10)
entry_x_axis.focus_set()

#### Create a LABEL for the entry widget DMPK study ID
label_studyID = tkinter.Label(main_win, text="Enter the study name (YYIV-0ABC):")
label_studyID.pack(padx=10, pady=10)

# Create an entry widget for DMPK study ID
entry_studyID = tkinter.Entry(main_win)
entry_studyID.pack(pady=10)
entry_studyID.focus_set()

### Create a LABEL for the entry widget study type
label_studytype = tkinter.Label(main_win, text=" Enter the type of LC-HRMS study (mice strain PK/applicability/in vitro stability):")
label_studytype.pack(padx=10, pady=10)

#Create an entry for the study type
entry_studytype = tkinter.Entry(main_win)
entry_studytype.pack(pady=10)
entry_studytype.focus_set()

#### Create a LABEL for the entry widget plasma matrix
label_plasma = tkinter.Label(main_win, text="Enter the matrix:")
label_plasma.pack(padx=10, pady=10)

#Create an entry widget for matrix plasma
entry_plasma = tkinter.Entry(main_win)
entry_plasma.pack(pady=10)
entry_plasma.focus_set()

### Create a LABEL for the entry widget for the desired output value
label_output_value = tkinter.Label(main_win, text=" Enter the desired output value to be calculated (DAR / SI):")
label_output_value.pack(padx=10, pady=10)

# Create an entry widget for the output value
entry_output_value = tkinter.Entry(main_win)
entry_output_value.pack(pady=10)
entry_output_value.focus_set()

# Function to retrieve the value from the entry widget
def get_entry_value():
    main_win.Mol_ID = entry_mol.get()  # Get the value from the entry widget
    print(f"Molecule ID: {main_win.Mol_ID}")
    main_win.x_axis = entry_x_axis.get()  # Get the value from the entry widget
    print(f"x axis definition: {main_win.x_axis}")
    main_win.studyID = entry_studyID.get()  # Get the value from the entry widget
    print(f"studyID definition: {main_win.studyID}")
    main_win.plasma = entry_plasma.get()  # Get the value from the entry widget
    print(f"matrix definition: {main_win.plasma}")
    main_win.studytype = entry_studytype.get() #Get the value from the entry widget
    print(f"study type definition : {main_win.studytype}")
    main_win.output_value = entry_output_value.get() #Get the value from the entry widget
    print(f"output value definition : {main_win.output_value}")

# Create a button to trigger the value retrieval
button = tkinter.Button(main_win, text="Submit", command=get_entry_value)
button.pack(pady=8)
button.focus_set()

def chooseFile():
    main_win.sourceFile = filedialog.askopenfilename(parent=main_win, initialdir="/", title='Please select a directory')

b_chooseFile = tkinter.Button(main_win, text="Browse input file", width=20, height=3, command=chooseFile)
b_chooseFile.pack(pady=8)
b_chooseFile.focus_set()

def close_window():
    main_win.destroy()

button = tkinter.Button(main_win, text="Quit the GUI and get plots", command=close_window)
button.pack()

main_win.mainloop()

print(main_win.sourceFile)
print(main_win.Mol_ID)
print(main_win.x_axis)
print(main_win.studyID)
print(main_win.plasma)
print(main_win.studytype)
print(main_win.output_value)

########################################## LOAD DATA ###########################################################################################################

# Load data
filepath = (main_win.sourceFile)
excluded_sheets = ["DAR summary", "Summary deglycosylated forms", "Summary", "Mouse_plasma01", "Mouse_plasma02", "LIst of glycans"]
# Load all sheets except "DAR summary"
xls = pd.ExcelFile(filepath)
sheets_to_read = [sheet for sheet in xls.sheet_names if sheet not in excluded_sheets]

# Read and concatenate all sheets
df_list = []
for sheet in sheets_to_read:
    temp_df = pd.read_excel(xls, sheet_name=sheet)
    temp_df[main_win.x_axis] = sheet  # Add sheet name as x_axis
    df_list.append(temp_df)

# Combine all sheets into a single DataFrame
df = pd.concat(df_list, ignore_index=True)

# Convert relevant columns to numeric, forcing errors to NaN
numeric_columns = ['Matched Mass Error (ppm)', 'Fractional Abundance', 'Relative Abundance', 'DAR', 'Sum Intensity']
for col in numeric_columns:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

# Filter the dataset based on your criteria
threshold_value = 85  # Example threshold for Matched Mass Error in ppm
fractional_abundance_threshold_value = 4 # Example threshold for Fractional Abundance in %

filtered_df = df[(df['Matched Mass Error (ppm)'] < threshold_value) &
                 (df['Relative Abundance'] > fractional_abundance_threshold_value) &
                 (df['DAR'].notna())].copy()

filtered_df_explained_area = df[(df['Matched Mass Error (ppm)'] < threshold_value) &
                 (df['Relative Abundance'] > fractional_abundance_threshold_value)
]

# Calculations for DAR Relative abundance
filtered_df['Product'] = filtered_df['DAR'] * filtered_df['Sum Intensity']

# Calculate the mean DAR for each x_axis unit
mean_dar_by_x_axis = (
    filtered_df
    .groupby(main_win.x_axis)
    .apply(lambda group: group['Product'].sum() / group['Sum Intensity'].sum() if group['Sum Intensity'].sum() > 0 else 0)
    .reset_index(name='Mean DAR')
)

# Calculations for DAR Explained AREA
def calculate_explained_area(group):
    total_identified_area = group.loc[group['DAR'].notna(), 'Sum Intensity'].sum()
    total_sum_intensity = group['Sum Intensity'].sum()
    explained_area_percentage = (total_identified_area / total_sum_intensity) * 100 if total_sum_intensity > 0 else 0
    return pd.Series({
        'Total Identified Area': total_identified_area,
        'Total Sum Intensity': total_sum_intensity,
        'Explained Area Percentage': explained_area_percentage
    })

# Apply the calculation to each group
explained_area_by_x_axis = filtered_df_explained_area.groupby(main_win.x_axis).apply(calculate_explained_area).reset_index()

# Function to extract numeric part from x_axis (handles decimals)
def extract_numeric_x_axis(tp):
    match = re.search(r"[-+]?\d*\.\d+|\d+", str(tp))  # Extract float or int
    return float(match.group()) if match else float('inf')  # Convert to float, default to inf for safety

# Ensure x_axis is sorted correctly
mean_dar_by_x_axis[main_win.x_axis] = pd.Categorical(
    mean_dar_by_x_axis[main_win.x_axis],
    categories=sorted(mean_dar_by_x_axis[main_win.x_axis].unique(), key=extract_numeric_x_axis),
    ordered=True
)

# Sort DataFrame to maintain order in the plot
mean_dar_by_x_axis = mean_dar_by_x_axis.sort_values(by=main_win.x_axis, key=lambda x: x.map(extract_numeric_x_axis))
explained_area_by_x_axis = explained_area_by_x_axis.sort_values(by=main_win.x_axis, key=lambda x: x.map(extract_numeric_x_axis))

# Print the mean DAR values for each x_axis value
print(mean_dar_by_x_axis)

################ Calculate DAR species

# Calculate the product of DAR and Sum Intensity
filtered_df['Product'] = filtered_df['DAR'] * filtered_df['Sum Intensity']

# Fill NaN in 'Modification' with empty string for concatenation
filtered_df['Modification'] = filtered_df['Modification'].fillna('')

# Create a unique identifier for each species
filtered_df['Species_ID'] = (
    filtered_df.groupby(['Chain', 'Modification']).ngroup() + 1
)

# Store mapping of Species_ID to Chain + Modification
species_mapping = filtered_df[['Species_ID', 'Chain', 'Modification']].drop_duplicates()

# Update species_grouped DataFrame using Species_ID
species_grouped = (
    filtered_df
    .groupby([main_win.x_axis, 'DAR', 'Species_ID'])
    ['Sum Intensity'].sum()
    .reset_index()
)

# Add total Sum Intensity per x_axis values
species_grouped[f'Total Intensity for {main_win.x_axis}'] = (
    species_grouped.groupby(main_win.x_axis)['Sum Intensity'].transform('sum')
)

# Calculate the relative percentage
species_grouped['Relative Percentage'] = (
    species_grouped['Sum Intensity'] / species_grouped[f'Total Intensity for {main_win.x_axis}'] * 100
)

# Map the relative percentage back to the filtered DataFrame
filtered_df = filtered_df.merge(
    species_grouped[[main_win.x_axis, 'DAR', 'Species_ID', 'Relative Percentage']],
    how='left',
    on=[main_win.x_axis, 'DAR', 'Species_ID']
)

# Function to extract numeric part from x_axis values (handles decimals)
def extract_numeric_x_axis(tp):
    match = re.search(r"[-+]?\d*\.\d+|\d+", str(tp))  # Extract float or int
    return float(match.group()) if match else float('inf')  # Convert to float, default to inf for safety

# Ensure x_axis value is sorted correctly
species_grouped[main_win.x_axis] = pd.Categorical(
    species_grouped[main_win.x_axis],
    categories=sorted(species_grouped[main_win.x_axis].unique(), key=extract_numeric_x_axis),
    ordered=True
)

# Define color palette for x_axis_values
x_axis_values = sorted(species_grouped[main_win.x_axis].unique())
colors = sns.color_palette("tab10", len(x_axis_values))
color_mapping = dict(zip(x_axis_values, colors))

############# PLOTTING ###########################################################################################################

# PLOTTING for DAR Relative abundance
plt.figure(figsize=(10, 6))
plt.plot(mean_dar_by_x_axis[main_win.x_axis], mean_dar_by_x_axis['Mean DAR'], marker='o', color='blue', linestyle='-')
plt.xlabel(main_win.x_axis, fontsize=20)
plt.ylabel(f'Mean {main_win.output_value} Value', fontsize=20)
plt.title(f'{main_win.studyID} {main_win.studytype} {main_win.Mol_ID} Mean {main_win.output_value} Value by {main_win.x_axis} {main_win.plasma}', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
# Set y-axis limits to max Mean DAR + 0.1
max_mean_dar = mean_dar_by_x_axis['Mean DAR'].max()
max_mean_dar_plus_20 = max_mean_dar + max_mean_dar * 0.20
max_mean_dar_minus_20 = max_mean_dar - max_mean_dar * 0.20
max_mean_dar_minus_50 = max_mean_dar - max_mean_dar * 0.50
max_mean_dar_minus_100 = max_mean_dar - max_mean_dar * 1.0
plt.axhline (y= max_mean_dar_plus_20, c='green')
plt.axhline (y=max_mean_dar_minus_20, c='green')
plt.fill_between(mean_dar_by_x_axis[main_win.x_axis], max_mean_dar_plus_20, max_mean_dar_minus_20, color='green', alpha=0.2, label='±20% Interval')
plt.fill_between(mean_dar_by_x_axis[main_win.x_axis], max_mean_dar_minus_20, max_mean_dar_minus_50, color='yellow', alpha=0.2)
plt.fill_between(mean_dar_by_x_axis[main_win.x_axis], max_mean_dar_minus_50, max_mean_dar_minus_100, color='red', alpha=0.2)
plt.ylim(0, max_mean_dar_plus_20 + 0.1)
plt.grid()
plt.xticks(rotation=0)  # Rotate x-axis labels if necessary
plt.tight_layout()  # Adjust layout to prevent clipping
plt.savefig(f'{main_win.Mol_ID}_mean_{main_win.output_value}_analysis_{main_win.plasma}_{main_win.studytype}_samples_append.png', dpi=300, bbox_inches='tight')
plt.show()

# PLOTTING for EXPLAINED AREA
plt.figure(figsize=(10, 6))
plt.plot(
    explained_area_by_x_axis[main_win.x_axis],
    explained_area_by_x_axis['Explained Area Percentage'],
    marker='o',
    color='blue',
    linestyle='-'
)
plt.xlabel(main_win.x_axis, fontsize=20)
plt.ylabel('Explained Area Percentage (%)', fontsize=20)
plt.title(f'{main_win.studyID} {main_win.studytype} {main_win.Mol_ID} Explained Area Percentage (%)', fontsize=20)
plt.ylim(0, 101)
plt.axhline(y=75, c='r', linestyle='dotted')
plt.grid()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig(f"Explained Area Percentage {main_win.studytype} {main_win.Mol_ID}.png", dpi=300)
plt.show()
print(explained_area_by_x_axis)

# PLOTTING for DAR SPECIES
for dar, group in species_grouped.groupby('DAR'):
    plt.figure(figsize=(12, 8))

    # Pivot table for plotting
    pivot_table = group.pivot(index='Species_ID', columns=main_win.x_axis, values='Relative Percentage').fillna(0)
    pivot_table.index = pivot_table.index.astype(str)  # Convert Species_ID to string for labeling

    # Generate bar plot
    pivot_table.plot(kind='bar', color=[color_mapping[tp] for tp in pivot_table.columns])

    plt.title(f'Relative Percentage of Species by {main_win.x_axis} for DAR {dar}')
    plt.xlabel('Species ID')
    plt.ylabel('Relative Percentage (%)')
    plt.legend(title=main_win.x_axis, bbox_to_anchor=(1.15, 0.8))
    plt.grid(axis='y')

    plt.xticks(rotation=0, fontsize=10)

    plt.savefig(f'{dar}_{main_win.Mol_ID}_{main_win.studytype}_{main_win.output_value}_analysis_%Relative_area_percentage.png', dpi=300, bbox_inches='tight')
    plt.show()

########################################################saving in excel file######################################################################

# Save the combined DataFrame to the original Excel file in "DAR summary" sheet
with pd.ExcelWriter(filepath, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    filtered_df.to_excel(writer, sheet_name="DAR summary", index=False)

# Write the filtered DataFrame to an Excel file (DAR explained and abundance)
filtered_df.to_excel(f'filtered_data_{main_win.studyID} {main_win.studytype} {main_win.Mol_ID}_{main_win.plasma}_samples_append.xlsx', index=False)

# Save Species_ID mapping for reference
species_mapping.to_excel(f'Species_ID_Mapping_{main_win.studyID} {main_win.studytype} {main_win.Mol_ID}_{main_win.plasma}.xlsx', index=False)

print("DAR summary updated and plots saved successfully!")
