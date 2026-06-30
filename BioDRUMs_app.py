import io
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st

st.set_page_config(page_title="DAR Analysis", layout="wide")

st.title("BioDRUMs(Biologics Drug Ratio and Unified intact Mass analysis)")
st.write("Upload your DAR Excel sheet, set parameters, and run the mean DAR analysis.")

# Sidebar inputs
with st.sidebar:
    st.header("Inputs")
    mol_id = st.text_input("Molecule ID (ADC name/chain)", "")
    x_axis_name = st.text_input("Define your x axis (e.g., conc or time)", "")
    study_id = st.text_input("Study name (YY-0ABC)", "")
    study_type = st.text_input("Type of LC-HRMS study (e.g., mice strain PK / applicability / in vitro stability)", "")
    plasma = st.text_input("Matrix", "")
    output_value = st.text_input("Desired output value label (e.g., DAR / SI)", "DAR")

    st.divider()
    st.subheader("Filters")
    threshold_value = st.number_input("Matched Mass Error threshold (ppm)", value=30, min_value=0, step=1)
    fractional_abundance_threshold_value = st.number_input("Relative Abundance threshold (%)", value=4, min_value=0, step=1)

    st.divider()
    excluded_sheets_default = ["DAR summary", "Summary deglycosylated forms", "Summary"]
    exclude_input = st.text_area(
        "Excluded sheet names (one per line)",
        value="\n".join(excluded_sheets_default),
        height=90
    )
    excluded_sheets = [s.strip() for s in exclude_input.splitlines() if s.strip()]

uploaded_file = st.file_uploader("Upload Excel file (.xlsx) exported from BiopharmaFinder", type=["xlsx"])
run = st.button("Run analysis")

def extract_numeric_x_axis(tp):
    match = re.search(r"[-+]?\d*\.\d+|\d+", str(tp))
    return float(match.group()) if match else float('inf')

def validate_inputs():
    missing = []
    if not uploaded_file: missing.append("Excel file")
    if not x_axis_name: missing.append("x axis name")
    if missing:
        st.warning("Please provide: " + ", ".join(missing))
        return False
    return True

if run:
    if not validate_inputs():
        st.stop()

    try:
        # Load all sheets except excluded ones
        xls = pd.ExcelFile(uploaded_file, engine="openpyxl")
        sheets_to_read = [s for s in xls.sheet_names if s not in excluded_sheets]

        if not sheets_to_read:
            st.error("No sheets to read after exclusions. Check your exclusions or file.")
            st.stop()

        df_list = []
        for sheet in sheets_to_read:
            temp_df = pd.read_excel(xls, sheet_name=sheet)
            # Add sheet name as x_axis value
            temp_df[x_axis_name] = sheet
            df_list.append(temp_df)

        df = pd.concat(df_list, ignore_index=True)

        # Convert numeric columns
        numeric_columns = [
            "Matched Mass Error (ppm)",
            "Fractional Abundance",
            "Relative Abundance",
            "DAR",
            "Sum Intensity",
        ]
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        # Validate columns needed for core calculations
        required_cols = ["Matched Mass Error (ppm)", "Relative Abundance", "DAR", "Sum Intensity", x_axis_name]
        missing_cols = [c for c in required_cols if c not in df.columns]
        if missing_cols:
            st.error(f"Missing required columns in the uploaded data: {missing_cols}")
            st.stop()

        # Filtered datasets
        filtered_df = df[
            (df["Matched Mass Error (ppm)"] < threshold_value) &
            (df["Relative Abundance"] > fractional_abundance_threshold_value) &
            (df["DAR"].notna())
        ].copy()

        filtered_df_explained_area = df[
            (df["Matched Mass Error (ppm)"] < threshold_value) &
            (df["Relative Abundance"] > fractional_abundance_threshold_value)
        ].copy()

        # Mean DAR by x_axis_name
        filtered_df["Product"] = filtered_df["DAR"] * filtered_df["Sum Intensity"]
        mean_dar_by_x_axis = (
            filtered_df
            .groupby(x_axis_name)
            .apply(lambda g: g["Product"].sum() / g["Sum Intensity"].sum() if g["Sum Intensity"].sum() > 0 else 0)
            .reset_index(name="Mean DAR")
        )

        # Explained Area calculation
        def calculate_explained_area(group):
            total_identified_area = group.loc[group["DAR"].notna(), "Sum Intensity"].sum()
            total_sum_intensity = group["Sum Intensity"].sum()
            explained_area_percentage = (total_identified_area / total_sum_intensity) * 100 if total_sum_intensity > 0 else 0
            return pd.Series({
                "Total Identified Area": total_identified_area,
                "Total Sum Intensity": total_sum_intensity,
                "Explained Area Percentage": explained_area_percentage
            })

        explained_area_by_x_axis = (
            filtered_df_explained_area.groupby(x_axis_name).apply(calculate_explained_area).reset_index()
        )

        # Sort x-axis values numerically (based on numbers inside labels)
        if not mean_dar_by_x_axis.empty:
            mean_dar_by_x_axis[x_axis_name] = pd.Categorical(
                mean_dar_by_x_axis[x_axis_name],
                categories=sorted(mean_dar_by_x_axis[x_axis_name].unique(), key=extract_numeric_x_axis),
                ordered=True
            )
            mean_dar_by_x_axis = mean_dar_by_x_axis.sort_values(
                by=x_axis_name, key=lambda x: x.map(extract_numeric_x_axis)
            )

        if not explained_area_by_x_axis.empty:
            explained_area_by_x_axis[x_axis_name] = pd.Categorical(
                explained_area_by_x_axis[x_axis_name],
                categories=sorted(explained_area_by_x_axis[x_axis_name].unique(), key=extract_numeric_x_axis),
                ordered=True
            )
            explained_area_by_x_axis = explained_area_by_x_axis.sort_values(
                by=x_axis_name, key=lambda x: x.map(extract_numeric_x_axis)
            )

        # Display computed tables
        st.subheader("Mean DAR by x axis")
        st.dataframe(mean_dar_by_x_axis)

        st.subheader("Explained Area by x axis")
        st.dataframe(explained_area_by_x_axis)

        # Plot 1: DAR Relative abundance (Mean DAR)
        figs_download = []
        if not mean_dar_by_x_axis.empty:
            fig1, ax1 = plt.subplots(figsize=(10, 6))
            ax1.plot(mean_dar_by_x_axis[x_axis_name], mean_dar_by_x_axis["Mean DAR"], marker='o', color='blue', linestyle='-')
            ax1.set_xlabel(x_axis_name, fontsize=14)
            ax1.set_ylabel(f"Mean {output_value} Value", fontsize=14)
            ax1.set_title(f"{study_id} {study_type} {mol_id} Mean {output_value} by {x_axis_name} {plasma}", fontsize=14)

            max_mean_dar = mean_dar_by_x_axis["Mean DAR"].max()
            max_mean_dar_plus_20 = max_mean_dar + max_mean_dar * 0.20
            max_mean_dar_minus_20 = max_mean_dar - max_mean_dar * 0.20
            max_mean_dar_minus_50 = max_mean_dar - max_mean_dar * 0.50
            max_mean_dar_minus_100 = max_mean_dar - max_mean_dar * 1.0

            ax1.axhline(y=max_mean_dar_plus_20, c='green')
            ax1.axhline(y=max_mean_dar_minus_20, c='green')
            ax1.fill_between(mean_dar_by_x_axis[x_axis_name], max_mean_dar_plus_20, max_mean_dar_minus_20, color='green', alpha=0.2, label='±20% Interval')
            ax1.fill_between(mean_dar_by_x_axis[x_axis_name], max_mean_dar_minus_20, max_mean_dar_minus_50, color='yellow', alpha=0.2)
            ax1.fill_between(mean_dar_by_x_axis[x_axis_name], max_mean_dar_minus_50, max_mean_dar_minus_100, color='red', alpha=0.2)
            ax1.set_ylim(0, max_mean_dar_plus_20 + 0.1)
            ax1.grid(True)
            st.pyplot(fig1)

            buf1 = io.BytesIO()
            fig1.savefig(buf1, format="png", dpi=300, bbox_inches="tight")
            buf1.seek(0)
            figs_download.append(("mean_dar_plot.png", buf1))

        # Plot 2: Explained Area
        if not explained_area_by_x_axis.empty and "Explained Area Percentage" in explained_area_by_x_axis.columns:
            fig2, ax2 = plt.subplots(figsize=(10, 6))
            ax2.plot(
                explained_area_by_x_axis[x_axis_name],
                explained_area_by_x_axis["Explained Area Percentage"],
                marker='o', color='blue', linestyle='-'
            )
            ax2.set_xlabel(x_axis_name, fontsize=14)
            ax2.set_ylabel("Explained Area Percentage (%)", fontsize=14)
            ax2.set_title(f"{study_id} {study_type} {mol_id} Explained Area Percentage (%)", fontsize=14)
            ax2.set_ylim(0, 101)
            ax2.axhline(y=75, c='r', linestyle='dotted')
            ax2.grid(True)
            st.pyplot(fig2)

            buf2 = io.BytesIO()
            fig2.savefig(buf2, format="png", dpi=300, bbox_inches="tight")
            buf2.seek(0)
            figs_download.append(("explained_area_plot.png", buf2))

        # DAR Species plot(s)
        species_mapping = None
        try:
            # Prep for species percentage per DAR
            filtered_df["Modification"] = filtered_df["Modification"].fillna('') if "Modification" in filtered_df.columns else ''
            if "Chain" in filtered_df.columns and "Modification" in filtered_df.columns:
                filtered_df["Species_ID"] = filtered_df.groupby(["Chain", "Modification"]).ngroup() + 1
                species_mapping = filtered_df[["Species_ID", "Chain", "Modification"]].drop_duplicates()

                species_grouped = (
                    filtered_df
                    .groupby([x_axis_name, "DAR", "Species_ID"])["Sum Intensity"].sum()
                    .reset_index()
                )
                species_grouped[f"Total Intensity for {x_axis_name}"] = (
                    species_grouped.groupby(x_axis_name)["Sum Intensity"].transform("sum")
                )
                species_grouped["Relative Percentage"] = (
                    species_grouped["Sum Intensity"] / species_grouped[f"Total Intensity for {x_axis_name}"] * 100
                )

                species_grouped[x_axis_name] = pd.Categorical(
                    species_grouped[x_axis_name],
                    categories=sorted(species_grouped[x_axis_name].unique(), key=extract_numeric_x_axis),
                    ordered=True
                )

                x_axis_values = sorted(species_grouped[x_axis_name].unique(), key=extract_numeric_x_axis)
                colors = sns.color_palette("tab10", len(x_axis_values))
                color_mapping = dict(zip(x_axis_values, colors))

                st.subheader("DAR Species Relative Percentages")
                for dar, group in species_grouped.groupby("DAR"):
                    pivot_table = group.pivot(index="Species_ID", columns=x_axis_name, values="Relative Percentage").fillna(0)
                    pivot_table.index = pivot_table.index.astype(str)

                    fig_s, ax_s = plt.subplots(figsize=(12, 8))
                    pivot_table.plot(kind="bar", color=[color_mapping[tp] for tp in pivot_table.columns], ax=ax_s)
                    ax_s.set_title(f"Relative Percentage of Species by {x_axis_name} for DAR {dar}")
                    ax_s.set_xlabel("Species ID")
                    ax_s.set_ylabel("Relative Percentage (%)")
                    ax_s.legend(title=x_axis_name, bbox_to_anchor=(1.15, 0.8))
                    ax_s.grid(axis='y')
                    ax_s.set_xticklabels(ax_s.get_xticklabels(), rotation=0, fontsize=10)

                    st.pyplot(fig_s)

                    buf_s = io.BytesIO()
                    fname = f"dar_{str(dar).replace('.', '_')}_{mol_id}_{study_type}_{output_value}_rel_area.png"
                    fig_s.savefig(buf_s, format="png", dpi=300, bbox_inches="tight")
                    buf_s.seek(0)
                    figs_download.append((fname, buf_s))
            else:
                st.info("Skipping species plots: required columns 'Chain' and/or 'Modification' not found.")
        except Exception as e:
            st.warning(f"Species plotting skipped due to error: {e}")

        # Downloads: images
        st.subheader("Download plots")
        for name, buf in figs_download:
            st.download_button(
                label=f"Download {name}",
                data=buf,
                file_name=name,
                mime="image/png"
            )

        # Downloads: Excel outputs
        st.subheader("Download data")
        # filtered data
        out1 = io.BytesIO()
        with pd.ExcelWriter(out1, engine="openpyxl") as writer:
            filtered_df.to_excel(writer, sheet_name="DAR summary", index=False)
        out1.seek(0)
        st.download_button(
            label="Download DAR_summary.xlsx",
            data=out1,
            file_name=f"DAR_summary_{study_id}_{study_type}_{mol_id}_{plasma}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )

        # species mapping (if available)
        if species_mapping is not None:
            out2 = io.BytesIO()
            with pd.ExcelWriter(out2, engine="openpyxl") as writer:
                species_mapping.to_excel(writer, sheet_name="Species_ID_Mapping", index=False)
            out2.seek(0)
            st.download_button(
                label="Download Species_ID_Mapping.xlsx",
                data=out2,
                file_name=f"Species_ID_Mapping_{study_id}_{study_type}_{mol_id}_{plasma}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )

        st.success("Analysis completed.")
    except Exception as e:
        st.error(f"An error occurred: {e}")
