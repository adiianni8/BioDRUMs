# app.py
import itertools
import io
import itertools

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import seaborn as sns
import streamlit as st
from networkx.algorithms import community


def calculate_proteoforms(
    antibody_masses,
    dar,
    linker_mass,
    delta_mass,
    degradation_masses,
    glycosylation_masses,
    accuracy_ppm,
    glycan_multiplier=2  # 2 for intact mAb, 1 for chain
):
    """
    Compute proteoforms, isobaric pairs, plots, and network visualization.

    glycan_multiplier controls how many glycans are summed in glycosylated forms:
      - 2 for intact antibodies (two chains -> two glycans)
      - 1 for single chain species
    """
    # List to store all results
    all_proteoforms = []

    # Apply delta mass correction to linker payload
    corrected_linker_mass = linker_mass - delta_mass

    # Function to compute intact ADC mass
    def calculate_intact_adc_mass(ab_mass_val, dar_val, corr_linker_mass):
        return ab_mass_val + (dar_val * corr_linker_mass)

    # Function to generate all proteoform combinations
    def generate_proteoforms(
        proteoform_id,
        antibody_label,
        ab_mass_val,
        max_dar,
        corr_linker_mass,
        deg_masses,
        glyco_masses,
        tol_ppm,
        include_glycosylation=True,
        glycan_mult=2
    ):
        proteoforms = []
        for current_dar in range(max_dar, -1, -1):  # max to 0
            intact_mass = calculate_intact_adc_mass(ab_mass_val, current_dar, corr_linker_mass)

            # Max 4 degradation forms relative to current DAR
            for num_degradations in range(0, 5 - current_dar):
                degradation_combos = itertools.combinations(deg_masses.items(), num_degradations)
                for deg_combo in degradation_combos:
                    deg_mass = sum(d[1] for d in deg_combo)
                    deg_names = "+".join([d[0] for d in deg_combo]) if deg_combo else "None"
                    total_mass = intact_mass + deg_mass

                    # Proteoform without glycosylation
                    proteoforms.append({
                        'Proteoform ID': proteoform_id + " (Deglycosylated)",
                        'Antibody Label': antibody_label,
                        'DAR': current_dar,
                        'Degradations': deg_names,
                        'Glycosylation': "None",
                        'Calculated Mass (Da)': total_mass,
                        'Accuracy, ppm': tol_ppm
                    })

                    # Include glycosylation combinations
                    if include_glycosylation and glycan_mult > 0:
                        glyco_combinations = list(
                            itertools.combinations_with_replacement(
                                glyko := glyco_masses.items(), glycan_mult
                            )
                        )
                        for glyco_combo in glyco_combinations:
                            glyco_mass = sum(g[1] for g in glyco_combo)
                            glyco_names = "+".join([g[0] for g in glyco_combo])
                            total_mass_with_glyco = total_mass + glyco_mass
                            proteoforms.append({
                                'Proteoform ID': proteoform_id + " (Glycosylated)",
                                'Antibody Label': antibody_label,
                                'DAR': current_dar,
                                'Degradations': deg_names,
                                'Glycosylation': glyco_names,
                                'Calculated Mass (Da)': total_mass_with_glyco,
                                'Accuracy, ppm': tol_ppm
                            })
        return proteoforms

    # Generate proteoform combinations for each antibody mass
    for idx, ab_mass in enumerate(antibody_masses):
        antibody_label = f"Antibody Mass {idx + 1}"
        proteoform_id = f"Proteoform_{idx + 1}"
        data_deglyco = generate_proteoforms(
            proteoform_id,
            antibody_label,
            ab_mass,
            dar,
            corrected_linker_mass,
            degradation_masses,
            glycosylation_masses,
            accuracy_ppm,
            include_glycosylation=False,
            glycan_mult=glycan_multiplier
        )
        data_glyco = generate_proteoforms(
            proteoform_id,
            antibody_label,
            ab_mass,
            dar,
            corrected_linker_mass,
            degradation_masses,
            glycosylation_masses,
            accuracy_ppm,
            include_glycosylation=True,
            glycan_mult=glycan_multiplier
        )
        all_proteoforms.extend(data_deglyco)
        all_proteoforms.extend(data_glyco)

    # Convert results to a DataFrame
    output_df = pd.DataFrame(all_proteoforms)

    # Identify isobaric species and include full details
    isobaric_species = []
    for (_, row1), (_, row2) in itertools.combinations(output_df.iterrows(), 2):
        mass1 = row1['Calculated Mass (Da)']
        mass2 = row2['Calculated Mass (Da)']
        if mass2 == 0:
            continue
        ppm_difference = abs((mass1 - mass2) / mass2 * 1e6)

        # Ensure IDs are different before considering as isobaric
        if ppm_difference <= row1['Accuracy, ppm'] and row1['Proteoform ID'] != row2['Proteoform ID']:
            isobaric_species.append({
                'Proteoform 1 ID': row1['Proteoform ID'],
                'Proteoform 1 Antibody Label': row1['Antibody Label'],
                'Proteoform 1 DAR': row1['DAR'],
                'Proteoform 1 Degradations': row1['Degradations'],
                'Proteoform 1 Glycosylation': row1['Glycosylation'],
                'Proteoform 1 Mass (Da)': mass1,
                'Proteoform 2 ID': row2['Proteoform ID'],
                'Proteoform 2 Antibody Label': row2['Antibody Label'],
                'Proteoform 2 DAR': row2['DAR'],
                'Proteoform 2 Degradations': row2['Degradations'],
                'Proteoform 2 Glycosylation': row2['Glycosylation'],
                'Proteoform 2 Mass (Da)': mass2,
                'PPM Difference': ppm_difference
            })

    # Convert isobaric species to DataFrame
    isobaric_df = pd.DataFrame(isobaric_species)

    # Merge isobaric species info into output DataFrame
    def find_confounding_species(proteoform_id):
        matches = isobaric_df[isobaric_df['Proteoform 1 ID'] == proteoform_id]['Proteoform 2 ID'].tolist()
        matches += isobaric_df[isobaric_df['Proteoform 2 ID'] == proteoform_id]['Proteoform 1 ID'].tolist()
        return ", ".join(matches) if matches else "None"

    if not output_df.empty:
        output_df['Confounding Proteoforms'] = output_df['Proteoform ID'].apply(find_confounding_species)

    # Save results to an Excel file (in-memory for Streamlit)
    excel_buffer = io.BytesIO()
    with pd.ExcelWriter(excel_buffer, engine="xlsxwriter") as writer:
        output_df.to_excel(writer, sheet_name="Proteoforms", index=False)
        isobaric_df.to_excel(writer, sheet_name="Isobaric Species", index=False)
    excel_bytes = excel_buffer.getvalue()

    # ================== GROUPED BAR PLOT split by Glycosylation Status ==================
    # Labels present (based on provided masses)
    if not output_df.empty:
        antibody_labels = sorted(output_df['Antibody Label'].unique(), key=lambda s: int(s.split()[-1]))
    else:
        antibody_labels = [f"Antibody Mass {i+1}" for i in range(len(antibody_masses))]

    # DAR order: DARmax, DARmax-1, ..., DAR0
    dar_levels = [f'DAR{i}' for i in range(dar, -1, -1)]

    # Initialize counts dict: DAR -> Antibody -> GlycoStatus -> count
    isobaric_counts = {
        dar_label: {
            mass_label: {'Deglycosylated': 0, 'Glycosylated': 0}
            for mass_label in antibody_labels
        }
        for dar_label in dar_levels
    }

    # Count occurrences for both sides of each isobaric pair, split by glyco status
    if not isobaric_df.empty:
        for _, r in isobaric_df.iterrows():
            # Proteoform 1
            dar1 = f"DAR{r['Proteoform 1 DAR']}"
            mass1_label = r['Proteoform 1 Antibody Label']
            status1 = 'Deglycosylated' if r['Proteoform 1 Glycosylation'] == "None" else 'Glycosylated'
            if dar1 in isobaric_counts and mass1_label in isobaric_counts[dar1]:
                isobaric_counts[dar1][mass1_label][status1] += 1
            # Proteoform 2
            dar2 = f"DAR{r['Proteoform 2 DAR']}"
            mass2_label = r['Proteoform 2 Antibody Label']
            status2 = 'Deglycosylated' if r['Proteoform 2 Glycosylation'] == "None" else 'Glycosylated'
            if dar2 in isobaric_counts and mass2_label in isobaric_counts[dar2]:
                isobaric_counts[dar2][mass2_label][status2] += 1

    # Convert counts to DataFrame
    plot_rows = []
    for dar_label in dar_levels:
        for mass_label in antibody_labels:
            for status_label in ['Deglycosylated', 'Glycosylated']:
                plot_rows.append({
                    'DAR Level': dar_label,
                    'Antibody Mass': mass_label,
                    'Glycosylation Status': status_label,
                    'Isobaric Species Count': isobaric_counts[dar_label][mass_label][status_label]
                })
    plot_df = pd.DataFrame(plot_rows)

    # Create two panels: Deglycosylated vs Glycosylated
    sns.set_theme(style="whitegrid")
    g = sns.catplot(
        data=plot_df,
        x='DAR Level',
        y='Isobaric Species Count',
        hue='Antibody Mass',
        col='Glycosylation Status',
        kind='bar',
        order=dar_levels,
        hue_order=antibody_labels,
        palette="Blues_r",
        height=5,
        aspect=1.2
    )
    g.set_axis_labels("DAR Level", "Isobaric Species Count")
    g.set_titles("{col_name}")
    # Ticks and grid
    for ax in g.axes.flat:
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.grid(axis="y", linestyle="--", alpha=0.6)
    fig_path = "Isobaric_counts_deglycosylated_vs_glycosylated.png"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()

    # ======================= NETWORK VISUALIZATION =======================
    def plot_isobaric_network(isobaric_data: pd.DataFrame, max_nodes=200):
        if isobaric_data.empty:
            return None
        G = nx.Graph()
        all_ids = set(isobaric_data["Proteoform 1 ID"]).union(isobaric_data["Proteoform 2 ID"])
        for node_id in all_ids:
            G.add_node(node_id)
        for _, r in isobaric_data.iterrows():
            p1 = r["Proteoform 1 ID"]
            p2 = r["Proteoform 2 ID"]
            ppm_diff = r["PPM Difference"]
            G.add_edge(p1, p2, weight=(1 / (1 + ppm_diff)))
        # Louvain communities (requires appropriate NetworkX version)
        try:
            coms = community.louvain_communities(G, weight='weight')
        except Exception:
            coms = community.greedy_modularity_communities(G, weight='weight')
        community_mapping = {}
        for i, comm_set in enumerate(coms):
            for node in comm_set:
                community_mapping[node] = i
        color_map = plt.get_cmap('tab10')
        plt.figure(figsize=(8, 8))
        pos = nx.spring_layout(G, seed=42, weight='weight', k=0.1)
        nx.draw_networkx_nodes(
            G, pos,
            node_size=500,
            node_color=[color_map(community_mapping.get(node, 0)) for node in G.nodes()]
        )
        nx.draw_networkx_edges(G, pos, width=1.5, edge_color="gray")
        nx.draw_networkx_labels(G, pos, font_size=8, font_color="black")
        plt.title("Isobaric Species Network")
        plt.axis("off")
        plt.tight_layout()
        return plt

    network_fig = plot_isobaric_network(isobaric_df)

    # Return artifacts
    return {
        "excel_bytes": excel_bytes,
        "excel_filename": "ADC_Proteoforms_Output.xlsx",
        "plot_path": fig_path,
        "network_plt": network_fig
    }


# ============================ STREAMLIT APP ============================

st.title("Proteoform Mass Finder for ADC")
st.write("This is a deconvoluted mass finder tool for ADC proteoforms.")

with st.expander("Inputs"):
    st.markdown("Enter antibody masses (Da). Leave blank entries empty.")
    col1, col2 = st.columns(2)
    with col2:
        antibody_mass1_v = st.number_input("Antibody Mass 1 (Da) (numeric)", value=0.0, step=0.1, format="%.4f")
        antibody_mass2_v = st.number_input("Antibody Mass 2 (Da) (numeric)", value=0.0, step=0.1, format="%.4f")
        antibody_mass3_v = st.number_input("Antibody Mass 3 (Da) (numeric)", value=0.0, step=0.1, format="%.4f")
        antibody_mass4_v = st.number_input("Antibody Mass 4 (Da) (numeric)", value=0.0, step=0.1, format="%.4f")
        antibody_mass5_v = st.number_input("Antibody Mass 5 (Da) (numeric)", value=0.0, step=0.1, format="%.4f")

    # Use numeric inputs to build masses (keep both for clarity)
    antibody_masses = [
        antibody_mass1_v,
        antibody_mass2_v,
        antibody_mass3_v,
        antibody_mass4_v,
        antibody_mass5_v,
    ]

    st.markdown("---")

    st.subheader("DAR and Linker")
    dar_val = st.number_input("DAR (max)", value=2, min_value=0, max_value=20, step=1)
    linker_mass_val = st.number_input("Linker Payload Mass (Da)", value=1500.0, step=0.1)
    delta_mass_val = st.number_input("Delta Mass (Da)", value=0.0, step=0.1)

    st.markdown("---")

    st.subheader("Degradation Masses (Da)")
    deg1 = st.number_input("Degradation Mass 1 (Da)", value=0.0, step=0.1)
    deg2 = st.number_input("Degradation Mass 2 (Da)", value=0.0, step=0.1)
    deg3 = st.number_input("Degradation Mass 3 (Da)", value=0.0, step=0.1)
    deg4 = st.number_input("Degradation Mass 4 (Da)", value=0.0, step=0.1)

    degradation_masses = {
        'deg1': deg1,
        'deg2': deg2,
        'deg3': deg3,
        'deg4': deg4,
    }

    st.markdown("---")

    st.subheader("Glycosylation Masses (Da)")
    gly_a2g0f = st.number_input("A2G0F (Da)", value=1444.35, step=0.1)
    gly_a2g0 = st.number_input("A2G0 (Da)", value=1299.18, step=0.1)
    gly_a2g1f = st.number_input("A2G1F (Da)", value=1607.46, step=0.1)
    gly_a2g2f = st.number_input("A2G2F (Da)", value=1769.60, step=0.1)

    glycosylation_masses = {
        'A2G0F': gly_a2g0f,
        'A2G0': gly_a2g0,
        'A2G1F': gly_a2g1f,
        'A2G2F': gly_a2g2f,
    }

    st.markdown("---")

    st.subheader("Accuracy (ppm)")
    accuracy_ppm_val = st.number_input("Accuracy (ppm)", value=5.0, step=0.1)

    st.markdown("---")

    st.subheader("Molecule Type, select Intact protein for an intact ADC, Chain for subunit analysis")
    molecule_type = st.selectbox("Molecule Type", ["Intact Protein", "Chain"])
    glycan_multiplier = 2 if molecule_type == "Intact Protein" else 1

    st.markdown("---")

    st.subheader("Number of Cysteines (Exclude Cys involved in conjugation), set 11 for the heavy chain and 5 for the light chain")
    num_cysteines = st.number_input("Number of Cysteines", value=32, min_value=0, step=1)

# Run calculation button
if st.button("Calculate Proteoforms"):
    # Build inputs similar to the original script
    # Ensure we only pass positive masses
    antibody_masses_input = [m for m in antibody_masses if m is not None]
    antibody_masses_input = [float(m) for m in antibody_masses_input if float(m) > -1e-9]

    # Adjust antibody masses for cysteines (subtract mass per cysteine)
    hydrogen_mass_per_cys = 1.0079
    adjusted_antibody_masses = [
        m - (num_cysteines * hydrogen_mass_per_cys) for m in antibody_masses_input
    ]

    results = calculate_proteoforms(
        adjusted_antibody_masses,
        int(dar_val),
        float(linker_mass_val),
        float(delta_mass_val),
        {k: float(v) for k, v in degradation_masses.items()},
        glycosylation_masses,
        float(accuracy_ppm_val),
        glycan_multiplier=glycan_multiplier
    )

    # Display and offer downloads
    if results is None:
        st.error("Calculation did not return results.")
    else:
        st.success("Proteoform calculations completed.")

        # Excel download
        st.download_button(
            label="Download Excel with Proteoforms",
            data=results["excel_bytes"],
            file_name=results["excel_filename"],
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )

        # Plot: Isobaric counts image
        st.image(results["plot_path"], caption="Isobaric counts (deglycosylated vs glycosylated)", use_column_width=True)

        # Network plot
        if results["network_plt"] is not None:
            # The network plot is produced via matplotlib; render it in Streamlit
            network_buf = io.BytesIO()
            results["network_plt"].tight_layout()
            results["network_plt"].savefig(network_buf, format="png", dpi=300)
            network_buf.seek(0)
            st.image(network_buf, caption="Isobaric Species Network", use_column_width=True)
