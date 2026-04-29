import itertools
import tkinter as tk
from tkinter import messagebox

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import seaborn as sns
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
      - 2 for intact antibodies (two chains → two glycans)
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
        for current_dar in range(max_dar, -1, -1):  # Iterate from max DAR down to DAR0
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
                        # glycan_mult = 2 → two glycans (intact mAb)
                        # glycan_mult = 1 → one glycan (chain)
                        glyco_combinations = list(
                            itertools.combinations_with_replacement(
                                glyco_masses.items(), glycan_mult
                            )
                        )
                        for glyco_combo in glyco_combinations:
                            glyco_mass = sum(g[1] for g in glyco_combo)
                            glyco_names = "+".join([g[0] for g in glyco_combo])
                            total_mass_with_glyco = intact_mass + deg_mass + glyco_mass
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

        # Ensure IDs are different before considering as isobaric and use tolerance from row1
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

    # Save results to an Excel file
    output_filename = "ADC_Proteoforms_Output.xlsx"
    with pd.ExcelWriter(output_filename) as writer:
        output_df.to_excel(writer, sheet_name="Proteoforms", index=False)
        isobaric_df.to_excel(writer, sheet_name="Isobaric Species", index=False)

    print(f"Proteoform calculations completed! Results saved in {output_filename}")

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
    plt.savefig("Isobaric_counts_deglycosylated_vs_glycosylated.tiff", dpi=300, bbox_inches="tight")
    plt.show()

    # ======================= NETWORK VISUALIZATION =======================
    def plot_isobaric_network(isobaric_data: pd.DataFrame):
        if isobaric_data.empty:
            return
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
            comms = community.louvain_communities(G, weight='weight')
        except Exception:
            comms = community.greedy_modularity_communities(G, weight='weight')
        community_mapping = {}
        for i, comm_set in enumerate(comms):
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
        plt.show()

    plot_isobaric_network(isobaric_df)


def on_calculate():
    try:
        # Read and validate inputs
        antibody_masses = [
            float(entry_antibody_mass1.get() or 0),
            float(entry_antibody_mass2.get() or 0),
            float(entry_antibody_mass3.get() or 0),
            float(entry_antibody_mass4.get() or 0),
            float(entry_antibody_mass5.get() or 0),
        ]
        dar_val = int(entry_dar.get())
        linker_mass_val = float(entry_linker_mass.get())
        delta_mass_val = float(entry_delta_mass.get())

        degradation_masses = {
            'deg1': float(entry_deg1.get() or 0),
            'deg2': float(entry_deg2.get() or 0),
            'deg3': float(entry_deg3.get() or 0),
            'deg4': float(entry_deg4.get() or 0),
        }
        glycosylation_masses = {
            'A2G0F': float(entry_glyco_a2g0f.get() or 0),
            'A2G0':  float(entry_glyco_a2g0.get() or 0),
            'A2G1F': float(entry_glyco_a2g1f.get() or 0),
            'A2G2F': float(entry_glyco_a2g2f.get() or 0),
        }
        accuracy_ppm_val = float(entry_accuracy_ppm.get())

        # Adjust based on cysteine count
        num_cys = int(entry_num_cysteines.get())
        hydrogen_mass_per_cys = 1.0079
        adjusted_antibody_masses = [m - (num_cys * hydrogen_mass_per_cys) for m in antibody_masses]

        # Determine glycan multiplier from dropdown
        molecule_type = molecule_type_var.get()
        glycan_multiplier = 2 if molecule_type == "Intact Protein" else 1

        calculate_proteoforms(
            adjusted_antibody_masses,
            dar_val,
            linker_mass_val,
            delta_mass_val,
            degradation_masses,
            glycosylation_masses,
            accuracy_ppm_val,
            glycan_multiplier=glycan_multiplier
        )
    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers in all fields.")


# ============================ GUI ============================

root = tk.Tk()
root.title("Proteoform Mass Calculator")
root.geometry("520x600")

# Antibody masses
tk.Label(root, text="Antibody Mass 1 (Da):").grid(row=0, column=0, sticky="w")
entry_antibody_mass1 = tk.Entry(root)
entry_antibody_mass1.grid(row=0, column=1)

tk.Label(root, text="Antibody Mass 2 (Da):").grid(row=1, column=0, sticky="w")
entry_antibody_mass2 = tk.Entry(root)
entry_antibody_mass2.grid(row=1, column=1)

tk.Label(root, text="Antibody Mass 3 (Da):").grid(row=2, column=0, sticky="w")
entry_antibody_mass3 = tk.Entry(root)
entry_antibody_mass3.grid(row=2, column=1)

tk.Label(root, text="Antibody Mass 4 (Da):").grid(row=3, column=0, sticky="w")
entry_antibody_mass4 = tk.Entry(root)
entry_antibody_mass4.grid(row=3, column=1)

tk.Label(root, text="Antibody Mass 5 (Da):").grid(row=4, column=0, sticky="w")
entry_antibody_mass5 = tk.Entry(root)
entry_antibody_mass5.grid(row=4, column=1)

# DAR and linker
tk.Label(root, text="DAR (max):").grid(row=5, column=0, sticky="w")
entry_dar = tk.Entry(root)
entry_dar.grid(row=5, column=1)

tk.Label(root, text="Linker Payload Mass (Da):").grid(row=6, column=0, sticky="w")
entry_linker_mass = tk.Entry(root)
entry_linker_mass.grid(row=6, column=1)

tk.Label(root, text="Delta Mass (Da):").grid(row=7, column=0, sticky="w")
entry_delta_mass = tk.Entry(root)
entry_delta_mass.grid(row=7, column=1)

# Degradation masses
tk.Label(root, text="Degradation Mass 1 (Da):").grid(row=8, column=0, sticky="w")
entry_deg1 = tk.Entry(root)
entry_deg1.grid(row=8, column=1)

tk.Label(root, text="Degradation Mass 2 (Da):").grid(row=9, column=0, sticky="w")
entry_deg2 = tk.Entry(root)
entry_deg2.grid(row=9, column=1)

tk.Label(root, text="Degradation Mass 3 (Da):").grid(row=10, column=0, sticky="w")
entry_deg3 = tk.Entry(root)
entry_deg3.grid(row=10, column=1)

tk.Label(root, text="Degradation Mass 4 (Da):").grid(row=11, column=0, sticky="w")
entry_deg4 = tk.Entry(root)
entry_deg4.grid(row=11, column=1)

# Glycan masses
tk.Label(root, text="A2G0F (Da):").grid(row=12, column=0, sticky="w")
entry_glyco_a2g0f = tk.Entry(root)
entry_glyco_a2g0f.grid(row=12, column=1)

tk.Label(root, text="A2G0 (Da):").grid(row=13, column=0, sticky="w")
entry_glyco_a2g0 = tk.Entry(root)
entry_glyco_a2g0.grid(row=13, column=1)

tk.Label(root, text="A2G1F (Da):").grid(row=14, column=0, sticky="w")
entry_glyco_a2g1f = tk.Entry(root)
entry_glyco_a2g1f.grid(row=14, column=1)

tk.Label(root, text="A2G2F (Da):").grid(row=15, column=0, sticky="w")
entry_glyco_a2g2f = tk.Entry(root)
entry_glyco_a2g2f.grid(row=15, column=1)

# Accuracy
tk.Label(root, text="Accuracy (ppm):").grid(row=16, column=0, sticky="w")
entry_accuracy_ppm = tk.Entry(root)
entry_accuracy_ppm.grid(row=16, column=1)

# Molecule type dropdown
tk.Label(root, text="Molecule Type:").grid(row=17, column=0, sticky="w")
molecule_type_var = tk.StringVar(value="Intact Protein")
molecule_type_menu = tk.OptionMenu(root, molecule_type_var, "Intact Protein", "Chain")
molecule_type_menu.grid(row=17, column=1, sticky="we")

# Number of cysteines
tk.Label(root, text="Number of Cysteines: (Exclude Cys involved in conjugation)").grid(row=18, column=0, sticky="w")
entry_num_cysteines = tk.Entry(root)
entry_num_cysteines.grid(row=18, column=1)

# Calculate button
button_calculate = tk.Button(root, text="Calculate Proteoforms", command=on_calculate)
button_calculate.grid(row=19, column=0, columnspan=2, sticky="we", pady=8)

root.mainloop()
