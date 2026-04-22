# BioDRUMs

**Getting started**

BioDRUMs (Biologics Drug Ratio & Unified intact Mass analysis) is a Python pipeline to streamline mean Drug‑to‑Antibody Ratio (DAR) and biologics structural integrity (SI) analysis from high‑resolution mass spectrometry (HRMS) data
<img width="900" height="600" alt="Figure 1" src="https://github.com/user-attachments/assets/3f820de9-8756-4976-ac62-582c1bb6d4bc" />

**Version: 1.0.0** **Language: Python 3** **GUI Interface: tkinter GUI**

**Overview**

BioDRUMs (Biologics Drug Ratio & Unified intact Mass analysis) is a Python pipeline with a graphical user interface (GUI) for automating:

**Mean drug‑to‑antibody ratio (DAR)** analysis of antibody-drug conjugates (ADCs) and other bioconjugates

**Structural integrity (SI)** analysis of unconjugated antibodies, bispecifics, and multispecific biologics

**Proteoform hypothesis generation** and assignment after intact MS deconvolution

**Quality control (QC)** on intact MS data (explained area %, isobaric species count)

It is designed to work on intact protein data from high‑resolution MS instruments and to standardize data processing across in vitro and in vivo studies (e.g., plasma stability, pharmacokinetics), including low‑volume small‑animal studies where replicate injections are often not feasible. It also works quite well for ADCs/biologics characterization in buffer, where high amount of sample is available.

If you use BioDRUMs in your work, please cite:

**Di Ianni et al., “BioDRUMs: a Python‑based pipeline to streamline mean drug‑to‑antibody ratio and biologics structural integrity analysis from high‑resolution mass spectrometry analyzers”, journal details to be added.**

**Features**

BioDRUMs provides:

**Automated proteoform hypothesis generation**

* Combines antibody backbone masses, theoretical DAR ladder, linker–payload and degradation products, and (optional) glycoforms

* Supports both glycosylated and deglycosylated intact proteins (glycan multiplier set to 2) or subunits(glycan multiplier set to 1).

**Unified DAR and SI analysis**

* Mean DAR over time or concentration for ADCs/bioconjugates
  
* Mean SI values for unconjugated biologics, where intact full‑domain species are scored SI = 1 and non‑intact species SI = 0

**Quality‑control metrics for intact MS experiments**

**Explained area %**: percentage of deconvoluted MS peaks that can be assigned to modeled species (typical threshold ≥ 75%)

**Isobaric species count** per DAR level or proteoform, to indicate assignment ambiguity (this is done on the whole set of generated DAR species, but you should look into your identified DAR species only)

**Flexible input compatibility**

Accepts table‑formatted deconvolution output from commercial and open‑access software (e.g., BioPharma Finder ReSpect™, UniDec, others), as long as the table includes:
a species mass / ID
**DAR value** for each species (when applicable)

**protein name**

**sum intensity**

Data filtering based on user‑tunable thresholds (default settings= 30 ppm mass accuracy, 4% relative abundance vs base peak)

**User‑friendly GUI** (No prior coding experience required)

Dedicated windows for:

* Proteoform mass calculation and hypothesis generation

* DAR / SI analysis, QC metrics, and plotting

**Installation**

**Requirements**

* Python ≥ 3.x

* Standard scientific Python stack (e.g., pandas, numpy, matplotlib)

* GUI library (e.g., tkinter in requirements.txt)

BioDRUMs operates on already deconvoluted intact MS data. You will also need:

* A deconvolution software (e.g., BioPharma Finder 4.0 with ReSpect™, or other tools that export tabular deconvolution results in .csv or .xlsx formats)
The corresponding **Excel template** for BioDRUMs input (see below: Input templates)

**Input templates**

BioDRUMs operates on a template Excel file that contains the deconvoluted species and metadata required for the analysis (see Tables S1–S3 in the paper).


Typical columns include (minimal example):

 * **Protein_Name** – identifier for the ADC or biologic

* **DAR_Value** – DAR of each individual species (or SI descriptor in unconjugated cases)

 * **Sum_Intensity** – deconvoluted peak intensity (used as absolute abundance measure)

 * **Mass** – experimental deconvoluted mass

* Additional metadata (timepoint, concentration, matrix type, study ID, etc.)

You can either:

Export a deconvolution table and adapt it manually following the template, or use a helper script (if provided) that converts vendor output into the BioDRUMs template format.

The repository includes:

* templates/Table S1_rev.xlsx/Table S2/Table S3

* examples/ – example input and corresponding outputs for the use cases in the paper (in vitro stability and PK studies)

**Modules & Workflow**

**1. Deconvoluted_mass_finder.py**

This module (with a GUI) is used to generate proteoform structure hypotheses and theoretical intact masses before DAR/SI analysis.


**Capabilities:**

**Correct the protein (or subunit) mass for:**

* Paired cysteines involved in disulfide bonds or conjugation (in case of cysteine-conjugated proteins, do not include this in the paired cysteines count)

* Engineered cysteines or hinge mutations (e.g., IgG4 stabilization, site‑specific conjugation)
  
Define:

**Theoretical DAR ladder (from DAR0 up to a specified maximum)**

* Linker–payload mass shifts (e.g., valine‑citrulline PABC MMAE; supports transglutaminase and interchain/engineered cysteine conjugation, but also click-chemistry or other conjugation types according to user's needs)

* Potential linker–payload degradation products (for in vivo experiments)

Include:

* Main glycosylated proteoforms (intact or at chain level), when deglycosylation is not possible ot residual glycoforms are still present in the sample after deglycosylation

* Deglycosylated forms (e.g., after PNGase‑F treatment)

**Output:**

A comprehensive Excel list of all predicted proteoforms, their DAR values, and the number of isobaric species per DAR level within the chosen mass tolerance

**Typical GUI steps**

Specify:

* Protein sequence / backbone mass

* Number and position of cysteines participating in disulfide bonds or conjugation

* Maximum DAR and linker–payload mass + degradation products

* Glycan masses / repertoire (if working with glycosylated samples)

Set:

* Mass accuracy tolerance (e.g., ±30 ppm, but it can be modified by the user)
* Relative abundance tolerance (e.g., set default as 4% of the base peak in MS)

Run hypothesis generation to produce:

* Proteoform list with theoretical masses

* Isobaric species statistics

* Exported Excel file for subsequent matching with experimental deconvoluted masses

**2. unified_DAR_analysis.py**

This module performs the core DAR/SI calculations, QC metrics, and plotting, operating on the template Excel file that links experimental deconvoluted peaks to predicted proteoforms.

Functions:

**Mean DAR / SI calculation**

Uses Equation 1 in the paper to compute mean DAR across all identified species at a given timepoint or concentration (intensity‑weighted average over DAR values). For SI, analogous calculations are performed using SI = 1 for intact full‑length species and SI = 0 for non‑intact species.


**Explained area percentage (%)**

Calculates the percentage of total deconvoluted peak area that is assigned to modeled species.
Typical acceptance threshold: ≥ 75% for good characterization (may drop near Limit Of Detection, especially in in vivo samples).

**DAR / proteoform species characterization**

For each timepoint/concentration, plots the relative intensity % of each identified species versus the total identified signal.
Assigns a unique Species ID to each proteoform and writes an Excel mapping (Species_ID_Mapping sheet).

**GUI workflow**

The GUI allows you to:

Select the input metadata Excel file (compiled template).

Set:

* Molecule ID

* Study type (e.g., in vitro stability, PK)

* Matrix (e.g. plasma)

* X‑axis variable (e.g., time, concentration)

* Indicate whether to compute DAR (conjugates) or SI (unconjugated biologics)

**Output files**

Run the analysis to generate:

* Plots of mean DAR or SI versus time/concentration (.png)

* Explained area % versus time/concentration (.png)

* Stacked plots of relative abundance of each species versus time/concentration (.png)

* Excel reports (DAR summary sheet from your template file), including Species ID mapping (.xlsx)

**Example use cases**

**Use case 1 – Mean DAR of a Trastuzumab ADC in an in vitro plasma stability study**

**Molecule:** Trastuzumab conjugated via a maleimide valine‑citrulline PABC MMAE linker on an engineered cysteine in the CH3 domain (maximum DAR 2).

**Samples:** Mouse plasma stability 0–168 h, deglycosylated with PNGase‑F to reduce isobaric ambiguity.

**Results (see paper Figures 4–5):**

* Mean DAR remains within 20% loss compared to t0, indicating good linker–payload stability.
* Explained area % stays above 75% for all timepoints.
* 21 ADC species identified and tracked over time.

This workflow is fully reproducible in BioDRUMs with the provided example dataset and template
