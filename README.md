# DeepPoW
[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Shiny](https://img.shields.io/badge/UI-Shiny-lightgrey)](https://shiny.posit.co/)

**DeepPoW: An interactive open-source R/Shiny platform for diagnosing and calibrating the time-dependent rate phenomenon (TDRP) in viral molecular evolution.**

Evolutionary rates estimated over short timescales are often higher than those inferred over longer timescales, giving rise to the **time-dependent rate phenomenon (TDRP)**. DeepPoW provides an integrated workflow for investigating this phenomenon by combining rapid rate-decay diagnostics with mechanistic **Prisoner of War (PoW) model** calibration across Bayesian MCMC posterior-tree ensembles.

## Key Features

### TDRP Diagnostics (Mode A)

Rapidly assess a time-scaled summary phylogenetic tree for evidence of time-dependent rate decay.

Mode A estimates:

- initial rate parameter (`μ0`);
- long-term substitution rate (`μ∞`); and
- decay coefficient(`k`).

It also provides model comparison, robustness assessment, and an intuitive **Traffic Light** classification to help determine whether further PoW calibration is warranted.

### Posterior Ensemble Calibration (Mode B)

Apply mechanistic PoW calibration across Bayesian posterior-tree ensembles.

Mode B is designed for large posterior datasets and uses efficient numerical inversion/interpolation and multi-core computation to transform posterior trees while preserving variation across posterior trees.

**Important:** Mode A and Mode B serve different analytical purposes. The `k` estimated from Mode A is **not directly passed to Mode B**.

### Structural Validation

Provide complementary evidence by comparing **PoW-calibrated evolutionary distances** with **protein structural divergence** when suitable structural data are available.

The workflow supports matrix-based comparisons and statistical assessment, including Mantel tests. Structural validation is intended as an independent, complementary assessment rather than direct proof of the PoW mechanism.

### Phylogenetic Utilities

DeepPoW also provides supporting tools for phylogenetic and structural workflows, including:

- **TreeTime TDRP Explorer**
- **RMSD Matrix Generator**
- **BEAST 2 XML Patcher**

---

## Workflow

A typical DeepPoW analysis follows this workflow:

```text
Bayesian phylogenetic inference
            │
            ▼
      Posterior trees
            │
            ▼
       Summary tree
            │
            ▼
 ┌──────────────────────┐
 │ Mode A: TDRP         │
 │ Diagnostics          │
 └──────────────────────┘
            │
            ▼
   Robustness assessment
            │
      ┌─────┴─────┐
      │           │
  supported    caution
      │           │
      ▼           ▼
   Mode B      further
 posterior    inspection
 calibration
      │
      ▼
 Calibrated posterior
      │
      ▼
 Structural validation
   (when applicable)
```
This diagram illustrates a recommended workflow; Mode B can also be used independently of Mode A when suitable posterior-tree data are available.

For Bayesian analyses, appropriate upstream quality control should be completed before posterior calibration, including burn-in assessment, MCMC convergence, and effective sample size evaluation.

---

## Installation

DeepPoW requires **R ≥ 4.2.0** and can be run on **macOS, Linux, or Windows**.

### Install Dependencies

Install the required R packages:

```r
cran_packages <- c(
  "shiny", "bslib", "plotly", "DT", "shinycssloaders",
  "xml2", "ape", "nleqslv", "ggplot2", "shinyFiles",
  "minpack.lm", "tidytree", "stringr", "dplyr", "tidyr", "vegan"
)

new_packages <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!("treeio" %in% installed.packages()[, "Package"])) {
  BiocManager::install("treeio")
}
```

### Launch DeepPoW

Clone or download this repository, then open the project in R/RStudio.

Launch the application with:

```r
setwd("path/to/DeepPoW")
shiny::runApp()
```

### External Software (For Structural Validation Only)

If you intend to use the **Structural Validation** module to calculate structural distances, DeepPoW relies on the following third-party tools. These must be installed independently and accessible in your system's PATH:

- **TM-align**: For protein structure alignment.

- **Foldseek**: For fast structural search and comparison.

*(Note: If you are only using Mode A or Mode B for rate diagnostics and calibration, these external tools are not required.)*

---

## Documentation

DeepPoW is accompanied by several complementary resources:

- **[Quick Guide](docs/DeepPoW_Quick_Guide_v1.0.pdf)** — a concise, step-by-step guide for running a typical DeepPoW analysis.
- **[User Manual](docs/DeepPoW_User_Manual_v1.0.pdf)** — detailed information on functions, parameters, interpretation, troubleshooting, and advanced usage.
- **[Research Article](#citation)** — theoretical background, computational framework, validation, and empirical applications.

> **Important:** DeepPoW is a diagnostic and calibration framework, not an automatic detector of biological truth. Results should be interpreted in the context of sampling design, evolutionary timescale, temporal signal, and phylogenetic uncertainty. See the **User Manual** for detailed limitations and interpretation guidance.

---

## Input Data

DeepPoW is primarily designed for viral phylogenetic analyses using time-scaled trees and Bayesian posterior-tree ensembles.

Depending on the analysis, typical inputs include:

- a time-scaled summary tree;
- Bayesian posterior trees;
- sampling dates or associated metadata; and
- protein structures or structural-distance matrices for structural validation.

Users should ensure that tree tip names, metadata, and posterior trees are correctly matched before analysis.

Detailed input requirements and preparation procedures are provided in the **Quick Guide** and **User Manual**.

---

## Repository Structure

```text
DeepPoW/
├── app.R                      # Main Shiny application
├── code/                      # Core R modules
│   ├── mod_treetime.R         # TreeTime-based TDRP exploration
│   ├── pow_modeling.R         # TDRP model fitting and robustness diagnostics
│   ├── pow_multitree.R        # PoW calibration of MCMC posterior-tree ensembles
│   ├── structure_validation.R # Structural validation and Mantel tests
│   ├── tree_parser.R          # Parsing annotated phylogenetic trees
│   ├── generate_rmsd_matrix.R # Protein structural distance matrix generation
│   └── Utilities.R            # Supporting utilities
├── demo/                      # Example datasets
├── docs/                      # User Manual and Quick Guide
├── README.md                  # Project overview and installation instructions
└── LICENSE                    # MIT License
```

The exact repository structure may change between releases.

---

## Interpretation and Limitations

DeepPoW is intended to support **diagnosis, robustness assessment, and mechanistic calibration** of TDRP patterns.

Users should be cautious when working with datasets characterized by:

- highly uneven temporal sampling;
- dense or strongly clustered sampling;
- weak temporal signal;
- substantial phylogenetic uncertainty;
- potential recombination;
- very short evolutionary timescales; or
- poorly converged Bayesian analyses.

A fitted PoW model does not by itself establish a biological mechanism. Similarly, a stable parameter estimate does not exclude all possible sampling or phylogenetic artifacts.

The **Traffic Light** classifications are practical robustness-screening criteria and should not be interpreted as universal statistical significance thresholds or biological laws.

A Red Light indicates a methodological calibration veto rather than a software-enforced lockout; exploratory Mode B analyses remain available at the user's discretion.

Recombination is not automatically resolved by DeepPoW and should be considered during upstream sequence and phylogenetic analyses when relevant.

---

## Reproducibility

For publication-quality analyses, we recommend retaining:

- original input sequences and phylogenetic trees;
- Bayesian posterior-tree files;
- sampling metadata;
- R and DeepPoW release or commit identifier;
- analysis parameters; and
- DeepPoW output files.

Recording the DeepPoW version used for each analysis is strongly recommended.

---

## Citation

If you use DeepPoW in your research, please cite:

> **DeepPoW: A rigorous and interactive computational framework for calibrating the time-dependent rate phenomenon in viral evolution.**

**Citation details:** *Zhang S, Shen J, Jin S, Du Z, Wu Z, Gao F. DeepPoW: A rigorous and interactive computational framework for calibrating the time-dependent rate phenomenon in viral evolution. 2026; MS in preparation.*

DeepPoW is based on the mechanistic PoW framework introduced by:

> Ghafari M, Simmonds P, Pybus OG, Katzourakis A. A mechanistic evolutionary model explains the time-dependent pattern of substitution rates in viruses. *Current Biology*. 2021;31:4689–4696.e5.

Please also cite the relevant upstream methods, databases, and software used in your analysis.

---

## Support & Issues

For questions, bug reports, feature requests, or suggestions, please open an issue on the **GitHub Issues** page.

When reporting a bug, please include:

- DeepPoW version;
- operating system;
- R version;
- relevant input-file information;
- complete error messages; and
- a minimal reproducible example, when possible.

---

## Contributing

We welcome contributions and feedback from the community! If you would like to contribute to DeepPoW:

- Check the **GitHub Issues** to see if your feature or bug is already being discussed.

- Fork the repository and create a new branch for your feature.

- Submit a Pull Request (PR) with a clear description of the changes and any relevant test cases.

---

## License

DeepPoW is released under the **MIT License**.

See the [LICENSE](LICENSE) file for details.

---

## Acknowledgements

DeepPoW builds on the mechanistic PoW framework developed by Ghafari and colleagues and incorporates methods and open-source software from the fields of viral evolution, Bayesian phylodynamics, phylogenetics, and structural bioinformatics.

We acknowledge the developers and maintainers of the open-source tools and databases that support these analyses.
