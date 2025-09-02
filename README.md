# Temperature-and-Age-predictions-of-young-stars
This repository contains the code and datasets used in the article:  

**Using a Neural Network approach and Starspots dependent models to predict Effective Temperatures and Ages of young stars**  
*Marco Tarantino, Loredana Prisinzano, Nicoletta D'Angelo, Francesco Damiani, Giada Adelfio*  

---
## Description  

This project develops a **statistical machine learning approach** to predict the **effective temperatures** of pre-main sequence (PMS) stars and to estimate their **stellar ages**, by combining:  
- A **Neural Network model** trained on high-quality spectroscopic temperatures from the **Gaia-ESO Survey**.  
- **Photometric data** from **Gaia DR3** and **2MASS** catalogs.  
- **Isochrone interpolation** with **starspot-dependent evolutionary models**.  

Main contributions:  
- Accurate prediction of stellar effective temperatures for large stellar populations with only photometric data.  
- Construction of Hertzsprung–Russell diagrams and stellar age predictions.  
- Improved age determinations using starspot evolutionary models.  
- Evidence of intrinsic age spreads in young clusters, suggesting multiple star formation events.  

---

## Table of Contents
1. [Repository Structure](#repository-structure)
2. [Usage](#usage)
   - [Python Code](#python-code)
   - [R Code](#r-code)
3. [Data](#datasets)
4. [License](#license)
5. [References](#references)

---

## Repository Structure
The repository is organized as follows:
```text
├── Usage/ # Scripts demonstrating the analysis
│   ├── PYTHON CODE TEMPERATURE PREDICTIONS.py # Python code for analysis
│   └── codes_repository_predict_age.R        # R code for analysis
├── Data/ # Datasets used in the analysis
│   ├── ISO_SPOTS_ph_id_complete.RData
│   ├── jackson_members_filt_binarie_final7000.RData
│   ├── test_set_df_stars_GG2M_model.csv
│   └── train_set_df_stars_GG2M_model.csv
├── LICENSE                # License information
├── README.md              # Project overview and instructions
├── requirements_R.txt     # R dependencies
└── requirements_py        # Python dependencies
```

**Description of folders and files:**
- **Usage/** – Contains Python and R scripts for performing the analysis.  
- **Data/** – Includes all datasets required for the analysis.  
- **LICENSE** – Specifies the terms under which the project can be used.  
- **README.md** – Provides an overview of the project, instructions, and documentation.  
- **requirements.txt** – Lists the Python packages needed to run the Python scripts.  
- **requirements.R** – Lists the R packages needed to run the R scripts.

---

## Usage
This section provides instructions on how to run the analysis scripts included in this repository. The examples are organized by programming language.

### Python Code
The Python script can be found in the `Usage/` folder: **PYTHON CODE TEMPERATURE PREDICTIONS.py**


