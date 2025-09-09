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
   - [Python Code for temperature prediction](#python-code)
   - [R Code for age prediction](#r-code)
   - [Python Code for extinction maps](#python-code-ext)
3. [Data](#datasets)
4. [Requirements](#requirements)

---

## Repository Structure
The repository is organized as follows:
```text
├── Usage/ # Scripts demonstrating the analysis
│   ├── PYTHON CODE TEMPERATURE PREDICTIONS.py # Python code for predicting the stars' temperature
│   ├── codes_repository_predict_age.R         # R code for predicting the stars' age
|   └── extinction_3dmap.py                    # Python code for interpolating the extinction maps
├── Data/ # Datasets used in the analysis
│   ├── ISO_SPOTS_ph_id_complete.RData
│   ├── jackson_members_filt_binarie_final7000.RData
│   ├── test_set_df_stars_GG2M_model.csv
│   └── train_set_df_stars_GG2M_model.csv
├── LICENSE                # License information
├── README.md              # Project overview and instructions
├── requirements_R.txt     # R dependencies
└── requirements_py.txt    # Python dependencies
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

### Python Code for temperature prediction
The Python script can be found in the `Usage/` folder: **PYTHON CODE TEMPERATURE PREDICTIONS.py**. These codes are necessary to obtain temperature predictions and their associated standard errors.
Contents:
- Import of the training and test set ("train_set_df_stars_GG2M_model.csv", "test_set_df_stars_GG2M_model.csv")
- Create the Neural Network model: 5 hidden layer with 128-64-32-16-8 neurons with relu activation functions and linear activation function for the output layer. Use of Adam optimizer and mae as loss function, 50 epochs and 8 as batch size. Use of k-fold cross validation with k = 10 to prevent from overfitting.
- Predict the effective temperature on the test set.
- Visualisation of the temperature predictions on the test set.
- Bootstrap procedure: Once the temperature predictions are obtained, it is necessary to derive the standard error associated with these predictions. To achieve this, we implemented a bootstrap procedure based on the neural network architecture described above. While the model assessment was carried out using 10-fold cross-validation, for the bootstrap procedure we adopted an early stopping criterion instead of cross-validation in order to reduce computational cost. Importantly, the architecture and hyperparameters of the network (5 layers, 50 epochs, batch size = 8, Adam optimizer, MAE loss) were kept fixed. Preliminary tests showed that the results obtained with early stopping were highly consistent with those from cross-validation, thus justifying the use of this more efficient procedure during the bootstrap phase.
- Visualisation of the predicted standard error values by bootstrap procedure.

--- 

### R Code for age prediction
The R script can be found in the `Usage/` folder: **codes_repository_predict_age.R**. These codes are necessary to predict the age of the stars and the associated standard errors. Furthermore, there are the codes to choose the best isochrone set for each cluster of the analysis, i.e. the starspot evolutionary model.
Contents:
- Import the two main dataset for the analysis: ISO_SPOTS_ph_id_complete.RData, which includes all the information related to the different isochrone sets, and jackson_members_filt_binarie_final7000.RData, which is the main dataset with the stars of interest for the analysis, which also includes the predicted temperature values.
- Definition of different functions in R to perform the interpolation of the position of the stars in the H-R diagram with respect to the position of the isochrones of the different sets.
- Prediction of the stars' age by applying the previous functions and summary of the results.
- Selection of the beta spot parameter from starspot evolutionary models, i.e. the different isochrone sets.
- Prediction of the standard errors associated to the age predictions. The estimation of the error spread was carried out using a Monte Carlo simulation approach. For each star, we randomly generated 100 values for Gaia’s 
G filter magnitude and 100 values for effective temperature, creating 100 fictitious realizations per star. These realizations account for the variability introduced by measurement errors.

---

### Python Code for extinction maps
The Python script is located in the root of the repository: **extinction_3dmap.py**.  
This code computes the line-of-sight extinction for a list of stars by applying the **3D extinction maps of Vergely et al. (2022)**. It also estimates extinction in several photometric bands commonly used in stellar astrophysics.  

#### Contents  

- **Input data**  
  The script requires a CSV file containing at least:  
  - `l`: galactic longitude (degrees)  
  - `b`: galactic latitude (degrees)  
  - `dist_par`: stellar distance in parsecs (e.g., from Gaia parallaxes or distance estimates).  

- **Extinction calculation**  
  - Galactic coordinates are converted into Cartesian coordinates.  
  - The line of sight is sampled in steps of 10 pc, from the Sun up to the star’s distance.  
  - The dust density is integrated along the line of sight using the **Vergely+22 3D dust maps** (three different FITS cubes at increasing distance ranges).  
  - The result of the integration is stored as `A0`, which corresponds to \(A_V\) in Vergely+22.  
  - The reddening value `E(B-V)` is also computed.  

- **Extinction in other bands**  
  - Using the **Gaia EDR3 extinction law**, the code derives extinction values in multiple filters:  
    - Gaia (G, BP, RP)  
    - 2MASS (J, H, Ks)  
    - Pan-STARRS1 (g, r, i, z, y, w)  
  - Gaia magnitudes are corrected for absorption to recover intrinsic colors, which are then used to refine the extinction estimates.  

- **Output**  
  - A CSV file (`extinction_3dmap.csv`) containing the original input data plus:  
    - Line-of-sight extinction `A0`  
    - Color excess `E(B-V)`  
    - Extinction values in Gaia, 2MASS, and Pan-STARRS1 bands  

- **Diagnostics**  
  - The script computes residuals between iterations (mean, RMS, median, MAD) to verify convergence.  
  - Since differences are generally small, a single iteration is sufficient for most cases.  

---

## Data

This folder provides the main datasets used in the analysis:

### ISO_SPOTS_ph_id_complete.RData
This table contains all the information related to the isochrones used in this study. It includes the following columns:
- **beta**: the beta spot value associated with each isochrone to identify the respective starspot evolutionary model.
- **ML**: mixing length parameter.
- **age_yr**: age identified by the isochrone.
- **G**: theoretical value of the G magnitude in the H-R diagram of the isochrone.
- **log_Te**: theoretical temperature in the H-R diagram of the isochrone (log10 scale).
- **phs**: phase of the isochrone.
- **isocrona**: ID of the isochrone.

### jackson_members_filt_binarie_final7000.RData
This is the main dataset with the stars of interest for the analysis. It includes the following columns:
- **ges_id_gaia**: Gaia ID of the star.
- **MG0_ML**: absolute G magnitude obtained from the predicted temperatures.
- **logTeff**: predicted temperature of the star by the Neural Network approach (log10 scale).
- **CLUSTER**: cluster of the star.

### train_set_df_stars_GG2M_model.csv and test_set_df_stars_GG2M_model.csv
These datasets include all the selected features to perform temperature prediction via the Neural Network approach. Columns include:
- **ges_id_gaia**: Gaia ID of the star.
- **TEFF**: spectroscopic temperature value obtained from the GES catalogue.
- **E_TEFF**: error associated with the spectroscopic effective temperature of GES.
- **MBP_corr_redd**: absolute magnitude in the Gaia *BP* filter, corrected for extinction.
- **mk_corr_redd**: apparent magnitude in the 2MASS *K* filter, corrected for extinction.
- **g_rp_corr_redd**: Gaia *G-Grp* color, corrected for extinction.
- **bp_rp_corr_redd**: Gaia *Gbp-Grp* color, corrected for extinction.
- **j_h_corr_redd**: 2MASS *J-H* color, corrected for extinction.
- **h_k_corr_redd**: 2MASS *H-K* color, corrected for extinction.
- **Teff_mont**: input photometric temperature.

**N.B.** There are additional datasets derived from the training and test sets, but they are too large to be uploaded on GitHub. They can be shared by the authors via Drive.

---

## Requirements
This project requires both Python and R to run. All necessary packages are listed in the corresponding requirements files.





