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
The Python script can be found in the `Usage/` folder: **PYTHON CODE TEMPERATURE PREDICTIONS.py**. These codes are necessary to obtain temperature predictions and their associated standard errors.
Contents:
- Import of the training and test set ("train_set_df_stars_GG2M_model.csv", "test_set_df_stars_GG2M_model.csv")
- Create the Neural Network model: 5 hidden layer with 128-64-32-16-8 neurons with relu activation functions and linear activation function for the output layer. Use of Adam optimizer and mae as loss function, 50 epochs and 8 as batch size. Use of k-fold cross validation with k = 10 to prevent from overfitting.
- Predict the effective temperature on the test set.
- Visualisation of the temperature predictions on the test set.
- Bootstrap procedure: Once the temperature predictions are obtained, it is necessary to derive the standard error associated with these predictions. To achieve this, we implemented a bootstrap procedure based on the neural network architecture described above. While the model assessment was carried out using 10-fold cross-validation, for the bootstrap procedure we adopted an early stopping criterion instead of cross-validation in order to reduce computational cost. Importantly, the architecture and hyperparameters of the network (5 layers, 50 epochs, batch size = 8, Adam optimizer, MAE loss) were kept fixed. Preliminary tests showed that the results obtained with early stopping were highly consistent with those from cross-validation, thus justifying the use of this more efficient procedure during the bootstrap phase.
- Visualisation of the predicted standard error values by bootstrap procedure.

--- 

### R Code
The R script can be found in the `Usage/` folder: **codes_repository_predict_age.R**. These codes are necessary to predict the age of the stars and the associated standard errors. Furthermore, there are the codes to choose the best isochrone set for each cluster of the analysis, i.e. the starspot evolutionary model.
Contents:
- Import the two main dataset for the analysis: ISO_SPOTS_ph_id_complete.RData, which includes all the information related to the different isochrone sets, and jackson_members_filt_binarie_final7000.RData, which is the main dataset with the stars of interest for the analysis, which also includes the predicted temperature values.
- Definition of different functions in R to perform the interpolation of the position of the stars in the H-R diagram with respect to the position of the isochrones of the different sets.
- Prediction of the stars' age by applying the previous functions and summary of the results.
- Selection of the beta spot parameter from starspot evolutionary models, i.e. the different isochrone sets.
- Prediction of the standard errors associated to the age predictions. The estimation of the error spread was carried out using a Monte Carlo simulation approach. For each star, we randomly generated 100 values for Gaia’s 
G filter magnitude and 100 values for effective temperature, creating 100 fictitious realizations per star. These realizations account for the variability introduced by measurement errors.

