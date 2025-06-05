# DBS fMRI Reproducibility Repository

This repository contains de-identified example data and resources for assessing the reproducibility of resting-state fMRI (rs-fMRI) scans in Parkinson's patients who are receiving deep brain stimulation (DBS) therapy. It includes code and patient-derived, de-identified data for analysis and processing. The repository is part of ongoing research in the Radiology-Morrison-lab-UCSF.

## Contents 

- **`testretest_ICC.m`**:
   - Code to calculate intraclass correlation coefficient (ICC) values between test and retest resting-state fMRI scans using functional connectivity and brain variability data, as well as corresponding statistical analyses in the lab's paper, "Test-retest reliability of resting-state functional magnetic resonance imaging during deep brain stimulation for Parkinson’s disease: an open-source dataset". This script works in conjunction with the ICC function (Salarian 2024): https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc .
- **`testretest_ICC_affectofmetalartifact.m`**:
   - Code to create a map of DBS electrode spatial trajectories and identify intersecting brain regions using T1-weighted MRI scans and rs-fMRI scans. Once intersecting regions are identified, the code calculates ICC values between test and retest resting-state fMRI scans of two subsets of brain regions: those affected by the electrode and those unaffected. Finally, the code conducts corresponding analyses in the lab's paper, "Test-retest reliability of resting-state functional magnetic resonance imaging during deep brain stimulation for Parkinson’s disease: an open-source dataset". 
- **`testretest_ICC_betweensessions.m`**:
   - Code to calculate ICC values between longitudinal test and retest resting-state fMRI scans, as well as conduct corresponding analyses from the lab's paper, "Test-retest reliability of resting-state functional magnetic resonance imaging during deep brain stimulation for Parkinson’s disease: an open-source dataset".
- **`testretest_makeconnectomes.m`**:
   - Code to calculate functional connectivity and brain variability values, as well as conduct corresponding analyses from the lab's paper, "Test-retest reliability of resting-state functional magnetic resonance imaging during deep brain stimulation for Parkinson’s disease: an open-source dataset".

- **`AdvMRIforDBS_InternalProtocol.pdf`**:
  - Advanced MRI for DBS: Internal UCSF Protocol for 3T Brain MRI with DBS implant.
- **`MDS-UPDRS_Scores_Template.xlsx`**:
   - ​An Excel template demonstrating MDS-UPDRS score calculations with dummy data. ​
- **`testretest_ROIs.xlsx`**:
   - A list of ROIs used for the whole brain, motor, limbic, and associative networks.
- **`testretest_demographics.xlsx`**:
    - An Excel sheet containing cross-sectional patient data, including demographics, DBS specifications, Movement Disorder Society-Unified Parkinson's Disease Rating Scale (MDS-UPDRS) sub-scores, and Levodopa Equivalent Daily Dose (LEDD) scores.-
- **`FC_BV_Data`**:
    - A folder containing all of our functional connectivity (FC) and brain variability (BV) data as Excel files.
- **`atlases`**:
  - A folder containing atlas files and their adjoining ROI text files for region identification.
- **`otsu-segmentation`**:
  - A folder containing code to run Otsu's method for image thresholding.



- **Data**:
  - All de-identified MRI scans and associated metadata will be made publicly available through [OpenNeuro](https://openneuro.org/) upon completion of ongoing analyses.

## License
- Copyright 2025 UCSF
- Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
- **`FC_BV_Data`** and **`testretest_demographics.xlsx`** were acquired at the University of California San Francisco and are licensed under the CC BY-NC 4.0 license

## Acknowledgments
This project is part of ongoing research in the Radiology-Morrison-lab-UCSF. Special thanks to all contributors and collaborators involved in this work.
