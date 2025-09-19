# About
This repository provides an official implementation of an algorithm for Activity cliff removal and a a set of machine learning and deep-learning models  develop over CHEMBL database considering  for classification:.

## Environment Setup
To set up the Python environment, use the provided `environment.yaml` file. Run the following command in Conda:
```
conda env create -f environment.yaml
```

## Removing Activity Cliffs
To remove activity cliffs from your dataset, use the `test.py` script. This script requires the following arguments:
1. A threshold value between 0.7 to 0.9.
2. The filename of the original dataset (".csv").
3. The filename for the output dataset.

### Input Data Requirements
The input dataset should contain:
- Active and inactive interactions (`CLASS`).
- Compound identifiers (`MOL_REGNO`).
- Target identifiers (`TARGET_ID`).

### Example Usage
Run the script with the following command:
```
python dd 0.8 data_with_ACs.csv data_without_ACs.csv
```

## Models
This repository also includes Target-Centric  (TCM) and Deep Learning models (DL) stored in the `MODELS` folder. These models were trained using ChEMBL release 34 and consider two bioactivity cutoffs: 10 µM and 1 µM. A web tool also expose an interface to make the predictions based on an smile input at "bioquimio.udla.edu.ec" version 2. 

For more information, please check our article at: 
