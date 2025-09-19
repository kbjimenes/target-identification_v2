# About

This repository provides the  implementation of an algorithm for Activity Cliff removal along with a collection of machine learning (ML) and deep learning (DL) models developed using the ChEMBL database. The models are designed for target identification and trained under different bioactivity thresholds of 10uM and 1uM.

## Environment Setup
To set up the Python environment, use the provided `environment.yaml` file. Run the following command in Conda:
```
conda env create -f environment.yaml
```

## Removing Activity Cliffs

The repository includes a script to automatically remove activity cliffs from datasets:

```
python remove_activity_cliff.py <threshold> <input_file> <output_file>

```

- Activity Cliff removal threshold value between 0.7–0.9.

- Filename of the original dataset (.txt).

- Filename for the output dataset (processed without activity cliffs).

### Input Data Requirements

The input dataset must include the following columns:

-TARGET_ID → target identifier.
-COMPOUND → SMILE compound .
-CLASS → activity label (active/inactive).

### Example Usage

```
python remove_activity_cliff.py 0.8 data_with_ACs.txt data_without_ACs.txt
```

## DATA

The DATA folder includes raw compound-target interaction datasets prepared with two different bioactivity thresholds(1 µM, 10 µM ) and three different cleanning strategies (MAD, only-activity, and large-only activity).

## Models

This repository also includes Target-Centric  (TCM) and Deep Learning models (DL) stored in the `MODELS` folder. These optimized models were trained on curated datasets derived from ChEMBL release 34, following activity cliff removal, and were evaluated under two bioactivity thresholds: 10 µM and 1 µM.

Additionally, a web-based tool is available at: "http://bioquimio.udla.edu.ec/tidentification02". 

## Citation

For more information, please check our article at:  Soon available. 