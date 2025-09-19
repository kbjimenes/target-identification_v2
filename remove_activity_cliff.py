#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script Name: remove_activity_cliff.py
Author: kbjimenes
Date: 01-11-2023
Description: This algorithm presents an activity cliff (ACs) removal algorithm for a compound-target interaction dataset. 
"""
import sys
import pandas as pd
import time
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import RDKFingerprint


RADIUS_FP = 8 

def validate_threshold(th):
    """
    Validates that the input value th is within the allowed range (0.1 to 1.0 in increments of 0.1).
    
    Parameters:
        th (float): The threshold value to validate.
    
    Returns:
        bool: True if valid, False otherwise.
    """
    valid_values = {round(i * 0.1, 1) for i in range(1, 10)}  # Set of valid values from 0.1 to 1.0
    
    if th in valid_values:
        return True
    else:
        print(f"Invalid value!. Please provide a value between 0.1 and 0.9 in increments of 0.1.")
        return False


def createTanimotoMatrix(pdDT, lst_pos, lst_neg):
    """
    Create tanimoto matrix

    Parameters:
    ------------------------------------------------------------------------
    pdDT : Input cleaned DataFrame containing at least the following columns:
                        - 'TARGET_ID': Identifier for the biological target.
                        - 'MOL_REGNO': Unique identifier for each molecule.
                        - 'CLASS': Classification label (e.g., 0 or 1).
    lst_pos : list of  active compounds in the dataset
    lst_neg : list of inactive compounds in the dataset   

    Returns :
    ------------------------------------------------------------------------
    pdMaTC : Computed tanimoto matrix
    """
    print("   - Computing tanimoto matrix - ")

    # Precompute fingerprints for all IDs in lst_pos and lst_neg
    id_to_fp = {}
    # Combine IDs to avoid duplicate fingerprint computations
    all_ids = set(lst_pos) | set(lst_neg)
    id_to_fp = {
        id_: RDKFingerprint(Chem.MolFromSmiles(pdDT.loc[pdDT.iloc[:, 3] == id_, pdDT.columns[1]].values[0]))
        for id_ in all_ids
    }

    Npos, Nneg = len(lst_pos), len(lst_neg)
    dist_ = np.zeros((Npos, Nneg))

    # Compute Tanimoto similarities
    for i, i_id in enumerate(lst_pos):
        fp1 = id_to_fp[i_id]
        for j, j_id in enumerate(lst_neg):
            fp2 = id_to_fp[j_id]
            dist_[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)

    pdMaTC = pd.DataFrame(dist_, index=lst_pos, columns=lst_neg)
    return pdMaTC

def recompute_cliff_positive_negative(dtPos): 
    """
    Compute and append the count of non-NA values for each row and its transpose. It refers to the activity cliff

    Parameters:
    ------------------------------------------------------------------------
        dtPos: The input DataFrame containing the original data.

    Returns:
    ------------------------------------------------------------------------
    - dtPos : The original DataFrame with an additional 'CLIFF' column representing the count of ACs
    - dtNeg : The transposed DataFrame with an additional 'CLIFF' column representing the count of ACs
    """
    #dtPos=pdData
    dtNeg = dtPos.transpose() 
    dtPos['CLIFF'] = dtPos.apply(lambda x: x.count(), axis=1)

    # return the transpse 
    dtNeg['CLIFF'] = dtNeg.apply(lambda x: x.count(), axis=1)
    return dtPos, dtNeg

def remove_column_CLIFF(pdData):
    """
    Remove the 'CLIFF' column if exist from a pandas DataFrame if it exists. 

    Parameters:
    ------------------------------------------------------------------------
        pdData: The DataFrame from which the 'CLIFF' column should be removed.

    Returns:
    ------------------------------------------------------------------------
    The DataFrame after removing the 'CLIFF' column if it existed; otherwise, the original DataFrame.
    """
    if 'CLIFF' in list(pdData.columns):
        pdData = pdData.drop('CLIFF', axis=1)
        return pdData
    return None


def find_compund_max_cliff(matPos, matNeg, tmp):
    """
    Identify the compound with the maximum 'CLIFF' value from an active/inactive matrices with the ACs description counts.

    Parameters:
    ------------------------------------------------------------------------
        matPos: DataFrame representing the active matrix with a 'CLIFF' column.
        matNeg: DataFrame representing the inactive matrix with a 'CLIFF' column.
        tmp : A control variable used to manage the iteration process in the calling function.

    Returns:
    ------------------------------------------------------------------------
        - compoundID: The identifier of the selected compound.
        - compoundType : The type of the compound, either 'P' for positive or 'N' for negative.
        - cliffValue (int or None): The 'CLIFF' value of the selected compound.
        - tmp (int): Updated control variable, typically returned unchanged.
    """
    print("     > max active Acs: ", np.max(matPos['CLIFF']), " & max inactive  ACs: ", np.max(matNeg['CLIFF']))
    if np.max(matPos['CLIFF'])==0 and np.max(matNeg['CLIFF'])==0:
        return None, None, None, 0

    # Exploring active compounds
    max_value = matPos['CLIFF'].max()
    pos_indices = matPos.index[matPos['CLIFF'] == max_value].tolist()
    id_compound_P = np.random.choice(pos_indices)

    # Exploring inactive compounds
    max_value = matNeg['CLIFF'].max()
    neg_indices = matNeg.index[matNeg['CLIFF'] == max_value].tolist()
    id_compound_N = np.random.choice(neg_indices)

    if(np.max(matPos['CLIFF'])==np.max(matNeg['CLIFF'])):
        #id_compound_P = np.random.choice(pos_indices); id_compound_N = matNeg['CLIFF'][random_neg]
        # print("+", id_compound_P, "-", id_compound_N)
        sel = np.random.choice([0,1]) 
        pX=matPos['CLIFF'].shape[0]; nX=matNeg['CLIFF'].shape[0]; 
            
        if pX>nX or sel==0:
            #print("remove from positive IDs")
            return id_compound_P, "P", matPos['CLIFF'].max(), tmp
        if nX>pX or sel==1:
            #print("remove from positive IDs")
            return id_compound_N, "N", matNeg['CLIFF'].max(), tmp
    elif(np.max(matPos['CLIFF'])>np.max(matNeg['CLIFF'])):
        #print("cliff with positive interaction")
        id_compound = id_compound_P #matPos['CLIFF'][random_pos]
        return id_compound, "P", matPos['CLIFF'].max(), tmp
    else:
        #print("cliff with negative interaction")
        id_compound = id_compound_N #matNeg['CLIFF'][random_neg]
        return id_compound, "N", matNeg['CLIFF'].max(), tmp

def remove_compound_in_matrix(pdData, cID, cType):
    """
    Remove a specific compound from a DataFrame based on its type active/inactive.

    Parameters:
    ------------------------------------------------------------------------
        pdData : The DataFrame from which the compound will be removed.
        cID : The label of the row or column to be removed.
        cType : The type of compound, determining whether to remove a row or a column.'N' refer to inactive and  'P' refer to active
    Returns:
    ------------------------------------------------------------------------
        The DataFrame after the specified row or column has been removed.
    """
    if cType=='N':
        pdData = pdData.drop(cID, axis=1)
    if cType=='P':
        pdData = pdData.drop(cID)
    return pdData 

def get_matrix_cliff(SM, lst_pos, lst_neg, th_ac):
    """
    Constructs a binary matrix indicating significant activity cliffs between
    positive and negative compound sets based on a similarity threshold.

    Parameters:
    ------------------------------------------------------------------------
    SM: Similarity tanimoto matrix containing structural similarity scores between compounds.
    lst_pos : List of compound identifiers with active interactions.
    lst_neg : List of compound identifiers  with inactive interactions.
    th_ac : Threshold value for activity cliff removal

    Returns:
    ------------------------------------------------------------------------
    dist_cliff: A binary DataFrame where '1' indicates a significant similarity based on th_ac threshold
    """

    # Asegúrate de que SM es un array de numpy
    matSM = np.array(SM)
    
    Npos = len(lst_pos);  Nneg = len(lst_neg)
    
    # Inicializar dist_cliff con ceros
    dist_cliff = np.zeros((Npos, Nneg))
    
    # Actualizar dist_cliff donde SM > th_cliff
    dist_cliff[matSM > th_ac] = 1
    dist_cliff[matSM <= th_ac] = None
    return pd.DataFrame(dist_cliff)
    


def check_byTH(SM, lst_pos, lst_neg, id_target, th_ac):
    """
    Filters compounds based on a similarity threshold, identifying and removing Activity Cliff

    Parameters:
    ------------------------------------------------------------------------
    SM : Similarity tanimoto matrix containing structural similarity scores between compounds.
    lst_pos: List of compound identifiers with active interactions.
    lst_neg: List of compound identifiers  with inactive interactions.
    id_target (int): Identifier for the target under investigation.
    th_ac:  Threshold value for activity cliff removal

    Returns:
    ------------------------------------------------------------------------
    list: List of selected compound identifiers after filtering.
    """

    thMat = get_matrix_cliff(SM, lst_pos, lst_neg, th_ac)
    # Asignar nombres a las filas y columnas del DataFrame
    thMat.index = lst_pos
    thMat.columns = lst_neg
    
    i=0
    tmp = 1
    while tmp !=0:
        [dtPos, dtNeg]  = recompute_cliff_positive_negative(thMat)
        [compoundID, compoundType, cliffValue, tmp]=find_compund_max_cliff(dtPos, dtNeg, tmp)
        dtPos = remove_column_CLIFF(dtPos)
        if compoundID != None:
            print(f"      in this step {i}, remove -- {compoundID} {compoundType} compound from a set of {cliffValue} ACs")
            #print(dtPos.shape)
            thMat=remove_compound_in_matrix(dtPos,compoundID, compoundType)
            #print(thMat.shape)
        if thMat.shape[0]<1 or thMat.empty:
            tmp=0       
            print("empty matrix ...")
        i=i+1

    lst_pos=[]; lst_neg = []; list_selected_ids = []
    if thMat.shape[1]>0:
        lst_pos= list(thMat.index)
        lst_neg = list(thMat.columns)
        list_selected_ids =lst_pos+lst_neg
        if "CLIFF" in list_selected_ids:
            list_selected_ids.remove("CLIFF")
    return list_selected_ids



def clean_dataset_cliff(pdData, tID, th_code):
    """
    Removing activity cliff in dataset based in a threshold.
    Parameters:
    ------------------------------------------------------------------------
    pdData : Input DataFrame containing at least the following columns:
                        - 'TARGET_ID': Identifier for the biological target.
                        - 'MOL_REGNO': Unique identifier for each molecule.
                        - 'CLASS': Classification label (e.g., 0 or 1).
    th      : Threhsold for removing activity cliff
    Assumes pdData has columns in order:
    - [0] Targets
    - [1] Compound smile
    - [2] class
    Returns :
    ------------------------------------------------------------------------
    None
    """
    try:

        # Group by TARGET_ID (col 0) and MOL_REGNO (col 1), count unique CLASS (col 2)
        pdA_filter = (
            pdData.groupby([pdData.columns[0], pdData.columns[1]], as_index=False)
                  .agg({pdData.columns[2]: 'nunique'})
                  .sort_values(by=pdData.columns[2], ascending=False)
        )

        # Keep only rows with more than one unique CLASS
        pdA_filter = pdA_filter[pdA_filter.iloc[:, 2] > 1].copy()
        
        # Remove the compounds from the original dataset that have inconsistent class tags. It means that the compound-target has two class labels 1 and 0.
        pdData = pdData[~pdData.iloc[:, 1].isin(set(pdA_filter.iloc[:, 1]))]
        pdData["idx"] = range(len(pdData))

        # Calculate the total number of compounds and the counts of positive and negative classifications.
        tot_mols = pdData.shape[0]
        pos_mols = pdData[pdData.iloc[:, 2] == 1].shape[0]
        neg_mols = pdData[pdData.iloc[:, 2] == 0].shape[0]
        print(f"FOR TARGET_ID {tID}")
        print("  Before removing ACs >  ", tot_mols, "P", pos_mols, "N", neg_mols)

        # Split the dataset into positive and negative class to create the tanimoto matrix
        # Separate by class
        pdMOL_0 = pdData[pdData.iloc[:, 2] == 0]
        pdMOL_1 = pdData[pdData.iloc[:, 2] == 1]

        # Create tanimoto matrix
        MT = createTanimotoMatrix(
            pdData,
            pdMOL_1.iloc[:, 3].tolist(),
            pdMOL_0.iloc[:, 3].tolist()
        )
        #print(pdMOL_0)
        #print(MT)
        #MT.to_csv(f"{tID}_v2_fast.csv", index=False, sep=",", encoding="utf-8")

 
         # Filter based on threshold
        list_selected_ids = check_byTH(
            MT,
            pdMOL_1.iloc[:, 3].tolist(),
            pdMOL_0.iloc[:, 3].tolist(),
            tID,
            th_code
        )

        print("Saving data: ", th_code, len(list_selected_ids))
 
        # Keep only filtered rows
        pdMOL_FILTER = pdData[pdData.iloc[:, 3].isin(list_selected_ids)].copy()

        # Stats after AC removal
        tot_mols = pdMOL_FILTER.shape[0]
        pos_mols = pdMOL_FILTER[pdMOL_FILTER.iloc[:, 2] == 1].shape[0]
        neg_mols = pdMOL_FILTER[pdMOL_FILTER.iloc[:, 2] == 0].shape[0]
        print("  After removing ACs >  ", tot_mols, "P", pos_mols, "N", neg_mols)

        return pdMOL_FILTER



    except Exception as err:

        print(f"{tID} Error: {err} ")
        return None


if __name__ == '__main__':
    #compute_ml_model(96)
    argumemts = sys.argv[1:]
    th_value = round(float(argumemts[0]), 1)
    in_filename = argumemts[1]
    out_filename = argumemts[2]


    if validate_threshold(th_value):
        print(th_value)

        pdCTI = pd.read_csv(in_filename, sep="\t") 
        print(pdCTI.head())
        print("Total interactions:", pdCTI.shape[0])
        print("Total targets :", len(set(pdCTI.iloc[:, 0])))
        print("Total  compounds:", len(set(pdCTI.iloc[:, 1])))
        print(f"Removing activity cliff using a threshold of {th_value}")

        lst_data=[]
        for tID in set(pdCTI.iloc[:, 0]):
            start_time = time.time()

            pdCTIFilter = pdCTI[pdCTI.iloc[:, 0] == tID]
            
            #print(pdCTIFilter.shape)
            cleanedDF = clean_dataset_cliff(pdCTIFilter, tID, th_value)
            del cleanedDF['idx']
            #print(cleanedDF.head())
            lst_data.append(cleanedDF)

            end_time = time.time()

            elapsed_seconds = end_time - start_time
            elapsed_minutes = elapsed_seconds / 60

            print(f"In {tID}: Elapsed time: {elapsed_minutes:.2f} minutes")
            #break
        pdClean=pd.concat(lst_data)
        pdClean.to_csv(out_filename, index=False, sep="\t", encoding="utf-8")
    else:
        print("Invalid threshold value!")
