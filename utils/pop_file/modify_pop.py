import pandas as pd
import numpy as np
import sys, os, math


def load_csv(csv_file, header=None):
    """
    Return a dataframe created from a csv file
    """
    return pd.read_csv(csv_file, sep=",", header=header)


def load_excel(file_name, sheet_name, header=0):
    """
    Return a dataframe from an excel file
    """
    return pd.read_excel(file_name, sheet_name, header=header)


def save_excel(dataframe, file_name, sheet_name, header=0):
    """
    Save a dataframe to excel. Default is without any header
    """
    dataframe.to_excel(file_name, sheet_name=sheet_name, header=header, index=False)


def select_pops(dataframe, pops):
    """
    Return dataframe of union of populations
    """
    return pd.DataFrame(np.union1d(dataframe[0].values, pops[0].values))



def remove_pops(dataframe, pops):
    """
    Return dataframe of disunion of populations
    """
    return pd.DataFrame(np.setdiff1d(dataframe[0].values, pops[0].values))
