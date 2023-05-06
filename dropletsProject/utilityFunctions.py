import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import Colormap
from typing import List, Tuple
from infix import shift_infix as infix

def k_mom_no_occ(df:pd.DataFrame, conc:str, conc_list:List[str], k:float)-> float:
    '''
    Function to calculate the kth moment of a given concentration from a dataframe df

    conc_list: list of concentrations present in the dataframe
    '''
    column_chosen_idx = conc_list.index(conc)
    matrix_sizes = np.array(df)
    column_chosen = matrix_sizes[:, column_chosen_idx]
    column_chosen = column_chosen[~np.isnan(column_chosen)]
    column_chosen = column_chosen[column_chosen >= 0]
    k_moment = np.mean(column_chosen ** k)
    return k_moment


def k_mom_no_occ_std(df:pd.DataFrame, conc:str, conc_list:List[str], k:float) -> float:
    '''
    A function equivalent to the above one but that returns the standard deviation of the kth moment
    '''
    column_chosen_idx = conc_list.index(conc)
    matrix_sizes = np.array(df)
    column_chosen = matrix_sizes[:, column_chosen_idx]
    column_chosen = column_chosen[~np.isnan(column_chosen)]
    column_chosen = column_chosen[column_chosen >= 0]
    k_moment_std = np.std(column_chosen ** k)
    return k_moment_std


def dict_moment_no_random(df:pd.DataFrame, conc_list:List[str], k:float) -> dict:
    '''
    Dictionary of kth moments for each concentration in conc_list
    '''
    return {conc: k_mom_no_occ(df, conc, conc_list, k) for conc in conc_list}

def dict_moment_no_random_std(df:pd.DataFrame, conc_list:List[str], k:float) -> dict:
    '''
    Dictionary of kth moment's std for each concentration in conc_list
    '''
    return {conc: k_mom_no_occ_std(df, conc, conc_list, k) for conc in conc_list}


def setFontMatplotlib(fontFamily:str = 'Helvetica') -> None:
    '''
    Function to set the font for matplotlib
    '''
    plt.rc('font', family=fontFamily)


def line(x, a, b):
    '''
    A simple straight line
    '''
    return a*x + b


def k_array_lines(k_to_try:List[float], df:pd.DataFrame, conc_list:List[str]) -> List[dict]:
    '''
    Rescales the k-th moment 
    '''
    k_array_lines:List[float] = []
    for k in k_to_try:
        k_mom:dict = dict_moment_no_random(df, conc_list, k)
        k_mom = dict(map(lambda kv: (kv[0], (kv[1])**(-1/k)), k_mom.items()))
        k_array_lines.append(k_mom)

    return k_array_lines

def createCanvas() -> Tuple[plt.Figure, plt.Axes]:
    '''
    Creates a canvas for plotting with the correct aspect ratio
    '''
    fig, ax = plt.subplots(figsize=(10, (10 * .65)))

    return fig, ax

def cmapGet(cmapName:str) -> Colormap:
    return matplotlib.cm.get_cmap(cmapName)

def strToFloat(listStr:List[str]) -> List[float]:
    return [float(i) for i in listStr]

def tail(listTail:list):
    return listTail[1:]

@infix
def functionComposition(f:callable, g:callable):
    return lambda x: f(g(x))

tailPlusStrToFloat = tail <<functionComposition>> strToFloat