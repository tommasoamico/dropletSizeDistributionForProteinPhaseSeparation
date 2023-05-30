import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import Colormap
from typing import List, Tuple, Union
from infix import shift_infix as infix
import re
from functools import reduce


def k_mom_no_occ(df: pd.DataFrame, conc: str, conc_list: List[str], k: float) -> float:
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


def k_mom_no_occ_std(df: pd.DataFrame, conc: str, conc_list: List[str], k: float) -> float:
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


def dict_moment_no_random(df: pd.DataFrame, conc_list: List[str], k: float) -> dict:
    '''
    Dictionary of kth moments for each concentration in conc_list
    '''
    return {conc: k_mom_no_occ(df, conc, conc_list, k) for conc in conc_list}


def dict_moment_no_random_std(df: pd.DataFrame, conc_list: List[str], k: float) -> dict:
    '''
    Dictionary of kth moment's std for each concentration in conc_list
    '''
    return {conc: k_mom_no_occ_std(df, conc, conc_list, k) for conc in conc_list}


def setFontMatplotlib(fontFamily: str = 'Helvetica') -> None:
    '''
    Function to set the font for matplotlib
    '''
    plt.rc('font', family=fontFamily)


def line(x, a, b):
    '''
    A simple straight line
    '''
    return a*x + b


def k_array_lines(k_to_try: List[float], df: pd.DataFrame, conc_list: List[str], fusUnfus: bool = False) -> List[dict]:
    '''
    Rescales the k-th moment 
    '''
    if fusUnfus:
        momentFunc: callable = dict_moment
    else:
        momentFunc: callable = dict_moment_no_random
    k_array_lines: List[float] = []
    for k in k_to_try:
        k_mom: dict = momentFunc(df, conc_list, k)
        k_mom = dict(map(lambda kv: (kv[0], (kv[1])**(-1/k)), k_mom.items()))
        k_array_lines.append(k_mom)

    return k_array_lines


def createCanvas() -> Tuple[plt.Figure, plt.Axes]:
    '''
    Creates a canvas for plotting with the correct aspect ratio
    '''
    fig, ax = plt.subplots(figsize=(10, (10 * .65)))

    return fig, ax


def cmapGet(cmapName: str) -> Colormap:
    '''
    Simple function to get a colormap from matplotlib
    '''
    return matplotlib.cm.get_cmap(cmapName)


def strToFloat(listStr: List[str]) -> List[float]:
    return [float(i) for i in listStr]


def tail(listTail: list):
    return listTail[1:]


@infix
def fC(f: callable, g: callable) -> callable:
    '''
    Infix function composition
    '''
    return lambda x: f(g(x))


def varPropagation(m: float, q: float, varm: float, varq: float):
    return (q**2 / m**4) * varm + (1 / q ** 2) * varq
# (popt[1]**2 / popt[0]**4) * pcov[0][0] + (1 / popt[1] **2) * pcov[1][1])


def weightedAverage(vars: float, values: float) -> Tuple[float, float]:
    k = np.sum(1 / np.array(vars))
    y_bar = np.sum(np.array(values) / np.array(vars)) / k
    var = np.sqrt(1/k)
    return y_bar, var


def tail2(x): return (tail << fC >> tail)(x)


def cumulative_data(array) -> Tuple[np.array]:
    '''
    Dunction to compute the survival of a given array/list
    '''
    array = np.array(array)
    array = np.sort(array)
    array = array[~np.isnan(array)]
    array = array[array >= 0]
    cumul = 1 - np.arange(0, len(array))/(len(array))
    return array, cumul


def cumulativeOccurances(sizes, array) -> Tuple[np.array]:
    array = np.array(array)
    array = array[~np.isnan(array)]
    array = array[array >= 0]
    cumul = 1 - np.cumsum(array) / np.sum(array)

    return sizes, cumul


def cumDict(concList: List[str], df: pd.DataFrame, fusUnfus: bool = False) -> dict:
    '''
    Create a dictionary of cumulatives
    '''
    if fusUnfus:
        return {conc: cumulativeOccurances(df['Size'], np.array(df)[:, concList.index(conc)]) for conc in concList}
    else:
        return {conc: cumulative_data(np.array(df)[:, concList.index(conc)]) for conc in concList}


def k_moment(sizes: Union[np.array, List[float]], occurances: int, k: float) -> float:
    '''
    compute k-th moment
    '''
    sizes: np.array = np.array(sizes)
    sizes: np.array = sizes[~np.isnan(sizes)]
    return np.sum(sizes**k * occurances) / np.sum(occurances)


############################
# Pappu specific functions #
############################
def k_moment_random(df: pd.DataFrame, concentration: str, pappu_conc: List[str], k: float) -> float:
    '''
    Function to calculate the kth moment of a given concentration from a dataframe df

    conc_list: list of concentrations present in the dataframe
    '''
    array_occurances: np.array = np.array(df)
    sizes: np.array = np.array(df['Size'])
    concentration_columns_idx: int = pappu_conc.index(concentration)
    column_occurrances: np.array = array_occurances[:, concentration_columns_idx*3 + 1:
                                                    concentration_columns_idx*3 + 4]

    column_chosen_idx: int = np.random.choice(
        [0, 1, 2], 1, p=[1/3, 1/3, 1/3])[0]
    column_chosen: np.array = column_occurrances[:, column_chosen_idx]

    k_th_moment: float = k_moment(sizes, column_chosen, k)

    return k_th_moment


def dict_moment_random(df, pappu_conc, k):
    return {conc: np.mean([k_moment_random(df, conc, pappu_conc, k) for _ in range(100)])
            for conc in pappu_conc}


def dict_moment(df, pappu_conc, k):
    '''
    Dictionary of kth moments for each concentration in pappu_conc
    '''
    sizes: np.array = np.array(df['Size'])
    return {conc: k_moment(sizes, np.array(df.iloc[:, 1 + pappu_conc.index(conc)]), k)
            for conc in pappu_conc}


def dict_moment_std(df, pappu_conc, k):
    '''
    Dictionary of kth moment's std for each concentration in conc_list
    '''
    return {conc: np.std([k_moment_random(df, conc, pappu_conc, k) for _ in range(100)])
            for conc in pappu_conc}


def k_array_lines_pappu(df: pd.DataFrame, concList: List[str], kToTry: List[float]) -> Tuple[List[float], List[float]]:
    k_array_lines_a: List[dict] = []
    k_array_std_a: List[dict] = []
    for k in kToTry:
        k_mom = dict_moment_random(df, concList, k)
        k_mom_std = dict_moment_std(df, concList, k)
        k_mom = dict(map(lambda kv: (kv[0], (kv[1])**(-1/k)), k_mom.items()))
        k_array_lines_a.append(k_mom)
        k_array_std_a.append(k_mom_std)
    return k_array_lines_a, k_array_std_a


def renameColumns(df: pd.DataFrame) -> List[str]:
    '''
    Function to rename column of a dtaFrame for tha Pappu paper data
    '''
    pappu_columns: List[str] = []
    for i, column in enumerate((tail << fC >> list)(df.columns)):
        if i % 3 == 0:
            pappu_columns += [re.findall('[0-9]\.[0-9]+|[0-9]', column)
                              [0] + f'_{j}' for j in range(3)]

    return pappu_columns


def column_random(df: pd.DataFrame, concentration, pappu_conc):
    array_occurances = np.array(df)
    concentration_columns_idx = pappu_conc.index(concentration)
    column_occurrances = array_occurances[:, concentration_columns_idx*3 + 1:
                                          concentration_columns_idx*3 + 4]

    column_chosen_idx = np.random.choice([0, 1, 2], 1, p=[1/3, 1/3, 1/3])[0]
    column_chosen = column_occurrances[:, column_chosen_idx]
    return column_chosen

##########################
# Class specific classes #
##########################


def tailDict(dictionary: dict) -> dict:
    firstKey = next(iter(dictionary))
    del dictionary[firstKey]
    return dictionary

#########################


def stackList(inputList: List[np.array]) -> np.array:
    return reduce(lambda x, y: np.vstack([x, y]), inputList)
