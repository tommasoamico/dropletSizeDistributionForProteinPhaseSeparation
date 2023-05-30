import matplotlib
from matplotlib.colors import Colormap
import numpy as np
from modules.utilityFunctions import cmapGet, tail
from typing import List
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

################
# DEAD protein #
################
pathDead: str = '/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/dead/capflex_point_corr.csv'
conc_dead: List[str] = ['3.33', '6.67', '13.33', '26.67', '33.33']
k_to_try_dead: List[float] = [0.5, 1, 1.5, 2]
markersDead: List[str] = ['^', 's', 'o', 'D']
cmapDead: Colormap = cmapGet('PuRd')
# Excluding values way undersampled and very far from the critical concentration
concDeadEffective: List[str] = tail(conc_dead)
conc_dead_phi: List[str] = tail(concDeadEffective)
rgbaValuesDead: List[float] = cmapDead(np.linspace(0, 1, len(conc_dead[1:])))
critical_c_dead: float = 40

###################
# pappuA/pappuB   #
###################
pathPappuA: str = '/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/pappu/Figure4A_B_C.xlsx'
pappu_a_conc: List[str] = ['0.125', '0.25', '0.5', '1', '2']
critical_c_a = 2
k_to_try_pappu = [0.5, 1, 1.5, 2]
cmapPappuA = matplotlib.cm.get_cmap('Greens')
num_tries: int = 50
rand_mean_dict: dict = {conc: pd.DataFrame() for conc in pappu_a_conc}

################
# pappuB       #
################
pathPappuB: str = '/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/pappu/Figure4A_B_C.xlsx'
pappu_b_conc = ['0.125', '0.25', '0.5', '1', '1.5', '2', '2.5', '3']
cmapPappuB = LinearSegmentedColormap.from_list(
    'rg', ['cornflowerblue', 'royalblue', 'orchid', 'blue', 'navy', 'black'], N=256)
critical_c_b = 3
cmapPappuBBlues = matplotlib.cm.get_cmap('Blues')
