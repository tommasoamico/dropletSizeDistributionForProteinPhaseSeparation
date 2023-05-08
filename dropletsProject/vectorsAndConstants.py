import matplotlib
from matplotlib.colors import Colormap
import numpy as np
from dropletsProject.utilityFunctions import cmapGet, tail
from typing import List
import pandas as pd

################
# DEAD protein # 
################
pathDead:str = './data/dead/capflex_point_corr.csv'
conc_dead:List[str] = ['3.33', '6.67', '13.33', '26.67', '33.33']
k_to_try_dead:List[float] = [0.5, 1, 1.5, 2]
markersDead:List[str] = ['^', 's', 'o', 'D']
cmapDead:Colormap = cmapGet('PuRd')
concDeadEffective:List[str] = tail(conc_dead)   #Excluding values way undersampled and very far from the critical concentration
conc_dead_phi:List[str] = tail(concDeadEffective)
rgbaValuesDead:List[float] = cmapDead(np.linspace(0, 1, len(conc_dead[1:]))) 
critical_c_dead:float = 40

################
# pappuA       # 
################
pathPappuA:str = './data/pappu/Figure4A_B_C.xlsx'
pappu_a_conc:List[str] = ['0.125', '0.25','0.5', '1', '2']
critical_c_a = 2
k_to_try_pappu_a = [0.5, 1, 1.5, 2]
cmapPappuA = matplotlib.cm.get_cmap('Greens') 
num_tries:int = 50 
rand_mean_dict:dict = {conc: pd.DataFrame() for conc in pappu_a_conc}