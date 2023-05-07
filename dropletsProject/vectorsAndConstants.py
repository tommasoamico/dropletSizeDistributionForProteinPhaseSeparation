import matplotlib
from matplotlib.colors import Colormap
import numpy as np
from dropletsProject.utilityFunctions import cmapGet, tail
from typing import List
pathDead:str = './data/dead/capflex_point_corr.csv'
conc_dead:list = ['3.33', '6.67', '13.33', '26.67', '33.33']
k_to_try_dead:list = [0.5, 1, 1.5, 2]
markersDead:list = ['^', 's', 'o', 'D']
cmapDead:Colormap = cmapGet('PuRd')
concDeadEffective:List[str] = tail(conc_dead)   #Excluding values way undersampled and very far from the critical concentration
conc_dead_phi:List[str] = tail(concDeadEffective)
rgbaValuesDead:List[float] = cmapDead(np.linspace(0, 1, len(conc_dead[1:]))) 
critical_c_dead:float = 40