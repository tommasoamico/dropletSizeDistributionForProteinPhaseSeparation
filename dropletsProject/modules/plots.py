import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from modules.vectorsAndConstants import coolwarm, k_to_try_dead
from matplotlib.colors import Colormap
import pandas as pd
from typing import List, Iterable
from scipy.stats import linregress

kToTry: List[float] = k_to_try_dead
plotWidth: float = 10
hightRatio: float = .65

markers: List[str] = ['o', 's', 'D', 'v', '>', '<']


def plotRhoEstimate(xData: np.ndarray, yData: np.ndarray, savePath: Path, cmap: Colormap, estimateDf: pd.DataFrame) -> None:
    colors = cmap(np.linspace(0, 1, len(kToTry)))
    plt.figure(figsize=(plotWidth, plotWidth * hightRatio))
    for i, k in enumerate(kToTry):
        xArray: np.ndarray = xData[i, :]
        yArray: np.ndarray = yData[i, :]
        plt.scatter(xArray, yArray, s=140,
                    color=colors[i], zorder=3, linewidth=.7, marker=markers[i])
        rhoEstimate = estimateDf['rho']
        errorEstimate = estimateDf['error']
        xAxis = np.linspace(0, rhoEstimate + 2 * errorEstimate, 100)
        slope, intercept, _, _, _ = linregress(xArray, yArray)
        plt.plot(xAxis, xAxis * slope + intercept, color=colors[i], zorder=2)
        plt.ylim(bottom=0)
        plt.scatter(rhoEstimate, 0, s=100, color='#FFD700',
                    label='$ \\rho_c$ estimate', zorder=5)
        plt.errorbar(rhoEstimate, 0, xerr=errorEstimate, color='#FFD700',
                     linewidth=3, capsize=5, capthick=3)
        plt.tick_params(axis='both', which='major', labelsize=17, length=6)
        plt.savefig(savePath, dpi=300, bbox_inches='tight')
