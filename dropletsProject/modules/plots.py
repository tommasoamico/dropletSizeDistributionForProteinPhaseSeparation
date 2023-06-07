import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from modules.vectorsAndConstants import coolwarm, k_to_try_dead
from matplotlib.colors import Colormap
import pandas as pd
from typing import List, Iterable, Union
from scipy.stats import linregress

kToTry: List[float] = k_to_try_dead
plotWidth: float = 10
hightRatio: float = .65

markers: List[str] = ['o', 's', 'D', 'v', '>', '<']


def plotRhoEstimate(xData: np.ndarray, yData: np.ndarray, savePath: Path, cmap: Colormap, estimateDf: pd.DataFrame, yError: Union[None, np.ndarray] = None) -> None:
    colors = cmap(np.linspace(0, 1, len(kToTry)))

    fig, ax = plt.subplots(1, 1, figsize=(plotWidth, plotWidth * hightRatio))
    if isinstance(yError, np.ndarray):
        inset = ax.inset_axes([.50, .58, .49, .4])
    for i, k in enumerate(kToTry):
        xArray: np.ndarray = xData[i, :]
        yArray: np.ndarray = yData[i, :]
        ax.scatter(xArray, yArray, s=140,
                   color=colors[i], zorder=3, linewidth=.7, marker=markers[i])
        rhoEstimate = estimateDf['rho']
        errorEstimate = estimateDf['error']
        xAxis = np.linspace(0, rhoEstimate + 2 * errorEstimate, 100)
        slope, intercept, _, _, _ = linregress(xArray, yArray)
        ax.plot(xAxis, xAxis * slope + intercept, color=colors[i], zorder=2)
        ax.set_ylim(bottom=0)
        ax.scatter(rhoEstimate, 0, s=100, color='#FFD700',
                   label='$ \\rho_c$ estimate', zorder=5)
        ax.errorbar(rhoEstimate, 0, xerr=errorEstimate, color='#FFD700',
                    linewidth=3, capsize=5, capthick=3)
        ax.tick_params(axis='both', which='major', labelsize=17, length=6)
        if isinstance(yError, np.ndarray):
            yErrorArray: np.ndarray = yError[i, :]
            inset.errorbar(xArray, yArray, yerr=yErrorArray,
                           fmt='none', capsize=3, color=colors[i])

    # plt.show()
    plt.savefig(savePath, dpi=300, bbox_inches='tight')
