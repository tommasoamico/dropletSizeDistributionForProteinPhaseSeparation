import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from modules.vectorsAndConstants import coolwarm, k_to_try_dead
from matplotlib.colors import Colormap
import pandas as pd
from typing import List, Iterable, Union
from scipy.stats import linregress
from scipy.optimize import curve_fit

kToTry: List[float] = k_to_try_dead
plotWidth: float = 10
hightRatio: float = .65

markers: List[str] = ['o', 's', 'D', 'v', '>', '<']


def plotRhoEstimate(xData: np.ndarray, yData: np.ndarray, savePath: Path, cmap: Colormap, estimateDf: pd.DataFrame, yError: np.ndarray) -> None:
    colors = cmap(np.linspace(0, 1, len(kToTry)))

    _, ax = plt.subplots(1, 1, figsize=(plotWidth, plotWidth * hightRatio))

    inset = ax.inset_axes([.50, .58, .49, .4])
    for i in range(len(kToTry)):
        xArray: np.ndarray = xData[i, :]
        yArray: np.ndarray = yData[i, :]
        ax.scatter(xArray, yArray, s=140,
                   color=colors[i], zorder=3, linewidth=.7, marker=markers[i], label=f'k = {kToTry[i]}')
        rhoEstimate = estimateDf['rho']
        errorEstimate = estimateDf['error']
        xAxis = np.linspace(0, rhoEstimate + 2 * errorEstimate, 100)
        slope, intercept, _, _, _ = linregress(xArray, yArray)
        ax.plot(xAxis, xAxis * slope + intercept, color=colors[i], zorder=2)
        ax.set_ylim(bottom=0)

        ax.errorbar(rhoEstimate, 0, xerr=errorEstimate, color='#14F70A',
                    linewidth=3, capsize=5, capthick=3)
        ax.tick_params(axis='both', which='major', labelsize=17, length=6)

        yErrorArray: np.ndarray = yError[i, :]
        inset.errorbar(xArray, yArray, yerr=yErrorArray,
                       fmt='none', capsize=3, color=colors[i])
    ax.scatter(rhoEstimate, 0, s=100, color='#14F70A',
               label='$ \\rho_c$ estimate', zorder=5, marker='X')
    # plt.show()
    plt.legend(loc=(0.2, 0.2), labelspacing=1)
    plt.savefig(savePath, dpi=300, bbox_inches='tight')


def plotPhiEstimation(xData: np.ndarray, yData: np.ndarray, cmap: Colormap, savePath: Path, yError: np.ndarray, dataset: str = 'FUS') -> None:
    colors: List[tuple] = [cmap(i) for i in np.linspace(0, 1, xData.shape[0])]
    intercepts = {'FUS': [4.49, 4.53, 4.58, 4.65],
                  'snapFUS': [4.59, 4.65, 4.75, 4.85]}

    xAxis: np.ndarray = np.linspace(
        np.min(xData[0, :]) - .1, np.max(xData[0, :]) + .1, 100)
    markers: List[str] = ['s', 'D', '^', 'o']
    _, ax = plt.subplots(1, 1, figsize=(plotWidth, plotWidth * hightRatio))
    axInset = ax.inset_axes([.58, .07, .415, .35])
    for i in range(xData.shape[0]):
        ax.scatter(xData[i, :], yData[i, :],
                   marker=markers[i], color=colors[i], zorder=2, s=100, label=f'k = {kToTry[i]}')

        def line(x, a, b):
            return a*x + b
        popt, pcov = curve_fit(line,
                               xData[i, 1:], yData[i, 1:])
        print(popt[0], pcov[0, 0])
        # plt.scatter(xAxis, xAxis * popt[0] + popt[1])
        ax.plot(xAxis, xAxis + intercepts[dataset][i],
                color=colors[i], zorder=1, linewidth=.7)

        axInset.errorbar(xData[i, :], yData[i, :], yerr=yError[i, :],
                         color=colors[i], fmt='none', capsize=3)

    ax.tick_params(axis='both', which='major', labelsize=17, length=6)

    plt.legend(loc=(.3, .8))
    # plt.savefig(savePath, dpi=300, bbox_inches='tight')


def plotAlphaEstimation(xData: np.ndarray, yData: np.ndarray, savePath: Path, yError: np.ndarray, dataset: str = 'FUS') -> None:
    _, ax = plt.subplots(1, 1, figsize=(plotWidth, plotWidth * hightRatio))
    blue: str = '#0F52BA'
    gray: str = '#36454F'
    xAxis: np.ndarray = np.linspace(
        np.min(xData) - .1, np.max(xData) + .1, 100)
    intercepts: dict = {'FUS': [4.40, 4.15],
                        'snapFUS': [4.48, 4.25]}

    ax.scatter(xData, yData, zorder=2,
               linewidth=.7, s=200, color=blue)
    axInset = ax.inset_axes([.58, .07, .415, .35])
    axInset.errorbar(xData, yData, yerr=yError,
                     zorder=2, fmt='none', capsize=3, color=blue)
    ax.plot(xAxis, xAxis + intercepts[dataset][0], color='black', zorder=1)
    ax.plot(xAxis, xAxis + intercepts[dataset][1], color=gray, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.savefig(savePath, dpi=300, bbox_inches='tight')


def plotLognormalParameters(dataPath: Path, savePath: Path, cmap: Colormap, concLists: List[List[str]], parameter: str = 'sNot') -> None:
    fusPath: Path = dataPath / 'FUS' / 'lognormalParameters'
    snapFusPath: Path = dataPath / 'snapFUS' / 'lognormalParameters'
    red: tuple = cmap(0.9)
    blue: str = '#0F52BA'
    _, ax = plt.subplots(1, 1, figsize=(plotWidth, plotWidth * hightRatio))
    if parameter == 'sNot':
        axInset = ax.inset_axes([.50, .05, .49, .4])
    else:
        ax.set_ylim(0.05, 0.3)
        axInset = ax.inset_axes([.50, .58, .49, .4])

    for path, color, concList in zip([fusPath, snapFusPath], [blue, red], concLists):
        floatConc: np.ndarray = np.array([float(conc) for conc in concList])
        yData: np.ndarray = np.loadtxt(
            path / f'{parameter}YJoined.txt')
        yError = np.loadtxt(
            path / f'{parameter}ErrorJoined.txt')
        ax.scatter(floatConc, yData, s=200, color=color)
        axInset.errorbar(floatConc, yData,
                         yerr=yError, zorder=2, fmt='none', capsize=3, color=color)

    ax.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.savefig(savePath, dpi=300, bbox_inches='tight')


def plotAllCollapsedPdf(dataPath: Path, savePath: Path, cmap: Colormap, concLists: List[List[str]], valuesToDiscard=[1, 2]) -> None:

    _, ax = plt.subplots(1, 1, figsize=(int(12 * 0.95), 12), tight_layout=True)
    ax.tick_params(axis='both', which='major', labelsize=17, length=6)
    xAxis = np.linspace(-4, 4, 1000)
    markers1 = ['o', 'D', '^', '>', 's']
    gray: str = '#36454F'

    def expFunc(x):
        return np.exp(-(x**2)/2)/np.sqrt(2*np.pi)
    markers2 = ['^', 'h', 'd', 'o', 'D', 'p', '8', 'v']
    for pathName, markers, bounds, concentations, values in zip(['FUS', 'snapFUS'], [markers1, markers2], [(0, .4), (.6, 1)], concLists, valuesToDiscard):
        floatConc = [float(conc) for conc in concentations]
        path: Path = dataPath / pathName / 'pdfCollapse'

        xData: np.ndarray = np.loadtxt(
            path / 'pdfCollapseXJoined.txt')[values:, :]
        yData: np.ndarray = np.loadtxt(
            path / 'pdfCollapseYJoined.txt')[values:, :]
        yError: np.ndarray = np.loadtxt(
            path / 'pdfCollapseYErrorJoined.txt')[values:, :]

        colors = [cmap(x) for x in np.linspace(*bounds, xData.shape[0])]
        for i in range(xData.shape[0]):
            ax.scatter(xData[i, :], yData[i, :], marker=markers[i], s=100,
                       color=colors[i], label=f'{pathName}, $\\rho = {floatConc[i]:.3f}$')
            ax.plot(xData[i, :], yData[i, :],
                    color=colors[i])
            ax.errorbar(xData[i, :], yData[i, :], yerr=yError[i,
                        :], zorder=2, fmt='none', capsize=5, color=colors[i])
    ax.set_xlim(-4, 4)
    ax.set_ylim(top=.5)
    ax.plot(xAxis, expFunc(xAxis), color=gray, linewidth=5, zorder=5,
            label='$\\dfrac{e^{-\\frac{x^2}{2}}}{\sqrt{(2 \pi)}}$', linestyle='--')
    ax.legend(loc=(1.1, .5))
    plt.savefig(savePath, dpi=300, bbox_inches='tight')


def plotPdf(xData: np.ndarray, yData: np.ndarray, yError: np.ndarray, savePath: Path, cmap: Colormap) -> None:
    _, ax = plt.subplots(1, 1, figsize=(plotWidth, plotWidth * hightRatio))
    ax.tick_params(axis='both', which='major', labelsize=17, length=6)
    colors = cmap(np.linspace(0, 1, xData.shape[0]))
    for i in range(xData.shape[0]):
        ax.plot(xData[i, :], yData[i, :], color=colors[i], linewidth=3)
        ax.errorbar(xData[i, :], yData[i, :], yerr=yError[i, :],
                    zorder=2, fmt='none', capsize=5, color=colors[i])

    plt.savefig(savePath, dpi=300, bbox_inches='tight')


def plotCriticalCollapse(xData: np.ndarray, yData: np.ndarray, yError: np.ndarray, savePath: Path, cmap: Colormap, rightLim: float) -> None:
    _, ax = plt.subplots(1, 1, figsize=(plotWidth, plotWidth * hightRatio))
    ax.tick_params(axis='both', which='major', labelsize=17, length=6)
    colors = cmap(np.linspace(0, 1, xData.shape[0]))
    for i in range(xData.shape[0]):
        ax.plot(xData[i, :], yData[i, :], color=colors[i], linewidth=3)
        ax.errorbar(xData[i, :], yData[i, :], yerr=yError[i, :],
                    zorder=2, fmt='none', capsize=5, color=colors[i])

    ax.set_xlim(left=0, right=rightLim)
    plt.savefig(savePath, dpi=300, bbox_inches='tight')
