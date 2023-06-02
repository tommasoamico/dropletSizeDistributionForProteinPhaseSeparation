from modules import plots
from pathlib import Path
import numpy as np
from pprint import pprint
from modules.vectorsAndConstants import coolwarm
import pandas as pd

saveDirectory: Path = Path.home() / 'Desktop' / 'newEstimates' / 'snapFUS'
basePath = Path.cwd().parent
dataPath = basePath / 'data' / 'appData' / \
    'fusUnfus' / 'snapFUS' / 'criticalRho'
experimentNumber: int = 3
# plots.plotRhoEstimate
xData: np.ndarray = np.loadtxt(
    dataPath / f'rhoEstimationXData{experimentNumber}.txt')
yData: np.ndarray = np.loadtxt(
    dataPath / f'rhoEstimationYData{experimentNumber}.txt')
estimateDf: pd.DataFrame = pd.read_csv(
    dataPath / f'estimateExperiment{experimentNumber}.csv')
plots.plotRhoEstimate(xData=xData, yData=yData,
                      savePath=saveDirectory / f'rhoEstimate{experimentNumber}.png', cmap=coolwarm, estimateDf=estimateDf)
