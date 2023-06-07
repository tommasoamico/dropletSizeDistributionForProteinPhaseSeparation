from modules import plots
from pathlib import Path
import numpy as np
from pprint import pprint
from modules.vectorsAndConstants import coolwarm
import pandas as pd

saveDirectory: Path = Path.home() / 'Desktop' / 'newEstimatesPlots' / 'snapFUS'
basePath = Path.cwd().parent
dataPath = basePath / 'data' / 'appData' / \
    'fusUnfus' / 'snapFUS' / 'criticalRho'
experimentNumber: int = 3
# plots.plotRhoEstimate
xData: np.ndarray = np.loadtxt(
    dataPath / f'rhoEstimationXData{experimentNumber}.txt')
yDataStacked: np.ndarray = np.stack([np.loadtxt(
    dataPath / f'rhoEstimationYData{number + 1}.txt') for number in range(experimentNumber)])
yData = np.mean(yDataStacked, axis=0)
yError = np.std(yDataStacked, axis=0)
estimateDf: pd.DataFrame = pd.read_csv(
    dataPath / f'estimateExperimentJoint.csv')
plots.plotRhoEstimate(xData=xData, yData=yData,
                      savePath=saveDirectory / f'rhoEstimateJoint.png', cmap=coolwarm, estimateDf=estimateDf, yError=yError)
