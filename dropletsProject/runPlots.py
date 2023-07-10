from modules import plots
from pathlib import Path
import numpy as np
from pprint import pprint
from modules.vectorsAndConstants import coolwarm  # , pappu_a_conc, pappu_b_conc
import pandas as pd


dataset: str = 'FUS'
dataStem: str = 'phiEstimation'

saveDirectory: Path = Path.home() / 'Desktop' / 'newEstimatePlots' / dataset

basePath = Path.cwd().parent
dataPath = basePath / 'data' / 'appData' / \
    'fusUnfus' / dataset / 'phiEstimation'


# experimentNumber: int = 3
# plots.plotRhoEstimate
xData: np.ndarray = np.loadtxt(
    dataPath / f'{dataStem}XJoined.txt')

yData: np.ndarray = np.loadtxt(
    dataPath / f'{dataStem}YJoined.txt')

yError: np.ndarray = np.loadtxt(
    dataPath / f'{dataStem}ErrorJoined.txt')
# estimateDf = pd.read_csv(dataPath / 'estimateExperimentJoint.csv')

# plots.plotRhoEstimate(xData=xData, yData=yData,
#                      savePath=saveDirectory / f'rhoEstimateJoint.png', cmap=coolwarm, estimateDf=estimateDf, yError=yError)
plots.plotPhiEstimation(xData=xData, yData=yData,
                        yError=yError, cmap=coolwarm, savePath=saveDirectory / f'fig1LeftPanels.png', dataset=dataset)
# plots.plotAllCollapsedPdf(dataPath=basePath / 'data' / 'appData' /
#                          'fusUnfus', savePath=saveDirectory / 'pdfCollapse', cmap=coolwarm, concLists=[pappu_a_conc, pappu_b_conc])
# plots.plotPdf(xData=xData, yData=yData, yError=yError,
#             savePath=saveDirectory / 'fig4LeftPanel', cmap=coolwarm)
# plots.plotCriticalCollapse(xData=xData, yData=yData,
#                           yError=yError, savePath=saveDirectory / 'fig4RightPanel.png', cmap=coolwarm, rightLim=200)
