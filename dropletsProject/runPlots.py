from modules import plots
from pathlib import Path
import numpy as np
from pprint import pprint
# , pappu_a_conc, pappu_b_conc
from modules.vectorsAndConstants import coolwarm, concOwnData
import pandas as pd


dataset: str = 'FUS'
dataStem: str = 'variances'

saveDirectory: Path = Path.home() / 'Desktop' / 'newEstimatePlots'

basePath = Path.cwd().parent
dataPath = basePath / 'data' / \
    'ownData' / 'plots' / dataStem


# experimentNumber: int = 3
# plots.plotRhoEstimate
# xData: np.ndarray = np.loadtxt(
#   dataPath / f'{dataStem}XJoined.txt', delimiter=',')

yData: np.ndarray = np.loadtxt(
    dataPath / f'{dataStem}YJoined.txt', delimiter=',')


# yError: np.ndarray = np.loadtxt(
#   dataPath / f'{dataStem}ErrorJoined.txt', delimiter=',')
# estimateDf = pd.read_csv(dataPath/ 'estimateExperimentJointOld.csv')

# plots.plotRhoEstimate(xData=xData, yData=yData,
#                     savePath=saveDirectory / f'rhoEstimateJointOld.png', cmap=coolwarm, estimateDf=estimateDf, yError=yError)
# plots.plotPhiEstimation(xData=xData, yData=yData,
#                       yError=yError, cmap=coolwarm, savePath=saveDirectory / f'fig1LeftPanels.png', dataset=dataset)
# plots.plotAlphaEstimation(xData=xData, yData=yData,
#                         savePath=saveDirectory / 'fig1RightPanel.png', yError=yError)

plots.plotLognormalParametersOwn(dataPath=Path(
    dataPath), savePath=saveDirectory / 'fig2LeftPanel.png', concList=concOwnData[:-1], parameter='variances')
# plots.plotAllCollapsedPdf(dataPath=basePath / 'data' / 'appData' /
#                          'fusUnfus', savePath=saveDirectory / 'pdfCollapse', cmap=coolwarm, concLists=[pappu_a_conc, pappu_b_conc])
# plots.plotPdf(xData=xData, yData=yData, yError=yError,
#             savePath=saveDirectory / 'fig4LeftPanel', cmap=coolwarm)
# plots.plotCriticalCollapse(xData=xData, yData=yData,
#                           yError=yError, savePath=saveDirectory / 'fig4RightPanel.png', cmap=coolwarm, rightLim=200)
