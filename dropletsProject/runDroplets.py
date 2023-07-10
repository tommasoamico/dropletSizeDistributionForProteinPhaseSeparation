import matplotlib.pyplot as plt
from modules.sizeMatrix import sizeInstance
from modules.vectorsAndConstants import pathPappuA, pappu_b_conc, pappu_a_conc, k_to_try_pappu, critical_c_dead, critical_c_b, criticalEstimationJoinedA, criticalEstimationJoinedB
from modules.utilityFunctions import cumulative_data, stackList
import numpy as np
from typing import Type
from pprint import pprint
import pandas as pd
from pathlib import Path

fusDict: dict = {'FUS': {'sheetName': 0, 'columnStart': 1, 'columnEnd': 17, 'rangeColumns': 5, 'concList': pappu_a_conc, 'criticalValue': criticalEstimationJoinedA},
                 'snapFUS': {'sheetName': 1, 'columnStart': 2, 'columnEnd': 27, 'rangeColumns': 8, 'concList': pappu_b_conc, 'criticalValue': criticalEstimationJoinedB}}
dataSet: str = 'snapFUS'

valuesToDiscard: int = 0
experimentNumber: int = 3
df: pd.DataFrame = pd.read_excel(
    pathPappuA, sheet_name=fusDict[dataSet]['sheetName'], header=[3])
# Selecting the columns of interest
df: pd.DataFrame = df.iloc[:, fusDict[dataSet]
                           ['columnStart']:fusDict[dataSet]['columnEnd']]
rangeColumns: int = fusDict[dataSet]['rangeColumns']
operation: str = 'variance'
saveDirectory: str = 'lognormalParameters'
savePath: str = f'/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/appData/fusUnfusFiniteSize/{dataSet}/{saveDirectory}/'
with sizeInstance.instantiateFusUnfus(df=df, concList=fusDict[dataSet]['concList'], kToTry=k_to_try_pappu, rangeColumns=rangeColumns, experimentNumber=experimentNumber, finiteSizePercentage=.5) as workingInstance:
    keys, values = workingInstance.getRhoEstimationLines(
        valuesToDiscard=valuesToDiscard, fusUnfus=True)

    rho, error = workingInstance.getRhoEstimate(keys, values)

    xDataPhi, yDataPhi = workingInstance.getPhiEstimationPoints(
        criticalValue=criticalEstimationJoinedA, valuesToDiscard=valuesToDiscard, fus=True)

    xDataAlpha, yDataAhlpha = workingInstance.getAlphaEstimationPoints(
        valuesToDiscard=valuesToDiscard, criticalValue=criticalEstimationJoinedA, fus=True)

    # xDataPdfCollapse, yDataPdfCollapse = workingInstance.getCollapsedPdf()

    allSnot = workingInstance.getAllSnot()
    allVariances = workingInstance.getAllSigmaLognormal()**2

    sizes: np.ndarray = np.array(workingInstance.df['Size'])

    allPdfs: np.ndarray = workingInstance.getPdf()

    xData, yData = workingInstance.getCriticalCollapse(
        criticalValue=fusDict[dataSet]['criticalValue'])


# cumulativeAnalisis: Type = workingInstance.criticalAnalysis(fusUnfus=True)

# cumulativeX, cumulativeY = cumulativeAnalisis.getCumulative(
#  valuesToDiscard=valuesToDiscard, criticalCollapse=critical_c_b)

'''
with open(savePath + f'estimateExperiment{experimentNumber}.csv', mode='w') as f:
    f.write(f'rho,error\n{rho},{error}')
'''
'''
np.savetxt(
    savePath + f'{operation}X{experimentNumber}.txt', xData, fmt='%.5e')
'''
np.savetxt(
    savePath + f'{operation}{experimentNumber}.txt', allVariances, fmt='%.5e')

# np.savetxt(
#    savePath + f'criticalCollapseX{experimentNumber}.txt', xData.T, fmt='%.5e')
# np.savetxt(
#    savePath + f'criticalCollapseY{experimentNumber}.txt', yData.T, fmt='%.5e')
