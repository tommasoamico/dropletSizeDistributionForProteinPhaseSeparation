import matplotlib.pyplot as plt
from modules.sizeMatrix import sizeInstance
from modules.vectorsAndConstants import concOwnData, k_to_try_dead, criticalOwn, pathPappuA, pappu_a_conc, pappu_b_conc, criticalEstimationJoinedA, criticalEstimationJoinedB
from modules.utilityFunctions import getHistogram, stackList, weightedAverage, line
import numpy as np
from typing import Type
from pprint import pprint
import pandas as pd
from pathlib import Path
from scipy.stats import linregress
from scipy.optimize import curve_fit
from functools import reduce

'''
fusDict: dict = {'FUS': {'sheetName': 0, 'columnStart': 1, 'columnEnd': 17, 'rangeColumns': 5, 'concList': pappu_a_conc, 'criticalValue': criticalEstimationJoinedA},
                 'snapFUS': {'sheetName': 1, 'columnStart': 2, 'columnEnd': 27, 'rangeColumns': 8, 'concList': pappu_b_conc, 'criticalValue': criticalEstimationJoinedB}}
dataSet: str = 'FUS'

experimentNumber: int = 0
df: pd.DataFrame = pd.read_excel(
    pathPappuA, sheet_name=fusDict[dataSet]['sheetName'], header=[3])
# Selecting the columns of interest
df: pd.DataFrame = df.iloc[:, fusDict[dataSet]
                           ['columnStart']:fusDict[dataSet]['columnEnd']]
rangeColumns: int = fusDict[dataSet]['rangeColumns']
operation: str = 'variance'
saveDirectory: str = 'lognormalParameters'
savePath: str = f'/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/appData/fusUnfusFiniteSize/{dataSet}/{saveDirectory}/'
'''
valuesToDiscard: int = 0
experimentNumber = 0
dirPath: str = '/Users/tommaso/Workspace/dropletSizeDistributionForProteinPhaseSeparation/data/ownData/separatedData/'
# writer = pd.ExcelWriter(
#    '/Users/tommaso/Desktop/newData.xlsx', engine='openpyxl')

allXDataRho = []
allYDataRho = []

'''
FINAL ESTIMATES
New Method: Rho = 124.50766078526506 7.559580692814351
Phi: 1.0568340353328873 0.18411994958235967
Alpha: 0.7898890631542947 0.15857685938591026
Old method: Rho = 137.29942023516128 10.142499101901835
Phi: 1.2643800866184773 0.21629771292328315
Alpha: 0.9443300731629096 0.18613121936116225
'''
allXDataPhi = []
allYDataPhi = []
allXDataAlpha = []
allYDataAlpha = []
allRhoInterceptsTotal = []
allSlopesPhiTotal = []
allSlopesAlphaTotal = []
allSnotTotal = []
allVariancesTotal = []
phi = 1
for experimentNumber in range(1, 6):
    with sizeInstance.instantiateOwnData(dirPath=dirPath, concList=concOwnData, kToTry=[0.25, 0.75, 1.25, 1.75], experimentNumber=experimentNumber, fusUnfus=False, exclude100=False, considerRestricetedConc=False) as workingInstance:
        # with sizeInstance.instantiateFusUnfus(df=df, rangeColumns=fusDict[dataSet]['rangeColumns'], experimentNumber=experimentNumber, concList=fusDict[dataSet]['concList'][:-1], kToTry=k_to_try_dead) as workingInstance:
        workingInstance: Type[sizeInstance]

        keys, values = workingInstance.getRhoEstimationLinesPhi(
            valuesToDiscard=valuesToDiscard, phi=1)

        '''
        pd.DataFrame(keys).to_excel(
            writer, sheet_name=f'xDataRhoRepetition{experimentNumber}')
        
        pd.DataFrame(values).to_excel(
            writer, sheet_name=f'yDataRhoRepetition{experimentNumber}')
        '''
        allXDataRho.append(keys)
        allYDataRho.append(values)

        rho, error = workingInstance.getRhoEstimate(keys, values)

        # Plot rho estimate
        # xAxis = np.linspace(20, 110, 1000)
        allRhoIntercepts = []
        allRhoVar = []

        for i in range(len(keys)):

            # plt.scatter(keys[i], values[i])
            # result = linregress(keys[i], values[i])
            # print(-result.intercept / result.slope)
            popt, pcov = curve_fit(line, keys[i], values[i])
            allRhoIntercepts.append(- popt[1] / popt[0])
            # allRhoIntercepts.append(-result.intercept / result.slope)
            # print(- popt[1] / popt[0])
            allRhoVar.append(
                (1 / popt[0])**2 * pcov[1, 1]**2 + (popt[1]*(1 / popt[0]**2))**2 * pcov[0, 0]**2)
            # plt.plot(xAxis, xAxis * result.slope + result.intercept)
        allRhoInterceptsTotal.append(weightedAverage(
            vars=allRhoVar, values=allRhoIntercepts)[0])
        # allRhoInterceptsTotal.append(np.mean(allRhoIntercepts))
        # plt.ylim(bottom=0)
        # plt.show()

        xDataPhi, yDataPhi = workingInstance.getPhiEstimationPoints(
            criticalValue=88, valuesToDiscard=valuesToDiscard)
        allXDataPhi.append(xDataPhi)
        allYDataPhi.append(yDataPhi)
        '''
        pd.DataFrame(xDataPhi).to_excel(
            writer, sheet_name=f'xDataPhiRepetition{experimentNumber}')
        allXDataPhi.append(xDataPhi)
        pd.DataFrame(yDataPhi).to_excel(
            writer, sheet_name=f'yDataPhiRepetition{experimentNumber}')
        allYDataPhi.append(yDataPhi)
        '''
        # Plot phi estimate

        allSlopesPhi = []
        allVarPhi = []
        for i in range(xDataPhi.shape[0]):
            # plt.scatter(xDataPhi[i], yDataPhi[i])
            result = linregress(xDataPhi[i], yDataPhi[i])
            popt, pcov = curve_fit(line, xDataPhi[i], yDataPhi[i])
            # allSlopesPhi.append(result.slope)
            print(popt[0], pcov[0, 0])
            allSlopesPhi.append(popt[0])
            allVarPhi.append(pcov[0, 0])

        allSlopesPhiTotal.append(weightedAverage(
            vars=allVarPhi, values=allSlopesPhi)[0])
        print(weightedAverage(
            vars=allVarPhi, values=allSlopesPhi))
        # allSlopesPhiTotal.append(allSlopesPhi)

        # plt.show()

        xDataAlpha, yDataAlpha = workingInstance.getAlphaEstimationPoints(
            valuesToDiscard=valuesToDiscard, criticalValue=88)
        allXDataAlpha.append(xDataAlpha)
        allYDataAlpha.append(yDataAlpha)
        '''
        pd.DataFrame(xDataAlpha).to_excel(
            writer, sheet_name=f'xDataAlphaRepetition{experimentNumber}')
        allXDataAlpha.append(xDataAlpha)
        pd.DataFrame(yDataAlpha).to_excel(
            writer, sheet_name=f'yDataAlphaRepetition{experimentNumber}')
        allYDataAlpha.append(yDataAlpha)
        '''

        # Plot alpha estimate

        result = linregress(xDataAlpha, yDataAlpha)
        allSlopesAlphaTotal.append(result.slope)
        # plt.scatter(xDataAlpha, yDataAlpha)
        # plt.show()

        allSnot = workingInstance.getAllSnot()
        allSnotTotal.append(allSnot)

        allVariances = workingInstance.getAllSigmaLognormal()**2
        allVariancesTotal.append(allVariances)


print('Rho:', np.mean(allRhoInterceptsTotal), np.std(allRhoInterceptsTotal))
print('Phi', np.mean(allSlopesPhiTotal), np.std(allSlopesPhiTotal))
print('Alpha', np.mean(allSlopesAlphaTotal), np.std(allSlopesAlphaTotal))

allXRho = np.mean(allXDataRho, axis=0)
allYRho = np.mean(allYDataRho, axis=0)
allYstdRho = np.std(allYDataRho, axis=0)
basePath: str = '/Users/tommaso/Workspace/dropletSizeDistributionForProteinPhaseSeparation/data/ownData/plots/'
# np.savetxt(fname=basePath + 'rhoEstimationOld/rhoEstimationOldXJoined.txt',
#           X=allXRho, delimiter=',')
# np.savetxt(fname=basePath + 'rhoEstimationOld/rhoEstimationOldYJoined.txt',
#           X=allYRho, delimiter=',')
# np.savetxt(fname=basePath + 'rhoEstimationOld/rhoEstimationOldErrorJoined.txt',
#           X=allYstdRho, delimiter=',')

allXPhi = np.mean(allXDataPhi, axis=0)

allYPhi = np.mean(allYDataPhi, axis=0)
allYstdPhi = np.std(allYDataPhi, axis=0)

# np.savetxt(fname=basePath + 'phiEstimation/phiEstimationXJoined.txt',
#           X=allXPhi, delimiter=',')
# np.savetxt(fname=basePath + 'phiEstimation/phiEstimationYJoined.txt',
#           X=allYPhi, delimiter=',')
# np.savetxt(fname=basePath + 'phiEstimation/phiEstimationErrorJoined.txt',
#           X=allYstdPhi, delimiter=',')

allPhi = []

for i in range(allXPhi.shape[0]):
    xAxis = np.linspace(np.min(allXPhi[i]), np.max(allXPhi[i]), 1000)
    plt.scatter(allXPhi[i], allYPhi[i])
    result = linregress(allXPhi[i], allYPhi[i])
    plt.plot(xAxis, xAxis * result.slope + result.intercept)

    allPhi.append(result.slope)

    # allIntercepts.append()

plt.show()

allXAlpha = np.mean(allXDataAlpha, axis=0)
allYAlpha = np.mean(allYDataAlpha, axis=0)
allYstdAlpha = np.std(allYDataAlpha, axis=0)

# np.savetxt(fname=basePath + 'alphaEstimation/alphaEstimationXJoined.txt',
#           X=allXAlpha, delimiter=',')
# np.savetxt(fname=basePath + 'alphaEstimation/alphaEstimationYJoined.txt',
# X=allYAlpha, delimiter=',')
# np.savetxt(fname=basePath + 'alphaEstimation/alphaEstimationErrorJoined.txt',
#           X=allYstdAlpha, delimiter=',')

xAxis = np.linspace(np.min(allXAlpha), np.max(allXAlpha), 1000)
plt.scatter(allXAlpha, allYAlpha)
result = linregress(allXAlpha, allYAlpha)
plt.plot(xAxis, xAxis * result.slope + result.intercept)

# print(result.slope)
# alpha = 1 - result.slope / phi

# allIntercepts.append()

plt.show()

allXIntercepts = []
xAxis = np.linspace(20, 150, 1000)
for i in range(allXRho.shape[0]):
    plt.scatter(allXRho[i], allYRho[i])
    result = linregress(allXRho[i], allYRho[i])
    plt.plot(xAxis, xAxis * result.slope + result.intercept)
    allXIntercepts.append(- result.intercept / result.slope)
    plt.ylabel(
        '$\\left (\\frac{\\langle s^{k} \\rangle}{\\langle s \\rangle} \\right )^{-1/\\varphi(k-1)} $', fontsize=15)
    plt.xlabel('$\\rho$', fontsize=15)

    # allIntercepts.append()
plt.ylim(bottom=0)
plt.show()

# print(np.mean(allXIntercepts), np.std(allXIntercepts), 'Phi:', phi)
allSnotY = np.mean(allSnotTotal, axis=0)
allYstdSnot = np.std(allSnotTotal, axis=0)
# np.savetxt(fname=basePath + 'sNot/sNotYJoined.txt',
#          X=allSnotY, delimiter=',')
# np.savetxt(fname=basePath + 'sNot/sNotErrorJoined.txt',
#          X=allYstdSnot, delimiter=',')


allVariancesY = np.mean(allVariancesTotal, axis=0)
allVariancesYstd = np.std(allVariancesTotal, axis=0)

# np.savetxt(fname=basePath + 'variances/variancesYJoined.txt',
#          X=allVariancesY, delimiter=',')
# np.savetxt(fname=basePath + 'variances/variancesErrorJoined.txt',
#          X=allVariancesYstd, delimiter=',')

'''
pd.DataFrame(allXRho).to_excel(
writer, sheet_name=f'xDataRhoMean')


pd.DataFrame(allYRho).to_excel(
writer, sheet_name=f'yDataRhoMean')


pd.DataFrame(allXPhi).to_excel(
writer, sheet_name=f'xDataPhiMean')

pd.DataFrame(allYPhi).to_excel(
writer, sheet_name=f'yDataPhiMean')

allXAlpha = (1/5) * reduce(lambda x, y: np.array(x) +
                        np.array(y), allXDataAlpha)
pd.DataFrame(allXAlpha).to_excel(
writer, sheet_name=f'xDataAlphaMean')

allYAlpha = (1/5) * reduce(lambda x, y: np.array(x) +
                        np.array(y), allYDataAlpha)
pd.DataFrame(allYAlpha).to_excel(
writer, sheet_name=f'yDataAlphaMean')
writer.close()
'''
# Plot alpha estimate
'''
for i in range(xDataAlpha.shape[0]):
plt.scatter(xDataAlpha[i], yDataAhlpha[i])
plt.show()
'''
# result = linregress(xDataAlpha, yDataAhlpha)

# print(result.slope)

#

# Plot sNot estimate
'''
for i in range(len(allSnot)):
plt.scatter(workingInstance.concArray[i], allSnot[i])
plt.show()
'''

# allVariances = workingInstance.getAllSigmaLognormal()**2

# Plot variance estimate
'''
for i in range(len(allVariances)):
plt.scatter(workingInstance.concArray[i], allVariances[i])
plt.show()
'''

# xDataCollapsePdf, yDataCollapsePdf = workingInstance.getCollapsedPdf(
#   bins=30)

# Plot pdf Collapse
'''
for i in range(xDataCollapsePdf.shape[0]):
plt.plot(xDataCollapsePdf[i, :], yDataCollapsePdf[i, :])
plt.show()
'''

# pdfX, pdfY = workingInstance.getPdf(bins=20)

# Plot pdf
'''
for i in range(pdfX.shape[0]):
plt.plot(pdfX[i, :], pdfY[i, :])
plt.show()
'''

# criticatCollapseX, criticatCollapseY = workingInstance.getCriticalCollapse(
#   criticalValue=criticalOwn, bins=20)

# Plot critical collapse
'''
for i in range(criticatCollapseX.shape[0]):
plt.plot(criticatCollapseX[i, :], criticatCollapseY[i, :])
plt.show()
'''

# Be aware of the different implementation between fiunfus and rp3 in critical collapse

# sizes: np.ndarray = np.array(workingInstance.df['Size'])

# allPdfs: np.ndarray = workingInstance.getPdf()

# xData, yData = workingInstance.getCriticalCollapse(
#   criticalValue=fusDict[dataSet]['criticalValue'])

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
# np.savetxt(
#    savePath + f'{operation}{experimentNumber}.txt', allVariances, fmt='%.5e')

# np.savetxt(
#    savePath + f'criticalCollapseX{experimentNumber}.txt', xData.T, fmt='%.5e')
# np.savetxt(
#    savePath + f'criticalCollapseY{experimentNumber}.txt', yData.T, fmt='%.5e')
