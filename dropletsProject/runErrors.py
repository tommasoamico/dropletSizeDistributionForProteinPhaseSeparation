import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint

dataset = 'FUS'
nameOperation: str = 'pdfCollapse'
dataPath: str = f'/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/appData/fusUnfus/{dataset}/pdfCollapse/'
nExperiments: int = 3


stackedArrayX: np.ndarray = np.stack([np.loadtxt(
    dataPath + f'{nameOperation}X{experiment + 1}.txt') for experiment in range(nExperiments)], axis=0)
stackedArrayY: np.ndarray = np.stack([np.loadtxt(
    dataPath + f'{nameOperation}Y{experiment + 1}.txt') for experiment in range(nExperiments)], axis=0)


joinedMeanX = np.mean(stackedArrayX, axis=0)
joinedMeanY = np.mean(stackedArrayY, axis=0)
joinedStdY = np.std(stackedArrayY, axis=0)/np.sqrt(3)

'''
for i in range(joinedMeanX.shape[0]):
    plt.errorbar(joinedMeanX[i, :], joinedMeanY[i, :], yerr=joinedStdY[i, :])
plt.show()
'''

np.savetxt(dataPath + f'{nameOperation}XJoined.txt', joinedMeanX, fmt='%.5e')
np.savetxt(dataPath + f'{nameOperation}YJoined.txt', joinedMeanY, fmt='%.5e')
np.savetxt(dataPath + f'{nameOperation}YErrorJoined.txt',
           joinedStdY, fmt='%.5e')
