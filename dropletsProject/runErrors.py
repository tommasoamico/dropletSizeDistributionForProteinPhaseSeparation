import numpy as np
import matplotlib.pyplot as plt

dataPath: str = '/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/appData/fusUnfus/snapFUS/pdfCollapse/'
nExperiments: int = 3

stackedArrayX: np.ndarray = np.stack([np.loadtxt(
    dataPath + f'pdfCollapseX{experiment + 1}.txt') for experiment in range(nExperiments)], axis=0)
stackedArrayY: np.ndarray = np.stack([np.loadtxt(
    dataPath + f'pdfCollapseY{experiment + 1}.txt') for experiment in range(nExperiments)], axis=0)

joinedMeanX = np.mean(stackedArrayX, axis=0)
joinedMeanY = np.mean(stackedArrayY, axis=0)
joinedStdY = np.std(stackedArrayY, axis=0)

'''
np.savetxt(dataPath + f'pdfCollapseXJoined.txt', joinedMeanX, fmt='%.5e')
np.savetxt(dataPath + f'pdfCollapseYJoined.txt', joinedMeanY, fmt='%.5e')
np.savetxt(dataPath + f'pdfCollapseYErrorJoined.txt',
           joinedStdY, fmt='%.5e')
'''
