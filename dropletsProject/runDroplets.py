import matplotlib.pyplot as plt
from modules.sizeMatrix import sizeInstance
from modules.vectorsAndConstants import pathPappuA, pappu_b_conc, k_to_try_pappu, critical_c_dead, critical_c_b
from modules.utilityFunctions import cumulative_data
import numpy as np
from typing import Type
from pprint import pprint
import pandas as pd

valuesToDiscard: int = 2
experimentNumber: int = 3
df: pd.DataFrame = pd.read_excel(pathPappuA, sheet_name=1, header=[3])
# Selecting the columns of interest
df: pd.DataFrame = df.iloc[:, 2:27]
rangeColumns: int = 8
savePath: str = '/Users/tommaso/Desktop/dropletSizeDistributionForProteinPhaseSeparation/data/appData/fusUnfus/snapFUS/criticalRho/'
with sizeInstance.instantiateFusUnfus(df=df, concList=pappu_b_conc, kToTry=k_to_try_pappu, rangeColumns=rangeColumns, experimentNumber=experimentNumber) as workingInstance:
    keys, values = workingInstance.getRhoEstimationLines(
        valuesToDiscard=valuesToDiscard, fusUnfus=True)

    rho, error = workingInstance.getRhoEstimate(keys, values)

    print(rho, error)
    cumulativeAnalisis: Type = workingInstance.criticalAnalysis(fusUnfus=True)

    cumulativeX, cumulativeY = cumulativeAnalisis.getCumulative(
       valuesToDiscard=valuesToDiscard, criticalCollapse=critical_c_b)

# with open(savePath + f'estimateExperiment{experimentNumber}.csv', mode='w') as f:
#    f.write(f'rho,error\n{rho},{error}')
# np.savetxt(
 #   savePath + f'cumulativeX{experimentNumber}Old.txt', cumulativeX, fmt='%.5e')
# np.savetxt(
 #   savePath + f'cumulativeY{experimentNumber}Old.txt', cumulativeY, fmt='%.5e')
#np.savetxt(savePath + 'concentrations.txt',
      #     np.array([float(conc) for conc in pappu_b_conc[valuesToDiscard:]]), fmt='%.5e')
#np.savetxt(
    #savePath + f'rhoEstimationXData{experimentNumber}.txt', keys, fmt='%.5e')
#np.savetxt(
    #savePath + f'rhoEstimationYData{experimentNumber}.txt', values, fmt='%.5e')
