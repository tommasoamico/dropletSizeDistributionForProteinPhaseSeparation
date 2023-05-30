from typing import Type, List, Tuple, Union
from modules.utilityFunctions import cumDict, stackList
import numpy as np
from pprint import pprint


class survival:
    def __init__(self, sizeInstance: Type, fusUnfus: bool = False) -> None:
        self.sizeInstance = sizeInstance
        self.fusUnfus = fusUnfus
        self.cumulativeDict = cumDict(
            concList=self.sizeInstance.concList, df=self.sizeInstance.df, fusUnfus=self.fusUnfus)

    def getCumulative(self, valuesToDiscard: int, criticalCollapse: Union[None, float] = None) -> Tuple[np.array]:
        allXData: List[np.array] = []
        allYData: List[np.array] = []

        for conc in self.sizeInstance.concList[valuesToDiscard:]:
            if criticalCollapse == None:
                collapseFactor = 1
            else:
                collapseFactor = np.abs(
                    ((float(conc) - criticalCollapse) / criticalCollapse))
            allXData.append(self.cumulativeDict[conc][0] * collapseFactor)
            allYData.append(self.cumulativeDict[conc][1])
        maxLenX = np.max([len(x) for x in allXData])
        maxLenY = np.max([len(y) for y in allYData])

        xMatrix = stackList(list(map(lambda x: np.pad(x.astype(float), pad_width=(
            0, maxLenX - len(x)), mode='constant', constant_values=(np.nan,)), allXData)))

        yMatrix = stackList(list(map(lambda y: np.pad(y.astype(float), pad_width=(
            0, maxLenY - len(y)), mode='constant', constant_values=(np.nan,)), allYData)))

        return xMatrix, yMatrix
