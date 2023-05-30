import numpy as np
import pandas as pd
from typing import List, Type, Tuple, Iterable
from modules.utilityFunctions import k_array_lines, strToFloat, varPropagation, line, weightedAverage, dict_moment_no_random, stackList
from scipy.optimize import curve_fit
from modules.survivalAnalysis import survival


class sizeInstance:
    allInstances: List[Type] = []

    def __init__(self, df, concList: List[str], kToTry: List[float]) -> None:

        self.concList: List[str] = concList
        self.kToTry: List[float] = kToTry
        self.numberOfK: int = len(self.kToTry)
        self.df: pd.DataFrame = df
        self.concArray: np.array = np.array(self.concList).astype(float)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        print('Exiting...')

    def getRhoEstimationLines(self, valuesToDiscard: int, fusUnfus: bool = False) -> None:
        linesData: List[dict] = k_array_lines(
            self.kToTry, self.df, self.concList, fusUnfus=fusUnfus)

        keysData: List[List[float]] = list(
            map(lambda x: strToFloat(list(x.keys()))[valuesToDiscard:], linesData))
        valuesData: List[List[float]] = list(
            map(lambda x: list(x.values())[valuesToDiscard:], linesData))

        return keysData, valuesData

    def getRhoEstimate(self, xData: Iterable, yData: Iterable) -> Tuple[float]:
        # slopes: np.array = np.zeros(self.numberOfK)
        # intercepts: np.array = np.zeros(self.numberOfK)
        variances: np.array = np.zeros(self.numberOfK)
        rhoEstimates: np.array = np.zeros(self.numberOfK)
        for i in range(self.numberOfK):
            popt, pcov = curve_fit(line, xData[i], yData[i])
            slope, intercept = popt
            variances[i] = varPropagation(*(list(popt) + list(np.diag(pcov))))
            # slopes[i] = slope
            # intercepts[i] = intercept
            rhoEstimates[i] = - intercept / slope

        finalRho, rhoError = weightedAverage(variances, rhoEstimates)

        return finalRho, rhoError

    @staticmethod
    def arrayValuesFromDict(dictionary: dict, valuesToDiscard: int = 0) -> np.array:
        return np.array(list(dictionary.values())[valuesToDiscard:])

    @staticmethod
    def arrayKeysFromDict(dictionary: dict, valuesToDiscard: int = 0) -> np.array:
        return np.array(list(dictionary.keys())[valuesToDiscard:])

    def getPhiEstimationPoints(self, criticalValue: float, valuesToDiscard: int) -> Tuple[np.array]:
        '''
        Thought for RP3 data
        '''
        xData: List[np.array] = []
        yData: List[np.array] = []
        for k in self.kToTry:
            momentRatio: np.array = self.arrayValuesFromDict(dict_moment_no_random(
                self.df, self.concList, k+1)) / self.arrayValuesFromDict(dict_moment_no_random(self.df, self.concList, k))

            x = np.array([np.abs((float(conc) - criticalValue)/criticalValue)
                         for conc in self.concList])[valuesToDiscard:]
            xData.append(-np.log(x))
            yData.append(np.log(momentRatio[valuesToDiscard:]))

        return stackList(xData), stackList(yData)

    def getAlphaEstimationPoints(self, valuesToDiscard: int, criticalValue: float) -> Tuple[np.array]:
        '''
        Thought for RP3 data
        '''
        yData: np.array = np.log(np.nanmean(
            self.df.iloc[:, valuesToDiscard:], axis=0))
        xData: np.array = -np.log(
            np.abs((self.concArray - criticalValue)[valuesToDiscard:] / criticalValue))

        return xData, yData

    @classmethod
    def instantiateFromPath(cls, pathData: str, concList: List[str], kToTry: List[float], kwargs: dict = {'header': 0}) -> Type:
        return sizeInstance(df=pd.read_csv(pathData, **kwargs), concList=concList, kToTry=kToTry)

    def criticalAnalysis(self, fusUnfus: bool = False) -> Type[survival]:
        return survival(self, fusUnfus=fusUnfus)

    @classmethod
    def instantiateFusUnfus(cls, df: pd.DataFrame, rangeColumns: int, experimentNumber: int, concList: List[str], kToTry: List[float]) -> Type:
        assert np.logical_and(experimentNumber > 0, experimentNumber <=
                              3), 'The number of experiment is has to be between 1 and 3'
        columnsToSelect: List[int] = [0] + \
            [i*3 + experimentNumber for i in range(rangeColumns)]

        dfColumns: List[str] = ['Size'] + list(
            filter(lambda x: 'Unnamed' not in x, list(df.columns)))[1:]
        dfToInsert = df.iloc[:, columnsToSelect]
        dfToInsert.columns = dfColumns
        return sizeInstance(df=dfToInsert, concList=concList, kToTry=kToTry)
