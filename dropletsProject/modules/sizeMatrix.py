import numpy as np
import pandas as pd
from typing import List, Type, Tuple, Iterable, Optional
from modules.utilityFunctions import k_array_lines, strToFloat, varPropagation, line, weightedAverage, dict_moment_no_random, stackList, dict_moment, k_moment, sNotLognormal, sigmaLognormal
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
        self.sizeArray: np.ndarray = np.array(self.df['Size'])
        self.occurrances: np.ndarray = np.array(self.df.iloc[:, 1:])

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

    def getPhiEstimationPoints(self, criticalValue: float, valuesToDiscard: int, fus: bool = False) -> Tuple[np.ndarray]:
        if fus:
            func: callable = dict_moment
        else:
            func: callable = dict_moment_no_random
        xData: List[np.array] = []
        yData: List[np.array] = []
        for k in self.kToTry:
            momentRatio: np.array = self.arrayValuesFromDict(func(
                self.df, self.concList, k+1)) / self.arrayValuesFromDict(func(self.df, self.concList, k))

            x: np.ndarray = np.array([np.abs((float(conc) - criticalValue)/criticalValue)
                                     for conc in self.concList])[valuesToDiscard:]
            xData.append(-np.log(x))
            yData.append(np.log(momentRatio[valuesToDiscard:]))

        return stackList(xData), stackList(yData)

    def getAlphaEstimationPoints(self, valuesToDiscard: int, criticalValue: float, fus: bool = False) -> Tuple[np.ndarray]:
        if fus:
            dfAlpha = self.df.iloc[:, 1 + valuesToDiscard:]
            nPoints = np.array(dfAlpha).shape[1]
            yData: np.array = stackList([np.log(k_moment(
                sizes=self.df.iloc[:, 0], occurances=self.df.iloc[:, 1 + i], k=1)) for i in range(nPoints)])
        else:
            yData: np.array = np.log(np.nanmean(
                self.df.iloc[:, valuesToDiscard:], axis=0))

        xData: np.array = -np.log(
            np.abs((self.concArray - criticalValue)[valuesToDiscard:] / criticalValue))

        return xData, yData

    def getAllSnot(self) -> np.ndarray:
        '''
        Thaught for FUS and snapFUS data
        '''

        allSnotLognormal: np.ndarray = np.apply_along_axis(
            func1d=lambda x: sNotLognormal(sizes=self.sizeArray, occurrances=x), axis=0, arr=self.occurrances)
        return allSnotLognormal

    def getAllSigmaLognormal(self) -> np.ndarray:

        allSigmaLognormal: np.ndarray = np.apply_along_axis(
            func1d=lambda x: sigmaLognormal(sizes=self.sizeArray, occurrances=x), axis=0, arr=self.occurrances)
        return allSigmaLognormal

    def getPdf(self) -> np.ndarray:

        sizeDiff: np.ndarray = np.diff(np.insert(self.sizeArray, 0, [0]))
        pS: np.ndarray = np.apply_along_axis(
            func1d=lambda x: x/(np.sum(x) * sizeDiff), axis=0, arr=self.occurrances)
        return pS

    def getCollapsedPdf(self) -> Tuple[np.ndarray]:
        '''
        Thaught for FUS and snapFUS data
        '''

        allSnotLognormal: np.ndarray = self.getAllSnot()

        allSigmaLognormal: np.ndarray = self.getAllSigmaLognormal()
        pS: np.ndarray = self.getPdf()
        sizeReplicated: np.ndarray = np.array(
            [self.sizeArray] * len(self.concList)).T
        yData: np.ndarray = ((pS * self.sizeArray.reshape((-1, 1))) *
                             allSigmaLognormal.reshape((1, -1)))
        xData: np.ndarray = (np.log(
            sizeReplicated) - np.log(allSnotLognormal).reshape((1, -1))) / allSigmaLognormal.reshape((1, -1))
        return xData.T, yData.T

    def getCriticalCollapse(self, criticalValue: float) -> Tuple[np.ndarray]:
        sizeReplicated: np.ndarray = np.array(
            [self.sizeArray] * len(self.concList)).T
        concFloat: np.ndarray = np.array(
            [float(conc) for conc in self.concList])
        xData: np.ndarray = sizeReplicated * (1 - (concFloat / criticalValue))

        pS: np.ndarray = self.getPdf()
        yData: np.ndarray = pS / (1 - (concFloat / criticalValue))

        return xData, yData

    @classmethod
    def instantiateFromPath(cls, pathData: str, concList: List[str], kToTry: List[float], kwargs: dict = {'header': 0}) -> Type:
        return sizeInstance(df=pd.read_csv(pathData, **kwargs), concList=concList, kToTry=kToTry)

    def criticalAnalysis(self, fusUnfus: bool = False) -> Type[survival]:
        return survival(self, fusUnfus=fusUnfus)

    @classmethod
    def instantiateFusUnfus(cls, df: pd.DataFrame, rangeColumns: int, experimentNumber: int, concList: List[str], kToTry: List[float], finiteSizePercentage: Optional[float] = None) -> Type:
        assert np.logical_and(experimentNumber > 0, experimentNumber <=
                              3), 'The number of experiment is has to be between 1 and 3'
        columnsToSelect: List[int] = [0] + \
            [i*3 + experimentNumber for i in range(rangeColumns)]

        dfColumns: List[str] = ['Size'] + list(
            filter(lambda x: 'Unnamed' not in x, list(df.columns)))[1:]
        dfToInsert = df.iloc[:, columnsToSelect]
        dfToInsert.columns = dfColumns
        if finiteSizePercentage is not None:
            dfToInsert.loc[:, 'Size':] = dfToInsert.loc[:,
                                                        'Size':] // (1 / finiteSizePercentage)
        return sizeInstance(df=dfToInsert, concList=concList, kToTry=kToTry)
