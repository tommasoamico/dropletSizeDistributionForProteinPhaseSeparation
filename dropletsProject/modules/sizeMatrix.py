import numpy as np
import pandas as pd
from typing import List, Type, Tuple, Iterable, Optional, Union
from modules.utilityFunctions import k_array_lines, strToFloat, varPropagation, line, weightedAverage, dict_moment_no_random, stackList, dict_moment, k_moment, sNotLognormal, sigmaLognormal, getHistogram, k_array_lines_phi
from scipy.optimize import curve_fit
from modules.survivalAnalysis import survival
from pathlib import Path
from typeguard import typechecked
from pprint import pprint


class sizeInstance:
    allInstances: List[Type] = []

    def __init__(self, df, concList: List[str], kToTry: List[float], fusUnfus=False) -> None:

        self.concList: List[str] = concList
        self.kToTry: List[float] = kToTry
        self.numberOfK: int = len(self.kToTry)
        self.df: pd.DataFrame = df
        self.concArray: np.array = np.array(self.concList).astype(float)
        self.fusUnfus = fusUnfus
        if self.fusUnfus:
            self.sizeArray: np.ndarray = np.array(self.df['Size'])
            self.occurrances: np.ndarray = np.array(self.df.iloc[:, 1:])

        else:
            self.sizeArray: np.ndarray = self.df.to_numpy()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        print('Exiting...')

    def getRhoEstimationLines(self, valuesToDiscard: int) -> None:
        linesData: List[dict] = k_array_lines(
            self.kToTry, self.df, self.concList, fusUnfus=self.fusUnfus)

        keysData: List[List[float]] = list(
            map(lambda x: strToFloat(list(x.keys()))[valuesToDiscard:], linesData))
        valuesData: List[List[float]] = list(
            map(lambda x: list(x.values())[valuesToDiscard:], linesData))

        return keysData, valuesData

    def getRhoEstimationLinesPhi(self, valuesToDiscard: int, phi: float) -> None:
        linesData: List[dict] = k_array_lines_phi(
            self.kToTry, self.df, self.concList, phi=phi)

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

    def getPhiEstimationPoints(self, criticalValue: float, valuesToDiscard: int) -> Tuple[np.ndarray]:
        if self.fusUnfus:
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

    def getAlphaEstimationPoints(self, valuesToDiscard: int, criticalValue: float) -> Tuple[np.ndarray]:
        if self.fusUnfus:
            dfAlpha = self.df.iloc[:, 1 + valuesToDiscard:]

            nPoints: np.ndarray = np.array(dfAlpha).shape[1]
            yData: np.ndarray = stackList([np.log(k_moment(
                sizes=self.df.iloc[:, 0], occurances=self.df.iloc[:, 1 + i], k=1)) for i in range(nPoints)])
        else:
            yData: np.ndarray = np.log(np.nanmean(
                self.df.iloc[:, valuesToDiscard:], axis=0))

        xData: np.ndarray = - (np.log(
            np.abs((self.concArray - criticalValue)[valuesToDiscard:] / criticalValue)))

        return xData, yData.flatten()

    def getAllSnot(self) -> np.ndarray:

        if self.fusUnfus:
            allSnotLognormal: np.ndarray = np.apply_along_axis(
                func1d=lambda x: sNotLognormal(sizes=self.sizeArray, occurrances=x), axis=0, arr=self.occurrances)
            return allSnotLognormal
        else:
            sizesArraySnot: np.ndarray = self.sizeArray.copy()
            sizesArraySnot[sizesArraySnot <= 0] = np.nan
            allSnotLognormal = np.exp(
                np.nanmean(np.log(sizesArraySnot), axis=0))
            return allSnotLognormal

    def getAllSigmaLognormal(self) -> np.ndarray:

        if self.fusUnfus:
            allSigmaLognormal: np.ndarray = np.apply_along_axis(
                func1d=lambda x: sigmaLognormal(sizes=self.sizeArray, occurrances=x), axis=0, arr=self.occurrances)
            return allSigmaLognormal
        else:
            sizesArraySnot: np.ndarray = self.sizeArray.copy()
            sizesArraySnot[sizesArraySnot <= 0] = np.nan

            allSigmaLognormal: np.ndarray = np.nanstd(
                np.log(sizesArraySnot), axis=0)
            return allSigmaLognormal

    def getPdf(self, bins: Union[str, int] = 'fd') -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        if self.fusUnfus:
            sizeDiff: np.ndarray = np.diff(np.insert(self.sizeArray, 0, [0]))
            pS: np.ndarray = np.apply_along_axis(
                func1d=lambda x: x/(np.sum(x) * sizeDiff), axis=0, arr=self.occurrances)
            return pS
        else:
            pS: np.ndarray = np.apply_along_axis(
                func1d=lambda x: getHistogram(x, bins=bins)[0], axis=0, arr=self.sizeArray)
            s: np.ndarray = np.apply_along_axis(
                func1d=lambda x: getHistogram(x, bins=bins)[1], axis=0, arr=self.sizeArray)
            return s.T, pS.T

    def getCollapsedPdf(self, bins: Union[str, int] = 'fd') -> Tuple[np.ndarray]:

        if self.fusUnfus:
            allSnotLognormal: np.ndarray = self.getAllSnot()

            allSigmaLognormal: np.ndarray = self.getAllSigmaLognormal(
            )
            pS: np.ndarray = self.getPdf()
            sizeReplicated: np.ndarray = np.array(
                [self.sizeArray] * len(self.concList)).T
            yData: np.ndarray = ((pS * self.sizeArray.reshape((-1, 1))) *
                                 allSigmaLognormal.reshape((1, -1)))
            xData: np.ndarray = (np.log(
                sizeReplicated) - np.log(allSnotLognormal).reshape((1, -1))) / allSigmaLognormal.reshape((1, -1))

            return xData.T, yData.T

        else:
            allSnotLognormal: np.ndarray = self.getAllSnot()
            allSigmaLognormal: np.ndarray = self.getAllSigmaLognormal(
            )

            pS: np.ndarray = np.apply_along_axis(
                func1d=lambda x: getHistogram(x, bins=bins)[0], axis=0, arr=self.sizeArray)
            s: np.ndarray = np.apply_along_axis(
                func1d=lambda x: getHistogram(x, bins=bins)[1], axis=0, arr=self.sizeArray)
            yData: np.ndarray = pS * s * \
                allSigmaLognormal.reshape((1, -1))
            xData: np.ndarray = (np.log(
                s) - np.log(allSnotLognormal).reshape((1, -1))) / allSigmaLognormal.reshape((1, -1))

            return xData.T, yData.T

    def getCriticalCollapse(self, criticalValue: float, bins: Union[str, int]) -> Tuple[np.ndarray]:
        if self.fusUnfus:
            sizeReplicated: np.ndarray = np.array(
                [self.sizeArray] * len(self.concList)).T
            concFloat: np.ndarray = np.array(
                [float(conc) for conc in self.concList])
            xData: np.ndarray = sizeReplicated * \
                (1 - (concFloat / criticalValue))

            pS: np.ndarray = self.getPdf()
            yData: np.ndarray = pS / (1 - (concFloat / criticalValue))

            return xData, yData
        else:
            pS: np.ndarray = np.apply_along_axis(
                func1d=lambda x: getHistogram(x, bins=bins)[0], axis=0, arr=self.sizeArray)
            s: np.ndarray = np.apply_along_axis(
                func1d=lambda x: getHistogram(x, bins=bins)[1], axis=0, arr=self.sizeArray)
            xData: np.ndarray = s * (1 - (self.concArray / criticalValue))

            return xData.T, pS.T

    @classmethod
    def instantiateFromPath(cls, pathData: str, concList: List[str], kToTry: List[float], kwargs: dict = {'header': 0}, fusUnfus: bool = False) -> "sizeInstance":
        return cls(df=pd.read_csv(pathData, **kwargs), concList=concList, kToTry=kToTry, fusUnfus=fusUnfus)

    @typechecked
    @classmethod
    def instantiateOwnData(cls, dirPath: str, concList: List[str], kToTry: List[float], experimentNumber: int, kwargs: dict = {'header': 0}, fusUnfus: bool = False, exclude100: bool = False, valuesToDiscard: int = 0, considerRestricetedConc=False) -> "sizeInstance":
        finalDf: pd.DataFrame = pd.DataFrame()
        if exclude100:
            concList = concList[:-1]
        concList = concList[valuesToDiscard:]
        if considerRestricetedConc and not exclude100:
            concList = ['20', '40', '50', '60', '75', '80']
        for conc in concList:
            df = pd.read_csv(
                dirPath + f'asyn{conc}muM.csv')
            finalDf = pd.concat(
                [finalDf, df.iloc[:, experimentNumber - 1]], axis=1)

        return cls(df=finalDf, concList=concList, kToTry=kToTry, fusUnfus=fusUnfus)

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
        return sizeInstance(df=dfToInsert, concList=concList, kToTry=kToTry, fusUnfus=True)
