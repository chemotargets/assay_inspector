#! /usr/bin/env python

__author__ = "Luca Menestrina"
__date__ = "20250130"
__copyright__ = "Copyright 2025, Chemotargets"
__license__ = ""
__credits__ = ["Data Science & Translational Research Group"]
__maintainer__ = "Luca Menestrina"
__version__ = "20250224"
__deprecated__ = False

import numpy as np

from DR_Utils import logging

from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
from sklearn.impute import SimpleImputer, KNNImputer

SCALER_TYPE = "minmax"  # TODO: move to configs
IMPUTER_STRATEGY = "mean"  # TODO: move to configs
MAX_NANS = 0.1  # TODO: find a better name and move to configs

class FeaturesPreprocessing():
    """
    Preprocess features by:
        - imputing missing values
        - scaling

    Methods:
    ----------
    fit_transform(dataObj, features_ids)
        Imputes missing values and scales the features, storing fitted values.
    transform()
        Transforms data using fitted values for imputation and scaling.
    """
    def __init__(self):
        self._is_fitted = False

        self.__imputer = None
        self.__scaler = None   

    def transform(self, dataObj):
        if not self._is_fitted:
            logging.error(f"Not fitted yet. Call 'fit_transform' before using 'transform'")
        if not dataObj.is_preprocessed:
            logging.debug("Feature preprocessing...")
            dataObj._deduplicate(self._endpoint2task)
            features_df = dataObj.DataFrame(features=self._features_ids)
            data_df = dataObj.DataFrame()

            features_df = self._impute(features_df)
            features_df = self._scale(features_df)
            # Add here new manipulations

            data_df[features_df.columns] = features_df
            dataObj = dataObj.copy(_dataframe=data_df, is_preprocessed=True)

        return dataObj

    def fit_transform(self, dataObj, features_ids, endpoint2task):
        self._features_ids = features_ids
        self._endpoint2task = endpoint2task

        if not self._is_fitted:
            if dataObj.is_preprocessed:
                logging.warning("Data already preprocessed, reprocessing...")
            else:
                logging.debug("Feature preprocessing...")
            dataObj._deduplicate(self._endpoint2task)
            features_df = dataObj.DataFrame(features=self._features_ids)
            self._features_columns = features_df.columns.tolist()
            data_df = dataObj.DataFrame()
            data_df = data_df.replace((np.inf, -np.inf), np.nan)

            # self._features_defaults = self.__get_feature_defaults(data_df)
            self.__imputer = self._get_imputer().fit(data_df.loc[:, self._features_columns].values)

            # Impute missing values assigning mean
            data_df = self._impute(data_df)

            # Scale feature values
            self.__scaler = self._get_scaler().fit(data_df.loc[:, self._features_columns].values)
            data_df = self._scale(data_df)

            # Add here new manipulations

            dataObj = dataObj.copy(_dataframe=data_df, is_preprocessed=True)
            self._is_fitted = True

        return dataObj

    def _impute(self, data_df):
        if self._get_imputer() is None:
            logging.error(f"Not fitted yet. Call 'fit_transform' before using 'transform'")
        else:
            data_df = data_df.replace((np.inf, -np.inf), np.nan)
            # Remove mols with more than 10% of not-finite features
            data_df = data_df.dropna(thresh=len(self._features_columns)*(1-MAX_NANS), axis=0)  # TODO: should it be here?
            data_df.loc[:, self._features_columns] = self._get_imputer().transform(data_df.loc[:, self._features_columns].values)

        return data_df

    def _scale(self, data_df):
        if self._get_scaler() is None:
            logging.error(f"Not fitted yet. Call 'fit_transform' before using 'transform'")
        else:
            data_df.loc[:, self._features_columns] = self._get_scaler().transform(data_df.loc[:, self._features_columns].values)

        return data_df

    def _get_imputer(self):
        if self.__imputer is None:
            if IMPUTER_STRATEGY == "knn":
                self.__imputer = KNNImputer()
            elif IMPUTER_STRATEGY in ["mean", "median", "most_frequent", "constant"]:
                self.__imputer = SimpleImputer(strategy=IMPUTER_STRATEGY)
            else:
                logging.error(f"Imputer strategy '{IMPUTER_STRATEGY}' not available")
        return self.__imputer

    def _get_scaler(self):
        if self.__scaler is None:
            match SCALER_TYPE:
                case "minmax":
                    self.__scaler = MinMaxScaler()
                case "standard":
                    self.__scaler = StandardScaler()
                case "robust":
                    self.__scaler = RobustScaler()
                case _:
                    logging.error(f"Scaler '{SCALER_TYPE}' not available")
        return self.__scaler