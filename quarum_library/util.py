# SPDX-FileCopyrightText: 2024 TU Dortmund Malte Mues / Constructor University Bremen gGmbH David Ernst
#
# SPDX-License-Identifier: Apache-2.0
# Copyright 2024 TU Dortmund Malte Mues / Constructor University Bremen gGmbH David Ernst
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os

import pandas as pd
import numpy as np
import json
from pandas import read_csv, read_excel
from pyrolite.geochem import REE

from .exceptions import ConfigurationError
from .constants import ANOMALIES
from scipy.cluster.hierarchy import dendrogram
import logging

logger = logging.getLogger(__name__)

def read_df(filename, **kwargs):
    if filename.endswith(".xlsx") or filename.endswith(".xls"):
        return read_excel(filename, **kwargs)
    else:
        kwargs.pop("sheet_name", None)
        return read_csv(filename, **kwargs)


def enrich_df_with_metadata(df, sid, filename, source_data_name):
    df.attrs["sid"] = sid
    df.attrs["filename"] = filename
    df.attrs["source_data_name"] = source_data_name
    df.attrs["isnormalized"] = False
    return df


def get_anomalies(input_data):
    expected_Eu = input_data["Sm"] * 0.67 + input_data["Tb"] * 0.33
    expected_Gd = input_data["Sm"] * 0.33 + input_data["Tb"] * 0.67
    expected_Ce = input_data["Pr"] * 2 - input_data["Nd"]
    expected_La = input_data["Pr"] * 3 - input_data["Nd"] * 2

    anEu = input_data["Eu"] / expected_Eu
    anGd = input_data["Gd"] / expected_Gd
    anCe = input_data["Ce"] / expected_Ce
    anLa = input_data["La"] / expected_La
    anomalies = []
    if not ((0.95 < anEu) and (anEu < 1.05)):
        anomalies.append("Eu")
    if not ((0.95 < anGd) and (anGd < 1.05)):
        anomalies.append("Gd")
    if not ((0.95 < anCe) and (anCe < 1.05)):
        anomalies.append("Ce")
    if not ((0.95 < anLa) and (anLa < 1.05)):
        anomalies.append("La")
    return pd.Series([anomalies], index=[ANOMALIES])

def remove_Gd_La(series):
    anomalies = series[ANOMALIES]
    if "La" in anomalies:
        anomalies.remove("La")
    if "Gd" in anomalies:
        anomalies.remove("Gd")
    series[ANOMALIES] = anomalies
    series["LTE"] = True
    return series

def filter_anomalies(series):
    anomalies = series[ANOMALIES]
    if "La" in anomalies and "Gd" in anomalies:
        anomalies.remove("La")
        anomalies.remove("Gd")
        series[ANOMALIES] = anomalies
        series["LTE"] = True
    else:
        series["LTE"] = False
    return series


def make_model_path(output, basename, model):
    filename = basename.split(".")
    filename = ".".join(filename[:-1]) + f"_{model}.csv"
    return os.path.join(output, filename)


def get_renaming(suffix="model", elements=REE()):
    renaming = {}
    for e in REE():
        renaming[e] = f"{e}_{suffix}"
    return renaming


## The plot_dendogram function is taken from the docs: https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
## It has a BSD License
def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

def filter_labels(df):
    df = df. replace("'","\"")
    labels = json.loads(df)
    try:
        return labels[0]
    except IndexError:
        return None


def load_mol_weights():
    df = read_df(os.path.join(os.path.dirname(__file__), "data", "mol_to_ppm.xlsx"))
    return df.loc[:, ["symbol", "g/mol"]].set_index("symbol").to_dict()["g/mol"]


def convert_mol_weight_to_ppm_weight(df, columns, source_unit):
    molmass = load_mol_weights()
    if isinstance(df, pd.Series):
        raise IndexError("This function only works on dataframes")
    for elem in columns:
        if source_unit == "pmol/kg":
            intermediate = df[elem] *10**-12
        elif source_unit == "nmol/kg":
            intermediate = df[elem] *10**-9
        elif source_unit == "mmol/kg":
            intermediate = df[elem] *10**-3
        elif source_unit == "fmol/kg":
            intermediate = df[elem] *10**-15
        else:
            raise ValueError(f"Unkonw unit: {source_unit}")
        intermediate = intermediate*molmass[elem]
        intermediate = intermediate * 1000
        df.loc[:, elem] = intermediate
    return df

def normalize(df, norm_to, units="ppm"):
    if "isnormalized" in df.attrs.keys() and df.attrs["isnormalized"]:
        logger.warning(f"Dataframe is already normalized to {df.attrs['normalization']}")
    df_normed = df.pyrochem.normalize_to(norm_to, units)
    df.attrs["normalization"] = norm_to
    df.attrs["isnormalized"] = True
    df.loc[:, df_normed.columns] = df_normed
    return df

def denormalize(df, units = "ppm"):
    if not df.attrs["isnormalized"]:
        logger.warning("Dataframe is not normalized")
        return None
    else:
        df_denormed = df.pyrochem.denormalize_from(df.attrs["normalization"], units)
        df.attrs["normalization"] = None
        df.attrs["isnormalized"] = False
        df.loc[:, df_denormed.columns] = df_denormed
        return df


def write_dataframe(df, outfile, write_index=False):
    if outfile.endswith(".csv"):
        df.to_csv(outfile, index=write_index)
    elif outfile.endswith(".xls") or outfile.endswith(".xlsx"):
        df.to_excel(outfile, index=write_index)
    else:
        raise ConfigurationError("Outfile does not end in .xls/.xlsx/.csv")
    return df