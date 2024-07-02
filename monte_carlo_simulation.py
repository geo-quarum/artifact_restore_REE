#!/usr/bin/env python3
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
import argparse
import json
import os
import sys
from configparser import ConfigParser
from time import time

import pandas
import numpy
import pyrolite.geochem
from pandarallel import pandarallel

from quarum_library.exceptions import ConfigurationError
from quarum_library.util import read_df
from quarum_library.simulations.monte_carlo import compute_monte_carlo_data, polynomial_fit, lambdatau_fit

def monte_carlo_series(series, outfolder, id_column, fitting, excluded_elements, error, repetitions,make_plots=True):
    outfile = os.path.join(outfolder, f"monte_carlo_simulation_{str(series[id_column])}.png")
    model_data= series.pyrochem.REE
    model_data = model_data.astype(numpy.float64)
    if "exclude" in series:
        excluded_elements = excluded_elements + json.loads(series["exclude"].replace("'", '"'))
    compute_monte_carlo_data(model_data, fitting, model_data.index, outfile, exclude=excluded_elements,
                             repetitions=repetitions, error=error,make_plots=make_plots)

def run_monte_carlo_simulation(df, norm_to, fitting_function, excluded_elements, outfolder, id_column, error, repetitions, make_plots=True):
    os.makedirs(outfolder, exist_ok=True)
    df_normed = df.pyrochem.normalize_to(norm_to)
    if fitting_function != lambdatau_fit:
        df_normed = numpy.log(df_normed)
    join_columns = [id_column]
    if "exclude" in df.columns:
        join_columns.append("exclude")
    if "LTE" in df.columns:
        join_columns.append("LTE")
    df[id_column] = df[id_column].astype(str)
    work_data = pandas.concat((df[join_columns], df_normed), axis=1)
    pandarallel.initialize(progress_bar=True)
    work_data.parallel_apply(lambda x: monte_carlo_series(x, outfolder, id_column, fitting_function, excluded_elements, error, repetitions, make_plots=make_plots), axis=1)

def add_anomalies(raw_data, model_file):
    model = read_df(model_file)
    raw_data["exclude"] = model["anomalies"]
    return raw_data

def _write_default_config():
    """
    Writes the default config for this file.
    """
    config = ConfigParser()
    config["sim"] = {
        "input_file": "The file that contains the dataframe",
        "outfolder": "The results of the simulation will be written to this folder",
        "model_excludes": "A model file from polyfitting or lambda fitting to copy excluded elements",
        "norm_to": "This describes the reference standard used for normalization",
        "sample_id": "The column with the sample ID.",
        "anomal_elements": "Elements that should be excluded during fitting",
        "nrows": "This is passed to the reader for processing only parts of the dataframe",
        "nsample": "This is passed to the reader for processing only a specific sample",
        "fitting": "Either 'poly' or 'lambda', default 'lambda'",
        "error": "Error applied in the monte carlo simulation, default 0.1",
        "repetitions": "Amount of repetitions for simulation, default 100"

    }
    conf_path = os.path.join(os.path.dirname(__file__), "configs", "default_monte_carlo_simulation.ini")
    print("Writing default config to: ", conf_path)
    with open(conf_path, "w") as conffile:
        config.write(conffile)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("configfile")
    args = parser.parse_args()
    if args.configfile == "write_default":
        _write_default_config()
        sys.exit(0)
    elif not os.path.isfile(args.configfile):
        print(f"{args.configfile} does not exist")
    config = ConfigParser()
    config.read(args.configfile)

    start = time()
    fconf = config["sim"]
    df_filename = fconf["input_file"]
    basename = os.path.basename(df_filename)
    sheet = fconf.get("sheet", 0)
    if "nrows" in fconf:
        df = read_df(df_filename, nrows=int(fconf["nrows"]), sheet_name=sheet)
    elif "nsample" in fconf:
        df = read_df(df_filename, sheet_name=sheet)
        df = df[df[fconf["sample_id"]].astype(str) == str(fconf["nsample"])]
    else:
        df = read_df(df_filename, sheet_name=sheet)
    fitting_function = lambdatau_fit
    if "fitting" in fconf:
        if fconf.get("fitting") == "poly":
            fitting_function = polynomial_fit
        elif fconf.get("fitting") == "lambda":
            fitting_function = lambdatau_fit
        else:
            raise ConfigurationError(f"Wrong config value  '{fconf['fitting']}' for fitting." 
            "Only 'lambda' or 'poly' are allowed.")
    plots = True
    if "make_plots" in fconf:
        if fconf.get("make_plots") == "False":
            plots = False
    excluded_elements = json.loads(fconf.get("anomal_elements", '[]'))
    sample_id = fconf.get("sample_id")
    error = float(fconf.get("error", 0.1))
    repetitions = int(fconf.get("repetitions", 100))

    if "model_excludes" in fconf:
        model_file = fconf["model_excludes"]
        df = add_anomalies(df, model_file)
        df.to_csv("/tmp/debug_mc.csv", index=False)

    run_monte_carlo_simulation(df, fconf["norm_to"], fitting_function,  excluded_elements, fconf["outfolder"], sample_id, error, repetitions,
                               make_plots=plots)
    print("Done", f"in {time() - start} seconds")


# On OSX with a M1, we need OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
# `OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES ./monte_carlo_simulation.py
if __name__ == "__main__":
    main()
#%%
