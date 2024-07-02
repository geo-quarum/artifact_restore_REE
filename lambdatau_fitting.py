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

import pandas as pd
from pandarallel import pandarallel

from quarum_library.constants import DEFAULT_UNIT
from quarum_library.fitting.lambdatau import optimize_model
from quarum_library.fitting.fitting_infra import OptimizationConfig
from quarum_library.util import read_df, enrich_df_with_metadata
from quarum_library.exceptions import NormalizationError



def fit_model(dataframe, norm, conf, outfile=None, MCsim_input_file=None, exclusion_threshold=0.1, units=DEFAULT_UNIT, **kwargs):
    """

    :param dataframe: The filename containing the data for which to compute fittings
    :param **kwargs: kwargs is passed to the reader
    :return: None, the script creates a folder next to the filename with the models in "output/models"

    We compute the model as follows:
        While excludable_elements:
            1. fit lambda_tau with current excludes and LTE=False
            2. compute difference measurement and model per element
            3. exclude Element with the highest difference that is in excludable_elements, iff threshold passed
        Fit model with LTE=True
    """

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    if norm is None:
        raise NormalizationError()
    pandarallel.initialize(progress_bar=True)
    dataframe_ree = dataframe.pyrochem.REE
    dataframe_ree = dataframe_ree.pyrochem.normalize_to(norm, units)
    if "LTE" in dataframe.columns:
        dataframe_ree["LTE"] = dataframe["LTE"]
    if "exclude" in dataframe.columns:
        dataframe_ree["exclude"] = dataframe["exclude"]
    optimized_model = optimize_model(dataframe_ree, norm, conf)
    pd.concat([dataframe[dataframe.attrs["sid"]],optimized_model], axis=1).to_csv(outfile, index=False)



def _write_default_config():
    """
    Writes the default config for this file.
    """
    config = ConfigParser()
    config["fit"] = {
        "input_file": "The file that contains the dataframe",
        "outfolder": "The folder to put the different models (use it with mode=4models instead of outfile)",
        "outfile": "Use this, if you expect a single file as output (use it with mode=optimize_fitting instead of outfolder)",
        "MCsim_input_file": "File location where the input file for the Monte Carlo Simulation is stored",
        "norm_to": "This describes the reference standard used for noramlization",
        "sample_id": "The column with the sample ID.",
        "lte": "LTE True or False during the fitting, default False",
        "anomal_elements": "Optional, default is [La, Ce, Eu, Gd, Lu]",
        "exclusion_threshold": "optional, default is 0.15",
        "nrows": "This is passed to the reader for processing only parts of the dataframe",
        "nsample": "This is passed to the reader for processing only a specific sample",
        "algorithm": "Choose 'oneil' or 'opt'. 'opt' is the default"
    }
    conf_path = os.path.join(os.path.dirname(__file__), "configs", "default_fitting_lambda_tau.ini")
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
    mode = "optimize_fitting"
    fconf = config["fit"]
    if "mode" in fconf:
        mode = fconf["mode"]
    df_filename = fconf["input_file"]
    basename = os.path.basename(df_filename)
    sheet = fconf.get("sheet", 0)
    if "outfile" in fconf:
        outpath = fconf["outfile"]
    else:
        script_dir = os.path.dirname(df_filename)
        outpath = os.path.join(script_dir, "optimized_model.csv")
    if "nrows" in fconf:
        df = read_df(df_filename, nrows=int(fconf["nrows"]), sheet_name=sheet)
    elif "nsample" in fconf:
        df = read_df(df_filename, sheet_name=sheet)
        df = df[df[fconf["sample_id"]] == fconf["nsample"]]
    else:
        df = read_df(df_filename, sheet_name=sheet)
    LTE = False
    if fconf.get("lte") =="True":
        LTE = True
    fit_conf = OptimizationConfig(
        float(fconf.get("exclusion_threshold", 0.15)),
        LTE,
        json.loads(fconf.get("anomal_elements", '["La", "Ce", "Eu", "Gd", "Lu"]')),
        algorithm=fconf.get("algorithm", None)
    )
    df = enrich_df_with_metadata(df, fconf["sample_id"], df_filename, "_".join(basename.split(".")[:-1]))
    fit_model(df, fconf["norm_to"], fit_conf, outpath)
    print("Done", f"in {time() - start} seconds")


# pandarallel does not work inside the DataSpell interactive shell.
# However, it can be run as script on a normal commandline.
if __name__ == "__main__":
    main()

#%%
