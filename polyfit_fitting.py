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

import numpy
import pandas
from pyrolite.geochem import REE

from quarum_library.exceptions import NormalizationError
from quarum_library.fitting.fitting_infra import OptimizationConfig
from quarum_library.fitting.polyfit import polyfit_model
from quarum_library.simulations.monte_carlo import polynomial_fit
from quarum_library.util import read_df
from pandarallel import pandarallel

def fit_model(dataframe, norm_to, sample_column, conf, units="ppm"):
    if norm_to is None:
        raise NormalizationError()
    pandarallel.initialize(progress_bar=True)
    optimized_model = polyfit_model(dataframe, norm_to, conf)
    return pandas.concat([dataframe[sample_column], optimized_model], axis=1)

def _write_default_config():
    """
    Writes the default config for this file.
    """
    config = ConfigParser()
    config["fit"] = {
        "input_file": "The file that contains the dataframe",
        "outfile": "This file concludes the final models constructed with polyfit",
        "norm_to": "This describes the reference standard used for normalization",
        "sample_id": "The column with the sample ID.",
        "anomal_elements": "Optional, default is [La, Ce, Eu, Gd, Lu]",
        "threshold": "Threshold for classifying something as anomaly, default is [0.1]",
        "nrows": "This is passed to the reader for processing only parts of the dataframe",
        "nsample": "This is passed to the reader for processing only a specific sample"
    }
    conf_path = os.path.join(os.path.dirname(__file__), "configs", "default_polyfit_fitting.ini")
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
    fconf = config["fit"]
    df_filename = fconf["input_file"]
    basename = os.path.basename(df_filename)
    sheet = fconf.get("sheet", 0)
    if "outfile" in fconf:
        outfile = fconf["outfile"]
    else:
        script_dir = os.path.dirname(df_filename)
        outfile = os.path.join(script_dir, "output", "polyfit_model.csv")
    if "nrows" in fconf:
        df = read_df(df_filename, nrows=int(fconf["nrows"]), sheet_name=sheet)
    elif "nsample" in fconf:
        df = read_df(df_filename, sheet_name=sheet)
        df = df[df[fconf["sample_id"]] == fconf["nsample"]]
    else:
        df = read_df(df_filename, sheet_name=sheet)
    excludedElements =  json.loads(fconf.get("anomal_elements", '["La", "Ce", "Eu", "Gd", "Lu"]'))
    threshold = fconf.getfloat("threshold", 0.1)
    opt_conf = OptimizationConfig(threshold, None, excludedElements)
    model = fit_model(df, fconf["norm_to"], fconf["sample_id"], opt_conf)
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    model.to_csv(outfile, index=False)
    print("Done", f"in {time() - start} seconds")


if __name__ == "__main__":
    main()

#%%
