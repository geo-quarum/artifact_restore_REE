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
from quarum_library.util import read_df
from pyrolite.geochem import REE
from textwrap import dedent
import os
import numpy


def make_columns_none(df, columns):
    for c in columns:
        df.loc[:, c] = numpy.nan
    return df


def keep_only_colmuns(df, columns):
    dropColumns = []
    for i in REE():
        if i not in columns:
            dropColumns.append(i)
    return make_columns_none(df, dropColumns)

original_data = read_df("./data/PetDB240122_mafic_ultramafic_cleared.csv")
inaa = original_data.copy()
inaa = make_columns_none(inaa, ["Pr", "Dy", "Ho", "Er"])
inaa.to_csv("./data/PetDB240122_mafic_ultramafic_cleared_naa.csv", index=False)

inaa_stepwise = original_data.copy()
inaa_stepwise = make_columns_none(inaa_stepwise,["Pr"])
inaa_stepwise.to_csv("./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr.csv", index=False)
inaa_stepwise = make_columns_none(inaa_stepwise,["Dy"])
inaa_stepwise.to_csv("./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy.csv", index=False)
inaa_stepwise = make_columns_none(inaa_stepwise,["Ho"])
inaa_stepwise.to_csv("./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy_Ho.csv", index=False)

id = original_data.copy()
id = make_columns_none(id, ["Pr", "Tb", "Ho", "Tm"])
id.to_csv("./data/PetDB240122_mafic_ultramafic_cleared_id.csv", index=False)

prefix = "./PetDB240122_mafic_ultramafic_cleared"
configs = "./data/configs"
os.makedirs(configs, exist_ok=True)

configs_required = [("./data/PetDB240122_mafic_ultramafic_cleared_id.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_id_model_lambda.csv"),
                    ("./data/PetDB240122_mafic_ultramafic_cleared_id.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_id_model_polyfit.csv"),
                    ("./data/PetDB240122_mafic_ultramafic_cleared_naa.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_naa_model_lambda.csv"),
                    ("./data/PetDB240122_mafic_ultramafic_cleared_naa.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_naa_model_polyfit.csv"),
                    ("./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr_model_lambda.csv"),
                    ("./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy_model_lambda.csv"),
                    ("./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy_Ho.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy_Ho_model_lambda.csv"),
                    ("./data/PetDB240122_mafic_ultramafic_cleared.csv",
                     "./data/PetDB240122_mafic_ultramafic_cleared_model_lambda.csv")]
for outname, modelfile in configs_required:
    config = dedent(f"""\
    [fit]
    mode = optimize_fitting
    input_file = {os.path.abspath(outname)}
    outfile = {os.path.abspath(modelfile)}
    norm_to = Chondrite_PON
    sample_id = sample
    lte = False
    anomal_elements=["Ce","Eu"]
    """)
    configfile = os.path.join(configs, os.path.basename(modelfile).replace(".csv", ".ini"))
    with open(configfile, "w") as outfile:
        print(config, file=outfile)
    #print(os.path.abspath(configfile))

#%%
