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

import math

import numpy as np
import pandas as pd
from pyrolite.geochem import REE
from pyrolite.geochem.ind import get_ionic_radii
from pyrolite.util.lambdas import _get_params
from pyrolite.util.lambdas.eval import get_function_components
from pyrolite.util.lambdas.opt import lambdas_optimize
from pyrolite.util.lambdas.oneill import lambdas_ONeill2016
from pyrolite.geochem.norm import update_database

from quarum_library.constants import LAMBDAS, TETRADEN, DEGREE, ANOMALIES, SAMPLE_ID, DEFAULT_UNIT
from quarum_library.fitting.fitting_infra import measured_model_diff, fitting_anomaly_detection
from quarum_library.util import get_renaming

def lambdatau_fitting_one_row(
        input_data,
        norm_to=None,
        LTE=False,
        log=True,
        charge=3,
        coordination=8,
        algorithm="opt"
):
    # update_database()
    REE_df = input_data.pyrochem.REE.astype(float)
    if norm_to:
        REE_df = REE_df.pyrochem.normalize_to(norm_to)
    if log:
        REE_df = np.log(REE_df)
    radii = get_ionic_radii(list(REE_df.index), charge, coordination)
    # quantify λτ parameter
    # use exclude_ele and LTE in the input to modify output
    if "LTE" in input_data.index:
        if str(input_data["LTE"]).lower() == "on":
            LTE = True
        else:
            LTE = False

    if algorithm.lower() == 'oneil' and not LTE:
        lambdatau_df = lambdas_ONeill2016(REE_df, radii,params=_get_params(None, DEGREE), add_X2=False, add_uncertainties=False)
    else:
        lambdatau_df = lambdas_optimize(
            REE_df,
            radii,
            params=_get_params(None, DEGREE),
            fit_tetrads=LTE,
            # We might want to set this true again in the future, but there is some problem with pyrolite
            add_uncertainties=False,
            add_X2=False
        )
    lambdatau_df["LTE"] = LTE
    return lambdatau_df


def lambdatau_model_one_row(series, normed_to, unit="ppm", charge=3, coordinates=8):
    lambdas = series.loc[LAMBDAS]
    tetraden = None
    if "τ0" in series.index:
        tetraden = series.loc[TETRADEN]

    radii = get_ionic_radii(REE(), coordination=8, charge=3)
    if "LTE" in series.index and series["LTE"]:
        all_data = pd.concat([lambdas, tetraden])
        _fit_tetrades = True
    elif "LTE" not in series.index and tetraden is not None:
        all_data = pd.concat([lambdas, tetraden])
        _fit_tetrades = True
    else:
        all_data = lambdas
        _fit_tetrades = False
    _, _, components = get_function_components(
        radii, fit_tetrads=_fit_tetrades
    )
    model = all_data.to_numpy() @ np.array(components)
    model = model.astype(float)
    pattern = pd.Series(np.exp(model), index=REE())
    # de-normalize the pattern_df
    return pattern.pyrochem.denormalize_from(normed_to, unit)


def compute_one_iteration(series, norm, units, opt_conf, anomalies):
    assert isinstance(series, pd.Series)
    try:
        model = lambdatau_fitting_one_row(series[[x for x in series.index if x not in anomalies]], LTE=opt_conf.LTE,
                                          algorithm=opt_conf.algorithm)
        model_generated = lambdatau_model_one_row(model, norm)
        model_gen_normed = model_generated.pyrochem.normalize_to(norm, units=units)
        dist, abs_dist, quotient, max_dist_element = fitting_anomaly_detection(series, model_gen_normed, anomalies,
                                                                               opt_conf.allowed_anomalies)
        return quotient, max_dist_element, model, model_generated, dist, abs_dist
    except ValueError as e:
        print(e)
        return 0.0, "", None, None, 4000, 4000

def _optimize_model_per_series(series, norm, opt_conf, units=DEFAULT_UNIT):
    anomalies = []
    if series.get("exclude", None) and str(series["exclude"]) != "nan":
        anomalies.extend(str(series.get("exclude")).replace(" ", "").split(","))
        _, _, model, model_generated, dist, abs_dist = compute_one_iteration(
            series, norm, units, opt_conf, anomalies)
        model[ANOMALIES] = anomalies
        return pd.concat([model, model_generated])
    else:
        quotient, max_dist_element, model, model_generated, dist, abs_dist = compute_one_iteration(series, norm, units,
                                                                                                   opt_conf, anomalies)
        while (quotient > opt_conf.threshold and len(anomalies) <= 6):
            anomalies.append(max_dist_element)
            quotient, max_dist_element, model, model_generated, dist, abs_dist = compute_one_iteration(
                series, norm, units, opt_conf, anomalies)
        _, _, model, model_generated, dist, abs_dist = compute_one_iteration(
            series, norm, units, opt_conf, anomalies)
        if model is not None and model_generated is not None:
            model[ANOMALIES] = anomalies
            return pd.concat([model, model_generated])
        else:
            return None


def optimize_model(real_data, norm, opt_conf, units=DEFAULT_UNIT):
    #Fitting does not work if 0 are in the data series 
    real_data.replace(0, np.nan, inplace=True)
    return real_data.parallel_apply(lambda x: _optimize_model_per_series(x, norm, opt_conf, units), axis=1)

