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
import pandas
import pandas as pd
from numpy.polynomial.polynomial import polyfit, polyval
import pyrolite.geochem
from pyrolite.geochem.ind import get_ionic_radii

from quarum_library.constants import ANOMALIES, DEFAULT_UNIT
from quarum_library.fitting.fitting_infra import fitting_anomaly_detection
from quarum_library.util import get_renaming

DEGREE = 3

def polynomial_fit(simulated_data, elements=pyrolite.geochem.REE(), allElements=pyrolite.geochem.REE()):
    radii = get_ionic_radii(list(allElements), 3, 8)
    if simulated_data.ndim == 1:
        if isinstance(simulated_data, pandas.Series):
            excluded_nan = [simulated_data.index[x] for x in np.argwhere(~np.isfinite(simulated_data))]
            filter = [True if s in elements and s not in excluded_nan else False for s in allElements]
        else:
            excluded_nan = [x for y in np.argwhere(~np.isfinite(simulated_data)) for x in y]
            filter = [True if s in elements and i not in excluded_nan else False for i, s in enumerate(allElements)]
        radii_np = np.array(radii)[filter]
        fitted = polyfit(radii_np, simulated_data[filter], deg=DEGREE)
        model = polyval(radii, fitted)
        return model
    else:
        return np.array([polynomial_fit(x, elements, allElements) for x in simulated_data.T])




def polyfit_one_row(
        input_data,
        norm_to=None,
        units="ppm",
        log=True,
        charge=3,
        coordination=8,
        anomalies=[]
):
    # update_database()
    REE_df = input_data.pyrochem.REE.astype(float)
    if norm_to:
        REE_df = REE_df.pyrochem.normalize_to(norm_to, units=units)
    if log:
        REE_df = np.log(REE_df)
    selected = [x for x in pyrolite.geochem.REE() if not x in anomalies]
    model = polynomial_fit(REE_df, selected)
    if log:
        model = np.exp(model)
    model = pd.Series(model, index=pyrolite.geochem.REE())
    if norm_to:
        model = model.pyrochem.denormalize_from(norm_to)
    return model




def compute_one_iteration(series, norm, units, opt_conf, anomalies, log):
    assert isinstance(series, pd.Series)
    try:
        nseries = series.pyrochem.normalize_to(norm, units=units)
        model_generated = polyfit_one_row(nseries, None, units, log, anomalies=anomalies)
        dist, abs_dist, quotient, max_dist_element = fitting_anomaly_detection(nseries, model_generated, anomalies,
                                                                               opt_conf.allowed_anomalies)
        model_generated = model_generated.pyrochem.denormalize_from(norm, units=units)
        return quotient, max_dist_element, model_generated, dist, abs_dist
    except ValueError as e:
        print(e)
        return 0.0, "", None, 4000, 4000


def _optimize_model_per_series(series, norm, opt_conf, units=DEFAULT_UNIT, log=True):
    anomalies = []
    if series.get("exclude", None) and str(series["exclude"]) != "nan":
        anomalies.extend(str(series.get("exclude")).replace(" ", "").split(","))
        _, _, model, dist, abs_dist = compute_one_iteration(
            series, norm, units, opt_conf, anomalies, log)
        model[ANOMALIES] = anomalies
        return model
    else:
        quotient, max_dist_element, model, dist, abs_dist = compute_one_iteration(series, norm, units,
                                                                                  opt_conf, anomalies, log)
        while (quotient > opt_conf.threshold and len(anomalies) <= 6):
            anomalies.append(max_dist_element)
            quotient, max_dist_element, model, dist, abs_dist = compute_one_iteration(
                series, norm, units, opt_conf, anomalies, log)
        _, _, model, dist, abs_dist = compute_one_iteration(
            series, norm, units, opt_conf, anomalies, log)
        if model is not None:
            model[ANOMALIES] = anomalies
            return model
        else:
            return None


def polyfit_model(real_data, norm, opt_conf, units=DEFAULT_UNIT, log=True):
    # Fitting does not work if 0 are in the data series
    real_data.replace(0, np.nan, inplace=True)
    return real_data.parallel_apply(lambda x: _optimize_model_per_series(x, norm, opt_conf, units, log), axis=1)

