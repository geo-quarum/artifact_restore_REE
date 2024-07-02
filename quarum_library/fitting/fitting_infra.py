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
from quarum_library.util import get_renaming
from pyrolite.geochem import REE
import pandas as pd

class OptimizationConfig(object):
    def __init__(self, threshold, LTE, allowed_anomalies, algorithm=None):
        self.threshold = threshold
        self.LTE = LTE
        self.allowed_anomalies = allowed_anomalies
        self.algorithm = algorithm if algorithm is not None else 'opt'

def measured_model_diff(series, model_prefix="", elements=REE()):
    renaming = get_renaming()
    abs_diff = 0
    diff = 0
    for i in elements:
        if not math.isinf(series.loc[renaming[i]]):
            dist = series.loc[i] - series.loc[renaming[i]]
            diff += dist
            abs_diff += abs(dist)
        else:
            diff = float("inf")
            abs_diff = float("inf")
    return pd.Series(
        [abs_diff, diff], index=["abs_diff_" + model_prefix, "diff_" + model_prefix]
    )

def fitting_anomaly_detection(series, model_gen, anomalies, allowed_anomalies):
    renaming = get_renaming()
    abs_diff = 0
    diff = 0
    highest_quotient = 0
    highest_abs_diff_element = None
    for i in model_gen.index:
        if not math.isinf(model_gen.loc[i]) and i in series.index:
            dist = series.loc[i] - model_gen.loc[i]
            diff += dist
            abs_diff += abs(dist)
            quotient = abs(1 - model_gen.loc[i] / series.loc[i])
            if quotient > highest_quotient and i not in anomalies and i in allowed_anomalies:
                highest_quotient = quotient
                highest_abs_diff_element = i
        else:
            diff = float("inf")
            abs_diff = float("inf")
    return diff, abs_diff, highest_quotient, highest_abs_diff_element