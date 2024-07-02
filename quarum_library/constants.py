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

DEGREE = 5
LAMBDAS = [f"λ{i}" for i in range(DEGREE)]
TETRADEN = ["τ0", "τ1", "τ2", "τ3"]

# df column names
ANOMALIES = "anomalies"
DEFAULT_UNIT="ppm"
SAMPLE_ID = "Sample_ID" # GEOTRACES data
# SAMPLE_ID = "sample_ID" # TRAM data

REE_FULL_NAME = {
"La": "Lanthanum",
"Ce": "Cerium",
"Pr": "Praseodymium",
"Nd":"Neodymium",
"Sm":"Samarium",
"Eu":"Europium",
"Gd":"Gadolinium",
"Tb":"Terbium",
"Dy":"Dysprosium",
"Ho":"Holmium",
"Er":"Erbium",
"Tm":"Thulium",
"Yb":"Ytterbium",
"Lu":"Lutetium",
"Y":"Yttrium"}