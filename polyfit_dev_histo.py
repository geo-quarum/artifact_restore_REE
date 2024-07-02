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
import os
import gc
import pandas as pd
from pyrolite.geochem.ind import REE, REY, get_ionic_radii
import matplotlib.pyplot as plt
import numpy as np
from seaborn import histplot, kdeplot

input_file = "/Users/dernst/Library/CloudStorage/OneDrive-ConstructorUniversity/author_paper/reconstruct_ree/isotope_dilution_ree_data/output_models/PetDB240122_mafic_ultramafic_cleared_POLYFIT_MEAS_MODEL.csv"
outfile = "/Users/dernst/Library/CloudStorage/OneDrive-ConstructorUniversity/author_paper/reconstruct_ree/isotope_dilution_ree_data/output_models/REE_comparison_plots/PetDB240122_polyfit"

input_df = pd.read_csv(input_file)

plot_data = input_df.loc[:, "La":"Lu"].copy()
plot_data = plot_data.to_numpy()
plot_data = plot_data.flatten()
plot_data = plot_data[~np.isnan(plot_data)]



# prepare the plot
fig1 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
gs1 = fig1.add_gridspec(1, 1, hspace=0.2, wspace=0.3)
f1_ax1 = fig1.add_subplot(gs1[0, 0])

# adjust the margine for the plots / subplots
fig1.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)

plt.style.use("default")
plt.rc("lines", markersize=2)
plt.rc("axes", linewidth=2, titlepad=20)
plt.rc("font", size=14)
props = dict(boxstyle="round", facecolor="white", linewidth=0)

f1_ax1.axes.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    width=2,
    length=8,
)
f1_ax1.axes.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    width=1,
    length=4,
)
# close and reload again so that the settings are also applied to the first plot
plt.close(fig1)
plt.style.use("default")
plt.rc("lines", markersize=2)
plt.rc("axes", linewidth=2, titlepad=20)
plt.rc("font", size=14)
props = dict(boxstyle="round", facecolor="white", linewidth=0)
fig1 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
gs1 = fig1.add_gridspec(1, 1, hspace=0.2, wspace=0.3)
f1_ax1 = fig1.add_subplot(gs1[0, 0])
# adjust the margine for the plots / subplots
fig1.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)
f1_ax1.axes.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    width=2,
    length=8,
)
f1_ax1.axes.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    width=1,
    length=4,
)

histplot(data=plot_data, ax=f1_ax1, stat="density", color="teal", zorder=1)
# kdeplot(data=plot_data, ax=f1_ax1, color="orange", zorder=2)

p05 = np.nanquantile(plot_data, 0.05)
p95 = np.nanquantile(plot_data, 0.95)
p50 = np.nanquantile(plot_data, 0.5)
num = np.count_nonzero(plot_data)
SD = "{:.1f}".format(np.nanstd(plot_data))
f1_ax1.text(0.7,
            0.05,
            "p0.05 = "
            + str("{:.1f}".format(p05))
            + "\np0.95 = "
            + str("{:.1f}".format(p95))
            + "\nmedian = "
            + str("{:.1f}".format(p50))
            + "\nSD = "
            + SD
            + "\nn = "
            + str(num),
            fontsize=10,
            transform=f1_ax1.transAxes
            )


f1_ax1.set_xlim(-38, 38)
f1_ax1.set_xlabel("Deviation from measured data in [%]")

fig1.savefig(os.path.join(outfile, "dev_hist_polyfit.png"))
plt.close(fig1)
del fig1
gc.collect()