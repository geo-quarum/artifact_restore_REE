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
import numpy
import os
import gc

from numpy import nanmin, nanmax
from quarum_library.fitting.fitting_infra import OptimizationConfig
from quarum_library.util import read_df
from quarum_library.fitting.lambdatau import optimize_model
from pandarallel import pandarallel
import matplotlib.pyplot as plt
from pyrolite.geochem.ind import REE

def make_REE_pattern_plot_series(
        real_data,
        model1,
        model2,
        source_data_name,
        PATCHES,
        PATCHES_LIST,
        COMBINED,
        COMBINED_ONLY,
        MEASURED_ONLY,
        out_folder=None,
        normto="Chondrite_PON",
        units="ppm",
        suffix=None,
        SAMPLE="sample",
        image_save_format=".png"
):
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

    REEpos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    REEelem = REE(dropPm=False)
    f1_ax1.set_xlim(0.1, 15.9)
    f1_ax1.set_xticks(range(1, 16))
    f1_ax1.set_xticklabels(REEelem)
    f1_ax1.set_yscale("log")
    f1_ax1.set_ylabel("REE$_{sample}$/REE$_{C1}$")

    s_name = real_data[SAMPLE]
    try:
        s_type = real_data.loc["sample_type"]
    except KeyError:
        s_type = "unknown"
    try:
        s_source = real_data.loc["source"]
    except KeyError:
        s_source = "unknown"

    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]

    # Plot measurment
    real_REE = real_data.pyrochem.REE
    real_REE = real_REE.pyrochem.normalize_to(normto, units)
    model_REE = model1.pyrochem.REE
    model_REE = model_REE.pyrochem.normalize_to(normto, units)
    model2_REE = model2.pyrochem.REE.pyrochem.normalize_to(normto, units)

    xplot = REEpos
    yplot_real = real_REE.loc["La":"Lu"]
    yplot_model = model_REE.loc["La":"Lu"]
    yplot_model2 = model2_REE.loc["La": "Lu"]

    if COMBINED_ONLY == True:
        real_model_combined = yplot_real.copy()
        for i in real_model_combined.index:
            if numpy.isnan(real_model_combined.loc[i]):
                real_model_combined.loc[i] = yplot_model.loc[i]
        f1_ax1.plot(xplot, real_model_combined, marker="o", markersize=7, color="black",
                    markerfacecolor="white", label="combined",
                    zorder=1)
        f1_ax1.scatter(xplot, yplot_real, s=10, c="black", label="measured", zorder=2)
        COMBINED = False
    elif MEASURED_ONLY == True:
        f1_ax1.plot(xplot, yplot_real, marker="o", markersize=7, color="black", label="measured")
    else:
        print(yplot_real)
        print(len(xplot), yplot_real.shape, yplot_model.shape)
        f1_ax1.scatter(xplot, yplot_real, s=20, label="measured", c="black", zorder=3)
        f1_ax1.plot(xplot, yplot_model, marker="o", markersize=7, alpha=0.5, color="teal",
                    label="model1", zorder=2)
        f1_ax1.plot(xplot, yplot_model2, marker="o", markersize=7, alpha=0.5, color="red",
                    label="model2", zorder=2)
        f1_ax1.scatter([2], [yplot_real["La"] + 0.33*(yplot_real["Nd"]-yplot_real["La"])], s=20, label="interpolated Ce", c="red", zorder=3)


    if COMBINED == True:
        real_model_combined = yplot_real.copy()
        for i in real_model_combined.index:
            if numpy.isnan(real_model_combined.loc[i]):
                real_model_combined.loc[i] = yplot_model.loc[i]
        f1_ax1.plot(xplot, real_model_combined, color="purple", label="combined")

    if PATCHES == True:
        patch_list = PATCHES_LIST
        for i in patch_list:
            xpos = model_REE.index.get_loc(i) + 1
            yplot = model_REE[i]
            f1_ax1.scatter(xpos, yplot, edgecolor="red", facecolor="white", marker="s", s=300, lw=2,
                           zorder=1)
        f1_ax1.scatter(0, 0, edgecolor="red", facecolor="white", marker="s", s=80, lw=2, zorder=1,
                       label="removed REE")

    # calculate deviation between model and measured data
    l1 = REE()
    # anomalies are removed from calculations
    l2 = model1["anomalies"]
    avg_dev_list = [x for x in l1 if x not in l2]
    sum_dev = 0
    # also measured data nan values are excluded from calculation
    number_nan = 0
    for i in avg_dev_list:
        if numpy.isfinite(yplot_real.loc[i]):
            aux = numpy.absolute((yplot_model.loc[i] / yplot_real.loc[i]) - 1)
            sum_dev = sum_dev + aux
        else:
            number_nan = number_nan + 1
    avg_dev = sum_dev / (len(avg_dev_list) - number_nan)

    anomalies_for_text = str(model1['anomalies']).replace("'", "")

    ymin = nanmin(real_REE) * 0.1
    ymax = nanmax(real_REE) * 10
    f1_ax1.set_ylim(ymin, ymax)

    plt.text(
        0.1,
        0.85,
        f"excluded anomalies: {anomalies_for_text}"
        + "\n\nsample: "
        + str(s_name)
        + "\nsample type: "
        + str(s_type)
        + "\ndata source: "
        + str(s_source)
        + "\n\naverage deviation: "
        + str("{:.4f}".format(avg_dev)),
        ha="left",
        va="center",
        fontsize=8,
        transform=f1_ax1.transAxes,
        )
    # f1_ax1.scatter(0, 0, marker="s", facecolor="white", edgecolor="red", s=75, lw=2, label="missing REE")
    legend = f1_ax1.legend(
        loc=2,
        bbox_to_anchor=(0.7, 0.95),
        edgecolor="inherit",
        scatteryoffsets=[0.5],
        fontsize=8,
    )

    if not out_folder:
        try:
            out_folder = os.path.dirname(__file__)
        except NameError:
            out_folder = os.getcwd()

    sample_file_name = (
        f"REE_model_measured_pattern_{str(source_data_name)}_{str(real_data[SAMPLE])}"
    )
    if suffix:
        sample_file_name += f"_{suffix}"
    sample_file_name += image_save_format
    os.makedirs(out_folder, exist_ok=True)
    plt.savefig(os.path.join(out_folder, sample_file_name))

    plt.close(fig1)
    del fig1
    gc.collect()

def make_fig10():
    pandarallel.initialize()
    raw_data = read_df("data/zindleretal79.csv")
    raw_data = raw_data[raw_data["sample"]=="RE 78"]
    raw_norm = raw_data.pyrochem.normalize_to("Chondrite_PON")
    print(raw_data)
    fit_conf = OptimizationConfig(
        0.1,
        False,
        ["Ce"],
        None
    )
    fit_conf2 = OptimizationConfig(
        0.1,
        False,
        [],
        None
    )
    optimized_model = optimize_model(raw_norm, "Chondrite_PON", fit_conf)
    raw_norm["Ce"] = raw_norm["La"] + 0.33*(raw_norm["Nd"]-raw_norm["La"])
    optimized_model2 = optimize_model(raw_norm, "Chondrite_PON", fit_conf2)
    print(optimized_model2)
    print(optimized_model)
    Pm_loc = optimized_model.columns.get_loc("Nd") + 1
    optimized_model.insert(loc=Pm_loc, column="Pm", value=numpy.nan)
    optimized_model2.insert(loc=Pm_loc, column="Pm", value=numpy.nan)
    optimized_model = optimized_model
    optimized_model2 = optimized_model2
    raw_data.insert(loc=Pm_loc, column="Pm", value=numpy.nan)

    make_REE_pattern_plot_series(raw_data.squeeze("rows"), optimized_model.squeeze("rows"), optimized_model2.squeeze("rows"), "Zindler et al. (1979)",
                                 False, None, False, False, False, "output/fig10")

if __name__ == "__main__":
    make_fig10()
