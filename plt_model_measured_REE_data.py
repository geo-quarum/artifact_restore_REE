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
import gc
import os
import re
import sys
import numpy
import json

import matplotlib.pyplot as plt
import pandas as pd
from pyrolite.geochem.ind import REE, REY
from pandarallel import pandarallel
from quarum_library.util import read_df, ANOMALIES
from time import time
from configparser import ConfigParser
from numpy import nanmin, nanmax, nanstd
from seaborn import histplot, kdeplot
from scipy.stats import skew, kurtosis
import statsmodels.stats.api as sms

def make_REE_pattern_plot_series(
        real_data,
        models,
        source_data_name,
        PATCHES,
        PATCHES_LIST,
        COMBINED,
        COMBINED_ONLY,
        MEASURED_ONLY,
        out_folder=None,
        # normto="EUS_Bau2018",
        normto="Chondrite_PON",
        units="ppm",
        suffix=None
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
    model_REE = models.pyrochem.REE
    model_REE = model_REE.pyrochem.normalize_to(normto, units)

    xplot = REEpos
    yplot_real = real_REE.loc["La":"Lu"]
    yplot_model = model_REE.loc["La":"Lu"]

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
        f1_ax1.scatter(xplot, yplot_real, s=20, label="measured", c="black", zorder=3)
        f1_ax1.plot(xplot, yplot_model, marker="o", markersize=7, alpha=0.7, color="teal",
                    label="model", zorder=2)

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
    l2 = models["anomalies"]
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

    anomalies_for_text = models['anomalies'].replace("'", "")

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

    # sort high deviation samples to extra folder
    # print(avg_dev, avg_dev > DEV_THRESHOLD, DEV_THRESHOLD)
    if avg_dev > DEV_THRESHOLD:
        results_dir = os.path.join(out_folder, source_data_name, f"high_deviation({DEV_THRESHOLD})")
    else:
        results_dir = os.path.join(out_folder, source_data_name,
                                   f"low_deviation({DEV_THRESHOLD})")
    os.makedirs(results_dir, exist_ok=True)
    plt.savefig(os.path.join(results_dir, sample_file_name))

    plt.close(fig1)
    del fig1
    gc.collect()


def make_histogram_plot(
        real_data,
        model_data,
        source_data_name,
        PATCHES,
        PATCHES_LIST,
        out_folder,
        xmin=-38,
        xmax=38,
        normto="Chondrite_PON",
        units="ppm"
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

    # open figure 2 for direct comparison plot
    fig2 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
    gs2 = fig2.add_gridspec(1, 1, hspace=0.2, wspace=0.3)
    f2_ax1 = fig2.add_subplot(gs2[0, 0])

    # adjust the margine for the plots / subplots
    fig2.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)

    f2_ax1.axes.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        width=2,
        length=8,
    )
    f2_ax1.axes.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        width=1,
        length=4,
    )
    f2_ax1.set_xlim(0.1, 15.9)
    f2_ax1.set_xticks(range(1, 16))
    f2_ax1.set_ylim(-28, 28)
    REEpos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    REEelem = REE()
    REEelem.insert(4, "Pm")
    f2_ax1.set_xticklabels(REEelem)

    # open figure 3 for kde of removed elements
    fig3 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
    gs1 = fig3.add_gridspec(1, 1, hspace=0.2, wspace=0.3)
    f3_ax1 = fig3.add_subplot(gs1[0, 0])

    # adjust the margine for the plots / subplots
    fig3.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)

    f3_ax1.axes.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        width=2,
        length=8,
    )
    f3_ax1.axes.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        width=1,
        length=4,
    )

    # prepare the data
    real_REE = real_data.pyrochem.REE
    real_REE = real_REE.pyrochem.normalize_to(normto, units)
    try:
        real_REE.insert(loc=4, column="Pm", value=numpy.nan)
    except ValueError:
        pass

    # remove anomaly data
    s_list = model_data[SAMPLE]
    histo_model_data = model_data.copy()
    anomalies_list = []
    anomalies_count = 0
    for s in s_list:
        aux = model_data.loc[model_data[SAMPLE] == s].copy()
        anomalies = aux["anomalies"]
        REElist = REE()
        for i in REElist:
            if i in str(anomalies):
                histo_model_data.loc[histo_model_data[SAMPLE] == s, i] = numpy.nan
                anomalies_list.append(i)
                anomalies_count = anomalies_count + 1

    model_REE = histo_model_data.pyrochem.REE
    model_REE = model_REE.pyrochem.normalize_to(normto, units)
    try:
        model_REE.insert(loc=4, column="Pm", value=numpy.nan)
    except ValueError:
        pass

    # calculate the final data
    plot_data_fig1 = ((model_REE / real_REE) - 1) * 100

    # extract all deviation data including sample names for filtering
    export_deviation = plot_data_fig1.copy()
    export_deviation["sample"] = histo_model_data["sample"]

    # copy data for figure 2
    plot_data_fig2 = plot_data_fig1.copy()
    # modify data for histogram in figure 1
    plot_data_fig1 = plot_data_fig1.to_numpy()
    plot_data_fig1 = plot_data_fig1[~numpy.isnan(plot_data_fig1)]
    plot_data_fig1 = plot_data_fig1.flatten()

    # plot data in figure 1
    histplot(data=plot_data_fig1, ax=f1_ax1, color="teal", stat="density", zorder=1)
    kdeplot(data=plot_data_fig1, ax=f1_ax1, color="orange", zorder=2)

    f1_ax1.set_xlim(xmin, xmax)
    f1_ax1.set_xlabel("Deviation from measured data in [%]")

    # plot data in figure 2
    marker_size = 200
    REE_col = plot_data_fig2.columns
    pos = 1
    for i in REE_col:
        yplot1 = plot_data_fig2[i]
        xplot = [pos] * len(yplot1)
        f2_ax1.scatter(xplot, yplot1, marker="_", c="teal", s=marker_size, alpha=0.6, zorder=2)
        pos = pos + 1

    if PATCHES == True:
        patch_list = PATCHES_LIST
        for i in patch_list:
            xpos = real_REE.columns.get_loc(i) + 0.5
            ymin = numpy.nanmin(plot_data_fig2[i]) - 1
            ymax = numpy.nanmax(plot_data_fig2[i]) + 1
            yheight = ymax - ymin
            REE_patch = plt.Rectangle((xpos, ymin), 1, yheight, facecolor=(1, 1, 1, 0),
                                      edgecolor="red", zorder=2)
            f2_ax1.add_patch(REE_patch)

    f2_ax1.hlines(0, 20, 0, color="black", lw=2, zorder=1)
    f2_ax1.fill_between((0, 16), -10, 10, alpha=0.3, color="gray", zorder=1)
    f2_ax1.fill_between((0, 16), -5, 5, alpha=0.3, linewidth=0, color="gray", zorder=1)
    f2_ax1.set_ylabel("Deviation from measured data in [%]")

    anomalies = set(anomalies_list)
    anomalies = str(anomalies).replace("{", "")
    anomalies = anomalies.replace("}", "")
    anomalies = anomalies.replace("'", "")
    num = numpy.count_nonzero(~numpy.isnan(plot_data_fig1))

    # skewness and kurtosis of KDE
    skew_val = "{:.2f}".format(skew(plot_data_fig1))
    kurtosis_val = "{:.2f}".format(kurtosis(plot_data_fig1))
    f1_ax1.text(0.6,
                0.5,
                f"skewness: {skew_val}"
                + f"\nkurtosis: {kurtosis_val}",
                fontsize=8,
                color="orange",
                transform=f1_ax1.transAxes
                )

    # percentiles, mean and SD
    p05 = numpy.nanquantile(plot_data_fig1, 0.05)
    p95 = numpy.nanquantile(plot_data_fig1, 0.95)
    p50 = numpy.nanquantile(plot_data_fig1, 0.5)
    avg = numpy.nanmean(plot_data_fig1)
    SD_val = "{:.1f}".format(nanstd(plot_data_fig1))
    f1_ax1.text(0.75,
                0.05,
                "p0.05 = "
                + str("{:.1f}".format(p05))
                + "\np0.95 = "
                + str("{:.1f}".format(p95))
                + "\nmean = "
                + str("{:.1f}".format(avg))
                + "\nSD = "
                + SD_val
                + "\nn = "
                + str(num),
                fontsize=10,
                transform=f1_ax1.transAxes
                )
    f2_ax1.text(0.7,
                0.05,
                "p0.05 = "
                + str("{:.1f}".format(p05))
                + "\np0.95 = "
                + str("{:.1f}".format(p95))
                + "\nmean = "
                + str("{:.1f}".format(avg))
                + "\nSD = "
                + SD_val
                + "\nn = "
                + str(num),
                fontsize=10,
                transform=f2_ax1.transAxes
                )
    # percentiles for the deleted elements
    if PATCHES == True:
        patch_list = PATCHES_LIST
        real_del_REE_data = real_REE[patch_list].copy()
        model_del_REE_data = model_REE[patch_list].copy()
        del_data_fig1 = ((model_del_REE_data / real_del_REE_data) - 1) * 100
        del_data_fig1_flat = del_data_fig1.copy().to_numpy()
        del_data_fig1_flat = del_data_fig1_flat.flatten()
        p05_del = numpy.nanquantile(del_data_fig1_flat, 0.05)
        p95_del = numpy.nanquantile(del_data_fig1_flat, 0.95)
        p50_del = numpy.nanquantile(del_data_fig1_flat, 0.5)
        avg_del = numpy.nanmean(del_data_fig1_flat)
        num_del = numpy.count_nonzero(~numpy.isnan(del_data_fig1_flat))
        SD_del = "{:.1f}".format(nanstd(del_data_fig1_flat))

        f2_ax1.text(0.1,
                    0.05,
                    "p0.05 = "
                    + str("{:.1f}".format(p05_del))
                    + "\np0.95 = "
                    + str("{:.1f}".format(p95_del))
                    + "\nmean = "
                    + str("{:.1f}".format(avg_del))
                    + "\nSD = "
                    + SD_del
                    + "\nn = "
                    + str(num_del),
                    fontsize=10,
                    color="red",
                    transform=f2_ax1.transAxes
                    )

    f1_ax1.text(
        0.05,
        0.9,
        f"excluded anomalies:\n{anomalies} ({anomalies_count} times)",
        ha="left",
        va="center",
        fontsize=8,
        transform=f1_ax1.transAxes,
    )
    f2_ax1.text(
        0.05,
        0.9,
        f"excluded anomalies:\n{anomalies} ({anomalies_count} times)",
        ha="left",
        va="center",
        fontsize=8,
        transform=f2_ax1.transAxes,
    )

    # plot data in figure 3
    if PATCHES == True:
        histplot(data=del_data_fig1_flat, ax=f3_ax1, color="teal", stat="density", zorder=1)
        kdeplot(data=del_data_fig1_flat, ax=f3_ax1, color="orange", zorder=2)
        # modify axes
        f3_ax1.set_xlim(xmin, xmax)
        f3_ax1.set_xlabel("Deviation from measured data in [%]")
        # skewness and kurtosis of KDE
        skew_val_del = "{:.2f}".format(skew(del_data_fig1_flat))
        kurtosis_val_del = "{:.2f}".format(kurtosis(del_data_fig1_flat))
        f3_ax1.text(0.6,
                    0.5,
                    f"skewness: {skew_val_del}"
                    + f"\nkurtosis: {kurtosis_val_del}",
                    fontsize=8,
                    color="orange",
                    transform=f3_ax1.transAxes
                    )

        f3_ax1.text(
            0.05,
            0.9,
            f"KDE for: {PATCHES_LIST}",
            ha="left",
            va="center",
            fontsize=8,
            fontweight="bold",
            transform=f3_ax1.transAxes,
        )
        f3_ax1.text(
            0.05,
            0.85,
            f"excluded anomalies:\n{anomalies} ({anomalies_count} times)",
            ha="left",
            va="center",
            fontsize=8,
            transform=f3_ax1.transAxes,
        )
        f3_ax1.text(0.7,
                    0.05,
                    "p0.05 = "
                    + str("{:.1f}".format(p05_del))
                    + "\np0.95 = "
                    + str("{:.1f}".format(p95_del))
                    + "\nmean = "
                    + str("{:.1f}".format(avg_del))
                    + "\nSD = "
                    + SD_del
                    + "\nn = "
                    + str(num_del),
                    fontsize=10,
                    color="black",
                    transform=f3_ax1.transAxes
                    )

    for individual_REE in REE():
        # for del_REE in PATCHES_LIST:
        fig4 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
        gs1 = fig4.add_gridspec(1, 1, hspace=0.2, wspace=0.3)
        f4_ax1 = fig4.add_subplot(gs1[0, 0])
        # adjust the margine for the plots / subplots
        fig4.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)
        f4_ax1.axes.tick_params(
            axis="both",
            which="major",
            direction="in",
            top=True,
            right=True,
            width=2,
            length=8,
        )
        f4_ax1.axes.tick_params(
            axis="both",
            which="minor",
            direction="in",
            top=True,
            right=True,
            width=1,
            length=4,
        )
        aux_real = real_REE[individual_REE].copy()
        aux_model = model_REE[individual_REE].copy()
        plot_data_fig4 = ((aux_model / aux_real) - 1) * 100
        plot_data_fig4 = plot_data_fig4.to_numpy()
        plot_data_fig4 = plot_data_fig4.flatten()
        p05_fig4 = numpy.nanquantile(plot_data_fig4, 0.05)
        p95_fig4 = numpy.nanquantile(plot_data_fig4, 0.95)
        p50_fig4 = numpy.nanquantile(plot_data_fig4, 0.5)
        avg_fig4 = numpy.nanmean(plot_data_fig4)
        num_fig4 = numpy.count_nonzero(~numpy.isnan(plot_data_fig4))
        SD_fig4 = "{:.1f}".format(nanstd(plot_data_fig4))
        skew_val_fig4 = "{:.2f}".format(skew(plot_data_fig4, nan_policy='omit'))
        kurtosis_val_fig4 = "{:.2f}".format(kurtosis(plot_data_fig4, nan_policy='omit'))
        # plot the data
        histplot(data=plot_data_fig4, ax=f4_ax1, color="teal", stat="density", zorder=1)
        kdeplot(data=plot_data_fig4, ax=f4_ax1, color="orange", zorder=2)
        f4_ax1.set_xlim(xmin, xmax)
        f4_ax1.set_xlabel("Deviation from measured data in [%]")
        f4_ax1.text(0.6,
                    0.5,
                    f"skewness: {skew_val_fig4}"
                    + f"\nkurtosis: {kurtosis_val_fig4}",
                    fontsize=8,
                    color="orange",
                    transform=f4_ax1.transAxes
                    )
        f4_ax1.text(
            0.05,
            0.9,
            f"KDE for: {individual_REE}",
            ha="left",
            va="center",
            fontsize=8,
            fontweight="bold",
            transform=f4_ax1.transAxes,
        )
        f4_ax1.text(
            0.05,
            0.85,
            f"excluded anomalies:\n{anomalies} ({anomalies_count} times)",
            ha="left",
            va="center",
            fontsize=8,
            transform=f4_ax1.transAxes,
        )
        f4_ax1.text(0.7,
                    0.05,
                    "p0.05 = "
                    + str("{:.1f}".format(p05_fig4))
                    + "\np0.95 = "
                    + str("{:.1f}".format(p95_fig4))
                    + "\nmean = "
                    + str("{:.1f}".format(avg_fig4))
                    + "\nSD = "
                    + SD_fig4
                    + "\nn = "
                    + str(num_fig4),
                    fontsize=10,
                    color="black",
                    transform=f4_ax1.transAxes
                    )
        # save figure
        if not out_folder:
            try:
                out_folder = os.path.dirname(__file__)
            except NameError:
                out_folder = os.getcwd()
        results_dir = os.path.join(out_folder, source_data_name)
        sample_file_name = (
            f"REE_model_measured_hist_{individual_REE}_{source_data_name}"
        )
        sample_file_name += image_save_format
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        fig4.savefig(os.path.join(results_dir, sample_file_name))

    # save the remaining figures above
    # save figure 1
    if not out_folder:
        try:
            out_folder = os.path.dirname(__file__)
        except NameError:
            out_folder = os.getcwd()
    results_dir = os.path.join(out_folder, source_data_name)
    sample_file_name = (
        f"REE_model_measured_hist_{source_data_name}"
    )
    sample_file_name += image_save_format
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    fig1.savefig(os.path.join(results_dir, sample_file_name))

    # save figure 2
    if not out_folder:
        try:
            out_folder = os.path.dirname(__file__)
        except NameError:
            out_folder = os.getcwd()
    results_dir = os.path.join(out_folder, source_data_name)
    sample_file_name = (
        f"REE_model_measured_comparison_{source_data_name}"
    )
    sample_file_name += image_save_format
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    fig2.savefig(os.path.join(results_dir, sample_file_name))

    # save figure 3
    if not out_folder:
        try:
            out_folder = os.path.dirname(__file__)
        except NameError:
            out_folder = os.getcwd()
    results_dir = os.path.join(out_folder, source_data_name)
    sample_file_name = (
        f"REE_model_measured_hist_removedREE_{source_data_name}"
    )
    sample_file_name += image_save_format
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    fig3.savefig(os.path.join(results_dir, sample_file_name))

    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)
    del fig1
    del fig2
    del fig3

    # save deviation data
    sample_file_name = (
        f"REE_model_measured_deviations_{source_data_name}"
    )
    export_deviation.to_csv(os.path.join(results_dir, (sample_file_name + ".csv")))


gc.collect()


def get_renaming(model, inverse=False, elements=REE(dropPm=False)):
    others = [ANOMALIES]
    renaming = {}
    for x in elements:
        if not inverse:
            renaming[x] = x + "_" + model
        else:
            renaming[x + "_" + model] = x
    for x in others:
        if not inverse:
            renaming[x] = x + "_" + model
        else:
            renaming[x + "_" + model] = x
    return renaming


def prepare_models(models):
    for k, v in models.items():
        columns = [SAMPLE, ANOMALIES]
        if "Y" in v.columns:
            columns.extend(REY(dropPm=False))
        else:
            columns.extend(REE(dropPm=False))
        v = v.loc[:, columns]
        v = v.rename(columns=get_renaming(k))
        models[k] = v
    for v in models.values():
        print(v.columns)
    return models


def merge_dataframe(real_data, models):
    all_data = real_data
    for v in models.values():
        all_data = all_data.join(v.set_index(SAMPLE), on=SAMPLE)
    return all_data


def extract_models(series, filter_models=None):
    model_names = set()
    for name in series.index:
        parts = name.split("_")
        if len(parts) > 1 and "model" in parts[1]:
            model_names.add(parts[1])
    models = {}
    for model in model_names:
        models[model] = series.loc[
            list(get_renaming(model).values()) + [SAMPLE]
            ].rename(index=get_renaming(model, True))
    if not filter_models:
        return models
    else:
        return {filter_models: models[filter_models]}


def make_plots_config(datafile, modelfile, outfolder, dataset_name, PATCHES, PATCHES_LIST, COMBINED,
                      COMBINED_ONLY, MEASURED_ONLY, make_REE_plots, make_histo_plots, sample_id= None, **kwargs):
    pandarallel.initialize(progress_bar=True, nb_workers=10)
    real_data = read_df(datafile, **kwargs)
    if sample_id is not None:
        real_data = real_data[real_data["sample"]==sample_id]
    if "Pm" not in real_data.columns:
        real_data = real_data.copy()
        Pm_loc = real_data.columns.get_loc("Nd") + 1
        real_data.insert(loc=Pm_loc, column="Pm", value=numpy.nan)
    else:
        real_data = real_data.copy()
    if "Unnamed: 0" in real_data.columns:
        real_data = real_data.drop(columns=["Unnamed: 0"])

    model_data = read_df(modelfile, **kwargs)
    model_data = model_data.copy()
    Pm_loc = model_data.columns.get_loc("Nd") + 1
    model_data.insert(loc=Pm_loc, column="Pm", value=numpy.nan)
    if "Unnamed: 0" in model_data.columns:
        model_data = model_data.drop(columns=["Unnamed: 0"])

    model = prepare_models({"model": model_data})
    all_data = merge_dataframe(real_data, model)

    if make_REE_plots:
        all_data.parallel_apply(
            lambda x: make_REE_pattern_plot_series(
                x.loc[real_data.columns],
                x.loc[model["model"].columns].rename(index=get_renaming("model", inverse=True)),
                dataset_name,
                PATCHES,
                PATCHES_LIST,
                COMBINED,
                COMBINED_ONLY,
                MEASURED_ONLY,
                out_folder=outfolder,
                # out_folder=os.path.join(outfolder, dataset_name)
            ),
            axis=1,
        )
    if make_histo_plots:
        make_histogram_plot(real_data, model_data, dataset_name, PATCHES, PATCHES_LIST,
                        out_folder=outfolder)


def _write_default_config():
    """
    Writes the default config for this file.
    """
    config = ConfigParser()
    config["plot"] = {
        "real_data_infile": "Real data",
        "model_data_infile": "Model data",
        "outfolder": "The folder to put the different plots",
        "dataset_name": "Used as prefix in the plots name"
    }
    conf_path = os.path.join(os.path.dirname(__file__), "configs",
                             "default_refill_REE_model_pattern_parallel.ini")
    with open(conf_path, "w") as conffile:
        config.write(conffile)


# On OSX with a M1, we need OBJ^C_DISABLE_INITIALIZE_FORK_SAFETY=YES
# `OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES python refill_REE_plot_model_real_pattern.py`
if __name__ == "__main__":
    start = time()
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
    fconfig = config["plot"]
    plot = config["plot"]
    os.makedirs(plot["outfolder"], exist_ok=True)
    if "DEV_THRESHOLD" in plot:
        DEV_THRESHOLD = float(plot["DEV_THRESHOLD"])
        print("DEV_THRESHOLD", DEV_THRESHOLD)
    else:
        DEV_THRESHOLD = 1.0
        print("DEV_THRESHOLD", DEV_THRESHOLD)
    if "SAMPLE" in plot:
        SAMPLE = plot["SAMPLE"]
    else:
        SAMPLE = "sample"
    if "IMAGE_SAVE_FORMAT" in plot:
        image_save_format = plot["IMAGE_SAVE_FORMAT"]
    else:
        image_save_format = ".png"

    make_plots_config(
        plot["real_data_infile"],
        plot["model_data_infile"],
        plot["outfolder"],
        plot["dataset_name"],
        bool(fconfig.get("PATCHES", False)),
        json.loads(fconfig.get("PATCHES_LIST", '["La"]')),
        bool(fconfig.get("COMBINED", False)),
        bool(fconfig.get("COMBINED_ONLY", False)),
        bool(fconfig.get("MEASURED_ONLY", False)),
        bool(fconfig.get("make_REE_plots", False)),
        bool(fconfig.get("make_histo_plots", False)),
        nrows=int(fconfig.get("nrows")) if "nrows" in fconfig else None,
        sample_id= plot["SAMPLE_ID"] if "SAMPLE_ID" in fconfig else None
    )
    print("Done in ", time() - start, "seconds")


