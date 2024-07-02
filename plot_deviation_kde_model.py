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
# import packages
import argparse
import gc
import os
import sys
from configparser import ConfigParser

import matplotlib.pyplot as plt
import numpy
import pandas as pd
from numpy import std
from pyrolite.geochem.ind import REE
from scipy.stats import skew, kurtosis
from seaborn import kdeplot

from quarum_library.constants import REE_FULL_NAME

SAMPLE = "sample"


def compute_histo(model):
    s_list = model[SAMPLE]
    histo_model = model.copy()
    for s in s_list:
        aux = model.loc[model[SAMPLE] == s].copy()
        anomalies = aux["anomalies"]
        REElist = REE()
        for i in REElist:
            if i in str(anomalies):
                histo_model.loc[histo_model[SAMPLE] == s, i] = numpy.nan
    return histo_model


def make_kde_plot(setup, data_path, model_path, model_prefix, model_suffix, output_path,
                  save_as_format="png", xmin=-23, xmax=23, normto="Chondrite_PON"):
    real_data = pd.read_csv(data_path)

    if setup == "ID":
        model_suffixes = ["Pr", "PrTb", "PrTbHo", "ID"]
        model0_text = "none"
        model1_text = "Pr"
        model2_text = "Pr, Tb"
        model3_text = "Pr, Tb, Ho"
        model4_text = "Pr, Tb, Ho, Tm (ID)"
    elif setup == "NAA":
        model_suffixes = ["nan_Pr", "nan_Pr_Dy", "nan_Pr_Dy_Ho", "naa"]
        model0_text = "none"
        model1_text = "Pr"
        model2_text = "Pr, Dy"
        model3_text = "Pr, Dy, Ho"
        model4_text = "Pr, Dy, Ho, Er (INNA)"
    else:
        print("setup not found")
        sys.exit()

    model_data0 = pd.read_csv(os.path.join(model_path, (f"{model_prefix}_{model_suffix}")))
    model_data1 = pd.read_csv(os.path.join(model_path, (f"{model_prefix}_{model_suffixes[0]}_{model_suffix}")))
    model_data2 = pd.read_csv(os.path.join(model_path, (f"{model_prefix}_{model_suffixes[1]}_{model_suffix}")))
    model_data3 = pd.read_csv(os.path.join(model_path, (f"{model_prefix}_{model_suffixes[2]}_{model_suffix}")))
    model_data4 = pd.read_csv(os.path.join(model_path, (f"{model_prefix}_{model_suffixes[3]}_{model_suffix}")))

    # open one axes to set outlay for first plot
    fig1 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
    gs1 = fig1.add_gridspec(1, 1, hspace=0, wspace=0)
    f1_ax1 = fig1.add_subplot(gs1[0, 0])
    fig1.subplots_adjust(left=0.15, right=0.98, bottom=0.1, top=0.98)
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

    # remove anomaly data
    # model0
    s_list = model_data0[SAMPLE]
    histo_model_data0 = compute_histo(model_data0)

    # model 1
    histo_model_data1 = compute_histo(model_data1)

    # model 2
    histo_model_data2 = compute_histo(model_data2)
    # model 3
    histo_model_data3 = compute_histo(model_data3)
    # model 4
    histo_model_data4 = compute_histo(model_data4)

    kde_parameter_table = pd.DataFrame(columns=["REE",
                                                "removed REE",
                                                "n",
                                                "p05",
                                                "p95",
                                                "mean",
                                                "SD",
                                                "skewness",
                                                "kurtosis"])

    for plotting_REE in REE():
        REE_full = REE_FULL_NAME[plotting_REE]

        plt.style.use("default")
        plt.rc("lines", markersize=2)
        plt.rc("axes", linewidth=2, titlepad=20)
        plt.rc("font", size=14)
        props = dict(boxstyle="round", facecolor="white", linewidth=0)

        fig1 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
        gs1 = fig1.add_gridspec(1, 1, hspace=0, wspace=0)
        f1_ax1 = fig1.add_subplot(gs1[0, 0])

        # adjust the margine for the plots / subplots
        fig1.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)

        pl = ([f1_ax1])
        pn = 0
        for x in pl:
            pl[pn].axes.tick_params(
                axis="both",
                which="major",
                direction="in",
                top=True,
                right=True,
                width=2,
                length=8,
            )
            pl[pn].axes.tick_params(
                axis="both",
                which="minor",
                direction="in",
                top=True,
                right=True,
                width=1,
                length=4,
            )
            pl[pn].set_xlim(xmin, xmax)
            # pl[pn].set_ylim(0, 0.42)
            if pn < (len(pl) - 1):
                pn = pn + 1
            else:
                StopIteration

        # individual adjustmenst for combined plot
        f1_ax1.set_xlabel("Deviation from measured data in [%]")

        # model 0 - no deleted data
        plot_data0 = ((histo_model_data0[plotting_REE] / real_data[plotting_REE]) - 1) * 100
        plot_data0 = plot_data0.to_numpy()
        plot_data0 = plot_data0[~numpy.isnan(plot_data0)]
        plot_data0 = plot_data0.flatten()
        # plot_data0 = plot_data0[(1e2 > plot_data0)]
        kdeplot(data=plot_data0, ax=f1_ax1, color="black", lw=3, label=f"removed: {model0_text}")

        f1_ax1.text(0.05,
                    0.93,
                    f"{REE_full}",
                    fontweight="bold",
                    va="top",
                    ha="left",
                    transform=f1_ax1.transAxes)

        # model 1
        plot_data1 = ((histo_model_data1[plotting_REE] / real_data[plotting_REE]) - 1) * 100
        plot_data1 = plot_data1.to_numpy()
        plot_data1 = plot_data1[~numpy.isnan(plot_data1)]
        plot_data1 = plot_data1.flatten()
        # plot_data1 = plot_data1[(1e2 > plot_data1)]
        kdeplot(data=plot_data1, ax=f1_ax1, alpha=0.6, lw=2, label=f"removed: {model1_text}")

        # model 2
        plot_data2 = ((histo_model_data2[plotting_REE] / real_data[plotting_REE]) - 1) * 100
        plot_data2 = plot_data2.to_numpy()
        plot_data2 = plot_data2[~numpy.isnan(plot_data2)]
        plot_data2 = plot_data2.flatten()
        # plot_data2 = plot_data2[(1e2 > plot_data2)]
        kdeplot(data=plot_data2, ax=f1_ax1, alpha=0.6, lw=2, label=f"removed: {model2_text}")

        # model 3
        plot_data3 = ((histo_model_data3[plotting_REE] / real_data[plotting_REE]) - 1) * 100
        plot_data3 = plot_data3.to_numpy()
        plot_data3 = plot_data3[~numpy.isnan(plot_data3)]
        plot_data3 = plot_data3.flatten()
        # plot_data3 = plot_data3[(1e2 > plot_data3)]
        kdeplot(data=plot_data3, ax=f1_ax1, alpha=0.6, lw=2, label=f"removed: {model3_text}")

        # model 4
        plot_data4 = ((histo_model_data4[plotting_REE] / real_data[plotting_REE]) - 1) * 100
        plot_data4 = plot_data4.to_numpy()
        plot_data4 = plot_data4[~numpy.isnan(plot_data4)]
        plot_data4 = plot_data4.flatten()
        # plot_data4 = plot_data4[(1e2 > plot_data4)]
        kdeplot(data=plot_data4, ax=f1_ax1, alpha=0.6, lw=2, label=f"removed: {model4_text}")

        # determine kde parameters for each model and print
        aux_kde_parameter_table = pd.DataFrame(columns=["REE",
                                                        "removed REE",
                                                        "n",
                                                        "p05",
                                                        "p95",
                                                        "mean",
                                                        "SD",
                                                        "skewness",
                                                        "kurtosis"])
        aux_counter = 0
        for n_plot_data in plot_data0, plot_data1, plot_data2, plot_data3, plot_data4:
            REE_name = str(plotting_REE)
            remove_REE_aux_list = [model0_text, model1_text, model2_text, model3_text, model4_text]
            removed_REE = str(remove_REE_aux_list[aux_counter])
            n_data_points = len(n_plot_data)
            p05 = numpy.nanquantile(n_plot_data, 0.05)
            p95 = numpy.nanquantile(n_plot_data, 0.95)
            p50 = numpy.nanquantile(n_plot_data, 0.5)
            avg = numpy.nanmean(n_plot_data)
            SD_val = std(n_plot_data)
            skew_val = skew(n_plot_data)
            kurtosis_val = kurtosis(n_plot_data)

            aux_df = pd.DataFrame(data={"REE": REE_name,
                                        "removed REE": removed_REE,
                                        "n": n_data_points,
                                        "p05": p05,
                                        "p95": p95,
                                        "mean": avg,
                                        "SD": SD_val,
                                        "skewness": skew_val,
                                        "kurtosis": kurtosis_val},
                                  index=[0])
            aux_kde_parameter_table = pd.concat([aux_kde_parameter_table, aux_df], ignore_index=True)
            aux_counter = aux_counter + 1

        # concat all data for one large kde parameter table
        kde_parameter_table = pd.concat([kde_parameter_table, aux_kde_parameter_table], ignore_index=True)

        # legend
        legend = f1_ax1.legend(
            loc=2,
            bbox_to_anchor=(0.55, 0.95),
            edgecolor="inherit",
            scatteryoffsets=[0.5],
            fontsize=7
        )

        # save figure 1
        fig1.savefig(output_path + f"{setup}_combined_kde_{plotting_REE}.{save_as_format}")

        plt.close(fig1)
        del fig1
        gc.collect()

    kde_parameter_table.to_csv(output_path + f"{setup}_combined_kde_parameters.csv")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("configfile")
    args = parser.parse_args()
    if not os.path.isfile(args.configfile):
        print(f"{args.configfile} does not exist")
    config = ConfigParser()
    config.read(args.configfile)
    pconf = config["plot"]

    setup = pconf["setup"]
    real_data = pconf["data_path"]
    model_path = pconf["model_path"]
    model_prefix = pconf["model_prefix"]
    model_suffix = pconf["model_suffix"]
    output_path = pconf["output"]
    save_as_format = pconf.get("save_format", ".png")
    norm = str(pconf["norm_to"])
    print("run make plot", setup, real_data, model_path, output_path, save_as_format, norm)
    os.makedirs(output_path, exist_ok=True)
    make_kde_plot(setup, real_data, model_path, model_prefix, model_suffix, output_path,
                  save_as_format=save_as_format, normto=norm)


if __name__ == "__main__":
    main()

# %%
