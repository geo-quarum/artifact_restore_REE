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
import gc
import matplotlib.pyplot as plt
import numpy
import os

import pandas
import seaborn as sns
import logging

import numpy as np
from numpy import nanmin, nanmax

from pyrolite.geochem import REE, REY
from quarum_library.exceptions import NormalizationError
from quarum_library.constants import ANOMALIES
from quarum_library.util import read_df, normalize

logger = logging.getLogger(__name__)

def setup_fig():
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
    plt.close(fig1)
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
    return fig1, f1_ax1
def make_plot_series(
        real_data,
        models,
        MCsim_path,
        MCsim_config,
        sample_column,
        source_data_name,
        out_folder,
        MCsim,
        norm_to,
        units="ppm",
        suffix=None
):
    fig1, f1_ax1 = setup_fig()
    if "Chondrite" in norm_to:
        f1_ax1.set_ylabel("REE$_{sample}$/REE$_{C1}$")
    elif "EUS" in norm_to:
        f1_ax1.set_ylabel("REE$_{sample}$/REE$_{EUS}$")
    else:
        raise NormalizationError(f"Unknown axis label for normalization: {norm_to}")

    s_name = real_data[sample_column]
    try:
        s_type = real_data.loc["sample_type"]
    except KeyError:
        s_type = "unknown"
    try:
        s_source = real_data.loc["source"]
    except KeyError:
        s_source = "unknown"
    try:
        s_location = (
            real_data.loc["Longitude [degrees_east]"],
            real_data.loc["Latitude [degrees_north]"],
        )
    except KeyError:
        try:
            s_location = real_data.loc["sample_location"]
        except KeyError:
            s_location = "unknown"

    # get the Monte-Carlo simulation data
    if MCsim:
        MCsim_file_name = s_name + MCsim_config + ".csv"
        MCsim_data = read_df(os.path.join(MCsim_path, MCsim_file_name))
        Pm_loc = MCsim_data.columns.get_loc("Nd") + 1
        MCsim_data.insert(loc=Pm_loc, column="Pm", value=numpy.nan)

    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]

    # plot the data
    # data preparation
    REEpos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    REEelem = REE(dropPm=False)

    real_REE = real_data.pyrochem.REE
    real_REE = real_REE.pyrochem.normalize_to(norm_to, units)
    model_REE = models.pyrochem.REE
    model_REE = model_REE.pyrochem.normalize_to(norm_to, units)

    xplot = REEpos
    yplot_real = real_REE.loc["La":"Lu"]
    yplot_model = model_REE.loc["La":"Lu"]

    # plot the model data and scatter the measured data
    f1_ax1.scatter(xplot, yplot_real, s=20, label="measured", c="black", zorder=3)
    f1_ax1.plot(xplot, yplot_model, marker="o", markersize=7, alpha=0.7, color="teal", label="model", zorder=2)

    # plot boxplots for Monte-Carlo simulation data if TRUE
    if MCsim:
        MCsim_REE = MCsim_data.pyrochem.REE
        MCsim_REE = MCsim_REE.pyrochem.normalize_to(norm_to, units)
        MCsim_columns = MCsim_REE.columns
        boxplot_pos = 1
        for col in MCsim_columns:
            y_boxplot = MCsim_REE[col]
            f1_ax1.boxplot(y_boxplot,
                           positions=[boxplot_pos],
                           zorder=1,
                           )
            boxplot_pos += 1


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
    if isinstance(models['anomalies'], list):
        anomalies_for_text = ",".join(models['anomalies']).replace("'", "")
    else:
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
        # + "\nsample location: "
        # + str(s_location)
        + "\ndata source: "
        + str(s_source)
        + "\n\naverage deviation: "
        + str("{:.4f}".format(avg_dev)),
        ha="left",
        va="center",
        fontsize=8,
        transform=f1_ax1.transAxes,
        )

    legend = f1_ax1.legend(
        loc=2,
        bbox_to_anchor=(0.7, 0.95),
        edgecolor="inherit",
        scatteryoffsets=[0.5],
        fontsize=8,
    )

    f1_ax1.set_xlim(0.1, 15.9)
    f1_ax1.set_xticks(range(1, 16))
    f1_ax1.set_xticklabels(REEelem)
    f1_ax1.set_yscale("log")

    if not out_folder:
        try:
            out_folder = os.path.dirname(__file__)
        except NameError:
            out_folder = os.getcwd()

    sample_file_name = (
        f"REE_model_measured_pattern_{str(source_data_name)}_{str(real_data[sample_column])}"
    )
    if suffix:
        sample_file_name += f"_{suffix}"
    sample_file_name += ".png"

    results_dir = os.path.join(out_folder, source_data_name)
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir, exist_ok=True)
    plt.savefig(os.path.join(results_dir, sample_file_name))
    plt.close(fig1)
    del fig1
    gc.collect()


def get_renaming(model, inverse=False, elements=REE(dropPm=False)):
    others = [ANOMALIES, "LTE"]
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


def prepare_models(models, sample_column):
    for k, v in models.items():
        columns = [sample_column, ANOMALIES, "LTE"]
        if "Y" in v.columns:
            columns.extend(REY())
        else:
            columns.extend(REE(dropPm=False))
        v = v.loc[:, columns]
        v = v.rename(columns=get_renaming(k))
        models[k] = v
    return models


def merge_dataframe(real_data, models, sample_column):
    all_data = real_data
    for v in models.values():
        all_data = all_data.join(v.set_index(sample_column), on=sample_column)
    return all_data

def insert_pm_to_df(df):
    Pm_loc = df.columns.get_loc("Nd") + 1
    if "Pm" not in df.columns:
        df_copy = df.copy()
        df_copy.insert(loc=Pm_loc, column="Pm", value=numpy.nan)
    return df_copy

def plot_rees(real,
              model_data,
              MCsim_path,
              MCsim_config,
              sample_column,
              dataset_name,
              outfolder,
              MCsim=False,
              norm_to="Chondrite_PON"):
    real = insert_pm_to_df(real)
    model_data = insert_pm_to_df(model_data)
    model = prepare_models({"model": model_data}, sample_column)
    all_data = merge_dataframe(real, model, sample_column)
    all_data.parallel_apply(
        lambda x: make_plot_series(
            x.loc[real.columns],
            x.loc[model["model"].columns].rename(index = get_renaming("model", inverse=True)),
            MCsim_path,
            MCsim_config,
            sample_column,
            dataset_name,
            out_folder=outfolder,
            MCsim=MCsim,
            norm_to=norm_to
        ),
        axis=1,
    )

def make_sinlge_REE_plot(rseries, outfolder, model=None, suffix=None):
    fig1, f1_ax1 = setup_fig()
    if rseries.attrs["isnormalized"]:
        if "Chondrite" in rseries.attrs.get("normalization", ""):
            f1_ax1.set_ylabel("REE$_{sample}$/REE$_{C1}$")
        elif "EUS" in rseries.attrs.get("normalization", ""):
            f1_ax1.set_ylabel("REE$_{sample}$/REE$_{EUS}$")
        else:
            raise NormalizationError(f"Unknown axis label for normalization: {rseries.attrs.get('normalization', '')}")
    else:
        logger.warning("Plotting data without normalization")
        f1_ax1.set_ylabel("REE$_{sample}$")
    try:
        s_type = rseries.loc["sample_type"]
    except KeyError:
        s_type = "unknown"
    try:
        s_source = rseries.loc["source"]
    except KeyError:
        s_source = "unknown"
    try:
        s_location = (
            rseries.loc["Longitude [degrees_east]"],
            rseries.loc["Latitude [degrees_north]"],
        )
    except KeyError:
        try:
            s_location = rseries.loc["sample_location"]
        except KeyError:
            s_location = "unknown"
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]

    # plot the data
    # data preparation
    if "Y" in rseries.index or hasattr(rseries,"columns") and "Y" in rseries.columns :
        REEpos = list(range(1,17))
        REEelem = REY(dropPm=False)
        real_REE = rseries.pyrochem.REY
    else:
        REEpos = list(range(1,16))
        REEelem = REE(dropPm=False)
        real_REE = rseries.pyrochem.REE

    sample_file_name = (
        f"REE_model_measured_pattern_{str(rseries.attrs['source_data_name'])}"
    )

    xplot = np.array(REEpos)
    if isinstance(rseries, pandas.Series):
        s_name = rseries[rseries.attrs["sid"]]
        sample_file_name += f"_{str(s_name)}"

        yplot_real = real_REE.loc["La":"Lu"].replace(0, np.nan).astype(np.float64)
        ymin = nanmin(yplot_real) * 0.1
        ymax = nanmax(real_REE) * 10
        f1_ax1.set_ylim(ymin, ymax)
        # plot the model data and scatter the measured data
        mask = np.isfinite(yplot_real)
        legend_text = (f"sample: "
                       + str(s_name)
                       + "\nsample type: "
                       + str(s_type)
                       # + "\nsample location: "
                       # + str(s_location)
                       + "\ndata source: "
                       + str(s_source))
        if model is not None:
            yplot_model = model.loc[:,"La":"Lu"].squeeze()
            legend_text += f"\nexcluded{model['anomalies']}"
            f1_ax1.scatter(xplot[mask], yplot_real[mask], s=20, label="measured", c="black", zorder=3)
            f1_ax1.plot(xplot, yplot_model, marker="o", markersize=7, alpha=0.7, color="teal", label="model", zorder=2)
        else:
            f1_ax1.plot(xplot[mask], yplot_real[mask], marker='+', markersize=8, linewidth=1, label="measured", c="orange", zorder=3)



    else:
        legend_text = f"samples from : {str(rseries.attrs['source_data_name'])}"
        for i in range(len(real_REE)):
            s_name = rseries.loc[i, rseries.attrs["sid"]]
            yplot_real = real_REE.loc[:,"La":"Lu"].replace(0, np.nan)
            ymin = nanmin(yplot_real) * 0.1
            ymax = nanmax(real_REE) * 10
            print(ymin, ymax, yplot_real)
            f1_ax1.set_ylim(ymin, ymax)
            # plot the model data and scatter the measured data
            mask = np.isfinite(yplot_real.loc[i, :])
            f1_ax1.plot(xplot[mask], yplot_real.loc[i,mask], marker='+', markersize=8, linewidth=1, label=s_name, c=colors[i], zorder=3)

    fig1.text(
        0.1,
        0.85,
        legend_text,
        ha="left",
        va="center",
        fontsize=8,
        transform=f1_ax1.transAxes,
        )

    legend = f1_ax1.legend(
        loc=2,
        bbox_to_anchor=(0.7, 0.95),
        edgecolor="inherit",
        scatteryoffsets=[0.5],
        fontsize=8,
    )

    f1_ax1.set_xlim(0.1, REEpos[-1]+0.9)
    f1_ax1.set_xticks(REEpos)
    f1_ax1.set_xticklabels(REEelem)
    f1_ax1.set_yscale("log")

    source_data_name = rseries.attrs['source_data_name']
    if not outfolder:
        try:
            outfolder = os.path.dirname(__file__)
            logger.warning(f"No output folder specified. Writing data to: {outfolder}")
        except NameError:
            outfolder = os.getcwd()

    if suffix:
        sample_file_name += f"_{suffix}"
    sample_file_name += ".png"
    results_dir = os.path.join(outfolder, source_data_name)
    os.makedirs(results_dir, exist_ok=True)
    plt.savefig(os.path.join(results_dir, sample_file_name))
    plt.close(fig1)
    del fig1
    gc.collect()


def plot_measured_REE_pattern(real, outfolder, suffix=None, aggregate=False):
    real = insert_pm_to_df(real)
    if aggregate:
        make_sinlge_REE_plot(real, outfolder, suffix)
    else:
        real.apply(lambda x: make_sinlge_REE_plot(x, outfolder, suffix), axis=1)
    return real

def plot_REE_pattern_with_model(real, model, outfolder, suffix=None):
    real_pm = insert_pm_to_df(real)
    model = insert_pm_to_df(model)
    if real_pm.attrs["isnormalized"] and not model.attrs["isnormalized"]:
        model = normalize(model, real_pm.attrs["normalization"])
    real_pm.apply(lambda x: make_sinlge_REE_plot(x, outfolder,
            model=model[model[model.attrs["sid"]] == x[x.attrs["sid"]]], suffix=suffix), axis=1)
    return real
