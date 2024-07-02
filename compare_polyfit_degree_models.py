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
import glob
import os
from configparser import ConfigParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from quarum_library.util import read_df


def get_file_name(file_path):
    file_path_components = file_path.split('/')
    file_name_and_extension = file_path_components[-1].rsplit('.', 1)
    return file_name_and_extension[0]


def compute_RMS_polyfit_various_degrees(measured_data,
                                        model_data_path,
                                        output_path,
                                        dataset_name):
    # load measured data
    measured_df = read_df(measured_data)
    # get list of all polyfit model data csv files
    model_list = glob.glob(model_data_path)

    # create empty dataframe that will save the output data
    RMS_df = pd.DataFrame()

    # iterate for all model data csv files
    for polyfit_model in model_list:
        # receive the polynomial degree for the model from the file name
        degree = get_file_name(polyfit_model)[-1]
        # load the model
        model_df = read_df(polyfit_model)
        # compute the deviation between model and measured data
        deviation_df = ((1 - (model_df.loc[:, "La":"Lu"] / measured_df.loc[:, "La":"Lu"])) * 100)
        # square the deviation
        dev_sqrd = deviation_df ** 2
        # compute the mean for the squared deviations for each REE
        dev_mean = dev_sqrd.apply(lambda x: np.sqrt(np.nanmean(x)), axis=0)
        # compute the mean for all REE for every degree
        dev_mean["all_REE_mean"] = np.mean(dev_mean[dev_mean.index.drop(["Ce", "Eu"])])

        # concat the data to the output dataframe
        RMS_df = pd.concat([RMS_df, dev_mean.rename(f"poly_deg_{degree}")], axis=1)

    # sort the columns to be in ascending order regarding the polynomial degree
    RMS_df = RMS_df.reindex(sorted(RMS_df.columns), axis=1)

    # write the output
    output_file = os.path.join(output_path, f"{dataset_name}polyfit_degrees_RMS.csv")
    RMS_df.to_csv(output_file)

    # return the output file
    return RMS_df


def plot_RMS_polyfit_various_degrees(RMS_df, dataset_name, output_path, image_format="png"):
    # plot style and setup
    plt.style.use("default")
    plt.rc("lines", markersize=2)
    plt.rc("axes", linewidth=2, titlepad=20)
    plt.rc("font", size=14)
    props = dict(boxstyle="round", facecolor="white", linewidth=0)
    fig1 = plt.figure(constrained_layout=False, figsize=(6, 6), dpi=300)
    gs1 = fig1.add_gridspec(1, 1, hspace=0.2, wspace=0.3)
    f1_ax1 = fig1.add_subplot(gs1[0, 0])
    fig1.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.9)

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

    xplot = np.arange(1, len(RMS_df.columns) + 1, 1)
    yplot_df = RMS_df.T
    cmap = sns.color_palette("tab20", 20)
    m = 0
    for elem in yplot_df.columns.drop(["all_REE_mean", "Ce", "Eu"]):
        yplot = yplot_df[elem]
        f1_ax1.plot(xplot, yplot, color=cmap[m], label=elem)
        m += 1

    # all REE mean
    yplot_mean = yplot_df["all_REE_mean"]
    f1_ax1.plot(xplot, yplot_mean, color="black", label="mean")

    legend = f1_ax1.legend(
        loc=2,
        bbox_to_anchor=(0.6, 0.95),
        ncols=2,
        edgecolor="inherit",
        scatteryoffsets=[0.5],
        fontsize=8,
    )

    f1_ax1.set_ylim(0, 21)
    f1_ax1.set_xlabel("Maximum polynomial degree", size=10)
    f1_ax1.set_ylabel("RMS deviation of model prediction from measured data", size=10)

    sample_file_name = f"{dataset_name}_polyfit_degrees_RMS.{image_format}"
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(os.path.join(output_path, sample_file_name))
    plt.close(fig1)
    del fig1
    gc.collect()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("configfile")
    args = parser.parse_args()
    if not os.path.isfile(args.configfile):
        print(f"{args.configfile} does not exist")
    config = ConfigParser()
    config.read(args.configfile)
    comp_config = config["comparison_config"]

    measured_data = comp_config["MEASURED_DATA"]
    model_data_path = comp_config["MODEL_DATA_PATH"]
    output_path = comp_config["OUTPUT_PATH"]
    dataset_name = comp_config["DATASET_NAME"]
    image_format = comp_config["IMAGE_SAVE_FORMAT"]

    os.makedirs(output_path, exist_ok=True)

    RMS_df = compute_RMS_polyfit_various_degrees(measured_data, model_data_path, output_path, dataset_name)
    plot_RMS_polyfit_various_degrees(RMS_df, dataset_name, output_path, image_format=image_format)


if __name__ == "__main__":
    main()
