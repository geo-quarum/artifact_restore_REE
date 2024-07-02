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
import glob
import gc

import argparse
from configparser import ConfigParser
import matplotlib.pyplot as plt
from seaborn import histplot, kdeplot
import numpy as np
from numpy import nanmin, nanmax, nanstd
import pandas as pd
from pyrolite.geochem.ind import REE, REY


def plot_MCsim_histogram(input_df,
                         output_path,
                         fitting_type,
                         error_scale,
                         repetitions,
                         dataset_name,
                         image_format = "png",
                         data_modification="none",
                         xmin=-110,
                         xmax=110
                         ):

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

    # prepare data for histogram
    lambda_collect_df_flat = input_df.copy().to_numpy()
    histo_plot_data = lambda_collect_df_flat.flatten()
    # plot data
    histplot(data=histo_plot_data,
             ax=f1_ax1,
             color="teal",
             stat="density",
             zorder=1
             )
    kdeplot(data=histo_plot_data,
            ax=f1_ax1,
            color="orange",
            zorder=2
            )

    # percentiles, mean and SD
    num = len(histo_plot_data)
    p05 = np.nanquantile(histo_plot_data, 0.05)
    p95 = np.nanquantile(histo_plot_data, 0.95)
    p50 = np.nanquantile(histo_plot_data, 0.5)
    avg = np.nanmean(histo_plot_data)
    print(histo_plot_data)
    SD_val = "{:.1f}".format(nanstd(histo_plot_data))
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
    fitting_type_title = "Standard Polynomial Modelling" if "polyfit" in fitting_type else "λ Polynomial Modeling"
    plt.title(f"Monte Carlo Simulation for\n{fitting_type_title}",
              fontsize=12)
    f1_ax1.set_xlim(xmin, xmax)
    f1_ax1.set_ylim(0, 0.068)
    f1_ax1.set_xlabel("Deviation from data mean [%]")
    f1_ax1.text(0.1, 0.85,
                f"error scale: {error_scale} %"
                f"\nrepetitions: {repetitions}"
                f"\ndata modification: {data_modification}",
                fontsize=8,
                transform=f1_ax1.transAxes
                )

    sample_file_name = (f"MCsim_{fitting_type}_"
                        f"{data_modification}mod_{dataset_name}_histogram{image_format}")
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(os.path.join(output_path, sample_file_name))
    plt.close(fig1)
    del fig1
    gc.collect()


def prepare_MCsim_histogram_data(input_path):
    csv_path = os.path.join(input_path, "*.csv")
    REEelem = REE()
    
    collect_df = pd.DataFrame(columns=REEelem, dtype="object")
    
    for fname in glob.glob(csv_path):
        MCsim_df = pd.read_csv(fname)
        MCsim_df = MCsim_df.drop(columns=["Unnamed: 0"])
        MCsim_df.columns = REEelem
        aux_df = MCsim_df.apply(lambda x: ((x / np.mean(x)) - 1) * 100)
        collect_df = pd.concat([collect_df, aux_df])

    return collect_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("configfile")
    args = parser.parse_args()
    if not os.path.isfile(args.configfile):
        print(f"{args.configfile} does not exist")
    config = ConfigParser()
    config.read(args.configfile)
    sim_conf = config["sim_config"]

    data_path = sim_conf["DATA_PATH"]
    output_path = sim_conf["OUTPUT_PATH"]
    image_format = sim_conf["IMAGE_SAVE_FORMAT"]
    data_set_name = sim_conf.get("DATASET", "PetDB240122")
    fitting = sim_conf["fitting"]
    data_df = prepare_MCsim_histogram_data(data_path)

    plot_MCsim_histogram(data_df,
                         output_path,
                         fitting_type=fitting,
                         error_scale=10,
                         repetitions=100,
                         dataset_name=data_set_name,
                         image_format=image_format,
                         data_modification="ID",
                        )

if __name__ == "__main__":
    main()