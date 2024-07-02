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

import matplotlib.pyplot as plt
import numpy as np
import pandas
import pandas as pd
from ..fitting.polyfit import polynomial_fit
from pyrolite.geochem.ind import get_ionic_radii, REE
from pyrolite.util.lambdas.eval import get_function_components
from scipy.stats import norm as norm_dist
from quarum_library.exceptions import ConfigurationError


def lambdatau_fit(simulated_data, elements=REE(), degree=5, allElements=REE()):
    sim_df = pd.DataFrame(simulated_data.T, columns=allElements)
    LTE = False
    lambdas = sim_df.pyrochem.lambda_lnREE(norm_to=None, exclude=[x for x in allElements if x not in elements],
                                           degree=degree, LTE=LTE)
    radii = get_ionic_radii(REE(), coordination=8, charge=3)
    _, _, components = get_function_components(radii, degree=degree, fit_tetrads=LTE)
    return lambdas @ np.array(components)


def make_monte_carlo_plot(measured, simulations, outfile):
    REE_PRD_QNT = np.quantile(np.exp(simulations), [0.5, 0.25, 0.75, 0.05, 0.95], axis=0)
    plt.style.use("default")
    plt.rc("lines", markersize=2)
    plt.rc("axes", linewidth=2, titlepad=20)
    plt.rc("font", size=14)
    props = dict(boxstyle="round", facecolor="white", linewidth=0)
    fig1 = plt.figure(figsize=(6, 6), dpi=300)
    plt.title('REE profile reconstruction: Monte Carlo ensemble', fontsize=12)
    # plt.xlabel('Ionic radii [Angstrom]')
    # plt.ylabel('Normalized REE concentrations')
    plt.ylabel("REE$_{sample}$/REE$_{C1}$")
    plt.yscale('log')
    plt.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.1)

    plt.gca().tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        width=2,
        length=8,
    )
    plt.gca().tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        width=1,
        length=4,
    )
    measured_min = np.nanmin(np.exp(measured))
    measured_max = np.nanmax(np.exp(measured))
    ymin = 10**(np.log10(measured_min) - 1)
    ymax = 10**(np.log10(measured_max) + 1)
    fig1.gca().set_ylim(ymin, ymax)
    fig1.gca().set_xlim(0.1, len(measured) + 1.9)
    fig1.gca().set_xticks(range(1, len(measured) + 2))
    REEpos = list(range(1, len(measured) + 2))
    # REEelem = REE(dropPm=False)
    REEelem = REE(dropPm=False)
    fig1.gca().set_xticklabels(REEelem)
    # .. line plot of the median
    # insert Pm at position 4 in all quantile sets
    REE_PRD_QNT_allREE = REE_PRD_QNT.copy()
    REE_PRD_QNT_allREE = np.insert(REE_PRD_QNT_allREE,
                                   4,
                                   np.nan, axis=1)
    fig1.gca().plot(REEpos,
                    REE_PRD_QNT_allREE[0, :],
                    color='teal',
                    linewidth=1,
                    linestyle='-',
                    label='Median estimate')
    # .. line plots of the 25% and 75% percentiles (quartiles)
    fig1.gca().plot(REEpos,
                    REE_PRD_QNT_allREE[1, :],
                    color='teal',
                    linewidth=1,
                    linestyle=':',
                    label='Percentiles (25,75)')
    fig1.gca().plot(REEpos,
                    REE_PRD_QNT_allREE[2, :],
                    color='teal',
                    linewidth=1,
                    linestyle=':')
    # .. shaded region between the 5% and 95% percentiles
    fig1.gca().fill_between(REEpos,
                            REE_PRD_QNT_allREE[3, :],
                            REE_PRD_QNT_allREE[4, :],
                            color='teal',
                            alpha=0.3,
                            linewidth=0,
                            label='Percentiles (5,95)')
    # .. measurements
    # insert Pm between Nd and Sm with some dirty transformations...
    measured_allREE = measured.copy()
    measured_allREE = measured_allREE.to_frame().transpose()
    Pm_loc = measured_allREE.columns.get_loc("Nd") + 1
    measured_allREE.insert(loc=Pm_loc, column="Pm", value=np.nan)
    measured_allREE = measured_allREE.squeeze(axis=0)
    fig1.gca().scatter(REEpos,
                       np.exp(measured_allREE),
                       color='black',
                       marker='+',
                       s=100,
                       label='Measurements')

    # fig1.gca().autoscale(enable=True, axis='y', tight=True)
    legend = fig1.gca().legend(
        loc=2,
        bbox_to_anchor=(0.65, 0.95),
        edgecolor="inherit",
        scatteryoffsets=[0.5],
        fontsize=8
    )
    plt.savefig(outfile)
    plt.close(fig1)

"""
Fit the Monte Carlo Simulation
Both method expect the data is normalised and already in log format.

Usage examples:
ln_df_norm.apply(lambda x: compute_monte_carlo_data(x, polynomial_fit,
                                                    x.index.drop(
                                                        excludedElems),
                                                    "./monte_carlo.png",
                                                    exclude=excludedElems,
                                                    repetitions=500,
                                                    error=0.2), axis=1)


ln_df_norm.apply(lambda x: compute_monte_carlo_data(x, lambdatau_fit,
                                                    x.index.drop(
                                                        excludedElems),
                                                    "./monte_carlo_lambda.png",
                                                    exclude=excludedElems,
                                                    repetitions=500,
                                                    error=0.2), axis=1)

ln_df_norm.apply(lambda x: compute_monte_carlo_individual_errors(x, lambdatau_fit,
                                                    x.index.drop(
                                                        excludedElems),
                                                    "./monte_carlo_lambdamixed_errors.png",
                                                    exclude=excludedElems,
                                                    repetitions=500,
                                                    individual_errors=[0.1, 0.01, 0.15, 0.07, 0.03, 0.03, 0.5, 0.12, 0.19, 0.23, 0.30, 0.40, 0.1]), axis=1)
"""
def compute_monte_carlo_data(measured_values, fitting_function, elements,
                             outfile, exclude=[], repetitions=100,
                             error=0.1, save_intermediate=True,make_plots=True):
    measured_fitting = measured_values.to_numpy()
    if fitting_function == lambdatau_fit:
        simulated_data = measured_fitting[:, None] + measured_fitting[:, None] * norm_dist.rvs(scale=error, size=(len(measured_fitting), repetitions))
    elif fitting_function == polynomial_fit:
        simulated_data = measured_fitting[:, None] + norm_dist.rvs(scale=error, size=(len(measured_fitting), repetitions))
    else:
        raise ConfigurationError(f"Unknown fitting function for monte carlo simulation: {str(fitting_function)}")
    fitted = fitting_function(simulated_data, elements=[x for x in elements if x not in exclude])
    if save_intermediate:
        f = fitted
        if not isinstance(f, pandas.DataFrame):
            f = pandas.DataFrame(fitted, columns=elements)
        f = np.exp(f)
        f.to_csv(outfile.replace(".png", "") + ".csv")
    if fitting_function == lambdatau_fit:
        measured_values = np.log(measured_values)
    if make_plots:
        make_monte_carlo_plot(measured_values, fitted, outfile)
