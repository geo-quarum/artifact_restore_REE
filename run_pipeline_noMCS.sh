#!/bin/bash
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
set -e
#./generate_test_data.py

#./lambdatau_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_id_model_lambda.ini
#./lambdatau_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_naa_model_lambda.ini
#./polyfit_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_id_model_polyfit.ini
#./polyfit_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_naa_model_polyfit.ini

mkdir -p figures
# Fig. 1
# The output of Fig. 1 can be found in the "low_deviation(1.0)".
./plt_model_measured_REE_data.py configs/fig1.ini
cp "output/REE_comparison_plots/PetDB240122_mafic_ultramafic_ID/low_deviation(1.0)/REE_model_measured_pattern_PetDB240122_mafic_ultramafic_ID_PETDB-1284-COCOS 36.png" figures/fig1.png

#  Fig. 2
# The output of Fig. 2 are (A) "REE_model_measured_hist..."
# and (B) "REE_model_measured_hist_removedREE..."
./plt_model_measured_REE_data.py configs/fig2.ini
cp "output/REE_comparison_plots_lambda/PetDB240122_mafic_ultramafic_ID/REE_model_measured_hist_PetDB240122_mafic_ultramafic_ID.png" figures/fig2a.png
cp "output/REE_comparison_plots_lambda/PetDB240122_mafic_ultramafic_ID/REE_model_measured_hist_removedREE_PetDB240122_mafic_ultramafic_ID.png" figures/fig2b.png

# Fig. 3
./compare_polyfit_degree_models.py configs/fig3.ini
cp "output/PetDB_cleared_polyfit_degrees_RMS..png" figures/fig3.png

# Fig. 4
#echo "Run monte carlo simulation with lambda fitting"
#./monte_carlo_simulation.py configs/lambda_monte_carlo.ini

#echo "Run monte carlo simulation with poly fitting"
#./monte_carlo_simulation.py configs/poly_monte_carlo.ini

#echo "Make Figure 4"
#./plt_MCsim_histogram_LPM_SPM.py configs/fig4a.ini
#./plt_MCsim_histogram_LPM_SPM.py configs/fig4b.ini

#cp "output/MCsim/MCsim_lambda_IDmod_PetDB240122_histogram.png" figures/fig4a.png
#cp "output/MCsim/MCsim_polyfit_IDmod_PetDB240122_histogram.png" figures/fig4b.png

#compute some data for Fig. 5 and Fig 6.
./lambdatau_fitting.py configs/pottscondie_lambda_model.ini
./lambdatau_fitting.py configs/stosch_lambda_model.ini
./lambdatau_fitting.py configs/zindler_lambda_model.ini

./monte_carlo_simulation.py configs/pottscondie_monte_carlo.ini
./monte_carlo_simulation.py configs/stosch_monte_carlo.ini
./monte_carlo_simulation.py configs/zindler_monte_carlo.ini


# Fig. 5
./plt_model_measured_REE_data.py configs/fig5.ini
cp "output/REE_comparison_plots/pottscondie71/low_deviation(1.0)/REE_model_measured_pattern_pottscondie71_91_ultramafite.png" figures/fig5.png

# Fig. 6
# The output of Fig. 6 A to C are "REE_model_measured_hist..."
# Fig. 6a (config is identical to Fig. 5)
./plt_model_measured_REE_data.py configs/fig6a.ini
# Fig. 6b
./plt_model_measured_REE_data.py configs/fig6b.ini
# Fig. 6c
./plt_model_measured_REE_data.py configs/fig6c.ini

cp "output/REE_comparison_plots/pottscondie71/REE_model_measured_hist_pottscondie71.png" figures/fig6a.png
cp "output/REE_comparison_plots/zindleretal79/REE_model_measured_hist_zindleretal79.png" figures/fig6b.png
cp "output/REE_comparison_plots/stoschetal87/REE_model_measured_hist_stoschetal87.png" figures/fig6c.png
# Fig. 6d
./plt_MCsim_histogram_LPM_SPM.py configs/fig6d.ini
# Fig. 6e
./plt_MCsim_histogram_LPM_SPM.py configs/fig6e.ini
# Fig. 6f
./plt_MCsim_histogram_LPM_SPM.py configs/fig6f.ini

cp "output/MCsim/MCsim_lambda_IDmod_pottscondie_histogram.png" figures/fig6d.png
cp "output/MCsim/MCsim_lambda_IDmod_zindler_histogram.png" figures/fig6e.png
cp "output/MCsim/MCsim_lambda_IDmod_stosch_histogram.png" figures/fig6f.png

# Fig. 7
# Fig. 7a
./plt_model_measured_REE_data.py configs/fig7a.ini
# Fig. 7b
./plt_model_measured_REE_data.py configs/fig7b.ini

cp "output/REE_comparison_plots/stoschetal87/low_deviation(1.0)/REE_model_measured_pattern_stoschetal87_S468.png" figures/fig7a.png
cp "output/REE_comparison_plots/stoschetal87_combined/low_deviation(1.0)/REE_model_measured_pattern_stoschetal87_combined_S468.png" figures/fig7b.png

# Fig. 8
./fig8.py

cp "output/fig10/REE_model_measured_pattern_Zindler et al. (1979)_RE 78.png" figures/fig8.png

# SUPPLEMENTARY FIGURES
# Fig. S1
#./lambdatau_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_model_lambda.ini
#./lambdatau_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_nan_Pr_model_lambda.ini
#./lambdatau_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy_model_lambda.ini
#./lambdatau_fitting.py data/configs/PetDB240122_mafic_ultramafic_cleared_nan_Pr_Dy_Ho_model_lambda.ini
./plot_deviation_kde_model.py configs/figS1.ini

# Fig. S2 and Fig. S3
# The output of Fig. S2 are (A) "REE_model_measured_hist..."
# and (B) "REE_model_measured_hist_removedREE..."
./plt_model_measured_REE_data.py configs/figS2.ini
# The outputs of Fig. S3 are generated together with Fig. S2

cp -r output/REE_comparison_plots_s2_s3 figures/supplementary

# Fig. S4
# The output of Fig. S4 are (A) "REE_model_measured_hist..."
# and (B) "REE_model_measured_hist_removedREE..."
# The difference to Fig. 2 is that the model_input is the ID mod optimised data via Standard
# Polynomial Modelling
./plt_model_measured_REE_data.py configs/figS4.ini
cp "output/REE_comparison_plots_polyfit/PetDB240122_mafic_ultramafic_ID/REE_model_measured_hist_PetDB240122_mafic_ultramafic_ID.png" figures/supplementary/figS4a.png
cp "output/REE_comparison_plots_polyfit/PetDB240122_mafic_ultramafic_ID/REE_model_measured_hist_removedREE_PetDB240122_mafic_ultramafic_ID.png" figures/supplementary/figS4b.png

# Fig. S5
# The output of Fig. S5 can be found in the "low_deviation(1.0)".
# The difference between (A) and (B) is that (A) uses LPM and (B) SPM data
# Fig. S5a
./plt_model_measured_REE_data.py configs/figS5a.ini
# Fig. S5b
./plt_model_measured_REE_data.py configs/figS5b.ini

cp "output/REE_comparison_plots/PetDB240122_mafic_ultramafic_ID_lpm/low_deviation(1.0)/REE_model_measured_pattern_PetDB240122_mafic_ultramafic_ID_lpm_IRVCAN-SOM-BAT-K11A18.png" figures/supplementary/figS5a.png
cp "output/REE_comparison_plots/PetDB240122_mafic_ultramafic_ID_spm/low_deviation(1.0)/REE_model_measured_pattern_PetDB240122_mafic_ultramafic_ID_spm_IRVCAN-SOM-BAT-K11A18.png" figures/supplementary/figS5b.png