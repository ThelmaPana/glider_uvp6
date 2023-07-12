# Glider UVP6

Plankton and particles distribution across a mesoscale front during the spring bloom, seen from a glider and a UVP6.

Thelma Panaïotis PhD Thesis

Laboratoire d'Océanographie de Villefranche (UMR 7093)

# Repo structure

## Folders

-   `lib`: useful functions and plots

-   `data`: raw and processed data

-   `plots`: exploratory plots

-   `figures`: figures for the paper

## Scripts

Scripts order is self-explanatory. Scripts `00` to `07` are for plankton data processing, scripts `08` to `12` are for data preparation, scripts `13` to `16` are for data analysis.

-   plankton data processing

    -   `00.morphocluster_data.R` Get list of samples for Morphocluster import.

    -   `01.download_data.py` Download data from training [7545] and test [7544] sets from Ecotaxa.

    -   `02.download_images.py`Download images from training [7545] and test [7544] sets from Ecotaxa.

    -   `03.cut_bottom_images.py` Remove the scale bar at the bottom of images.

    -   `04.prepare_datasets_rf.R` Prepare datasets for RF  training: learning and test set, and new data to predict.

    -   `05.train_rf.py` Train a RF model on learning data, predict test data and new data.

    -   `06.calibrate_rf.Rmd` Calibrate the RF model based on cumulative precision on the test set.

    -   `07.evaluate_rf.Rmd` Evaluate RF performance on test set after calibration.

-   data preparation

    -   `08.download_data_science.py` Download data from Ecotaxa project for science.

    -   `09.regroup_plankton_taxo.R` Apply calibration threshold and regroup plankton taxonomy at higher level.

    -   `10.process_ctd.R` Download and process CTD and particle data.

    -   `11.process_objects.R` Compute plankton and particles concentration per bin.

    -   `12.interpolate.R` Interpolate environmental, particles and plankton data on transects.

-   data analysis

    -   `13.plots_for_paper.R` Generate plots for paper.

    -   `14.exploratory_plots.R` Generate exploratory plots for CTD and particles.

    -   `15.multivariate_analysis.R` Perform multivariate analyses and generate associated plots.

    -   `16.additional_analyses.R` Perform additional analyses, mostly to check biases and produce supplementary figures.
# glider_uvp6
