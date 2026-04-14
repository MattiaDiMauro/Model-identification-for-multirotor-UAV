# Model Identification for Multirotor UAV

This repository contains the project developed for the **Estimation and Learning in Aerospace** exam (A.Y. 2025/2026) at **Politecnico di Milano**.

The project investigates the system identification and statistical analysis of the longitudinal closed-loop dynamics of the **ANT-X**, a small-scale quadrotor UAV (270g), using frequency-domain techniques.

## Authors
Marianna Boschini
Filippo Cavalli
Mattia Di Mauro

**Supervisor:** Prof. Marco Lovera

---

## Project Overview

The study follows a structured identification pipeline to estimate aerodynamic derivatives and evaluate model uncertainty:

### 1. Grey-Box Identification
A linearised single-axis longitudinal model was developed to represent the drone's behavior around a hovering condition. 
- **Framework:** Grey-box state-space model with states $x = [u, q, \theta]^T$.
- **Methodology:** Frequency-domain identification using the MATLAB `greyest` function applied to sweep-excited simulated data.
- **Preprocessing:** Data windowing and filtering to isolate meaningful maneuvers and remove simulator transients.

### 2. Sensitivity Analysis
The study evaluates the impact of measurement noise on the identification performance. By perturbing the standard deviations of the accelerometer and gyroscope, we identified:
- The optimal sensor configuration for the most accurate model.
- The resulting parameter variance under different noise-to-signal ratios.
- The model's "FIT" and "VAF" (Variance Actually Fraction) metrics.

### 3. Monte Carlo Study
A statistical analysis was performed on the selected best configuration using **300 independent runs** to quantify:
- Statistical dispersion of estimated parameters (histograms).
- Distribution of system poles and zeros in the complex plane.
- Variability of Frequency Response Functions (FRFs) across noise realizations.

---

## Repository Structure

The project is organized as follows:

### Core
* main_ELAproject26.m: Main entry point. Run this script to execute the entire identification pipeline and analysis.
* common.m: Utility script containing core functions and shared routines required for the project execution.

### Data and Excitation
* ExcitationM.m: Script for the generation and management of the sine-sweep excitation signals.
* Excitation_3211.m and u3211.mat: Files containing the 3211 maneuver used for model validation.

### Parameters and Control
* parameters.m: Defines the physical parameters of the UAV and the linearised model matrices.
* controller.m: Contains the controller settings and gains for the closed-loop system.

### Simulation Models (Simulink)
* Simulator_single_axis.slx: The primary non-linear simulator of the quadrotor longitudinal dynamics.
* Validation.slx: Simulink model dedicated to running validation tests on the identified models.

---

## Getting Started

1. Clone the repository:
   git clone https://github.com/MattiaDiMauro/Model-identification-for-multirotor-UAV.git
   
2. Open MATLAB and set the repository folder as the current working directory.

3. Run the analysis:
   Simply run the main script in the Command Window:
   run('main_ELAproject26.m')

---

## Requirements
To run the simulations and identification scripts, the following are required:
* **MATLAB** (R2023b or later recommended)
* **Simulink**
* **System Identification Toolbox**
* **Control System Toolbox**
