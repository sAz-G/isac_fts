# Implementation of an optimization algorithm for Integrated Sensing and Communication (ISAC)  

This repository contains a custom implementation of **_"ISAC from the Sky: UAV Trajectory Design for Joint Communication and Target Localization"_**. The implementation in this repository is based on [this paper](https://arxiv.org/abs/2207.02904). This repository is independent of the original paper.

## Installation

To run the code, you will need to install the [cvx](https://cvxr.com/cvx/) toolbox.

## Running the Code

1. **Set Parameters**: Open `parameters.m` and follow the instructions in the comments to set the parameters according to your requirements.

2. **Multi-Stage Algorithm**:
   - To run the multi-stage algorithm, execute `multi_stage_script.m` located in `src/main_algorithm/multi_stage/`.

3. **Monte Carlo Simulations**:
   - If you want to run Monte Carlo simulations, navigate to `src/main_algorithm/monte_carlo/`.

