# DNA Mutation Kinetics Simulations in Rust
A Rusty reimplementation of series of DNA mutation kinetics simulations first carried out in [this work](https://github.com/SahakyanLab/GenomicPR2Simulations), available [here](https://www.biorxiv.org/content/10.1101/2022.12.23.521832v1). Please refer to the original paper for scientific background on the justification for kinetics implemented here. 

>This is a project to practice using Rust, and in it's current state does not reproduce the full paper results.

This simulation is an abbreviated version of the original simulations, all contained in a single executable parameterised in the `main.rs` script. I have included the original mutation rate data from the [Trek paper]<>, and the DNA base-pair contents from the original work in compressed `parquet` format.


## Overall Pipeline
To emulate the simulations in the original work, this script: 
1. Loads the base-pair statistics for the target organism
2. Loads the base-pair mutation rates for sampling. 
    - Currently only supports the truncated normal distribution.
3. Runs the ODE system with these original base-pair mutations and sampled rates
4. Plots the AT/GC skews in a simple 2d scatter-plot using plotters.


## Notes
As a learning tool it's been very useful to get used to Rust's syntax and behaviours, and it's clear that a lot of the work I've been doing in Python can be moved to Rust. 

The speed of this implementation has massively exceeded my expectations, but I can't yet guarentee that it's a stable solution to the original problem. 
