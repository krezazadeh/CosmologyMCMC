# Cosmology MCMC (MPI + LAPACK/BLAS) - Fortran 90

## Overview

This program implements a parallel Markov Chain Monte Carlo (MCMC) sampler for cosmological parameter estimation using the Metropolis-Hastings algorithm. It runs multiple chains in parallel using MPI and checks convergence dynamically using the Gelman-Rubin R̂ statistic.

The code models cosmological observables based on:
- Type Ia Supernovae (SN)
- Baryon Acoustic Oscillations (BAO)
- Cosmic Microwave Background (CMB)

It is written in Fortran 90, optimized with LAPACK/BLAS linear algebra libraries, and parallelized via MPI.

## Features

- Parallel MCMC sampling using MPI (one process per chain)
- Uses SN, BAO, and CMB data to compute likelihoods
- Supports flat and curved ΛCDM cosmological models
- Convergence check using Gelman-Rubin statistic (R̂)
- Output of parameter chains and best-fit estimates
- Modular design with easy-to-edit priors and parameters

## Requirements

- Intel Fortran compiler with MPI support (`mpiifort`)
- LAPACK and BLAS libraries
- Observational data files (SN, BAO, CMB) in suitable format

## Compilation

To compile the program, run:

```bash

mpiifort -O2 CosmologyMCMC.f90 -o CosmologyMCMC -llapack -lblas

## Execution

mpirun -np 8 ./CosmologyMCMC
