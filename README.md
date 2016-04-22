# ledabayes

MCMC model fitting for the global HI spectrum.

## Overview

- Uses a fully bayesian framework to fit models to the measured global HI spectrum
- Fully bayesian = parameter estimation AND model selection via the evidence (Occam's razor quantified)
- This is NOT maximum likelihood or expensive least squares etc. - rather, use sampler to explore full posterior probability distribution of parameters and unmask degeneracies, multimodalities, correlations, wings, skirts, etc. between parameters
- Sampler is MultiNEST (Feroz et al.) which is geared towards calculating the evidence, but gives the posterior samples 'for free' (higher density of samples <-> higher probability)
- Calculating evidence via nested sampling is expensive compared to vanilla MCMC, but perfectly doable for < 30 parameters, say
- Code is in python + MPI and takes O(minutes) to run on laptop depending on model complexity
- Output for given model run is a bayesian evidence value plus error and a 'chain' of posterior samples
- Use evidence to select winning model, then examine corresponding triangle plot and derive reconstructed model
- Current models are polynomial foreground + gaussian HI; eventually connect to physical; any possible
- Noise supplied by Danny rather than using equation in slides; assumed gaussian (-> gaussian likelihood); errors propagated automatically by the sampling process
- Of course there's an inescapable choice of priors on parameters (and models!), but evidence quantifies this
- Easy to add a joint likelihood over multiple data sets/expts, or incorporate models for telescope systematics

## Install

## Usage

