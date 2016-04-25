# ledabayes

MCMC model fitting for the global HI spectrum.

We provide a python framework for a fully Bayesian analysis of the
global HI spectrum from the cosmic dawn. As written the code uses
multinest for sampling in order to obtain the evidence as well as the
posterior, but feel free to plug in your own sampler, MCMC or
otherwise.

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

## Software requirements

- python 2.7.x
- MPI (optional)
- [multinest](http://ccpforge.cse.rl.ac.uk/gf/project/multinest)
- [pymultinest](http://johannesbuchner.github.io/PyMultiNest)

## Install

1. Fetch and install multinest, enabling MPI for an optional ~ ```NPROCS```
speed-up (see below). Check the multinest examples run.

2. Install pymultinest:

```pip install pymultinest```

Don't forget to set the ```(DY)LD_LIBRARY_PATH``` environment! No output
from the following command indicates success:

```python -c 'import pymultinest'```

Then check the pymultinest examples run.

3. Don't forget to

```chmod +x *py```


## Usage

From the project root directory,

```mpiexec -n NPROCS ./hi_multinest.py examples/1_simulated_data/config.ini```

where ```NPROCS``` is the number of cores you wish to use (execution with
MPI typically takes a few minutes on a laptop). Without MPI, just run

```./hi_multinest.py examples/1_simulated_data/config.ini```

The Multinest output goes to ```examples/1_simulated_data/output/```,
roughly as follows (see Multinest README for more info):

- ```1-stats.dat``` is written out every so often; the top line contains the
  evidence plus uncertainty, with the (optional)
  importance-nested-sampling (INS) evidence on line 2. Mean, ML and
  MAP parameter estimates follow. Caution - be wary of
  overinterpreting these averages and point estimates without
  eyeballing the full posteriors!
- ```1-post_equal_weights.dat```, populated once the multinest run is
  completed, contains the equally-weighted posterior samples, the last
  column being the value of the posterior for each sample. This file
  is used for plotting and reconstruction.
- ```1-summary.txt``` is used for reconstruction.
- ```1-ev.dat``` and ```1-phys_live.points``` are written out as sampling proceeds.
- ```1-IS.*``` are the corresponding files for INS.
- ```1-resume.txt``` for checkpointing.
- ```1-.txt``` for analysis in [cosmomc](http://cosmologist.info/cosmomc).

Now create a triangle plot (PDF):

```./hi_plot.py examples/1_simulated_data/config.ini```

And generate a text file containing a MAP-centred reconstruction with
error bars (see code for details):

```./hi_recon.py examples/1_simulated_data/config.ini```


