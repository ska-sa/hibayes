# ledabayes

MCMC model fitting for the global HI spectrum.

We provide a python framework for a fully Bayesian analysis of the
global HI spectrum from the cosmic dawn. As written the code uses
[(py)multinest](http://ccpforge.cse.rl.ac.uk/gf/project/multinest) for
sampling in order to obtain the evidence as well as the posterior, but
feel free to plug in your own sampler, MCMC or otherwise.

## Overview

- Uses a fully bayesian framework to fit models to the measured global
  HI spectrum. Fully bayesian here means parameter estimation **and**
  model selection via the bayesian evidence (Occam's razor
  quantified). See
  e.g. [Mackay 2003](http://www.inference.phy.cam.ac.uk/mackay/itila/book.html).
- Note that the algorithm is **not** maximum likelihood or an
  expensive least squares (though the likelihood function could be
  used for either of these). Rather, the sampler is used to explore
  the full posterior probability distribution of the model parameters
  and so unmask any degeneracies, multimodalities, correlations,
  wings, skirts, etc. between parameters
- For sampling we have used ```multinest```, though any nested or MCMC
  sampler could be used. ```multinest``` has the advantage of offering
  the evidence as well as the posterior samples (a higher density of
  samples corresponding to a higher probability). Calculating the
  evidence is intrinsically expensive compared to 'vanilla' MCMC, but
  is perfectly doable for problems, such as this one, that have a small
  (< 30, say) number of parameters.
- The code is in python plus MPI and takes a few minutes to run on a
  laptop, runtime depending mainly on the model complexity
  (i.e. number of parameters).
- The output for a given model is a bayesian evidence (plus
  uncertainty) and a 'chain' of posterior samples.
- The *modus operandi* is to use the evidence to select the winning
  model from a field of single-multinest-run competitors, then
  examining the corresponding triangle plot (the 'final answer') and
  deriving reconstructed and residual spectra.
- Currently implemented models are polynomial foregrounds and a
  gaussian empirical ```hi``` decrement, but any other parametric
  model can be coded in straightforwardly (e.g. simulated/physical
  parameters).
- Note that the form of the likelihood function assumes gaussian
  uncertainties on the input data, the uncertainties being propagated
  automatically by the sampling process.
- Note that there is an inescapable choice of priors on parameters
  (and models), but the evidence quantifies this.
- It is easy to add a joint likelihood over multiple data
  sets/experiments or incorporate models for telescope systematics.

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


## Supplied examples (see Bernardi et al. 2016

1. Simulation

2. LEDA data

## Citation

Use of this code should be cited as Zwart et al. 2016:

```@misc{ascl_hibayes,
author = {{Zwart}, J.~T.~L. and others},
title = "{HIBAYES: Global 21-cm Bayesian Monte-Carlo Model Fitting}",
howpublished = {Astrophysics Source Code Library},
year = 2016,
month = mar,
archivePrefix = "ascl",
eprint = {xxxx.xxx},
adsurl = {http://adsabs.harvard.edu/abs/xxxx.xxx},
adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

The algorithm, applied to both simulated and LEDA data, is described
in Bernardi et al. 2016:

(bibtex here)




