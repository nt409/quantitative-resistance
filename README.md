# Polygenic resistance / quantitative resistance model

## Model

Model of quantitative fungicide resistance, applied to septoria of winter
wheat. The model is fitted to field data.

## Code

The model is implemented in python.

In the work 'Modelling quantitative fungicide resistance and breakdown of resistant cultivars: designing integrated disease management strategies for Septoria of winter wheat', we describe the model found in `src\polymodel`.

The most important file in this folder is `simulator.py`, which contains the
classes `RunSingleTactic` and `RunGrid` which are the main ways we test tactic
performance. These have docstrings which describe their use and how to
configure them.

### Scans

We ran three scans over different parameter values/scenarios in the work:

- Robustness to variation in parameter values; see `src\param_scan`
- Effect of between-season sexual reproduction; see `src\sr_scan`
- Mixtures vs Alternations (with between-season sexual reproduction); see
  `src\alternation_scan`

### Plotting

Some plotting functions are found in `src\plotting`, although you may prefer to
write your own custom plotting functions. Most of these are called in various
scripts in `src\figs` or in the notebooks.
