# Polygenic resistance / quantitative resistance model

## Model

Model of quantitative fungicide resistance and quantitative cultivar
protection, applied to septoria of winter wheat. The model is fitted to field
data.

## Code

The model is implemented in python.

In the work 'Modelling quantitative fungicide resistance and breakdown of resistant cultivars: designing integrated disease management strategies for Septoria of winter wheat', we describe the model found in `src/polymodel`.

The most important file in this folder is `simulator.py`, which contains the
classes `SimulatorOneTrait` and `SimulatorBothTraits` which are the fungicide
only and full model respectively. These have docstrings describing their use
and configuration.

Also see `fitting.py`, which contains functions and classes involved in model
fitting.

Much of the model analysis and fitting was done using the jupyter notebooks in
the `notebooks` folder.

### Scans

We ran various larger ensemble runs where the infection rate (beta) varied over
years and between runs. We used these to assess the effect of environmental
variability on the optimal strategies. The code that generated the outputs of
the scans is found in `src/cluster`.

### Plotting

Some plotting functions are found in `src/plots`, although you may prefer to
write your own custom plotting functions. Much of the plotting for the actual
paper was done in the notebooks in the `notebooks` folder.
