# emergentPredictability
MATLAB code to reproduce figures in Moran &amp; Tikhonov (2024), "Emergent predictability in microbial ecosystems"

Code by Jacob Moran and Mikhail Tikhonov

To use, run RunMe with no argument. The code will generate all data-dependent figures in the paper.

The precomputed data files will be used, if available. Remove (or rename) a data file to recompute it from scratch.

**Note:** The GitHub version does not include the three largest precomputed data files. Without them, the de novo simualtions will take several hours. For faster plotting, download the three additional files using [this link](https://wustl.box.com/s/q8h6fqmwy9aimn5ta6kqnxq8stb7jlz4) (data repository hosted at Washington University) and place them into the /figureData fodler alongside the .mat files that folder already contains.

Code uses "parfor" and will use the parallel computing toolbox with the default local configuration, if available.
