# hashpy2

hashpy2 is an interactive HASH [1] wrapper to meassure earthquake P wave first
motion polarities and S/P amplitude ratios and to compute focal mechanisms.

To install add the subfolders hashpy2 and hash to your PATH variable and make
sure that the files:

* hashpy2/hashpy2.py
* hashpy2/strdiprake2ptnaxes.py
* hashpy2/stereonet.py
* hashpy2/plot_mechanism.sh
* hash/hash_hashpy1D

are executeable. Then execute the hashpy2.py in your working directory. Make
sure that a config.yaml is present.

---

[1] HASH is a collection of fortran routines written by Jeane Hardebeck and
Peter M. Shearer. The original source code appears to be in the public domain
and can be found at:

https://earthquake.usgs.gov/research/software/index.php#HASH

When using this software in an academic context, please cite:

Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
2264-2276, 2002.

Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
Seismological Society of America, 93, 2434-2444, 2003.
