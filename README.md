# Selection of fungicide insensitivity in *Phytophthora infestans* on potato.

This repository contains a model that simulates the selection for fungicide insensitivity within a *Phytophthora infestans* pathogen population on potato.

Once compiled, the program should be run with a command line argument between 0 and 25.

Each simulation produces two files:

1) Day*.csv - a file that gives the daily leaf area index of healthy, latent and infectious canopy, for each strain. Only simulations 0:10 produce this file.
2) Year*.csv - the yearly severity, AUDPC, HAD and gene frequencies of each simulation.

The following table details the simulations

|Number|Description|
|--|------------|
|0|Run the simulation with no pathogen|
|1-5|The simulation is run with a resistant cultivar (R3-R7 respectively), but no fungicide|
|6-10|The cultivar resistance is varied from R3-R7. A single fungicide to which resistance is developing is included. The resistant pathogen strain is completely insensitive to the fungicide.|
|11-15|As for simulations 6:10, but with the addition of a mixing partner|
|16-20|As simulations 6:10, but the resistant pathogen strain is only partially insensitive to the fungicide|
|21-25|As simulations 16:20, but with the addition of a mixing aprtner|

# Graphs

The repository also contains a file Graphs.r, which specifies how the graphs for the paper were put together.
