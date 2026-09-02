This repository contains the data and code necessary to generate the results of the manuscript "Constant individuals, generalist species: floral constancy sustains plant–pollinator network functioning".

The original dataset of visitation sequences and associated plant data are included in ./Data folder. To generate the results and figures, run the code of each of the files as indicated below.

Figures 3, 5, and S2-S4 run in Python, using the libraries: numpy, pandas, seaborn, matplotlib, json, natsort, scipy, statsmodels, and hexalattice. The figure files (Figure_x.py), and helper functions (Functions.py) are included here.

To generate Figure 3 and Figures S2-S4, run: 

>> python Figures_3_S2_S3_S4.py 

Note: Because these figures rely on null model simulations, and the spatial null model takes time for computation, we have included a file containing simulated visits ("./Output/Sim_visitation_sequences*). If you want to create your own, simply remove or rename the original file and the script will automatically generate a a new sampling of the null models. 

To generate Figure 5 run:

>> python Figure_5.py

Figures 2, 4, and S1 run in R version 4.3.1 (2023-06-16). To generate these figures, run:

>> r Figures_2_4_S1.R

In all cases the output figures will appear in ./Output/Figure_*.pdf
