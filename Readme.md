This code generates the figures in the manuscript "Floral constancy can sustain plant–pollinator network functioning despite species-level generalization"

The original databases of visitation sequences are included in ./Data folder
To generate the asociated figures, just run the code of each of the files like indicated below.

Figures 3,5, S2-S4 run in python, using the libraries: numpy, pandas, seaborn, matplotlib, json, natsort, scipy, statsmodels, and hexalattice. 
The figure files (Figure_x.py), and helper functions (Functions.py) are included here.

To generate Figure 3 and Figures S2-S4 run: 

>> python Figures_3_S2_S3_S4.py 

Note: Because these figures rely on null model simulations, and the spatial null model takes time for computation, we have included a file containing simulated visits ("./Outoput/Sim_visitation_sequences*). If you want to create your own, simply remove or renamos the original file and the script will automatically generate a a new sampling of the null models. 

To generate Figure 5 run:

>> python Figure_5.py

In both cases the output figures will appear in ./Output/Figure_*.pdf