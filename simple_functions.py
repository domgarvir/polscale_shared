import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from scipy import stats
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

def random_visitation_switches(visitation_array, abundance_df,randomization=True,null_model=2):
    """
    Simulates random sampling of plant species based on abundance and calculates species switches.

    Parameters:
    - visitation_array (list): A list of plant species names representing visitation.
    - abundance_df (pd.DataFrame): A dataframe where the index contains plant species names
      and the column contains their respective abundances.
    - type of null model to use

    Returns:
    - int: The number of species switches in the randomly sampled visitation sequence.
    """
    if (randomization):
        #Null model 1: We fix the species in the visitation array so the randomized visits can only be fone in these species
        if (null_model==1):
            # Step 1: Extract species in visitation array and their corresponding abundances
            species_in_visits = list(set(visitation_array))  # Unique species in visitation array
            available_abundances = abundance_df.loc[species_in_visits].values.flatten()  # Get abundances
            species_labels = np.repeat(species_in_visits, available_abundances).astype(object)  # Expand species list by abundance
            # Step 2: Sample randomly with replacement to match visitation array length
            sampled_sequence = list(np.random.choice(species_labels, size=len(visitation_array), replace=True))
        #Null model 2: We allow the randomized visits in all the abundances in a given site and week.
        if (null_model==2):
            # Step 1: Now consider all species with their corresponding abundances
            available_abundances=abundance_df.values.flatten()
            species_labels = np.repeat(list(abundance_df.index), available_abundances).astype(object)
            # Step 2: Sample randomly with replacement to match visitation array length
            sampled_sequence = list(np.random.choice(species_labels, size=len(visitation_array), replace=True))
        
        # Step 2: Sample randomly with replacement to match visitation array length
        #sampled_sequence = list(np.random.choice(species_labels, size=len(visitation_array), replace=True))

    else:
        sampled_sequence=visitation_array

    # Step 3: Count the number of species switches
    switches = sum(sampled_sequence[i] != sampled_sequence[i - 1] for i in range(1, len(sampled_sequence)))


    return sampled_sequence, switches

def power_law(x, a, b):
    return a * x**b

#get psuedo R2
def get_pseudo_R2(model, Dvariable):
    y = model.data[Dvariable]
    y_hat = model.fits
    # Total variance (SST)
    sst = np.var(y, ddof=1)
    # Explained variance (SSR)
    ssr = np.var(y_hat, ddof=1)
    # Pseudo R²
    r2_pseudo = ssr / sst
    return r2_pseudo