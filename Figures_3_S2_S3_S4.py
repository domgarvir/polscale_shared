import numpy as np
import pandas as pd
import seaborn as sns
from Functions import *
import json
from natsort import natsorted
import matplotlib.colors as mcolors
from scipy.stats import spearmanr


# #load all visitation sequences
filename="../Data/pollinator_sequences_all_weeks_new.csv"#"./Data/pollinator_sequences_week2.csv"
pol_sequences_df=pd.read_csv(filename)
sequences=pol_sequences_df["Obs_id"].unique() #these are all possible sequences

#load plant abundances per plot and week
filename2="../Data/plant_abundances_all_weeks_new.csv"#"./Data/plant_abundances_week2.csv"
abundances_df=pd.read_csv(filename2, index_col=["Plot","Week"])
abundances_df=abundances_df.sort_index()

#load dataframe with positions in case we need to generate the simulation sequences with the spatial null model
filename3="../Data/plant_species_coordinates.csv"
coord_df=pd.read_csv(filename3, index_col=0)

#create the species maps and colors for later, such that they are common
species_names = pol_sequences_df["Plant_sp"].unique()#abundances_df['Plant_sp'].unique()
species_map = {name: i + 1 for i, name in enumerate(species_names)}
index_map={v: k for k, v in species_map.items()}
species_map["EMPTY"]=0
index_map[0]="EMPTY"

cmap = mpl.colormaps['tab20'].resampled(len(species_names))
color_list = [mcolors.to_hex(cmap(i)) for i in range(len(species_names))]
species_colors = {0: '#FFFFFF'}
for name, color in zip(species_names, color_list):
    species_colors[species_map[name]] = color 


#full names and short names of pollinators
pol_short_name_dcit={'Lasioglossum_sp':'Las','Eristalis_sp':'Eri','Megachile_sicula':'Msi',
'Anthophora_sp':'Adi','Xylocopa_cantabrita':'Xca', 'Apis_mellifera':'Ame',
'Bombus_terrestris':'Bte','Eucera_sp':'Eru','Bombylius_sp':'Bom','Andrena_hispania':'Ahi',
'Andrena_sp':'Afl','Vanessa_atalanta':'Vat','Leptotes_pirithous':'Lpi',
'Plebejus_argus':'Par','Amegilla_sp':'Am','Cyaniris_semiargus':'Cse','Dasypoda_sp':'Dci'}
pol_long_name_dcit={'Lasioglossum_sp':'Lasioglossum sp. (Las)','Eristalis_sp':'Eri',
'Megachile_sicula':'Megachile sicula (Msi)', 'Anthophora_sp':'Anthophora dispar (Adi)',
'Xylocopa_cantabrita':'Xylocopa cantabrita (Xca)', 'Apis_mellifera':'Apis mellifera (Ame)',
'Bombus_terrestris':'Bte','Eucera_sp':'Eucera rufa (Eru)','Bombylius_sp':'Bombilius sp. (Bom)',
'Andrena_hispania':'Andrena hispania (Ahi)', 'Andrena_sp':'Andrena flavipes (Afl)',
'Vanessa_atalanta':'Vat','Leptotes_pirithous':'Lpi', 'Plebejus_argus':'Par','Amegilla_sp':'Am',
'Cyaniris_semiargus':'Cse','Dasypoda_sp':'Dasypoda cingulata (Dci)'}
pol_long_name_simple_dcit={'Lasioglossum_sp':'Lasioglossum sp.','Eristalis_sp':'Eri',
'Megachile_sicula':'Megachile sicula', 'Anthophora_sp':'Anthophora dispar',
'Xylocopa_cantabrita':'Xylocopa cantabrita', 'Apis_mellifera':'Apis mellifera',
'Bombus_terrestris':'Bte','Eucera_sp':'Eucera rufa','Bombylius_sp':'Bombilius sp.',
'Andrena_hispania':'Andrena hispania', 'Andrena_sp':'Andrena flavipes',
'Vanessa_atalanta':'Vat','Leptotes_pirithous':'Lpi', 'Plebejus_argus':'Par','Amegilla_sp':'Am',
'Cyaniris_semiargus':'Cse','Dasypoda_sp':'Dasypoda cingulata'}

# SPATIALLY EXPLICIT NULL MODEL (NM2)
#load the simulated visitation sequences: These are the options for the spatial null model simulation
# those with NM1 are generated on the fly, as the model has a very low computation demand
start="species" # species  (in a hex of the given species) or random start      
mytype="random" # random or linear,  type of sampling sequende
her_arb_ratio="real" # real or fixed
M=10 # number of different simulated maps
S=50 # number of simulated sequences in each simulated map
filename="./Output/Sim_visitation_sequences_%s_%s_M%s_S%s_R%s.json" % (start,mytype,M,S,her_arb_ratio)

try:
    with open(filename, "r") as f:
        All_visitation= json.load(f)
except:
    #generate simulations of these sequences
    All_visitation={}
    for seq in sequences:
        print(seq)
        seq_df=pol_sequences_df[pol_sequences_df["Obs_id"]==seq]
        plot=seq_df["Plot"].unique()[0]
        week=seq_df["Week"].unique()[0]
        sequence=list(seq_df["Plant_sp"]) 
        if (len(sequence)>1): 
            sim_visits=get_sim_visitation(sequence,plot,week,abundances_df,coord_df,species_map,index_map,species_colors,M=5,S=20,type="random",start="species")
            All_visitation[seq]=sim_visits

    #store so don't have to simulate each time
    with open(filename, "w") as f:
        json.dump(All_visitation, f, indent=2)

#with open("nested_sequences.json", "r") as f:
#    nested_dict = json.load(f)

#Build a dataframe with the simulated visitation sequences
records = []

for seq_id, maps in All_visitation.items():
    for map_id, reps in maps.items():
        for rep_id, visit_array in reps.items():
            records.append((seq_id, map_id, rep_id, visit_array))

# Create the DataFrame
simVisitation_df = pd.DataFrame(records, columns=["seq", "map", "rep", "visits"])
simVisitation_df.set_index(["seq", "map", "rep"], inplace=True)
simVisitation_df["jumps"] = simVisitation_df["visits"].apply(calc_switches)

#we now have the real visitation sequences and the simualted sequences in the spatial null model
# Get sorted seqs
sorted_seqs = natsorted(simVisitation_df.index.get_level_values("seq").unique())
# Group and reindex for output
mean_jumps = simVisitation_df.groupby(level="seq")["jumps"].mean()
mean_jumps = mean_jumps.loc[sorted_seqs]
mean_jumps.name="mean"
print(mean_jumps)
std_jumps = simVisitation_df.groupby(level="seq")["jumps"].std()
std_jumps = std_jumps.loc[sorted_seqs]
std_jumps.name="std"
print(std_jumps)


# ABUNDACE BASED NULL MODEL (NM1)
## Now get the same simVisitation_df pero para el modelo nulo más simple en que las visitas se repartes de forma prop a la abundancia
Nrand=M*S
Vrand1_df=pd.DataFrame(index=range(len(sequences)*Nrand),columns=["seq","Plot","Week","rnd","secuence","steps","jumps","Distance"])

#Now I just need to combine this with the original dataframe of th empirical visitation frequencies
#create a database to store these visitation sequences
V_df=pd.DataFrame(index=sequences, columns=["Plot","Week","Pollinator_sp","secuence","steps","Distance","jumps_EMP"])

idx=0
for seq in sequences:
    print(seq)
    seq_df=pol_sequences_df[pol_sequences_df["Obs_id"]==seq]
    plot=seq_df["Plot"].unique()[0]
    week=seq_df["Week"].unique()[0]
    pol_sp=seq_df["Pollinator_sp"].unique()[0]
    #filter the abundances by site and week!!
    try:
        abn_df=abundances_df.loc[(plot,week)].drop("Date",axis=1)
        abn_df=abn_df.set_index("Plant_sp")
    except:
        abn_df=pd.DataFrame()
        print("%s no abundances" % seq)
    
    #visitation sequence length
    steps=seq_df.max()["Step"]

    #distance travelled
    distance=seq_df.sum(numeric_only=True)["Distance"]
    
    #get visitation secuence
    sequence=list(seq_df["Plant_sp"])

    #get number of jumps
    V,jumps=random_visitation_switches(sequence,abn_df, randomization=False)
    long=len(V)>1
    #print(V)
    V_df.loc[seq]=[plot,week,pol_sp,V,steps,distance,jumps]

    #NM
    #befor comparing with null lets check if all plants in the visitation sequence have abundances
    p_in_seq=set(sequence)
    p_in_abn=set(abn_df.index)
    my_subset=p_in_seq.issubset(p_in_abn)
    
    if (long & my_subset):
        #print("doing NM")
        for i in range(Nrand):
            #print('\r\033[K', end='') 
            #print(i, end='\r', flush=True)
            V_i1,jumps_i1=random_visitation_switches(V,abn_df,randomization=True,null_model=2)
            Vrand1_df.loc[idx]=[seq,plot,week,i,V_i1,len(V_i1),jumps_i1,0]
            idx += 1

#get mean values of null model
Vrand1_df.dropna(inplace=True)
cols_to_convert = ['steps', 'Distance', 'jumps']
Vrand1_df[cols_to_convert] = Vrand1_df[cols_to_convert].apply(pd.to_numeric, errors='coerce')
grouped_stats = Vrand1_df.groupby("seq", sort=False)["jumps"].agg(['mean', 'std'])
grouped_stats.columns=["mean_NM1", "std_NM1"]


Vall_df=pd.concat([V_df,mean_jumps,std_jumps,grouped_stats],axis=1)
#Vall_df=pd.concat([V_df,grouped_stats],axis=1)

Vall_df=Vall_df[Vall_df["steps"]>1]
Vall_df=Vall_df.dropna()
# Create a mask for non-zero standard deviation
non_zero_std = Vall_df["std"] != 0

# Initialize the column with NaN
Vall_df["z_score"] = np.nan
Vall_df["z_score_NM1"] = np.nan
# Compute z-scores only where std ≠ 0
Vall_df.loc[non_zero_std, "z_score"] = (
    (Vall_df.loc[non_zero_std, "jumps_EMP"] - Vall_df.loc[non_zero_std, "mean"]) /
    Vall_df.loc[non_zero_std, "std"]
)
Vall_df.loc[non_zero_std, "z_score_NM1"] = (
    (Vall_df.loc[non_zero_std, "jumps_EMP"] - Vall_df.loc[non_zero_std, "mean_NM1"]) /
    Vall_df.loc[non_zero_std, "std_NM1"]
)
#convert some columns to numeric
cols_to_convert = ['steps','Distance','jumps_EMP','z_score', 'z_score_NM1']
Vall_df[cols_to_convert] = Vall_df[cols_to_convert].apply(pd.to_numeric, errors='coerce')

Vall_df_copy=Vall_df.copy(deep=True)#copy priginal just in case

#we are going to retain only those pollinator that have at least 4 sequuences -----
# Count observations per pollinator and sort descending
order = Vall_df['Pollinator_sp'].value_counts().index.tolist()

#if we want to filter for the pollinators with more than 1 entry
order_trimmed = order[:-4]  
color_order=species_alphabetical = sorted(order_trimmed)

#consistent palette for pollinators
# Build a consistent palette (dictionary: species → color)
palette = dict(zip(color_order, sns.color_palette("tab10", n_colors=len(order_trimmed))))

#Para medir cuantas secuencias tenemos de los polinizadores seleccionados
V_df[(V_df["Pollinator_sp"].isin(order_trimmed)) & (V_df["steps"]>1)]

#Now we are going to retain only the sequences of the pollinators that have mroe than 3 individuals
Vselected_df = Vall_df[Vall_df["Pollinator_sp"].isin(order_trimmed)]
Nseq=Vselected_df.shape[0] #219 sequences
##### OBTAIN PERCENTAGES of SEQUENCES ABOVE and BELOW z=0 
percentage_below_zero = (Vselected_df["z_score"] < 0).mean() * 100 #78%
percentage_not_significant= percentage_within = ((Vselected_df["z_score"] >= -1.96) & (Vselected_df["z_score"] <= 1.96)).mean() * 100 #76%
percentage_below_zero_NM1 = (Vselected_df["z_score_NM1"] < 0).mean() * 100 # 88%
percentage_not_significant_NM1= percentage_within = ((Vselected_df["z_score_NM1"] >= -1.96) & (Vselected_df["z_score_NM1"] <= 1.96)).mean() * 100 #51%

L=Vselected_df["z_score"]<-1.96
H=Vselected_df["z_score"]>1.96
((L|H).sum())/Nseq #percentage of sequences not explained by NM2 (space+abundance)
((L|H).sum())
LNM1=Vselected_df["z_score_NM1"]<-1.96
HNM1=Vselected_df["z_score_NM1"]>1.96
((LNM1|HNM1).sum())/Nseq #percentage of sequences not explained by NM1 (abundance)
(LNM1|HNM1).sum() 

#LETS START FIGURES 
##### FIGURE 2 ######## Z-SCORE in Spatail null model by species
sns.set_context("paper", font_scale=0.8)  
#violin plot of z-score by pollinator type:
#plt.figure(figsize=(6,4))   # smaller, compact panel
fig, ax = plt.subplots(figsize=(8,5))  # create and capture explicitly

# horizontal reference lines
plt.axhline(0, color="grey", ls="-", lw=1,zorder=0)
plt.axhline(-1.96, color="grey", ls="--", lw=1,zorder=0)
plt.axhline(1.96, color="grey", ls="--", lw=1,zorder=0)
#plt.axhline( 1.96, color="black", ls="--", lw=1)

# Swarm plot (points)
sns.swarmplot( 
    data=Vselected_df, 
    x="Pollinator_sp",
    y="z_score", 
    ax=ax,
    hue="Pollinator_sp", 
    order=order_trimmed,
    palette=palette, 
    dodge=False, 
    legend=False,
    size=3.5,
    alpha=1.,
    linewidth=0.,          # 🔹 thickness of marker edge
    edgecolor="black"       # 🔹 color of marker edge
    )

""" sns.violinplot(
    data=Vall_df,
    x="Pollinator_sp",
    y="z_score",
    inner=None,
    cut=0,
    order=order,   # apply order
    color="lightgray"
)
 """
# Boxplot background
sns.boxplot(
    data=Vselected_df,
    x="Pollinator_sp",
    y="z_score",
    ax=ax,
    color="white",
    order=order_trimmed,
    showfliers=False, 
)
# Labels and ticks
plt.xlabel("Pollinator sp.", fontsize=9)
plt.ylabel("Z-score(J)", fontsize=9)
plt.yticks(fontsize=9)
# Replace tick labels using dictionary
ax.set_xticklabels([pol_short_name_dcit.get(label.get_text(), label.get_text()) 
                    for label in ax.get_xticklabels()], rotation=90, fontsize=9,ha='right')
#plt.xticks(rotation=90, fontsize=8)  # rotate labels if needed
""" plt.legend(
    title="Pollinator",
    bbox_to_anchor=(1.05, 1),   # (x, y) position outside plot
    loc="upper left",
    borderaxespad=0,
    fontsize=8,
    title_fontsize=9
) """

# shaded band between -1.96 and +1.96
plt.axhspan(-1.96, 1.96, color="#E0E0E0", alpha=0.2, zorder=0)
#plt.tight_layout()

# Legend outside at the bottom, two rows
handles = [mpatches.Patch(color=palette[sp], label=pol_long_name_dcit[sp]) for sp in species_alphabetical]

fig.legend( handles=handles, 
            loc="lower center", 
            bbox_to_anchor=(0.5, -0.15),
            ncol=5, 
            fontsize=8, 
            frameon=False)

filename="./Output/Fig_2.pdf"

# --- Adjust spacing so legend isn’t cut off ---
plt.subplots_adjust(bottom=0.1)  # increase bottom margin
#plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(filename, bbox_inches="tight")
plt.show()

##### FIGURE S2 - Comparative of Z-scores in NM1 and NM2 in all these sequences
##
ncol=1
nrow=1
plt.rcParams.update({'font.size': 16})
fig, axs = plt.subplots(nrow, ncol,figsize=(7*ncol,5*nrow),sharey=True)#row for different Nspecies, columns for plants and animals

#fixed bin 
bins = np.arange(-13, 5, 0.5)  # Bins: [0.0, 0.1, ..., 1.0]

sns.histplot(Vselected_df["z_score"], kde=True,stat='probability', bins=bins,  edgecolor="white",label="Spatially-explicit",ax=axs,alpha=0.8)  # Plot histogram with default fill beh

sns.histplot(Vselected_df["z_score_NM1"], kde=True,stat='probability', bins=bins,color='#C1CDCD',edgecolor="white",label="Abudance-weighted",ax=axs,alpha=0.8)  # Plot histogram with default fill beh
axs.axvline(x=-1.96, color='k', linestyle='dashed',linewidth=1.5)
axs.axvline(x=+1.96, color='k', linestyle='dashed',linewidth=1.5)
axs.axvline(x=0, color='k',linewidth=1.5)
#axs.set_title("NM 1")
axs.set_ylabel("Porbability Z-score(J)")
axs.set_xlabel("J")
axs.legend()

filename="./Output/Figure_S2.pdf"
plt.tight_layout()
plt.savefig(filename)
plt.show()  



########### FIGURE S3 - ZSCORE VS L 
plt.figure(figsize=(7.5,5.5))   # smaller, compact panel
# Compute Spearman correlation
rho, pval = spearmanr(Vselected_df["steps"], Vselected_df["z_score"])

ax=sns.scatterplot(
    data=Vselected_df,
    x="steps",
    y="z_score",
    hue="Pollinator_sp",
    hue_order=species_alphabetical,
    alpha=0.7,
    palette=palette
)
ax.set_xscale("log")
plt.axhline(-1.96, color="grey", ls="--", lw=1,zorder=0)
plt.axhline(1.96, color="grey", ls="--", lw=1,zorder=0)
plt.axhline(0, color="grey", ls="-", lw=1,zorder=0)

# Add correlation text
plt.text(
    0.65, 0.9,
    f"$\\rho$= {rho:.2f}({pval:.2f})",
    transform=plt.gca().transAxes,   # relative to axes
    fontsize=15,
)
# Legend outside at the bottom, two rows
handles = [mpatches.Patch(color=palette[sp], label=pol_long_name_simple_dcit[sp]) for sp in species_alphabetical]
plt.legend(
    #title="Pollinator",
    handles=handles,
    bbox_to_anchor=(0.5, -0.4),   # move legend outside below
    loc="lower center",
    ncol=5,                        # split into 2 columns (rows of entries)
    fontsize=15,
    frameon=False
)
ax.tick_params(axis='both', which='major', labelsize=12)
plt.xlabel("Total sequence length (L)", fontsize=15)
plt.ylabel("Z-score(J)", fontsize=15)
plt.subplots_adjust(bottom=0.25)  # 🔹 leave space for legend
#plt.tight_layout()
filename="./Output/Figure_S3.pdf"
plt.savefig(filename,bbox_inches="tight")
plt.show()


######## FIGURE S4
# Count distinct weeks per pollinator
V_long_df=Vselected_df[Vselected_df["steps"]>1]

weeks_present = (
    V_long_df.groupby("Pollinator_sp")["Week"]
      .nunique()
      .reset_index(name="Weeks_present")
).set_index("Pollinator_sp")
plots_present = (
    V_long_df.groupby("Pollinator_sp")["Plot"]
      .nunique()
      .reset_index(name="Plots_present")
).set_index("Pollinator_sp")
pollinator_plots = (
    V_long_df.groupby("Pollinator_sp")["Plot"]
      .unique()
      .reset_index(name="Plots_sampled")
).set_index("Pollinator_sp")
pollinator_weeks = (
     V_long_df.groupby("Pollinator_sp")["Week"]
      .unique()
      .reset_index(name="Weeks_sampled")
).set_index("Pollinator_sp")
samples=V_long_df.groupby("Pollinator_sp")["Pollinator_sp"].count()
samples.name="Samples"
z_mean=V_long_df.groupby("Pollinator_sp")["z_score"].mean()
z_std=V_long_df.groupby("Pollinator_sp")["z_score"].std()
z_mean.name="z_mean"
z_std.name="z_std"

Polis_df=pd.concat([samples,plots_present,pollinator_plots,weeks_present,pollinator_weeks],axis=1)

zscore_range = (
    Vselected_df.groupby("Pollinator_sp")["z_score"]
      .agg(z_min="min", z_max="max")
      .reset_index()
).set_index("Pollinator_sp")
# Add a column with the actual range
zscore_range["z_range"] = zscore_range["z_max"] - zscore_range["z_min"]

Polis_df=pd.concat([Polis_df,zscore_range,z_mean,z_std],axis=1)

#Ya tengo una base de datos con  para cada polinizador su máximo y minimo de z-score, en cuantos plots está, y cuantas semanas está
#plot de Samples, plots y weeks present vs z-score range, z-score mean
ncol=3
nrow=2

xs=["Samples","Plots_present","Weeks_present"]
ys=["z_mean","z_range"]

label_dict={"z_mean":"Mean Z-score(J)","z_range": "Range of Z-score(J)","Samples":"Samples","Plots_present": "Plots present","Weeks_present": "Weeks present"}

fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(ncol*5, nrow*3), sharey='row',sharex='col')

ix=0
for i, y in enumerate(ys):             # rows
    for j, x in enumerate(xs):         # columns
        ax = axes[i, j]
        
        
        
         # Add error bars only for the first row (y = Z_range)
        if y == "z_mean":
            colors = Polis_df.reset_index()["Pollinator_sp"].map(palette)

            ax.axhline(-1.96, color="lightgrey", ls="--", lw=1,zorder=0)
            ax.axhline(1.96, color="lightgrey", ls="--", lw=1,zorder=0)
            ax.axhline(0, color="lightgrey", ls="-", lw=1,zorder=0)

            for xi, yi, err, c in zip(Polis_df[x], Polis_df[y], Polis_df["z_std"], colors):
                ax.errorbar(
                xi, yi, yerr=err,
                fmt="o",
                ecolor=c,
                elinewidth=1.2, alpha=0.8, capsize=2, zorder=1
                )
            


        else:
            # Scatterplot
            sns.scatterplot(
                data=Polis_df.reset_index(), x=x, y=y,
                hue="Pollinator_sp", 
                s=50,  
                alpha=0.8,
                palette=palette,           # 🔹 same palette dictionary as swarmplot
                hue_order=species_alphabetical,  # optional: alphabetical order
                legend=False,
                ax=ax
            )

        # Only label leftmost y-axis
        if j == 0:
            ax.set_ylabel(label_dict[y], fontsize=12)
        else:
            ax.set_ylabel("")
        
        if i== 1:
            ax.set_xlabel(label_dict[x], fontsize=12)
        else:
            ax.set_xlabel("", fontsize=12)

        rho, pval = spearmanr(Polis_df[x], Polis_df[y])
        ax.text(
            0.65, 0.1,
            f"$\\rho$= {rho:.2f}({pval:.2f})",
            transform=ax.transAxes,   # relative to axes
            fontsize=12,
            )
        ax.tick_params(axis='both', which='major', labelsize=12)       

        # Add annotation at the upper left corner
        ax.text(0.05, 0.98, chr(65 + ix), transform=ax.transAxes, fontsize=15, va='top', ha='right')
        ix += 1
# Put one legend outside for the whole figure
# Build handles manually from your palette
handles = [mpatches.Patch(color=palette[sp], label=pol_long_name_simple_dcit[sp]) for sp in species_alphabetical]

fig.legend(
    handles=handles,
    loc="lower center",
    #bbox_to_anchor=(0.5, -0.12),   # 🔹 move lower (adjust -0.12 as needed)
    ncol=5,                        # 🔹 5 columns
    fontsize=12,
    frameon=False
)
plt.subplots_adjust(bottom=0.22)  # 🔹 leave space for legend
plt.savefig(filename, bbox_inches="tight")  # ensures legend is included

filename="./Output/Figure_S4.pdf"
plt.savefig(filename)
plt.show()