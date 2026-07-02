import numpy as np
import pandas as pd
import seaborn as sns
from Functions import *
from sklearn.preprocessing import StandardScaler
import statsmodels.formula.api as smf
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.optimize import curve_fit

# #load visitation sequences
filename="./Data/pollinator_sequences_all_weeks_new.csv"#"./Data/pollinator_sequences_week2.csv"
pol_sequences_df=pd.read_csv(filename)

filename2="./Data/plant_abundances_all_weeks_new.csv"#"./Data/plant_abundances_week2.csv"
abundances_df=pd.read_csv(filename2, index_col=["Plot","Week"])
abundances_df=abundances_df.sort_index()

sequences=pol_sequences_df["Obs_id"].unique() #all possible sequences

#nombres completos y acronimos de polinizadores
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
    print(V)
    V_df.loc[seq]=[plot,week,pol_sp,V,steps,distance,jumps]

V_df_clean=V_df[V_df["steps"]>1]

# Count observations per pollinator and sort descending
n_seq_pol=V_df_clean['Pollinator_sp'].value_counts()
#if we want to filter for the pollinators with more than 3 entries
order_trimmed=n_seq_pol[n_seq_pol>3].index.tolist()
#consistent palette for pollinators
color_order=sorted(order_trimmed)
# Build a consistent palette (dictionary: species → color)
palette = dict(zip(color_order, sns.color_palette("tab10", n_colors=len(order_trimmed))))

#only retain information about specie with more than 3 sequences
V_df_clean=V_df_clean[V_df_clean["Pollinator_sp"].isin(order_trimmed)]

#lets build the network
# ---------------------------------------------
# 1 Number of distinct plants across all sequences per pollinator
# ---------------------------------------------
V_df_clean.loc[:,'unique_plants_in_seq'] = V_df_clean['secuence'].apply(lambda seq: set(seq))  # convert each to a set

# Combine across all sequences of the same pollinator
total_unique_plants = (
    V_df_clean.groupby('Pollinator_sp')['unique_plants_in_seq']
      .apply(lambda sets: len(set().union(*sets)))  # union of all sets per pollinator
)

# ---------------------------------------------
# 2️ Average number of distinct plants per sequence per pollinator
# ---------------------------------------------
avg_unique_per_sequence = (
    V_df_clean.groupby('Pollinator_sp')['unique_plants_in_seq']
      .apply(lambda sets: sum(len(s) for s in sets) / len(sets))
)
# Standard deviation of number of distinct plants per sequence per pollinator
std_unique_per_sequence = (
    V_df_clean.groupby('Pollinator_sp')['unique_plants_in_seq']
      .apply(lambda sets: np.std([len(s) for s in sets], ddof=1))  # ddof=1 for sample SD
)
# Plots and weeks the pollinators are present
pollinator_plots = (
    V_df_clean.groupby("Pollinator_sp")["Plot"]
      .unique()
      .reset_index(name="Plots_sampled")
).set_index("Pollinator_sp")
plots_present = (
    V_df_clean.groupby("Pollinator_sp")["Plot"]
      .nunique()
      #.reset_index(name="Plots_present")
)#.set_index("Pollinator_sp")
pollinator_weeks = (
     V_df_clean.groupby("Pollinator_sp")["Week"]
      .unique()
      .reset_index(name="Weeks_sampled")
).set_index("Pollinator_sp")
weeks_present = (
    V_df_clean.groupby("Pollinator_sp")["Week"]
      .nunique()
      #.reset_index(name="Weeks_present")
)

#total number of sequences for each polliantor sp
total_counts=V_df_clean["Pollinator_sp"].value_counts()


# ---------------------------------------------
# Combine into one DataFrame
# ---------------------------------------------
degree_df = pd.DataFrame({
    'total_unique_plants': total_unique_plants,
    'avg_unique_per_sequence': avg_unique_per_sequence,
    'std_unique_per_sequence':std_unique_per_sequence,
    'nweeks':weeks_present,
    #'weeks':pollinator_weeks,
    #'plots':pollinator_plots,
    'nplots':plots_present,
    'total_counts':total_counts
})

# --- Apply to all rows and collect results ---
rows = []
for _, row in V_df_clean.iterrows():
    rows.extend(split_sequence(row["Pollinator_sp"], row["secuence"]))

# --- Build the new dataframe ---
interactions_df = pd.DataFrame(rows, columns=["Pollinator_sp", "Plant_sp", "sequence_segment"])
print(interactions_df)

# Lets try to obtein the degree as a function of the number of sites and weeks the pollinator is in the field.
scaler = StandardScaler()
degree_df_norm = degree_df.copy()
normalized_values = scaler.fit_transform(degree_df)
cols=list(degree_df)
degree_df_norm.loc[:,cols]=normalized_values

model = smf.ols("total_unique_plants ~ nweeks + nplots", data=degree_df_norm).fit()
print(model.summary())

my_columns=["Model","nweeks", "nplots", "R2","AIC","BIC"]
my_results_df=pd.DataFrame(index=[1],columns=my_columns)
mod=1

# Extract and format results
nweeks_val = f"{model.params['nweeks']:.2f} ({model.pvalues['nweeks']:.2f})"
nplots_val = f"{model.params['nplots']:.2f} ({model.pvalues['nplots']:.2f})"

# Fill the DataFrame
my_results_df = pd.DataFrame([{
    "Model": "degree ~ nweeks + nplots",
    "nweeks": nweeks_val,
    "nplots": nplots_val,
    "R2": round(model.rsquared, 2),
    "AIC": round(model.aic, 2),
    "BIC": round(model.bic, 2)
}], columns=my_columns)

print(my_results_df)

#export table
table_latex=my_results_df.to_latex(float_format="%.3f", index=False)
filename="./Output/Table_S2.tex"
text_file = open(filename, "w")
text_file.write(table_latex)
text_file.close()

#now lets build the network ##### FIGURE 5 ###############
#quantity interaction
interactions_df["n_visits"] = interactions_df["sequence_segment"].apply(len)

##
quality_pol=interactions_df.groupby("Pollinator_sp").mean(numeric_only=True).squeeze() 
qualitystd_pol=interactions_df.groupby("Pollinator_sp").std(numeric_only=True).squeeze() 
quantity_pol=interactions_df.groupby("Pollinator_sp").sum(numeric_only=True).squeeze() 
degree_df = pd.DataFrame({
    'total_unique_plants': total_unique_plants,
    'avg_unique_per_sequence': avg_unique_per_sequence,
    'std_unique_per_sequence':std_unique_per_sequence,
    'nweeks':weeks_present,
    #'weeks':pollinator_weeks,
    #'plots':pollinator_plots,
    'nplots':plots_present,
    'total_counts':total_counts,
    'quality_avg':quality_pol,
    'quality_std':qualitystd_pol,
    'quantity':quantity_pol
})

# Set jitter amplitude (adjust 0.1 as needed)
jitter_strength = 0.1
nrow=1
ncol=2
fig, axs = plt.subplots(nrow, ncol,figsize=(5*ncol,4*nrow),sharey=False) 
# Default edge color: light gray
plt.rcParams['axes.axisbelow'] = True
axs[0].grid(color='gainsboro',linewidth = 0.5)
axs[1].grid(color='gainsboro',linewidth = 0.5)

x = degree_df["total_unique_plants"]
y = degree_df["quality_avg"]
yerr = degree_df["quality_std"]

# regresson line (lineal fit)
sns.regplot(
    x=x, y=y,scatter=False, color="grey", ci=False,ax=axs[0],line_kws={"lw":1,"ls":"--"})

# scatter points with errorbars
for xi, yi, err in zip(x, y, yerr):
                axs[0].errorbar(
                xi+np.random.uniform(-jitter_strength, jitter_strength), yi, yerr=err,
                fmt="o",
                #ecolor=c,
                elinewidth=1.2, alpha=0.8, capsize=0, zorder=2
                )

# 3️ correlation text
# --- compute correlation ---
r, p = pearsonr(x, y)
print(f"r = {r:.2f}, p = {p:.3f}")
axs[0].text(
    0.65, 0.95, f"$r={r:.2f}$({p:.2f})",
    transform=axs[0].transAxes, fontsize=10, va="top", ha="left",
    bbox=dict(facecolor="white", alpha=0.0, edgecolor="none")
)

axs[0].set_xlabel("Pollinator sp. degree")
axs[0].set_ylabel("Pollinator species'\nmean interaction quality")
#second figure ######################################

xdata = degree_df["total_counts"].values
ydata = degree_df["total_unique_plants"].values

# --- Compute correlation between original variables ---
r, p = spearmanr(xdata, ydata)
# Fit the parameters a, b
popt, pcov = curve_fit(power_law, xdata, ydata, p0=[1, 0.5])  # p0 gives initial guesses
a_fit, b_fit = popt
print(f"Best-fit parameters: a = {a_fit:.3f}, b = {b_fit:.3f}")
# Fitted curve
x_fit = np.linspace(min(xdata), max(xdata), 200)
y_fit = power_law(x_fit, *popt)

axs[1].plot(x_fit, y_fit, color="black", lw=1, ls="-" ,label=f"Fit: y = {a_fit:.2f}·x^{b_fit:.2f}",zorder=1)
# --- Add fit equation text ---
fit_text = (
    f"$y = {a_fit:.2f} \\, x^{{{b_fit:.2f}}}$\n"
    f"$\\rho = {r:.2f}({p:.2f})$"
)
axs[1].text(0.7, 0.15, fit_text, transform=axs[1].transAxes,
        fontsize=10, va='top', ha='left',
        bbox=dict(facecolor=None, alpha=0.0, edgecolor='none'))

#
sns.scatterplot(
    data=degree_df.loc[order_trimmed,:].reset_index(),
    x="total_counts",
    y="total_unique_plants",
    ax=axs[1],
    hue="Pollinator_sp",
    hue_order=color_order,
    alpha=0.9,

    legend=False,
    s=60,
    palette=palette
)

axs[1].set_xlabel("Foraging sequences")
axs[1].set_ylabel("Pollinator sp. degree")

# Add annotation at the upper left corner
axs[0].text(-0.05, 1.1, chr(65), transform=axs[0].transAxes, fontsize=15, va='top', ha='right')
axs[1].text(-0.05, 1.1, chr(66), transform=axs[1].transAxes, fontsize=15, va='top', ha='right')

handles = [mpatches.Patch(color=palette[sp], label=pol_long_name_simple_dcit[sp]) for sp in color_order]

fig.legend(
    handles=handles,
    loc="lower center",
    #bbox_to_anchor=(0.5, -0.12),   # 🔹 move lower (adjust -0.12 as needed)
    ncol=5,                        # 🔹 5 columns
    fontsize=9,
    frameon=False
)

plt.subplots_adjust(bottom=0.25)  # 🔹 leave space for legend
filename="./Output/Figure_5.pdf"
plt.savefig(filename,bbox_inches="tight")
plt.show()
