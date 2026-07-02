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
from hexalattice.hexalattice import create_hex_grid
from scipy.spatial import cKDTree

def calc_switches(sampled_sequence):
    switches = sum(sampled_sequence[i] != sampled_sequence[i - 1] for i in range(1, len(sampled_sequence)))
    return switches

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

def get_sim_visitation(sequence,plot,week,abundances_df,coord_df,species_map,index_map,species_colors,M=1,S=10,type="random",start="random",herb_arb_ratio="fixed",farbust=0.5,fherb=0.8):
    L=len(sequence)
    sim_sequences={}
    initial=sequence[0]#initial species
    initial_id=species_map[initial]
    herb_arb_ratio=herb_arb_ratio
    start=start
    try:
        abn_df=abundances_df.loc[(plot,week)].drop("Date",axis=1)
        abn_df=abn_df.set_index("Plant_sp")

        my_coord_df=coord_df[coord_df["Plot"]==plot].reset_index().drop(["Plant_id","Plot"],axis=1)
        my_df=my_coord_df[my_coord_df["Plant_sp"].isin(list(abn_df.index))]
        p_in_seq=set(sequence)
        p_in_abn=set(abn_df.index)
        my_subset=p_in_seq.issubset(p_in_abn)

    except:
        abn_df=pd.DataFrame()
        my_coord_df=pd.DataFrame()
        my_df=pd.DataFrame()
        my_subset=False

    if(my_subset): #si tenemos las abunancias podemos seguir
        #determine herbaceous and non-herbaceous plants
        for m in range(M):#for each simulated map
            sim_sequences[m]={}
            final_hex_grid=create_nullmodel_map(plot,week,abundances_df,coord_df,species_map,index_map,species_colors,herb_arb_ratio=herb_arb_ratio)
            s=0
            while(s<S): #for each repetititon generate a apth
                if (start=="random"):
                    current = np.random.randint(len(final_hex_grid))
                else:
                    indices = np.where(final_hex_grid == initial_id)[0]
                    current = np.random.choice(indices)
                if (type=="random"):
                    my_path,my_visits=create_random_path(L,start=current,species_grid=final_hex_grid)
                elif (type=="linear"):
                    my_path,my_visits=create_straight_path(L,start=current,species_grid=final_hex_grid)
                if (len(my_path)==L):
                    species_visited=[index_map[int(final_hex_grid[x])] for x in my_path]
                    sim_sequences[m][s]=species_visited
                    s += 1
        return sim_sequences

    else:#si no, no podemos aplicar el modelo, salir
        return sim_sequences

def corr_func_annotate(x, y, axi=None, method="pearson", color="black", fontsize=13,square=False,pvalue=True, **kws):

    axi = axi or plt.gca()

    nan_policy='propagate'
    try:
        nan_policy=kws["nan_policy"]
    except:
        pass    

    if (method=="pearson"):
        try:
            my_df=pd.DataFrame({'x': x, 'y': y})
            my_df=my_df.dropna()
            r,p=stats.pearsonr(my_df["x"], my_df["y"])
        except:
            r,p=(np.nan,np.nan)

    elif(method =="spearman"):
        try:
            r, p = stats.spearmanr(x, y,nan_policy=nan_policy)
        except:
            r,p=(np.nan,np.nan)  
    print(r)
    print(p)
    x = 0.2
    y = 0.2
    try:
        x=kws["xy"][0]
        y=kws["xy"][1]
    except:
        pass

    try:
         color=kws["color"]
    except:
        pass

    if (square==True):
        axi.annotate(r'$r^{2} = {%.2f}$' % (pow(r,2)), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)

    else:
        if (pvalue==True):
            if (method=="spearman"):
                axi.annotate(r'$\rho = {%.2f}(%.2f)$'% (r,p), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
            else:
                axi.annotate(r'$r = {%.2f}(%.2f)$'% (r,p), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
        else:
            if (method=="spearman"):
                axi.annotate(r'$\rho = {%.2f}$'% (r), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
            else:
                axi.annotate(r'$r = {%.2f}$'% (r), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
            
    return

#helper for plotting bipartite networks
def calculate_positions(nodes, node_widths, buffer):
    positions = {}
    current_x = 0
    for node in nodes:
        width = node_widths[node]
        positions[node] = (current_x + width / 2, 0)  # Position at the center of the width
        current_x += width + buffer  # Move to the next position with buffer space
    return positions, current_x

def plot_bipartite_network_to_ax(Net,draw_labels=False,node_cathegories=None,height = 2,w=2,edgecolor="None",h=0.4,buffer=4,lw=1.2,ax=None,colored=True):
    node_color_dict={'common':'#009455','y1':'#73d2de','y2':"#ffd23f"}
    beta_link_dict={1.0:'y1',2.0:'y2',3.0:'common'}
    counts = [Net[u][v]['count'] for u, v in Net.edges()]
    norm = mcolors.Normalize(vmin=min(counts), vmax=max(counts))
    original_cmap = cm.Blues
    cmap = truncate_colormap(original_cmap, minval=0.2, maxval=1.0)

    plants = [n for n, d in Net.nodes(data=True) if d['bipartite'] == 'Plant']
    pollinators = [n for n, d in Net.nodes(data=True) if d['bipartite'] == 'Pollinator']
# Calculate degree of nodes
    degree_dict = dict(Net.degree(weight='weight'))

    sorted_plants = sorted(plants, key=lambda x: (degree_dict[x]), reverse=True)
    sorted_pollinators = sorted(pollinators, key=lambda x: (degree_dict[x]), reverse=True)

    # Normalize the degree for node width scaling
    max_degree = max(degree_dict.values())
    w=w
    node_widths = {node: (degree / max_degree) * w for node, degree in degree_dict.items()}  # Adjust scale as needed

    #sorted edges
    edge_categories = set([d['weight'] for u, v, d in Net.edges(data=True)])
    # Sort edge categories so higher categories are plotted last
    sorted_edge_categories = sorted(edge_categories, reverse=False)

    # Create a layout for the nodes with buffer space
    pos = {}
    height = 2
    buffer = w/buffer

    # Position plants at the bottom, sorted by degree
    plant_positions, max_x_plant = calculate_positions(sorted_plants, node_widths, buffer)
    for node, (x, y) in plant_positions.items():
        pos[node] = (x, 0)

    # Position pollinators at the top, sorted by degree
    pollinator_positions, max_x_pollinator = calculate_positions(sorted_pollinators, node_widths, buffer)
    for node, (x, y) in pollinator_positions.items():
        pos[node] = (x, height)

    # Determine the maximum x value for setting plot limits
    max_x = max(max_x_plant, max_x_pollinator)

    # Draw edges
    for category in sorted_edge_categories:
        edges = [(u, v, d) for u, v, d in Net.edges(data=True) if d['weight'] == category]
        for edge in edges:
            x1, y1 = pos[edge[0]]
            x2, y2 = pos[edge[1]]
            cosa=edge[2]['weight']
            if (node_cathegories):
                category = beta_link_dict[cosa]
                link_color = node_color_dict[category]
                if (cosa>2):
                    lwi=lw*1.8
                else:
                    lwi=lw
                ax.plot([x1, x2], [y1, y2], 'k-', lw=lwi, color=link_color,zorder=1)
            else:
                count=edge[2]["count"]
                color = cmap(norm(count))
                ax.plot([x1, x2], [y1, y2], '-',color=color, lw=category, zorder=1)   
    
    #edgecolor=edgecolor
    h=h
    # Draw nodes as rectangles after edges
    for node in sorted_plants:
        x, y = pos[node]
        width = node_widths[node]
        if (node_cathegories):
            category = node_cathegories[node]
            color = node_color_dict[category]
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor=color, zorder=2)
        else:
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor='grey', zorder=2)
        ax.add_patch(rect)

    for node in sorted_pollinators:
        x, y = pos[node]
        width = node_widths[node]
        if (node_cathegories):
            category = node_cathegories[node]
            color = node_color_dict[category]
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor=color, zorder=2)
        else:
            rect = Rectangle((x - width / 2, y - 0.1), width, h, linewidth=1, edgecolor=edgecolor, facecolor='black', zorder=2)
        ax.add_patch(rect)

    if (draw_labels):
    # Draw labels for better readability
        for node, (x, y) in pos.items():
            if node in sorted_plants:
                ax.text(x, y - 0.25, node, ha='center', va='top', fontsize=10, rotation=90, zorder=3)  # Below the rectangles for plants
            else:
                ax.text(x, y + 0.25, node, ha='center', va='bottom', fontsize=10, rotation=90, zorder=3)  # Above the rectangles for pollinators

    
    # Set limits and display the graph
    ax.set_xlim(-1, max_x + buffer)  # Adjust x limits based on cumulative width and buffer
    ax.set_ylim(-0.5, height+0.5)
    ax.set_aspect('equal')
    #plt.title('Bipartite Network of Plants and Pollinators')
    ax.axis('off')
    #plt.tight_layout()
    #filename="../Output/Network.pdf"
    #plt.savefig(filename)
    #plt.show()

    return ax

def truncate_colormap(cmap, minval=0.2, maxval=1.0, n=256):
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap
# Creating grid ###################
# HEXAGONAL
def get_hexa_neighbors(idx, nx, ny):
    """
    Get neighbors in a pointy-topped hex grid using odd-r layout (hexalattice style).
    Indexing is row-major: left to right, bottom to top.
    """
    row = idx // nx
    col = idx % nx

    if row % 2 == 0:  # even row
        directions = [(-1, -1), (-1, 0), (0, -1), (0, +1), (+1, -1), (+1, 0)]
    else:  # odd row
        directions = [(-1, 0), (-1, +1), (0, -1), (0, +1), (+1, 0), (+1, +1)]

    neighbors = []
    for dr, dc in directions:
        r, c = row + dr, col + dc
        if 0 <= r < ny and 0 <= c < nx:
            neighbors.append(r * nx + c)
    return neighbors

def are_hexa_neighbors(i, j, nx, ny):
    return j in get_hexa_neighbors(i, nx, ny)
#create initial hexagoanal grid
def create_initial_hexa_grid(my_df,my_abn_df,species_map,index_map,species_colors,herbaceous=False):
    # La funcion tiene los datos de abundancias y el mapa de las plantas arbustivas. 
    # Si herbaceous = False, no introduce las plantas herbaceas en el mapa inicial

    #Create the hexagonal lattice and place the plants
    nx, ny = 60, 60 # 2500 cells
    n_cells=nx*ny
    radius = 5.0
    min_diam = 2 * radius

    # Create hex grid (centered near 0,0)
    hex_centers, _ = create_hex_grid(
        nx=nx,
        ny=ny,
        min_diam=min_diam,
        align_to_origin=True,
        do_plot=False  # if you want to visualize it
    )

    #Bounding box of hex grid
    hx_min, hy_min = hex_centers.min(axis=0)
    hx_max, hy_max = hex_centers.max(axis=0)

    # Bounding box of plant coordinates (lets include a buffer to seaprate the outermost plants from the border)
    buffer=1
    px_min, py_min = my_df['X'].min()-buffer, my_df['Y'].min()-buffer
    px_max, py_max = my_df['X'].max()+buffer, my_df['Y'].max()+buffer

    if (herbaceous):
        #determine herbaceous plants without positions
        V_plants=set(my_df["Plant_sp"].unique())
        H_plants=list((set(my_abn_df["Plant_sp"].unique())).difference(V_plants))
        V_seeds=my_df["Plant_sp"].value_counts()
        mean_plants = V_seeds.mean()
        std_plants =  V_seeds.std() 
    
        new_rows = []
        for hp in H_plants:
            n = max(1, int(np.random.normal(mean_plants, std_plants)))  # at least 1
            xs = np.random.uniform(px_min, px_max, n)
            ys = np.random.uniform(py_min, py_max, n)
            for x, y in zip(xs, ys):
                new_rows.append({'Plant_sp': hp, 'X': x, 'Y': y})

        new_df = pd.DataFrame(new_rows)
        my_df = pd.concat([my_df, new_df], ignore_index=True)

    # Normalize plant coordinates to hex grid space
    my_df.loc[:,'X_norm'] = ((my_df['X'] - px_min) / (px_max - px_min)) * (hx_max - hx_min) + hx_min
    my_df.loc[:,'Y_norm'] = ((my_df['Y'] - py_min) / (py_max - py_min)) * (hy_max - hy_min) + hy_min

    tree = cKDTree(hex_centers)
    my_df.loc[:,'hex_index'] = my_df[['X_norm', 'Y_norm']].apply(lambda row: tree.query([row['X_norm'], row['Y_norm']])[1], axis=1)

    hex_assignments = my_df.groupby('hex_index').sample(n=1, random_state=42)#this avoids having more than one plant by cell
    #hex_assignments

    # --- Map species to IDs and colors ---
    #species_names = hex_assignments['Plant_sp'].unique()
    #species_map = {name: i + 1 for i, name in enumerate(species_names)}
    hex_assignments['species_id'] = hex_assignments['Plant_sp'].map(species_map)
    #print(species_map)
    #hex_assignments
    #cmap = mpl.colormaps['tab20'].resampled(len(species_names))
    #color_list = [mcolors.to_hex(cmap(i)) for i in range(len(species_names))]
    #species_colors = {0: '#FFFFFF'}
    #for name, color in zip(species_names, color_list):
        #species_colors[species_map[name]] = color

    

    # --- Initialize species grid and colors ---
    species_grid = np.zeros(n_cells, dtype=int)
    for _, row in hex_assignments.iterrows():
        species_grid[row['hex_index']] = row['species_id']
    colors = grid_to_colors(species_grid, species_colors)

    return species_grid, colors
#simulate colonization in hexa grid
def simulate_colonization(initial_species_grid,P_colonize,index_map,species_colors,N_STEPS=10,FILL=0.5,stop_condition="steps"):
    # La funcion simula la colonizacion de celdas vecinas por las celdas ocupadas. La condición de estop puede ser:
    # "steps" (la simulación corre un numero detemrinado de pasos), 0.7 (o cualquier otro numero entre 0 y 1, la simulación corre hasta
    # que haya tantas celdas ocuopadas como se pide)

    # --- Run the simulation ---
    nx, ny = 60, 60 # 2500 cells
    n_cells=nx*ny
    species_grid = initial_species_grid.copy()
    #p_colonize =0.3
    step=0
    Fill_percent=0

    #while np.any(species_grid == 0):
    if (stop_condition=="steps"):
        
        while (step<N_STEPS):
            new_grid = species_grid.copy()
            empty_cells = np.where(species_grid == 0)[0]
            np.random.shuffle(empty_cells)

            for idx in empty_cells:
                neighbors = get_hexa_neighbors(idx, nx, ny)
                neighbor_species = [species_grid[n] for n in neighbors if species_grid[n] > 0]
                if neighbor_species:
                    chosen = np.random.choice(neighbor_species)
                    if np.random.rand() < P_colonize[index_map[chosen]]:
                        new_grid[idx] = chosen
                        #colors[idx] = species_colors[chosen]

            species_grid = new_grid
            step += 1
    
        colors = grid_to_colors(species_grid, species_colors)
    
    else:
       
        #while np.any(species_grid == 0):
        while ((Fill_percent<FILL)& (step<1000)):
            #print(Fill_percent)
            new_grid = species_grid.copy()
            empty_cells = np.where(species_grid == 0)[0]
            Fill_percent=1-(len(empty_cells)/n_cells)
            np.random.shuffle(empty_cells)

            for idx in empty_cells:
                neighbors = get_hexa_neighbors(idx, nx, ny)
                neighbor_species = [species_grid[n] for n in neighbors if species_grid[n] > 0]
                if neighbor_species:
                    chosen = np.random.choice(neighbor_species)
                    if np.random.rand() < P_colonize[index_map[chosen]]:
                        new_grid[idx] = chosen
                        #colors[idx] = species_colors[chosen]

            species_grid = new_grid
            step += 1
        
        colors = grid_to_colors(species_grid, species_colors)


    #arbustive_grid=species_grid

    return species_grid, colors

# Creating PATHS ####################
def get_direction_offsets(row): #helper for straigth path function 
    """Return direction offsets for a pointy-topped odd-r layout depending on row parity."""
    if row % 2 == 0:  # even row
        return [(-1, -1), (-1, 0), (0, -1), (0, +1), (+1, -1), (+1, 0)]
    else:  # odd row
        return [(-1, 0), (-1, +1), (0, -1), (0, +1), (+1, 0), (+1, +1)]

def create_straight_path(N,start=None,species_grid=None,nx=60, ny=60):
    """
    Build a straight path of length N (or shorter if border is hit),
    starting at a random cell and going in a consistent direction based on hex parity.
    """
    total_cells = nx * ny
    
    # Choose a starting point
    if start is not None:
        start = start
    else:
        start = np.random.randint(total_cells)
    
    path = []
    visited = []
    current = start

    direction_index = np.random.randint(6)  # Save fixed direction
    steps = 0
    pvisits=0
    while pvisits < N:
        visited.append(current)
        if species_grid is None or species_grid[current] > 0:
            path.append(current)
            pvisits +=1

        row = current // nx
        col = current % nx

        # Get direction for current row
        offsets = get_direction_offsets(row)
        dr, dc = offsets[direction_index]

        new_row = row + dr
        new_col = col + dc

        if 0 <= new_row < ny and 0 <= new_col < nx:
            next_idx = new_row * nx + new_col
            current = next_idx
            steps += 1
        else:
            break  # Out of bounds

    return path, visited

#nx,ny,L_seq,index_map,start=current,species_grid=final_hex_grid
def create_random_path(N,start=None,species_grid=None,nx=60,ny=60):
#this function returns a random path with N visited plants
    total_cells = nx * ny

    if (start):
        current=start
    else:
        current = np.random.randint(total_cells)
    
    if (species_grid[current]>0):
        path = [current]
    else:
        path = []

    visited = set([current]) 

    while len(path) < N:
        neighbors = get_hexa_neighbors(current, nx, ny)
        np.random.shuffle(neighbors)
        for n in neighbors:
            if n not in visited:# si no se ha visitado ya
                visited.add(n)
                if species_grid[n]>0: #solo la añade al path si hay una planta a visitar
                    path.append(n)
                current = n
                break
        else:
            break  # no unvisited neighbors found



    return path,visited

#SQUARE LATTICE
def get_square_neighbors(row, col, nx, ny):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]  # up, down, left, right
    neighbors = []
    for dr, dc in directions:
        r, c = row + dr, col + dc
        if 0 <= r < ny and 0 <= c < nx:
            neighbors.append((r, c))
    return neighbors



def grid_to_colors(grid, color_map):
    """Convert a 1D grid of species IDs to a NumPy array of color hex codes."""
    return np.array([color_map[sp] for sp in grid], dtype='<U7')  # or just dtype=str

def get_sample_sequences_simple(plot,week,abundances_df,coord_df,species_map,index_map,species_colors,L_seq,N_seq=1,type_seq="random"):
    #esta función toma un plot, una semana, las abundancias de las plantas en esa semana en ese plot, las coordenadas espaciales de las arbustivas
    # , los códigos de las plantas, los códigos de color, la longitud de la secuencia a simular, el número de secuencias a devolver y el tipo de recorrido
 
    sequences=[]
    paths=[]

    #species_names = abundances_df['Plant_sp'].unique()
    #filter by plot and week
    my_coord_df=coord_df[coord_df["Plot"]==plot].reset_index().drop(["Plant_id","Plot"],axis=1)
    my_abn_df=abundances_df.loc[(plot,week)].drop("Date",axis=1)
    
    #get positions of non-herbaceous plants from selected plot and week
    my_df=my_coord_df[my_coord_df["Plant_sp"].isin(list(my_abn_df["Plant_sp"].values))]
    my_df["Plant_sp"].value_counts()

    #get noramlized abundances of each plant species
    P_colonize=my_abn_df.groupby("Plant_sp")["N_ind_flowers"].sum().pipe(lambda x: x / x.sum())

    #create hexagonal grid and simulations #(my_df,my_abn_df,species_map,index_map,species_colors)
    initial_hex_grid,initial_species_colors=create_initial_hexa_grid(my_df,my_abn_df,species_map,index_map,species_colors)

    #plant_colonization #(initial_species_grid,P_colonize,index_map,species_colors,N_STEPS=10)
    final_hex_grid,final_colors=simulate_colonization(initial_hex_grid,P_colonize,index_map,species_colors)

    nx, ny = 60, 60 # 2500 cells
    n_cells=nx*ny
    radius = 5.0
    min_diam = 2 * radius
    print_grid=False
    if (print_grid):
        
    
        nrow=1
        ncol=2
        fig, axs = plt.subplots(nrow, ncol,figsize=(5*ncol,4*nrow),sharey=True) 
        initial_colors=grid_to_colors(initial_hex_grid, species_colors)
        create_hex_grid(nx=60, ny=60, min_diam=min_diam, h_ax=axs[0],align_to_origin=True, do_plot=True, face_color=initial_colors)
        axs[0].set_title("Initial Grid from Empirical Plant Positions")

        final_colors=grid_to_colors(final_hex_grid, species_colors)
        create_hex_grid(nx=60, ny=60, min_diam=min_diam, h_ax=axs[1],align_to_origin=True, do_plot=True, face_color=final_colors)
        axs[1].set_title(f"Final Grid After Colonization (Steps: {10})")
        filename="./Output/Initial_conf.png"
        #plt.savefig(filename)
        plt.show()

    #generate the N_seq sequences, same floral map
    #for s in range(N_seq):
    n_sampled_seq=0
    while (n_sampled_seq < N_seq):
        #start from different place each time
        current = np.random.randint(n_cells)
        if (type_seq=="linear"):
            my_path,my_visits=create_straight_path(nx, ny,L_seq,start=current,species_grid=final_hex_grid)
        else: #random path
            my_path,my_visits=create_random_path(nx,ny,L_seq,start=current,species_grid=final_hex_grid)

        species_visited=[index_map[int(final_hex_grid[x])] for x in my_path]
        species_visited_clean=[x for x in species_visited if x != "EMPTY"]
        #print("L:%d myL:%d" % (L_seq,len(species_visited_clean)))
        if (len(species_visited_clean)==L_seq):
            sequences.append(species_visited_clean)
            paths.append(my_visits)
            n_sampled_seq += 1

        else:
            a=1


    return sequences, paths, final_hex_grid

def plant_herbaceous(initial_species_grid,my_df,H_plants,species_map,species_colors):
    
    species_grid = initial_species_grid.copy()
    V_seeds=my_df["Plant_sp"].value_counts()
    mean_plants = V_seeds.mean()
    std_plants =  V_seeds.std()
     
   
    for hp in H_plants:
        empty_cells = np.where(species_grid == 0)[0]
        n = max(1, int(np.random.normal(100, std_plants)))  # at least 1
        #now select n random cells
        to_fill = subset = np.random.choice(empty_cells, size=n, replace=False)

        for cell in to_fill:
            species_grid[cell]=species_map[hp]

        colors = grid_to_colors(species_grid, species_colors)
    
    return species_grid,colors

def get_sample_sequences(plot,week,abundances_df,coord_df,species_map,index_map,species_colors,L_seq,N_seq=1,type_seq="random", herbaceas="random",farbust=0.5,fherb=0.8):
    #esta función toma un plot, una semana, las abundancias de las plantas en esa semana en ese plot, las coordenadas espaciales de las arbustivas
    # , los códigos de las plantas, los códigos de color, la longitud de la secuencia a simular, el número de secuencias a devolver y el tipo de recorrido
    #    y la distribución de las plantas herbaceas (rellenando los huecos que quedan de las arbustivas) 
    sequences=[]
    paths=[]
   
    #species_names = abundances_df['Plant_sp'].unique()
    #filter by plot and week
    my_coord_df=coord_df[coord_df["Plot"]==plot].reset_index().drop(["Plant_id","Plot"],axis=1)
    my_abn_df=abundances_df.loc[(plot,week)].drop("Date",axis=1)
    
    #get positions of non-herbaceous plants from selected plot and week
    my_df=my_coord_df[my_coord_df["Plant_sp"].isin(list(my_abn_df["Plant_sp"].values))]
    my_df["Plant_sp"].value_counts()

    #determine herbaceous and non-herbaceous plants
    V_plants=set(my_df["Plant_sp"].unique())
    H_plants=list((set(my_abn_df["Plant_sp"].unique())).difference(V_plants))

    #get noramlized abundances of each plant species
    P_colonize=my_abn_df.groupby("Plant_sp")["N_ind_flowers"].sum().pipe(lambda x: x / x.sum()) 

    #create intial hexagonal grid with seeds of arbustive plants #(my_df,my_abn_df,species_map,index_map,species_colors)
    initial_hex_grid,initial_species_colors=create_initial_hexa_grid(my_df,my_abn_df,species_map,index_map,species_colors)

    if (len(H_plants)>0): #if herbaceous plants
        #run first colonization considering only arbustives untille 50% filled # (initial_species_grid,P_colonize,index_map,species_colors,N_STEPS=10,FILL=0.5,stop_condition="steps"))
        arbustive_hex_grid,arbustive_colors=simulate_colonization(initial_hex_grid,P_colonize,index_map,species_colors,FILL=farbust,stop_condition="fill")

        #randomply place herbaceus species in empty cells ((species_grid,my_df,H_plants,species_map)
        combined_hex_grid,combined_colors=plant_herbaceous(arbustive_hex_grid,my_df,H_plants,species_map,species_colors)

        #herbaceus_colonization #(initial_species_grid,P_colonize,index_map,species_colors,N_STEPS=10,FILL=0.5,stop_condition="steps")
        H_colonize=P_colonize.copy()
        H_colonize.loc[list(V_plants)]=0
        final_hex_grid,final_colors=simulate_colonization(combined_hex_grid,H_colonize,index_map,species_colors,FILL=fherb, stop_condition="fill")
    else:
        #run first colonization considering only arbustives untille 50% filled # (initial_species_grid,P_colonize,index_map,species_colors,N_STEPS=10,FILL=0.5,stop_condition="steps"))
        arbustive_hex_grid,arbustive_colors=simulate_colonization(initial_hex_grid,P_colonize,index_map,species_colors,FILL=fherb,stop_condition="fill")
        final_hex_grid=arbustive_hex_grid

    nx, ny = 60, 60 # 2500 cells
    n_cells=nx*ny
    radius = 5.0
    min_diam = 2 * radius
    print_grid=True

    if (print_grid):

        legend_handles = []
        for species_id, rgb_color in species_colors.items():
            species_name = index_map[species_id]
            patch = mpatches.Patch(color=rgb_color, label=species_name)
            legend_handles.append(patch)
    
        nrow=1
        ncol=4
        fig, axs = plt.subplots(nrow, ncol,figsize=(5*ncol,4*nrow),sharey=True) 
        #initial_colors=grid_to_colors(initial_hex_grid, species_colors)
        create_hex_grid(nx=60, ny=60, min_diam=min_diam, h_ax=axs[0],align_to_origin=True, do_plot=True, face_color=initial_species_colors)
        axs[0].set_title("Initial Grid from Empirical APlant Positions")

        #final_colors=grid_to_colors(final_hex_grid, species_colors)
        create_hex_grid(nx=60, ny=60, min_diam=min_diam, h_ax=axs[1],align_to_origin=True, do_plot=True, face_color=arbustive_colors)
        axs[1].set_title(f"Grid After Arbustive Colonization (Fill: {farbust})")

        create_hex_grid(nx=60, ny=60, min_diam=min_diam, h_ax=axs[2],align_to_origin=True, do_plot=True, face_color=combined_colors)
        axs[2].set_title(f"Grid with Herbaceous random placement")

        create_hex_grid(nx=60, ny=60, min_diam=min_diam, h_ax=axs[3],align_to_origin=True, do_plot=True, face_color=final_colors)
        axs[3].set_title(f"Grid with Herbaceous random placement")
        
        fig.legend(
        handles=legend_handles,
        loc="center right",    # You can change to 'lower center', 'upper right', etc.
        ncol=1,  # Number of columns in the legend
        frameon=False,
        #bbox_to_anchor=(0.5, 1.05)  # Positioning above the subplots
        )
        plt.tight_layout
        #filename="./Output/Initial_conf.png"
        #plt.savefig(filename)
        plt.show()

    #generate the N_seq sequences, same floral map
    #for s in range(N_seq):
    n_sampled_seq=0
    while (n_sampled_seq < N_seq):
        #start from different place each time
        current = np.random.randint(n_cells)
        if (type_seq=="linear"):
            my_path,my_visits=create_straight_path(nx, ny,L_seq,start=current,species_grid=final_hex_grid)
        else: #random path
            my_path,my_visits=create_random_path(nx,ny,L_seq,start=current,species_grid=final_hex_grid)

        species_visited=[index_map[int(final_hex_grid[x])] for x in my_path]
        species_visited_clean=[x for x in species_visited if x != "EMPTY"]
        #print("L:%d myL:%d" % (L_seq,len(species_visited_clean)))
        if (len(species_visited_clean)==L_seq):
            sequences.append(species_visited_clean)
            paths.append(my_visits)
            n_sampled_seq += 1

        else:
            a=1


    return sequences, paths, final_hex_grid

def create_nullmodel_map(plot,week,abundances_df,coord_df,species_map,index_map,species_colors,farbust=0.5,fherb=0.8,herb_arb_ratio="fixed"):
    #filter by plot and week
    my_coord_df=coord_df[coord_df["Plot"]==plot].reset_index().drop(["Plant_id","Plot"],axis=1).copy()
    my_abn_df=abundances_df.loc[(plot,week)].drop("Date",axis=1).copy()
    
    #get positions of non-herbaceous plants from selected plot and week
    my_df=my_coord_df[my_coord_df["Plant_sp"].isin(list(my_abn_df["Plant_sp"].values))].copy()
    #my_df["Plant_sp"].value_counts()

    #determine herbaceous and non-herbaceous plants
    V_plants=set(my_df["Plant_sp"].unique())
    H_plants=list((set(my_abn_df["Plant_sp"].unique())).difference(V_plants))

    #get noramlized abundances of each plant species
    P_colonize=my_abn_df.groupby("Plant_sp")["N_ind_flowers"].sum().pipe(lambda x: x / x.sum()) 
    
    if (herb_arb_ratio=="real"): #si quiero el ratio real lo calculo aqui
        farbust=(P_colonize.loc[list(V_plants)].sum())*0.8
        fherb=0.8

    #create intial hexagonal grid with seeds of arbustive plants #(my_df,my_abn_df,species_map,index_map,species_colors)
    initial_hex_grid,initial_species_colors=create_initial_hexa_grid(my_df,my_abn_df,species_map,index_map,species_colors)  
    #print("initial hex grid created")
    if (len(H_plants)>0):
        arbustive_hex_grid,arbustive_colors=simulate_colonization(initial_hex_grid,P_colonize,index_map,species_colors,FILL=farbust,stop_condition="fill")

    # if (len(H_plants)>0):
    #     #print("including herbaceous plants")
    #     #run first colonization considering only arbustives untille 50% filled # (initial_species_grid,P_colonize,index_map,species_colors,N_STEPS=10,FILL=0.5,stop_condition="steps"))
    #     if (herb_arb_ratio=="fixed"):
    #         arbustive_hex_grid,arbustive_colors=simulate_colonization(initial_hex_grid,P_colonize,index_map,species_colors,FILL=farbust,stop_condition="fill")
    #     else:
    #         #see the % of abundance that corresponds to arbustive and to herbaceous
    #         total_flowers = my_abn_df['N_ind_flowers'].sum()
    #         # Sum for each group
    #         arbustive_flowers = my_abn_df[my_abn_df['Plant_sp'].isin(V_plants)]['N_ind_flowers'].sum()
    #         herbaceous_flowers = my_abn_df[my_abn_df['Plant_sp'].isin(H_plants)]['N_ind_flowers'].sum()
    #         arbustive_pct = arbustive_flowers / total_flowers
    #         farbust=bounded_ratio_transform(arbustive_pct)*0.8
    #         fherb=0.8
    #         arbustive_hex_grid,arbustive_colors=simulate_colonization(initial_hex_grid,P_colonize,index_map,species_colors,FILL=farbust,stop_condition="fill")

            
        #print("arbustive grid ready!")

        NotComplete=True
        while (NotComplete):
            #randomply place herbaceus species in empty cells ((species_grid,my_df,H_plants,species_map)
            combined_hex_grid,combined_colors=plant_herbaceous(arbustive_hex_grid,my_df,H_plants,species_map,species_colors)
            #print("planting herbs!")
            #herbaceus_colonization #(initial_species_grid,P_colonize,index_map,species_colors,N_STEPS=10,FILL=0.5,stop_condition="steps")
            H_colonize=P_colonize.copy()
            H_colonize.loc[list(V_plants)]=0
            final_hex_grid,final_colors=simulate_colonization(combined_hex_grid,H_colonize,index_map,species_colors,FILL=fherb, stop_condition="fill")
            #check that total fill is fherb, if not fill with more
            # Count non-zero elements
            nonzero_count = np.count_nonzero(arbustive_hex_grid)
            # Total number of elements
            total_count = arbustive_hex_grid.size
            # Calculate percentage of non-zero values
            percentage_filled = (nonzero_count / total_count) * 100
            if (percentage_filled>=fherb):
                NotComplete=False
        #print("final grid ready!")
    else:
        #print("only arbustives here")
        arbustive_hex_grid,arbustive_colors=simulate_colonization(initial_hex_grid,P_colonize,index_map,species_colors,FILL=fherb,stop_condition="fill")
        final_hex_grid=arbustive_hex_grid
        #print("final grid ready!")
    
    return final_hex_grid

def plot_map(species_grid,species_colors,index_map):
    
    grid_colors=grid_to_colors(species_grid, species_colors)
    nx, ny = 60, 60 # 2500 cells
    n_cells=nx*ny
    radius = 5.0
    min_diam = 2 * radius

    # Create a list of legend handles
    legend_handles = []
    for species_id, rgb_color in species_colors.items():
        species_name = index_map[species_id]
        patch = mpatches.Patch(color=rgb_color, label=species_name)
        legend_handles.append(patch)

    nrow=1
    ncol=1
    fig, axs = plt.subplots(nrow, ncol,figsize=(5*ncol,4*nrow),sharey=True) 
    
    create_hex_grid(nx=nx, ny=ny, min_diam=min_diam, h_ax=axs,align_to_origin=True, do_plot=True, face_color=grid_colors)
    axs.set_title("Initial Grid from Empirical Plant Positions")

    fig.legend(
    handles=legend_handles,
    loc="center right",    # You can change to 'lower center', 'upper right', etc.
    ncol=1,  # Number of columns in the legend
    frameon=False,
    #bbox_to_anchor=(0.5, 1.05)  # Positioning above the subplots
    )
    plt.tight_layout
    filename="./Output/Initial_conf.png"
    #plt.savefig(filename)
    plt.show()

    return


def bounded_ratio_transform(p, min_ratio=0.1, max_ratio=0.9):
    # convert p ∈ [0,1] to log-odds
    odds = np.log(p / (1 - p))
    
    # Rescale odds to fit into bounded proportion
    # Cap log-odds at some threshold (say, ±2.2) so that extreme values become 0.1 / 0.9
    cap = 2.2
    odds = np.clip(odds, -cap, cap)
    
    # Rescale back to [min_ratio, max_ratio]
    scaled = 1 / (1 + np.exp(-odds))
    
    # Linearly map [1 / (1 + exp(cap)), 1 / (1 + exp(-cap))] → [min_ratio, max_ratio]
    min_input = 1 / (1 + np.exp(cap))   # ~0.1
    max_input = 1 / (1 + np.exp(-cap))  # ~0.9
    
    return min_ratio + (scaled - min_input) * (max_ratio - min_ratio) / (max_input - min_input)

def power_law(x, a, b):
    return a * x**b

def split_sequence(pollinator, seq):
    runs = []
    current_run = [seq[0]]
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            current_run.append(seq[i])
        else:
            runs.append((pollinator, seq[i-1], current_run.copy()))
            current_run = [seq[i]]
    runs.append((pollinator, seq[-1], current_run.copy()))
    return runs