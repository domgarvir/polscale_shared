# 03-07-2026
# Code to create Fig. 2, Fig. 4 and Fig. S1

# Load packages
library(plyr)
library(sand)
library(robustbase)
library(ergm) # Will load package 'network' as well.
library(igraph)
library(ggplot2)
library(dplyr)
library(maditr)
library(textshape)
library(bipartite)
library(sjstats)
library(Hmisc)
library(network)
library(nlme)
library(ergm.count)
library(tibble)
library(tidyverse)
library(scales)
library(patchwork)
library(ggridges)
library(ggpubr)



## Characterizing sequences
pol.data <- read.csv("Data/pollinator_sequences_all_weeks_new.csv", sep=",")

sum.data <- pol.data %>% dplyr::select(Obs_id, Pollinator_sp) %>% unique()
pol.yes <- names(which(table(sum.data$Pollinator_sp) > 3))

pol.data %<>% filter(Pollinator_sp %in% pol.yes)


# number of plant species visited (S)
dat2 <- pol.data %>% group_by(Obs_id, Pollinator_sp, Week) %>% 
  summarise(S= n_distinct(Plant_sp)) 
# Total sequence length (L)
dat3 <- pol.data %>% group_by(Obs_id, Pollinator_sp, Week) %>% 
  summarise(L=n())
# Conspecific sequence length (L_C)
dat5 <- pol.data %>% group_by(Obs_id, Pollinator_sp, Week) %>% 
  reframe(L_C=sequence(rle(as.character(Plant_sp))$lengths)) %>%
  #filter(consp_visits>1) %>% 
  group_by(Obs_id, Pollinator_sp) %>% filter(L_C == max(L_C))
# number of heterospecific jumps (J)
dat6 <- pol.data %>%
  arrange(Obs_id, Pollinator_sp, Step) %>% 
  group_by(Obs_id, Pollinator_sp) %>%
  mutate(
    Plant_switch = Plant_sp != lag(Plant_sp, default = first(Plant_sp))
  ) %>%
  summarise(
    J = sum(Plant_switch, na.rm = TRUE),
    .groups = "drop"
  )

test <- left_join(dat2, dat3) %>% left_join(dat5) %>% left_join(dat6) %>% 
  mutate(Obs_id=as.numeric(Obs_id)) %>% filter(L>1)
glimpse(test)

head(test)
length(which(test$L_C>=2)) / nrow(test)
sd(test$L_C)

palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)

ggplot(test, aes(x=J, y=L_C, colour=Pollinator_sp)) +
  geom_point(size=2) +
  scale_colour_manual(values = palette) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()


mean(test$L)
sum(test$J < 5)/nrow(test)


test.long <- gather(test, variable, value, S:J, factor_key=TRUE)


labs <- c("Plant species \nrichness (S)", 
          "Total sequence \nlength (L)",
          "Conspecific \nsequence length (L_C)",
          "Heterospecific \njumps (J)")
names(labs) <- c("S", "L", "L_C", "J")



sp.labs <- rev(c("Xca", "Msi", "Las", "Eru", "Dci", "Bom", "Ame", "Adi", "Afl", "Ahi"))


# Fig. 2

fig.2 <- ggplot(test.long,
                  aes(x = value, y = Pollinator_sp,
                      fill = Pollinator_sp, colour = Pollinator_sp)) +
  geom_density_ridges(scale = 1, alpha = 0.8, stat = "binline", bins = 40) +
  facet_wrap(~variable, ncol = 4, scales = "free_x",
             labeller = labeller(variable = labs),
             strip.position = "bottom") +
  scale_fill_manual(values = palette) +
  scale_colour_manual(values = palette) +
  theme_bw() +
  #labs(tag="A") +
  guides(fill = "none", colour = "none") +
  theme(
    strip.background = element_rect(
      color = "white", fill = "white", size = 1.5, linetype = "solid"
    ),
    strip.text.x = element_text(size = 12, color = "black"),
    axis.text=element_text(size=12),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    plot.tag = element_text(face="bold", size=18)) + 
  scale_y_discrete(labels= sp.labs)

ggsave(
  filename = "Figure_2.pdf",
  plot = fig.2,
  path = "Output",
  width = 8,
  height = 6,
  units = "in"
)

## Fig. S1

# Perform PCA (center & scale is usually recommended)
pca_res <- prcomp(test[, c("S","L","L_C","J")], center = TRUE, scale. = TRUE)

# Extract scores
scores <- as.data.frame(pca_res$x)
scores$Pollinator_sp <- test$Pollinator_sp

# Extract loadings for arrows
loadings <- as.data.frame(pca_res$rotation)
loadings$var <- rownames(loadings)

# Scale arrows to fit nicely on plot
arrow_multiplier <- min(
  (max(scores$PC1) - min(scores$PC1)) / (max(loadings$PC1) - min(loadings$PC1)),
  (max(scores$PC2) - min(scores$PC2)) / (max(loadings$PC2) - min(loadings$PC2))
) * 0.7

loadings <- loadings %>%
  mutate(PC1 = PC1 * arrow_multiplier,
         PC2 = PC2 * arrow_multiplier)

summary(pca_res)
# Plot


fig.s1 <- ggplot(scores, aes(x = PC1, y = PC2, color = Pollinator_sp)) +
  geom_point(size = 3, alpha=0.7, position = position_dodge(2.5)) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black") +
  geom_text(data = loadings,
            aes(x = PC1, y = PC2, label = var),
            vjust = 1.2, hjust=-0.5, color = "black") +
  theme_bw() +
  labs(x = "PC1 (59.6%)", y = "PC2 (33.5%)"
       #, tag = "B"
  ) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.tag = element_text(face="bold", size=18)) +
  #guides(color="none") + 
  scale_color_manual(values=palette,
                     labels = rev(c("Xylocopa cantabrita (Xca)", 
                                    "Megachile sicula (Msi)", 
                                    "Lasioglossum sp. (Las)", 
                                    "Eucera rufa (Eru)", 
                                    "Dasypoda cingulata (Dci)", 
                                    "Bombilius sp. (Bom)", 
                                    "Apis mellifera (Ame)", 
                                    "Anthophora dispar (Adi)", 
                                    "Andrena flavipes (Afl)", 
                                    "Andrena hispania (Ahi)")))

ggsave(
  filename = "Figure_S1.pdf",
  plot = fig.s1,
  path = "Output",
  width = 8,
  height = 6,
  units = "in"
)

# Compute quality and quantity interaction components

pol.data <- read.csv("Data/pollinator_sequences_all_weeks.csv", sep=",")

pol.data %>% dplyr::select(Obs_id, Pollinator_sp) %>% 
  unique() %>% group_by(Pollinator_sp) %>% summarise(n=n())
#

sum.data <- pol.data %>% dplyr::select(Obs_id, Pollinator_sp) %>% unique()
pol.yes <- names(which(table(sum.data$Pollinator_sp) > 3))

pol.data %<>% filter(Pollinator_sp %in% pol.yes)

# Step 1: Order by Obs_id and Step to respect the sequence
pol.data <- pol.data %>%
  arrange(Obs_id, Step)

# Step 2: Identify contiguous runs of Plant_sp within each Obs_id
pol.data <- pol.data %>%
  group_by(Obs_id) %>%
  mutate(run_group = with(rle(Plant_sp), rep(seq_along(lengths), lengths)))

# Step 3: Summarize the runs
run_summary <- pol.data %>%
  group_by(Obs_id, run_group, Plant_sp) %>%
  summarise(run_length = n(), .groups = "drop")

# Step 4: Keep only runs of length >= 2 (i.e. sequences that count toward quality)
filtered_runs <- run_summary %>%
  filter(run_length >= 1)

# Step 5: Compute average run length per Obs_id and Plant_sp → this is "quality"
quality_summary <- filtered_runs %>%
  group_by(Obs_id, Plant_sp) %>%
  summarise(quality = mean(run_length), .groups = "drop")

# Step 6: Compute quantity as total number of times the species appears per Obs_id
quantity_summary <- pol.data %>%
  group_by(Obs_id, Plant_sp) %>%
  summarise(quantity = n(), .groups = "drop")

# Step 7: Merge both and join back to original data
pol.data <- pol.data %>%
  left_join(quantity_summary, by = c("Obs_id", "Plant_sp")) %>%
  left_join(quality_summary, by = c("Obs_id", "Plant_sp")) %>%
  select(-run_group) %>%
  mutate(quality = replace_na(quality, 0)) %>%
  unique()


pol.data %<>% dplyr::select(Obs_id, Pollinator_sp, Plant_sp, quantity, quality) %>% unique()


summary_df <- pol.data %>%
  group_by(Plant_sp, Pollinator_sp) %>%
  summarise(
    total_quantity = sum(quantity),
    mean_quantity = mean(quantity, na.rm = TRUE),
    sd_quantity = sd(quantity, na.rm = TRUE),
    mean_quality = mean(quality, na.rm = TRUE),
    sd_quality = sd(quality, na.rm = TRUE),
    .groups = "drop"
  )


# Plot network of interactions with quality and quantity components

interaction_matrix <- reshape2::acast(summary_df, 
                                      Plant_sp ~ Pollinator_sp, 
                                      value.var = "total_quantity",
                                      fill = 0)

# Create a matching matrix for mean_quality (same shape as interaction_matrix)
quality_matrix <- reshape2::acast(summary_df, 
                                  Plant_sp ~ Pollinator_sp, 
                                  value.var = "mean_quality", fill = NA)

# Normalize mean_quality for color mapping (0 to 1)
norm_quality <- (quality_matrix - min(quality_matrix, na.rm = TRUE)) / 
  (max(quality_matrix, na.rm = TRUE) - min(quality_matrix, na.rm = TRUE))

norm_quality <- quality_matrix / max(quality_matrix, na.rm = TRUE)

# 1. Build edge list explicitly
idx <- which(interaction_matrix > 0, arr.ind = TRUE)

edges <- data.frame(
  from   = rownames(interaction_matrix)[idx[, 1]],
  to     = colnames(interaction_matrix)[idx[, 2]],
  weight = interaction_matrix[idx],
  quality = norm_quality[idx]
)

# 2. Build graph from edge list
g <- graph_from_data_frame(edges, directed = FALSE)

# 3. Degree and node attributes
deg <- igraph::degree(g)

V(g)$size  <- 10 + deg * 2
V(g)$shape <- "circle"
V(g)$color <- c(rep("grey90", 19), palette)
V(g)$frame <- c(rep("grey40", 19), rep("white", 10))

# 4. Vertex labels
V(g)$name <- c(
  rep("", 19),
  rev(c("Xca", "Msi", "Las", "Eru", "Dci",
        "Bom", "Ame", "Adi", "Afl", "Ahi"))
)

# 5. Edge width = interaction strength
w <- E(g)$weight

w_scaled <- sqrt(w)   # log(1 + w) avoids log(0)

E(g)$width <- w_scaled / max(w_scaled) * 10


# 6. Edge color = quality
q <- E(g)$quality

q_scaled <- (q - min(q, na.rm = TRUE)) /
  (max(q, na.rm = TRUE) - min(q, na.rm = TRUE))

col_fun <- colorRamp(c("lightblue", "blue4"))

cols <- col_fun(q_scaled)

E(g)$color <- rgb(cols[,1], cols[,2], cols[,3],
                  maxColorValue = 255)

# Optional transparency (recommended for dense networks)
E(g)$color <- adjustcolor(E(g)$color, alpha.f = 0.7)

# 7. Layout
layout_circ <- layout_in_circle(g)

# 8. Plot (Fig.4B)
plot(
  g,
  layout = layout_circ,
  vertex.label.cex = 1.2,
  vertex.label.color = "white",
  vertex.frame.color = V(g)$frame,
  vertex.frame.width = 2,
  vertex.label.font = 2,
  edge.color = E(g)$color
)



# 1) Export to GML
# GML will store vertex/edge attributes if present.
# write_graph(g, file = "graph.gml", format = "gml")

# 2) Export to GraphML
# GraphML preserves attributes and is good for interoperability with many tools.
# write_graph(g, file = "graph.graphml", format = "graphml")

# 3) Export to Pajek (.net)
# write_graph(g, file = "graph.net", format = "pajek")

# 4) Export to an edge list (two ways)
# (a) igraph native: simple edgelist format
# write_graph(g, file = "graph.edgelist", format = "edgelist")



# Compute degree of pollinator species

pollinator_plant_counts <- summary_df  %>%  # Make sure each plant-pollinator pair is counted once
  group_by(Pollinator_sp) %>%
  summarise(degree = n(), 
            mean_quality_all=mean(mean_quality, na.rm = TRUE),
            sd_quality_all=sd(mean_quality, na.rm = TRUE),
            mean_quantity_all=mean(total_quantity, na.rm = TRUE)) %>%
  arrange(desc(degree))  # optional: to sort by number of plant species


degree.quality <- ggplot(pollinator_plant_counts, aes(degree, mean_quality_all,
                                                      color=Pollinator_sp)) +
  geom_smooth(method="lm", color="grey40") +
  geom_errorbar(aes(ymin = mean_quality_all - sd_quality_all,
                    ymax = mean_quality_all + sd_quality_all), width = 0.,
                position=position_dodge(width=0.5)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  theme_bw() + theme(text = element_text(size=16)) +
  ylab("Pollinator species' \nmean interaction quality") + 
  xlab("Pollinator \nspecies' degree") +
  guides(color="none", size="none") + 
  scale_color_manual(values=palette)

### quality versus quantity

subset.pol.data <-pol.data
summary_df <- subset.pol.data %>%
  group_by(Plant_sp, Pollinator_sp) %>%
  summarise(
    total_quantity = sum(quantity),
    mean_quantity = mean(quantity, na.rm = TRUE),
    sd_quantity = sd(quantity, na.rm = TRUE) / sqrt(n()),
    mean_quality = mean(quality, na.rm = TRUE),
    sd_quality = sd(quality, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Create a grid of quantity and quality values
quantity_range <- seq(min(summary_df$total_quantity) * 0.8, max(summary_df$total_quantity) * 1.2, length.out = 1500)
quality_range  <- seq(min(summary_df$total_quantity) * 0.8, max(summary_df$total_quantity) * 1.2, length.out = 1500)

grid <- expand.grid(
  quantity = quantity_range,
  quality = quality_range
) %>%
  mutate(effectiveness = quantity * quality)


# Fig. 4A

fig.4a <- ggplot(summary_df, aes(x = total_quantity, y = mean_quality, color=Pollinator_sp)) +
  # Isoclines (contours of constant effectiveness)
  #geom_contour(data = grid, aes(x = quantity, y = quality, z = effectiveness), 
  #             breaks = c(5, 10, 20, 40, 80, 160),  # choose your own meaningful levels
  #             colour = "gray60", linetype = "solid") +
  
  # vertical error bars (quality)
  geom_errorbar(aes(ymin = mean_quality - sd_quality,
                    ymax = mean_quality + sd_quality), width = 0.1,
                position=position_dodge(width=0.5)) +
  geom_point(size = 3, alpha = 0.8, position=position_dodge(width=0.5)) +
  labs(x = "Pairwise species \ninteraction quantity (QT)", 
       y = "Pairwise species \ninteraction quality (QL)") +
  theme_bw() +
  theme(text=element_text(size=17),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_x_log10() +
  scale_color_manual(values=palette,
                     labels = rev(c("Xylocopa cantabrita (Xca)", 
                                    "Megachile sicula (Msi)", 
                                    "Lasioglossum sp. (Las)", 
                                    "Eucera rufa (Eru)", 
                                    "Dasypoda cingulata (Dci)", 
                                    "Bombilius sp. (Bom)", 
                                    "Apis mellifera (Ame)", 
                                    "Anthophora dispar (Adi)", 
                                    "Andrena flavipes (Afl)", 
                                    "Andrena hispania (Ahi)")))

cor.test(summary_df$total_quantity, summary_df$mean_quality, method="spearman")

ggsave(
  filename = "Figure_4.pdf",
  plot = fig.4a,
  path = "Output",
  width = 8,
  height = 6,
  units = "in"
)

