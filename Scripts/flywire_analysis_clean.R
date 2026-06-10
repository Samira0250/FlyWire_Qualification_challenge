# ============================================================
# FlyWire Connectome Challenge
# Shared Circuit Identification: FAFB + BANC + MCNS
# ============================================================

setwd("~/Downloads/wire")
library(tidyverse)
library(igraph)

# ============================================================
# SECTION 1: Load & Validate Edge Lists
# ============================================================

fafb <- read_csv("fafb_783_edge_list.csv")   %>% rename(from = `source neuron id`, to = `target neuron id`)
banc <- read_csv("banc_626_edge_list.csv")   %>% rename(from = `source neuron id`, to = `target neuron id`)
manc <- read_csv("manc_1.2.1_edge_list.csv") %>% rename(from = `source neuron id`, to = `target neuron id`)
maol <- read_csv("maol_1.1_edge_list.csv")   %>% rename(from = `source neuron id`, to = `target neuron id`)
mcns <- read_csv("mcns_0.9_edge_list.csv")   %>% rename(from = `source neuron id`, to = `target neuron id`)

# Sanity check table
sanity <- tibble(
  dataset        = c("fafb", "banc", "manc", "maol", "mcns"),
  edges          = c(nrow(fafb), nrow(banc), nrow(manc), nrow(maol), nrow(mcns)),
  unique_neurons = c(
    n_distinct(c(fafb$from, fafb$to)),
    n_distinct(c(banc$from, banc$to)),
    n_distinct(c(manc$from, manc$to)),
    n_distinct(c(maol$from, maol$to)),
    n_distinct(c(mcns$from, mcns$to))
  )
)
print(sanity)
write_csv(sanity, "output_01_dataset_summary.csv")

# ============================================================
# SECTION 2: Build Directed Graphs
# ============================================================

g_fafb <- graph_from_data_frame(fafb, directed = TRUE)
g_banc <- graph_from_data_frame(banc, directed = TRUE)
g_manc <- graph_from_data_frame(manc, directed = TRUE)
g_maol <- graph_from_data_frame(maol, directed = TRUE)
g_mcns <- graph_from_data_frame(mcns, directed = TRUE)

graph_summary <- data.frame(
  dataset  = c("fafb", "banc", "manc", "maol", "mcns"),
  vertices = c(vcount(g_fafb), vcount(g_banc), vcount(g_manc), vcount(g_maol), vcount(g_mcns)),
  edges    = c(ecount(g_fafb), ecount(g_banc), ecount(g_manc), ecount(g_maol), ecount(g_mcns))
)
print(graph_summary)
write_csv(graph_summary, "output_02_graph_summary.csv")

# ============================================================
# SECTION 3: Degree Distribution
# ============================================================

get_degrees <- function(g, name) {
  tibble(
    dataset = name,
    neuron  = V(g)$name,
    in_deg  = degree(g, mode = "in"),
    out_deg = degree(g, mode = "out"),
    total   = degree(g, mode = "all")
  )
}

deg_fafb <- get_degrees(g_fafb, "fafb")
deg_banc <- get_degrees(g_banc, "banc")
deg_manc <- get_degrees(g_manc, "manc")
deg_maol <- get_degrees(g_maol, "maol")
deg_mcns <- get_degrees(g_mcns, "mcns")

# Summary table
degree_summary <- bind_rows(deg_fafb, deg_banc, deg_manc, deg_maol, deg_mcns) %>%
  group_by(dataset) %>%
  summarise(
    median_degree    = median(total),
    mean_degree      = round(mean(total), 1),
    max_degree       = max(total),
    neurons_deg1     = sum(total == 1),
    neurons_deg5plus = sum(total >= 5),
    .groups = "drop"
  )
print(degree_summary)
write_csv(degree_summary, "output_03_degree_summary.csv")

# Degree distribution plot — all 5 datasets
bind_rows(deg_fafb, deg_banc, deg_manc, deg_maol, deg_mcns) %>%
  ggplot(aes(x = total, fill = dataset)) +
  geom_histogram(bins = 60, show.legend = FALSE) +
  facet_wrap(~dataset, scales = "free") +
  scale_x_log10() +
  labs(
    title   = "Degree Distribution per Dataset",
    subtitle = "Log scale; each panel independent y-axis",
    x       = "Total degree (log scale)",
    y       = "Neuron count"
  ) +
  theme_minimal()
ggsave("output_04_degree_distribution.png", width = 12, height = 8, dpi = 300)

# ============================================================
# SECTION 4: Dataset Selection Rationale
# Note: FAFB + BANC + MCNS selected based on overlapping
# degree distributions (medians 40, 24, 46). MANC (median 335)
# and MAOL (median 183) are nerve cord / optic lobe datasets
# with incompatible degree profiles for direct matching.
# ============================================================

# ============================================================
# SECTION 5: Degree Signature Filtering
# A signature is the (in_deg, out_deg) pair per neuron.
# Neurons with identical signatures are candidate homologs.
# Filter range 3-200 captures the overlapping degree zone.
# ============================================================

sig_fafb <- deg_fafb %>% filter(total >= 3, total <= 200) %>% mutate(sig = paste(in_deg, out_deg, sep = "_"))
sig_banc <- deg_banc %>% filter(total >= 3, total <= 200) %>% mutate(sig = paste(in_deg, out_deg, sep = "_"))
sig_mcns <- deg_mcns %>% filter(total >= 3, total <= 200) %>% mutate(sig = paste(in_deg, out_deg, sep = "_"))

# Find shared signatures without joining neuron IDs (avoids memory explosion)
shared_sigs <- list(
  unique(sig_fafb$sig),
  unique(sig_banc$sig),
  unique(sig_mcns$sig)
) %>% Reduce(intersect, .)

cat("Shared degree signatures across FAFB + BANC + MCNS:", length(shared_sigs), "\n")
write_csv(tibble(sig = shared_sigs), "output_05_shared_signatures.csv")

# Filter each dataset to shared signatures only
sig_fafb <- sig_fafb %>% filter(sig %in% shared_sigs)
sig_banc <- sig_banc %>% filter(sig %in% shared_sigs)
sig_mcns <- sig_mcns %>% filter(sig %in% shared_sigs)

cat("Candidate neurons after filtering:\n",
    "FAFB:", nrow(sig_fafb), "\n",
    "BANC:", nrow(sig_banc), "\n",
    "MCNS:", nrow(sig_mcns), "\n")

# ============================================================
# SECTION 6: Identify Rare Signatures (max 1 neuron/dataset)
# These are the most tractable for isomorphism testing.
# ============================================================

counts <- tibble(sig = shared_sigs) %>%
  left_join(sig_fafb %>% count(sig, name = "n_fafb"), by = "sig") %>%
  left_join(sig_banc %>% count(sig, name = "n_banc"), by = "sig") %>%
  left_join(sig_mcns %>% count(sig, name = "n_mcns"), by = "sig") %>%
  mutate(max_count = pmax(n_fafb, n_banc, n_mcns)) %>%
  arrange(max_count)

write_csv(counts, "output_06_signature_counts.csv")
print(head(counts, 30))

# Perfect signatures: exactly 1 neuron per dataset
perfect_sigs <- counts %>% filter(max_count == 1)
cat("Perfect signatures (1 neuron per dataset):", nrow(perfect_sigs), "\n")

perfect_neurons <- perfect_sigs %>%
  left_join(sig_fafb %>% select(sig, neuron), by = "sig") %>% rename(neuron_fafb = neuron) %>%
  left_join(sig_banc %>% select(sig, neuron), by = "sig") %>% rename(neuron_banc = neuron) %>%
  left_join(sig_mcns %>% select(sig, neuron), by = "sig") %>% rename(neuron_mcns = neuron)

write_csv(perfect_neurons, "output_07_perfect_neurons.csv")
print(head(perfect_neurons, 10))

# ============================================================
# SECTION 7: Induced Subgraph from All Shared-Signature Neurons
# ============================================================

sub2_fafb <- induced_subgraph(g_fafb, vids = which(V(g_fafb)$name %in% as.character(sig_fafb$neuron)))
sub2_banc <- induced_subgraph(g_banc, vids = which(V(g_banc)$name %in% as.character(sig_banc$neuron)))
sub2_mcns <- induced_subgraph(g_mcns, vids = which(V(g_mcns)$name %in% as.character(sig_mcns$neuron)))

cat("Induced subgraphs (filtered candidates):\n",
    "FAFB:", vcount(sub2_fafb), "neurons,", ecount(sub2_fafb), "edges\n",
    "BANC:", vcount(sub2_banc), "neurons,", ecount(sub2_banc), "edges\n",
    "MCNS:", vcount(sub2_mcns), "neurons,", ecount(sub2_mcns), "edges\n")

# ============================================================
# SECTION 8: Translate Edges to Signature Space
# Replace neuron IDs with degree signatures so edges are
# comparable across datasets.
# ============================================================

lkp_fafb <- sig_fafb %>% select(neuron, sig) %>% distinct()
lkp_banc <- sig_banc %>% select(neuron, sig) %>% distinct()
lkp_mcns <- sig_mcns %>% select(neuron, sig) %>% distinct()

translate_edges <- function(g, lkp) {
  as_data_frame(g, what = "edges") %>%
    left_join(lkp, by = c("from" = "neuron")) %>% rename(sig_from = sig) %>%
    left_join(lkp, by = c("to"   = "neuron")) %>% rename(sig_to   = sig) %>%
    filter(!is.na(sig_from), !is.na(sig_to)) %>%
    mutate(sig_edge = paste(sig_from, sig_to, sep = "->"))
}

esig_fafb <- translate_edges(sub2_fafb, lkp_fafb)
esig_banc <- translate_edges(sub2_banc, lkp_banc)
esig_mcns <- translate_edges(sub2_mcns, lkp_mcns)

# Shared signature-level edges across all three
common_sig_edges <- Reduce(intersect, list(
  esig_fafb$sig_edge,
  esig_banc$sig_edge,
  esig_mcns$sig_edge
))

cat("Signature-level edges shared across all three datasets:", length(common_sig_edges), "\n")
write_csv(tibble(sig_edge = common_sig_edges), "output_08_common_sig_edges.csv")

# ============================================================
# SECTION 9: Build Shared Circuit Graph & Extract Largest
# Connected Component
# ============================================================

shared_edge_df <- tibble(
  from = sub("->.*", "", common_sig_edges),
  to   = sub(".*->", "", common_sig_edges)
)

g_shared <- graph_from_data_frame(shared_edge_df, directed = TRUE)

comps <- components(g_shared, mode = "weak")
cat("Connected components:", comps$no,
    "| Largest:", max(comps$csize), "nodes\n")

largest_nodes <- names(which(comps$membership == which.max(comps$csize)))
g_circuit     <- induced_subgraph(g_shared, vids = largest_nodes)

cat("Circuit (largest component):",
    vcount(g_circuit), "nodes,",
    ecount(g_circuit), "edges\n")

# ============================================================
# SECTION 10: Map Circuit Back to Neuron IDs
# One representative neuron per signature per dataset.
# This is valid: any neuron sharing a signature satisfies
# the degree-isomorphism condition.
# ============================================================

one_fafb <- sig_fafb %>% select(sig, neuron) %>% group_by(sig) %>% slice(1) %>% ungroup()
one_banc <- sig_banc %>% select(sig, neuron) %>% group_by(sig) %>% slice(1) %>% ungroup()
one_mcns <- sig_mcns %>% select(sig, neuron) %>% group_by(sig) %>% slice(1) %>% ungroup()

circuit_neurons <- tibble(sig = V(g_circuit)$name) %>%
  left_join(one_fafb, by = "sig") %>% rename(neuron_fafb = neuron) %>%
  left_join(one_banc, by = "sig") %>% rename(neuron_banc = neuron) %>%
  left_join(one_mcns, by = "sig") %>% rename(neuron_mcns = neuron)

cat("Final circuit neuron mapping:\n",
    "N =", nrow(circuit_neurons), "neurons\n",
    "FAFB mapped:", sum(!is.na(circuit_neurons$neuron_fafb)), "\n",
    "BANC mapped:", sum(!is.na(circuit_neurons$neuron_banc)), "\n",
    "MCNS mapped:", sum(!is.na(circuit_neurons$neuron_mcns)), "\n")

# ============================================================
# SECTION 11: Verify Isomorphism
# ============================================================

gf_fafb <- induced_subgraph(g_fafb, which(V(g_fafb)$name %in% circuit_neurons$neuron_fafb))
gf_banc <- induced_subgraph(g_banc, which(V(g_banc)$name %in% circuit_neurons$neuron_banc))
gf_mcns <- induced_subgraph(g_mcns, which(V(g_mcns)$name %in% circuit_neurons$neuron_mcns))

iso_check <- tibble(
  comparison   = c("FAFB vs BANC", "FAFB vs MCNS", "BANC vs MCNS"),
  n_nodes_A    = c(vcount(gf_fafb), vcount(gf_fafb), vcount(gf_banc)),
  n_nodes_B    = c(vcount(gf_banc), vcount(gf_mcns), vcount(gf_mcns)),
  edges_A      = c(ecount(gf_fafb), ecount(gf_fafb), ecount(gf_banc)),
  edges_B      = c(ecount(gf_banc), ecount(gf_mcns), ecount(gf_mcns)),
  isomorphic   = c(
    isomorphic(gf_fafb, gf_banc),
    isomorphic(gf_fafb, gf_mcns),
    isomorphic(gf_banc, gf_mcns)
  )
)
print(iso_check)
write_csv(iso_check, "output_09_isomorphism_check.csv")

# ============================================================
# SECTION 12: Save Submission CSV
# Format: 3 columns (one per dataset), N rows (matched neurons)
# ============================================================

submission <- circuit_neurons %>% select(neuron_fafb, neuron_banc, neuron_mcns)
write_csv(submission, "flywire_circuit_submission.csv")
cat("Submission saved:", nrow(submission), "rows x 3 columns\n")

# ============================================================
# SECTION 13: Circuit Visualization
# ============================================================

# Degree of each node in the circuit graph
V(g_circuit)$degree <- degree(g_circuit, mode = "all")

circuit_edge_df  <- as_data_frame(g_circuit, what = "edges")
circuit_node_df  <- as_data_frame(g_circuit, what = "vertices") %>%
  rename(sig = name) %>%
  mutate(
    in_deg  = as.integer(sub("_.*", "", sig)),
    out_deg = as.integer(sub(".*_", "", sig)),
    total   = in_deg + out_deg
  )

# Node color by total degree
ggplot() +
  geom_segment(
    data = circuit_edge_df,
    aes(
      x    = match(from, circuit_node_df$sig),
      xend = match(to,   circuit_node_df$sig),
      y    = 0, yend = 0
    ),
    alpha = 0.1, color = "gray60"
  ) +
  geom_point(
    data = circuit_node_df,
    aes(x = seq_len(nrow(circuit_node_df)), y = 0, color = total),
    size = 2
  ) +
  scale_color_viridis_c(name = "Total degree") +
  labs(
    title    = paste0("Shared Circuit: ", vcount(g_circuit), " neurons, ", ecount(g_circuit), " edges"),
    subtitle = "FAFB + BANC + MCNS | Nodes colored by degree",
    x = "Neuron index", y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("output_10_circuit_linear.png", width = 14, height = 4, dpi = 300)

# Circular layout using igraph
png("output_11_circuit_graph.png", width = 2400, height = 2400, res = 300)
plot(
  g_circuit,
  layout         = layout_with_fr(g_circuit),
  vertex.size    = 3,
  vertex.label   = NA,
  vertex.color   = adjustcolor("steelblue", alpha.f = 0.7),
  edge.arrow.size = 0.2,
  edge.color     = adjustcolor("gray50", alpha.f = 0.4),
  main           = paste0("Shared Circuit (N=", vcount(g_circuit), ")")
)
dev.off()

cat("\nAll outputs saved to ~/Downloads/wire/\n")
cat("Files generated:\n")
cat(" output_01_dataset_summary.csv\n")
cat(" output_02_graph_summary.csv\n")
cat(" output_03_degree_summary.csv\n")
cat(" output_04_degree_distribution.png\n")
cat(" output_05_shared_signatures.csv\n")
cat(" output_06_signature_counts.csv\n")
cat(" output_07_perfect_neurons.csv\n")
cat(" output_08_common_sig_edges.csv\n")
cat(" output_09_isomorphism_check.csv\n")
cat(" output_10_circuit_linear.png\n")
cat(" output_11_circuit_graph.png\n")
cat(" flywire_circuit_submission.csv\n")
