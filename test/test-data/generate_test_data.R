library(ape)

seed = 6868
taxa_sizes <- c(10, 50, 100, 200)

# Function to generate random tree with metadata
generate <- function(n_taxa, seed = 6868) {
  set.seed(seed)
  # tree
  raw_tree <- rtree(n = n_taxa, rooted = TRUE)
  add_random_bootstrap <- function(tree, min = 50, max = 100) {
      n_internal <- tree$Nnode  # number of internal nodes
      ufboot = sample(min:max, n_internal, replace = TRUE)
      shalrt = sample(min:max, n_internal, replace = TRUE)
      tree$node.label <- paste(ufboot, shalrt, sep = "/")
      return(tree)
  }
  tree = add_random_bootstrap(tree = raw_tree)


  # Define consistent mappings: each year gets one color and one shape
  years <- 2020:2024
  year_colors <- sample(colors(), length(years))
  year_shapes <- sample(c(21:25), length(years))
  names(year_colors) <- years
  names(year_shapes) <- years

  # Sample years for each taxon
  assigned_years <- sample(years, n_taxa, replace = TRUE)

  metadata <- data.frame(
    taxon = tree$tip.label,
    country = sample(c("USA", "UK", "Germany", "France", "Japan", "Australia", 
                       "Canada", "Brazil", "India", "South Africa"), 
                   n_taxa, replace = TRUE),
    vaccination_status = sample(c("vaccinated", "unvaccinated", "partial", "unknown"), 
                               n_taxa, replace = TRUE, 
                               prob = c(0.4, 0.3, 0.2, 0.1)),
    collection_year = assigned_years,
    collection_year_col = year_colors[as.character(assigned_years)],
    collection_year_shape = year_shapes[as.character(assigned_years)],
    epicluster = sample(paste0("cluster_", LETTERS[1:8]), n_taxa, replace = TRUE),
    epicluster_col = sample(colors(), n_taxa, replace = TRUE)
  )

  list(tree = tree, metadata = metadata)
}

# write out trees
for (n in taxa_sizes) {
    tm = generate(n_taxa = n)
    
    write.table(x = tm$metadata, file = paste0("test.", n, ".tsv"), sep = "\t", quote = T, row.names = F)
    write.tree(phy = tm$tree, file = paste0("test.", n, ".nwk"))
}

# write out a single nexus file for testing purposes
write.nexus(phy = tm$tree, file = paste0("test.", n, ".nexus"))
