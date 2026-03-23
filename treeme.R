#!/usr/bin/env Rscript

version = "0.0.1"
set.seed(1987)
options(warn = -1)
rlang::global_handle()

suppressPackageStartupMessages(
  invisible(
    lapply(
      c("ggnewscale", "optparse", "ggtree", "ggplot2", "treeio", "cowplot", "dplyr", "Cairo", "viridisLite"),
      require,
      character.only = TRUE
    )
  )
)

###############################
########## Functions ##########
###############################
start_message = paste0("Treeme.R - v", version, " Ammar Aziz \n")

is_var_empty = function(var) {
  # returns TRUE if the input is empty string, null
  
  if (is.null(var)) {
    return(TRUE)
  }
  if (is.character(var) && var == '') {
    return(TRUE)
  }
  if (length(var) > 0 & !is.null(var)) {
    return(FALSE)
  } else {
    logger(
      "What the hell did you pass to? This is a bug, go to github and submit an issue",
      "critical"
    )
  }
}

print_inputs = function(arguments) {
  logger("Input arguments:" ,simple = TRUE)
  arg_names = names(arguments)

  for (n in seq_along(arg_names)) {
    if (!is_var_empty(arguments[n]) & !(arg_names[n] == "help")) {
      message = paste0("--", arg_names[n], ": ", arguments[n])
      logger(message, "info", TRUE)
    }
  }
  logger(" ", simple = TRUE)

}

logger = function(text, level = "info", simple = FALSE) {
  start = format(Sys.time(), "%H:%M:%S")
  off = "\033[0m\n"
  codes = c(
    black = "\033[0;30m", # Black
    red = "\033[0;31m", # Red
    green = "\033[0;32m", # Green
    yellow = "\033[0;33m", # Yellow
    blue = "\033[0;34m", # Blue
    purple = "\033[0;35m", # Purple
    cyan = "\033[0;36m", # Cyan
    white = "\033[0;37m" # White
  )
  levels = c(
    info = codes[["blue"]],
    success = codes[["green"]],
    warning = codes[["yellow"]],
    critical = codes[["red"]]
  )

  if (simple) {
    cat(levels[level], text, off, file = stderr())
  } else {
    cat(levels[level], start, "|", toupper(level), "|", text, off, file = stderr())
  }
}

read_tree_auto = function(infile) {
  ext = strsplit(infile, ".", fixed = T)[[1]][-1]
  # newick
  if (length(ext) > 1) {
    ext = tail(ext, 1)
  }
  if (ext == "nwk" | ext == "newick") {
    tree_format <<- "newick"
    return(read.newick(infile))
  }
  if (ext == "nexus") {
    tree_format <<- "nexus"
    return(read.nexus(infile))
  }
  if (ext == "beast") {
    tree_format <<- "beast"
    return(read.beast(infile))
  }
  if (ext == "nhx") {
    tree_format <<- "nhx"
    return(read.nhx(infile))
  }
  if (ext == "json") {
    logger("Reading in nextstrain json, this can take a while!", "warning")
    tree_format <<- "json"
    return(read.nextstrain.json())
  }
}

reader = function(infile, type) {
  # safely read in files
  tryCatch(
    expr = {
      if (type == "tree") {
        t = read_tree_auto(infile)
        logger(paste0(basename(infile), " successfully read in."))
        return(t)
      }
      if (type == "tsv") {
        m = read.delim(
          infile,
          sep = "\t",
          na.strings = "",
          stringsAsFactors = F
        )
        logger(paste0(basename(infile), " successfully read in."))
        return(m)
      }
    },
    error = function(e) {
      logger(
        text = paste0(
          "Unable to read in file - unknown file type: ",
          infile,
          "\n",
          "\t",
          e
        ),
        level = "critical"
      )
      quit(status = 1)
    }
  )
}

check_var_in_meta = function(metafile, variable) {
  if (!(variable %in% names(metafile))) {
    logger(
      text = paste0(
        "Column '",
        variable,
        "' does not exist in metafile. Check inputs."
      ),
      level = "critical"
    )
    quit(status = 1)
  } else {
    return(TRUE)
  }
}

set_device = function(outfile) {
  ext = strsplit(outfile, ".", fixed = T)[[1]][-1]
  if (ext == "pdf") {
    return(cairo_pdf)
  }
  if (ext == "svg") {
    return("svg")
  }
}

get_page_size = function(size) {
  paper = list(
    "A2p" = c(594, 840),
    "A2l" = c(840, 594),
    "A3p" = c(297, 420),
    "A3l" = c(420, 297),
    "A4p" = c(210, 297),
    "A4l" = c(297, 210),
    "A5p" = c(148, 210),
    "A5l" = c(210, 148),
    "A6p" = c(105, 148),
    "A6l" = c(148, 105)
  )
  if (!any(size %in% names(paper))) {
    logger(
      paste0("Page size of: ", size, " is not accepted. Choose from: ", paste0(names(paper), collapse = ", ")),
      "critical"
    )
    quit(status = 1)
  }
  return(paper[[size]])
}

tree_title = function(title) {
  today = format(Sys.time(), format = "%d %b %Y")
  out = paste0(as.character(title), "\n", today)
}

###############################################
########### Tree helper functions #############
###############################################

calc_taxa_label_size = function(phylo, page_size, type = "logistic") {
  # phylo is is the output of func read_tree_auto
  # type controls how text size is calculated
  # area - page size aware font scaling
  # logistic - for scaling font size based on the number of taxa using sigmoid
  # inkdensity - total ink area of all tip labels constant relative to the usable plotting area

  if (class(phylo) == "treedata") {
    ntaxa = phylo@phylo$Nnode
  }
  if (class(phylo) == "phylo") {
    ntaxa = phylo$Nnode
  } else {
    logger("Unable to get num of tips in input tree. Check tree type is nwk or nhx", "critical")
    quit(status = 1)
  }

  width = page_size[1]
  height = page_size[2]
  min_font = 1
  max_font = 5

  if (type == "area") {
    mid_vir = max_font / 2
    font_size = min_font + (max_font - min_font) / (1 + exp(0.1 * (mid_vir - ntaxa)))
  } else if (type == "logistic") {
    # Kimi K2 hallucination
    # Heuristic: font size decreases with more labels and increases with page area
    base_size = 0.5
    area = width * height
    font_size = base_size * sqrt(area) / (20 * log(ntaxa + 1))
    font_size = max(min_font, min(font_size, max_font))
  } else if (type == "inkdensity") {
    # Kimi K2 hallucination
    # Heuristic: scale font size linearly so that the total ink used by all tip labels stays roughly constant.
    usable_frac = 0.80 # amount of plotting area that is useable - this is about 0.8 assuming plotting with a heatmap
    usable_area <- usable_frac * width * height
    char_area <- 0.65 # target aggregate label area (empirical constant 0.35 mm² per character)
    avg_label_len <- 20 # default average label length in characters
    font_size <- sqrt(usable_area / (ntaxa * char_area * avg_label_len))
    font_size <- max(min_font, min(max_font, font_size)) # clamp to readable range
  } else {
    logger("Internal bug - report this issue on github", "critical")
  }
  return(font_size)
}

calc_shape_size = function(phylo, page_size, k = 0.01) {
  if (class(phylo) == "treedata") {
    ntaxa = phylo@phylo$Nnode
  }
  if (class(phylo) == "phylo") {
    ntaxa = phylo$Nnode
  } else {
    logger("Unable to get num of tips in input tree. Check tree type is nwk or nhx", "critical")
    quit(status = 1)
  }
  min_size = 0.1
  max_size = 3
  width = page_size[1]
  height = page_size[2]

  area = width * height
  size = k * (area / ntaxa)
  return(max(min_size, min(size, max_size)))
}

calc_offset = function(phylo, page_size) {
  return(0.01)
}

check_taxa_names = function(tree, meta_desig) {
  # tree - must be of object type X and Y
  # meta_desig :  vector list, input should be like metadata[, 'designation']

  # via treeio::read.beast or ape::read.nexus
  if (class(tree) == "treedata") {
    taxa_labels = tree@phylo$tip.label
  }
  # read in via treeio::read.newick, treeio::read.iqtree
  if (class(tree) == "phylo") {
    taxa_labels = tree$tiplabs
  } else {
    logger("Unable to extract tip labels from input tree. Check tree type is nwk or nhx", "critical")
    quit(status = 1)
  }

  mismatch = taxa_labels[!taxa_labels %in% meta_desig]

  if (!length(mismatch)) {
    logger(text = "All tip labels present in metafile")
  } else {
    logger(
      text = "The following tip labels not present in metafile:",
      level = "warning"
    )
    cat(paste("\t", mismatch, collapse = "\n"))
    cat("\n")
  }
}

add_clades = function(cladesFile, tree_data, plot_dim_x) {
  # clade file must be a dataframe with vars:
  #   node_number
  #   clade_name
  # tree_data  :  tree@data from read.beast

  create_clade = function(node_number, clade_name, offset) {
    geom_cladelabel(
      node_number,
      label = clade_name,
      align = FALSE,
      angle = 270,
      hjust = "center",
      offset = offset,
      offset.text = plot_dim_x * 0.01,
      barsize = 1,
      fontsize = 5,
      extend = 0.2
    )
  }

  get_clade_node = function(muts, tree_data) {
    # returns node number given a mutation

    ind = which(tree_data$aa_muts %in% muts)
    if (length(ind) == 0) {
      return(NA_integer_)
    }
    node = tree_data[[ind, "node"]]
    return(node)
  }

  cladesFile$node_number = apply(
    cladesFile %>% select(mutations),
    1,
    get_clade_node,
    tree_data
  )
  cladesFile = cladesFile[complete.cases(cladesFile), ]

  num_clade = nrow(cladesFile)
  min_offset = plot_dim_x * 0.02
  from = plot_dim_x - (plot_dim_x * 0.80)
  to = plot_dim_x + (plot_dim_x * 0.30)
  offset_values = seq(from, to, by = min_offset)[1:3]
  offset = rep(offset_values, ceiling(num_clade / 3))[1:num_clade]

  return(
    Map(
      create_clade,
      as.integer(cladesFile$node_number),
      cladesFile$clade_name,
      offset
    )
    # pmap(
    #   list(
    #     as.integer(cladesFile$node_number),
    #     cladesFile$clade_name,
    #     offset
    #   ),
    #   create_clade
    # )
  )
}

builder_tiplab = function(tplot, taxa_size, offset, ucolors, color_var = NULL) {
  if (!is.null(color_var)) {
    logger(paste0("Adding geom_tiplab, coloring by ", color_var), "info")
    tplot = tplot +
      #new_scale_color() +
      geom_tiplab(
        aes(color = !!sym(color_var)),
        size = taxa_size,
        offset = offset,
        family = "Arial"
      ) +
      scale_color_manual(
        values = ucolors,
        na.translate = FALSE
      ) +
      guides(
        color = guide_legend(
          override.aes = list(
            label = "\u25A0",
            size = 10,
            linetype = 3
          )
        )
      )
  } else {
    logger("Adding taxa labels, no color specified", "info")
    tplot = tplot +
      geom_tiplab(
        geom = "text",
        size = taxa_size,
        offset = offset,
        family = "Arial"
      )
  }
  return(tplot)
}

builder_tippoint = function(tplot, fill_by, shape_by, shape_colors, shape_maps, shape_size) {
  logger(paste0("Adding tippoint shapes, variable: ", shape_by), "info")
  tplot = tplot +
    #new_scale_color() +
    geom_tippoint(
      size = shape_size,
      stroke = 0.2,
      aes(
        # convert to factor protects against error
        # "Continuous value supplied to a discrete scale."
        fill = as.factor(!!sym(fill_by)), 
        shape = as.factor(!!sym(shape_by))
      ),
    ) +
    scale_fill_manual(
      fill_by,
      values = shape_colors,
      na.value = "grey50"
    ) +
    scale_shape_manual(
      "Disabled",
      values = shape_maps
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          size = 5,
          label = "",
          shape = unname(shape_maps)
        )
      ),
      shape = "none"
    )
  return(tplot)
}

builder_bootstrap = function(tplot, format, size) {
  logger("Adding bootstrap values to branches", "info")

  if (format == "nhx") {
    bootstrap_label = "B"
  }
  if (format == "newick") {
    bootstrap_label = "label"
  } else {
    logger(paste0("Plotting bootstrap values are not implemented for tree format: ", bootstrap_label), "warning")
    quit(status = 1)
  }

  tplot = tplot +
    geom_text(
      data = td_filter(!isTip),
      aes(
        x = branch,
        label = !!sym(bootstrap_label)
      ),
      size = size * 0.75,
      vjust = -0.3
    )
  return(tplot)
}

builder_heatmap = function(tplot, meta, tstart, tend, ufills, offset) {
  logger("Adding heatmap to phylotree", "info")
  #tplot = tplot + new_scale_fill()
  tdata = meta[, c(tstart, tend)]
  rownames(tdata) = meta[, 1]

  # tplot = gheatmap(p = tplot, data = table, offset = offset, colnames = FALSE, legend_title = "") +
  #   scale_x_ggtree() +
  #   scale_y_continuous(expand = c(0, 0.3)) +
  #   scale_fill_manual()

  return(ggtree::gheatmap(tplot, tdata))
}

master_builder = function(
  phylo,
  meta,
  args,
  t_color,
  p_shape,
  p_color,
  taxa_size,
  shape_size,
  taxa_offset,
  tree_format
) {
  tplot = ggtree(tree, linewidth = 0.1) %<+% meta

  tplot = builder_tiplab(
    tplot = tplot,
    color_var = arguments$`color-by`,
    taxa_size = taxa_size,
    offset = taxa_offset,
    ucolors = t_color
  )

  if (!is_var_empty(arguments$`shape-by`)) {
    tplot = builder_tippoint(
      tplot = tplot,
      fill_by = arguments$`shape-by`,
      shape_by = arguments$`shape-by`,
      shape_colors = p_color,
      shape_map = p_shape,
      shape_size = shape_size
    )
  }

  if (isTRUE(arguments$bootstrap)) {
    tplot = builder_bootstrap(
      tplot = tplot,
      format = tree_format,
      size = taxa_size
    )
  }

  if (!anyNA(ufill_map)) {
    tplot = builder_heatmap(tplot = tplot, meta = meta, offset = 0.01, tstart = tstart, tend = tend)
  }

  tplot = tplot +
    theme(
      legend.position = c(0.1, 0.85),
      legend.key.size = unit(5, "mm"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0, "mm"),
      legend.spacing.y = unit(0, "mm"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 15),
      plot.title = element_text(hjust = 0.06, vjust = -15, size = 20),
      plot.subtitle = element_text(hjust = 0.02, vjust = -12, size = 20)
    )

  # Fix tip label clipping
  plot_dim_x = ggplot_build(tplot)$layout$panel_scales_x[[1]]$range$range[2]
  tplot = tplot +
    coord_cartesian(clip = "off", expand = FALSE) +
    xlim(NA, ((0.40 * plot_dim_x) + plot_dim_x))

  return(tplot)
}

create_umaps = function() {}

#########################################
############# CLI Parser ################
#########################################

option_list = list(
  make_option(
    c("-t", "--tree"),
    help = "Required: input tree - auto detect extension.",
    action = "store",
    type = "character",
    default = NA
  ),
  make_option(
    c("-o", "--output"),
    help = "Required: output name - must end in pdf or svg",
    action = "store",
    type = "character",
    default = NA
  ),
  make_option(
    c("-m", "--meta"),
    help = "Required: tsv file containing metafile.",
    action = "store",
    type = "character",
    default = NA
  ),
  make_option(
    c("-c", "--color-by"),
    help = "Optional: Column name in metafile to control the color of taxa labels. If not provided, tiplabs are black.",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-s", "--shape-by"),
    help = "Optional: Column name in metafile to control the shape of tip points. If not provided, tips are blank (no shape).",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-v", "--heatmap-value-range"),
    help = "Optional: Specify the columns in the metadata containing the VALUES for plotting the heatmap. Format required: 'X:Y'.",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-f", "--heatmap-fill-range"),
    help = "Optional: Specify the columns in the metadata containing the FILLS for plotting the heatmap.  Format required: 'X:Y'",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-b", "--bootstrap"),
    help = "Optional: Show boostrap values on branches. Note: Treeme does not show terminal bootstrap values.",
    action = "store_true"
  ),
  make_option(
    c("--clades-file"),
    help = "Optional: TSV file - 'clade\\tnode number' - used for vertical bar lines beside tree to identify clades.",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-i", "--title"),
    help = "Optional: Optional: The title of the final output",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-p", "--paper-size"),
    help = "Optional: Output size - options: A3p, A3l, A4l, [A4p] - p/l is the orientation. [default %default]",
    action = "store",
    type = "character",
    default = "A4p"
  ),
  make_option(
    c("--taxa-font-size"),
    help = "Optional: Specify the taxa labels font size. Set to 0 to turn off taxa labels. Treeme will try to calculate the best value.",
    action = "store",
    type = "numeric"
  ),
  make_option(
    c("--tippoint-shape-size"),
    help = "Optional: Specify the shape size. Set to 0 to turn off tippoint shapes. Treeme will try to calculate the best value.",
    action = "store",
    type = "numeric"
  )
)

parser = OptionParser(
  epilogue = c(
    "Metafile must be a tab separated file, with the first column containing the sample id matching taxa labels on tree.",
    "All subsequent columns are optional - except those specified with --color-by and --shape-by arguments.",
    "\n",
    "Common Issues and solutions:",
    "\t- Text is squished together; Increase output size from A4 to A3.",
    "\t- A warning appears about taxa/metafile labels; open the tree file in a text editor and check the label names. These much match exactly, no spaces, no underscores."
  ),
  option_list = option_list,
  usage = c(
    "Basic; treeme.R -t tree.nwk -m metadata.tsv -o tree.pdf",
    "With taxa colored; treeme.R -t tree.nwk -m metadata.tsv -o tree.pdf -c state",
    "With tip point shapes; treeme.R -t tree.nwk -m metadata.tsv -o tree.pdf -c {var} -s {var}"
  )
)

tryCatch(
  expr = {
    arguments = parse_args(
      object = parser,
      positional_arguments = TRUE
    )$options
  },
  finally = {
    if (any(is.na(arguments))) {
      print_help(parser)
      cat("\n")
      message = c(
        "Missing arguments:",
        paste0("--", names(which(is.na(arguments))))
      )
      logger(message, "critical", TRUE)
      cat("\n")
      quit(status = 1)
    }
  }
)

# for manual testing
if (interactive()) {
  arguments = list(
    metadata = "test/test-data/test.10.tsv",
    tree = "test/test-data/test.10.nwk",
    `color-by` = "epicluster",
    `shape-by` = "collection_year",
    `paper-size` = "A3p"
  )
}

#########################################
############# Input checks ##############
#########################################
logger(start_message, "success", simple=TRUE)

tree = reader(arguments$tree, "tree")
metadata = reader(arguments$meta, "tsv")

# set taxa label colors
if (!is_var_empty(arguments$`color-by`)) {
  if (arguments$`color-by` %in% colnames(metadata)) {
    # check if a column called X_col exists
    taxa_col_name = paste0(arguments$`color-by`, "_col")
    if (taxa_col_name %in% colnames(metadata)) {
      t_color_map = setNames(
        metadata[, taxa_col_name],
        metadata[, arguments$`color-by`]
        )
    } else {
      # set to viridis colors
      t_color_map = setNames(
        viridis(length(metadata[, arguments$`color-by`]), option = "D"), 
        unique(metadata[, arguments$`color-by`])
        )
    }
  } else {
    logger(paste0("The column ", arguments["color-by"] ," specified by `--color-by` does not exist in metafile"), "critical")
    quit(status = 1)
  }
} else {
  # here we set t_color_map to empty vector
  # which results in default color (black)
  t_color_map = c()
}

# heatmap
heatmap_args = sum(!is_var_empty(arguments$`heatmap-fill-range`), !is_var_empty(arguments$`heatmap-value-range`))
if (heatmap_args == 2) {
  heatmap_fill_range = as.numeric(
    unlist(
      strsplit(x = arguments$`heatmap-fill-range`, split = ":")
      )
    )
  if ((length(heatmap_fill_range) != 2) & any(is.numeric(heatmap_fill_range))) {
    logger("Can not parse the heatmap fill range. It must be in this format: X:Y eg 5:7 (inclusive)", "critical")
    quit(status = 1)
  }

  heatmap_value_range = as.numeric(unlist(strsplit(x = arguments$`heatmap-value-range`, split = ":")))
  if ((length(heatmap_value_range) != 2) & any(is.numeric(heatmap_value_range))) {
    logger("Can not parse the heatmap category range. It must be in this format: X:Y eg 5:7 (inclusive)", "critical")
    quit(status = 1)
  }
  tstart = heatmap_value_range[1]
  tend = heatmap_value_range[2]
  fstart = heatmap_fill_range[1]
  fend = heatmap_fill_range[2]
  ufill_map = setNames(c(metadata[, fstart], metadata[, fend]), c(metadata[, tstart], metadata[, tend]))
} else if (heatmap_args == 1) {
  logger("Both --heatmap-value-range and --heatmap-fill-range are required.", "critical")
  quit(status = 1)
} else {
  ufill_map = NA
}

# set geom point shape and color
if (!is_var_empty(arguments$`shape-by`)) {
  if (arguments$`shape-by` %in% colnames(metadata)) {
    
    # shapes
    p_shape_name = paste0(arguments$`shape-by`, "_shape")
    if (p_shape_name %in% colnames(metadata)) {
      p_shape_map = setNames(
        metadata[, p_shape_name],
        metadata[, arguments$`shape-by`]
        )
    } else {
      # set to circle - default when no shape column is specified
      p_shape_map = setNames(
        sample(c(15:25), length(unique(metadata[, arguments$`shape-by`]))),
        unique(metadata[, arguments$`shape-by`])
        )
    }
    # fill
    p_col_name = paste0(arguments$`shape-by`, "_col")
    if (p_col_name %in% colnames(metadata)) {
      p_color_map = setNames(
        metadata[, p_col_name],
        metadata[, arguments$`shape-by`]
        )
    } else {
      #  set to black - default when no _col column is exists
      p_color_map = setNames(
        rep(c("black"), nrow(metadata)),
        unique(metadata[, arguments$`shape-by`])
        )
    }
  } else {
    logger(paste0("The column ", arguments["shape-by"] ," specified by `--shape-by` does not exist in metafile"), "critical")
    quit(status = 1)
  }
} else {
  # here we set ucolor_map to black/circle if nothing is specified
  p_shape_map = c(
    rep(c(21), nrow(metadata)),
    nrow(metadata)
  )
  p_color_map = c(
    rep(c("black"), nrow(metadata)),
    nrow(metadata)
  )
}

# clades file
if (!is_var_empty(arguments$clade_var)) {
  clades = reader(arguments$clade_var, "tsv")
}

check_taxa_names(tree, metadata[, 1])
output_size = get_page_size(arguments$`paper-size`)

# set the taxa label size
if (!is_var_empty(arguments$`taxa-font-size`)) {
  if (arguments$`taxa-font-size` == 0) {
    size_taxa_labal = 0
  } else {
    size_taxa_labal = arguments$`taxa-font-size`
  }
} else {
  size_taxa_labal = calc_taxa_label_size(phylo = tree, page_size = output_size)
}

# set the tippoint shape size
if (!is_var_empty(arguments$`tippoint-shape-size`)) {
  if (arguments$`tippoint-shape-size` == 0) {
    tippoint_size = 0
  } else {
    tippoint_size = arguments$`tippoint-shape-size`
  }
} else {
  tippoint_size = calc_shape_size(phylo = tree, page_size = output_size)
}

taxa_offset = calc_offset(phylo = tree, page_size = output_size)

logger(paste0("Using font size: ", round(size_taxa_labal, 3)))
logger(paste0("Using shape size: ", round(tippoint_size, 3)))

##############################################
################ Tree plotting ###############
##############################################

tplot = master_builder(
  phylo = tree,
  meta = metadata,
  args = arguments,
  t_color = t_color_map,
  p_shape = p_shape_map,
  p_color = p_color_map,
  taxa_offset = taxa_offset,
  taxa_size = size_taxa_labal,
  shape_size = tippoint_size,
  tree_format = tree_format
)

logger(paste0("Saving plot to ", arguments$output), "info")

ggsave(
  filename = basename(arguments$output),
  path = dirname(arguments$output),
  plot = tplot,
  device = set_device(arguments$output),
  width = output_size[1],
  height = output_size[2],
  units = "mm"
)
logger("---- Processing complete. I hope it was a pleasent experience -----", "success", TRUE)

# output to svg. this is done to perserve text objects that seem to fail in inkscape
# current issue is that the font is not registered with svglite. this needs to be done with
# register_font(). see https://www.tidyverse.org/blog/2021/02/svglite-2-0-0/
# svglite::svglite(gsub(".pdf", ".svg", arguments$output),
#        width = 8.3,
#        height = 11.7)
# treeplot
# dev.off()
