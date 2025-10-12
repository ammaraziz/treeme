#!/usr/bin/env Rscript
options(warn = -1)

pacman::p_load(
  optparse,
  ggtree,
  ggplot2,
  treeio,
  cowplot,
  dplyr,
  Cario
)

###############################
########## Functions ##########
###############################

logger = function(text, level = "info", simple = FALSE) {
  time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
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
    warning = codes[["yellow"]],
    critical = codes[["red"]]
  )

  if (simple) {
    cat(levels[level], text, off, file = stderr())
  } else {
    cat(levels[level], time, "|", toupper(level), "|", text, off, file = stderr())
  }
}

read_tree_auto = function(infile) {
  ext = strsplit(infile, ".", fixed = T)[[1]][-1]
  # newick
  if (length(ext) > 1) {
    ext = tail(ext, 1)
  }
  if (ext == "nwk" | ext == "newick") {
    return(read.newick(infile))
  }
  if (ext == "nexus") {
    return(read.nexus(infile))
  }
  if (ext == "nhx" | ext == "beast") {
    return(read.beast(infile))
  }
  if (ext == "json") {
    logger("Reading in nextstrain json, this can take a while!", "warning")
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
          "Unable to read infile - unknown file type: ",
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
      paste0("Page size of: ", size, " is not accepted. Options: A3p A3l A4p A4l. See help for further info"),
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
########### Functions for Tree ################
###############################################

calc_text_size = function(phylo, page_size, type = "logistic") {
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
  logger(paste0("Number of taxa on tree: ", ntaxa, ", using font size: ", font_size))
  return(font_size)
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

check_empty = function(var) {
  if (is.null(var)) {
    return(FALSE)
  }
  if (is.character(var) && var == '') {
    return(FALSE)
  }
  if (length(var) > 0 & !is.null(var)) {
    return(TRUE)
  } else {
    logger(
      "What the hell did you pass to? This is a bug, go to github and submit an issue",
      "critical"
    )
  }
}

print_inputs = function(arguments) {
  logger("CLI Inputs:", simple = TRUE)
  arg_names = names(arguments)
  for (n in seq_along(arg_names)) {
    if (check_empty(arguments[n]) & !(arg_names[n] == "help")) {
      message = paste0("  ", arg_names[n], ": ", arguments[n])
      logger(message, "info", TRUE)
    }
  }
}

builder_tiplab = function(tplot, color_var, size) {
  logger("Adding geom_tiplab", "info")

  offset = 0.001
  if (check_empty(color_var)) {
    tplot = tplot +
      geom_tiplab(
        aes(color = !!sym(color_var)),
        size = size,
        offset = offset,
        family = "Arial",
        key_glyph = rectangle_key_glyph(
          fill = color,
          padding = margin(0, 0, 0, 0),
          color = "black",
          linetype = 3
        )
      ) +
      scale_color_manual(
        values = ucolors_maps,
        na.value = "#000000"
      )
  } else {
    tplot = tplot +
      geom_tiplab(
        geom = "text",
        size = size,
        offset = offset,
        family = "Arial",
        key_glyph = rectangle_key_glyph(
          padding = margin(0, 0, 0, 0),
          color = "black",
          linetype = 3
        )
      )
  }
  return(tplot)
}

builder_tippoint = function(tplot, tip_var, ushape_df) {
  logger("Adding geom_tippoint", "info")

  if (check_empty(tip_var)) {
    tplot = tplot +
      geom_tippoint(
        size = text_size,
        aes(fill = !!sym(tip_var), shape = !!sym(tip_var)),
      ) +

      scale_fill_manual(
        tip_var,
        values = ushape_col_maps,
        na.value = "#000000"
      ) +

      scale_shape_manual(
        "Legend",
        values = ushape_maps,
      ) +

      guides(
        fill = guide_legend(
          override.aes = list(
            size = text_size * 1.5,
            label = "",
            shape = ushape_maps
          )
        )
      )
    return(tplot)
  } else {
    return(tplot + geom_tippoint(size = text_size))
  }
}

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
    help = "Column name in metafile to control the color of taxa labels. If not provided, tiplabs are black.",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-s", "--shape-by"),
    help = "Column name in metafile to control the color AND shape of tip points. If not provided, tips are blank (no shape).",
    action = "store",
    type = "character"
  ),
  make_option(
    c("--clades-file"),
    help = "TSV file - 'clade\\tnode number' - used for vertical bar lines beside tree to identify clades.",
    action = "store",
    type = "character"
  ),
  make_option(
    c("-i", "--title"),
    help = "Optional: The title of the final output",
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
    c("-f", "--font-size"),
    help = "Optional: Specify the taxa labels font size. Set to 0 to turn off taxa labels. Treeme will try to auto calculate the best font size for you.",
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
    "\t- Text is squished together; Increase output side from A4 to A3.",
    "\t- A warning appears about taxa/metafile labels; open the tree file in a text editor and check the label names. These much match exactly, no spaces, no underscores."
  ),
  option_list = option_list,
  usage = "probably out of date: treeme.r -t {tree} -o {pdf} -c {file} -m {meta} -l {variable} -p {variable} -g {title} -s {size}"
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

#########################################
############# Input checks ##############
#########################################
print_inputs((arguments))

tree = reader(arguments$tree, "tree")
metadata = reader(arguments$meta, "tsv")

# colors
has_cols_file <- check_empty(arguments$`colors-file`)
has_cols_var <- check_empty(arguments$`color-by-var`)

if (xor(has_cols_file, has_cols_var)) {
  logger("Both --colors-file and --colors-by-var needed.", "critical")
  quit(status = 1)
}
if (has_cols_file && has_cols_var) {
  check_var_in_meta(metadata, arguments$`color-by-var`)
  ucolors = reader(arguments$`colors-file`, "tsv")
  ucolors_maps = setNames(ucolors$color, ucolors$category)
}

# shapes
has_shapes_file <- check_empty(arguments$`shapes-file`)
has_shapes_var <- check_empty(arguments$`shape-by-var`)

if (xor(has_shapes_file, has_shapes_var)) {
  logger("Both --shapes-file and --shape-by-var needed.", "critical")
  quit(status = 1)
}

if (has_shapes_file && has_shapes_var) {
  check_var_in_meta(metadata, arguments$`shape-by-var`)
  ushapes = reader(arguments$`shapes-file`, "tsv")
  ushape_maps = setNames(ushapes$shape, ushapes$category)
  ushape_col_maps = setNames(ushapes$color, ushapes$category)
}

# clades file
if (check_empty(arguments$clade_var)) {
  clades = reader(arguments$clade_var, "tsv")
}

check_taxa_names(tree, metadata[, 1])
output_size = get_page_size(arguments$`paper-size`)

# set the tip lab size
if (check_empty(arguments$`font-size`)) {
  if (arguments$`font-size` == 0) {
    taxa_text_size = 0
  } else {
    taxa_text_size = arguments$`font-size`
  }
} else {
  taxa_text_size = calc_text_size(phylo = tree, page_size = output_size)
}

text_size = calc_text_size(phylo = tree, page_size = output_size)

logger("----- All Checks Okay - Plotting tree -----", "info", TRUE)

##############################################
################ Tree plotting ###############
##############################################

tplot = ggtree(tree) %<+% metadata

tplot = builder_tiplab(
  tplot = tplot,
  color_var = arguments$`color-by-var`,
  size = taxa_text_size
)

tplot = builder_tippoint(
  tplot = tplot,
  tip_var = arguments$`shape-by-var`,
  ushape_df = ushapes
)

# guides(
#   color = guide_legend(
#     override.aes = list(
#       size = text_size * 5,
#       label = "\u25A0",
#       linetype = 3
#     )
#   )
# ) +

# ggtitle(tree_title(arguments$title))

tplot = tplot +
  theme(
    legend.position = c(0.1, 0.65),
    legend.key.size = unit(text_size * 2, "mm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    plot.title = element_text(hjust = 0.06, vjust = -15, size = 20),
    plot.subtitle = element_text(hjust = 0.02, vjust = -12, size = 20)
  )

# # mutations on branches
# nudge = ggplot_build(tplot)$layout$panel_scales_y[[1]]$range$range[1] / 3

# tplot = tplot +
#   geom_text(
#     aes(x = branch, label = aa_muts),
#     size = text_size - 0.5,
#     nudge_y = nudge
#   )

# Fix tip label clipping
# plot_dim_x = ggplot_build(tplot)$layout$panel_scales_x[[1]]$range$range[2]

# tplot = tplot +
#   coord_cartesian(clip = "off", expand = FALSE) +
#   xlim(NA, ((0.40 * plot_dim_x) + plot_dim_x))

# # add clades
# if (check_empty(arguments$`clades-file`)) {
#   treeplot = treeplot + add_clades(cladesFile, tree@data, plot_dim_x)
# }

ggsave(
  filename = basename(arguments$output),
  path = dirname(arguments$output),
  plot = tplot,
  device = set_device(arguments$output),
  width = output_size[1],
  height = output_size[2],
  units = "mm"
)

# output to svg. this is done to perserve text objects that seem to fail in inkscape
# current issue is that the font is not registered with svglite. this needs to be done with
# register_font(). see https://www.tidyverse.org/blog/2021/02/svglite-2-0-0/
# svglite::svglite(gsub(".pdf", ".svg", arguments$output),
#        width = 8.3,
#        height = 11.7)
# treeplot
# dev.off()
