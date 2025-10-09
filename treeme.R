#!/usr/bin/env Rscript
options(warn = -1)

pacman::p_load(
  optparse,
  ggtree,
  ggplot2,
  treeio,
  cowplot,
  dplyr,
  Cairo
)

###############################
########## Functions ##########
###############################

logger = function(text, level, simple = FALSE) {
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
  if (!simple) {
    cat(levels[level], time, "|", toupper(level), "|", text, off)
  } else {
    cat(levels[level], toupper(level), "|", text, off)
  }
}

read_tree_auto = function(infile) {
  ext = strsplit(infile, ".", fixed = T)[[1]][-1]
  # newick
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
        logger(
          paste0(basename(infile), " successfully read in."),
          "info"
        )
        return(t)
      }
      if (type == "tsv") {
        m = read.delim(
          infile,
          sep = "\t",
          na.strings = "",
          stringsAsFactors = F
        )
        logger(paste0(basename(infile), " successfully read in."), "info")
        return(m)
      }
    },
    error = function(e) {
      logger(
        text = paste0(
          "Unable to read infile - unknown file type: ",
          infile,
          "\n Printing error: ",
          e
        ),
        level = "critical"
      )
      quit()
    }
  )
}

check_var_in_meta = function(metafile, variable) {
  if (!(variable %in% names(metafile))) {
    logger(
      text = paste0(
        "Color column'",
        variable,
        "' does not exist in metafile. Check inputs."
      ),
      level = "warning"
    )
    quit()
  }
}

set_device = function(outfile) {
  ext = strsplit(outfile, ".", fixed = T)[[1]][-1]
  if (ext == "pdf") {
    device = "cairo_pdf"
  }
  if (ext == "svg") {
    device = "svg"
  }
  return(device)
}

get_paper_size = function(size) {
  paper = list(
    "A2p" = c(297 * 2, 420 * 2),
    "A2l" = c(420 * 2, 297 * 2),
    "A3p" = c(297, 420),
    "A3l" = c(420, 297),
    "A4p" = c(210, 297),
    "A4l" = c(297, 210)
  )
  if (!any(size %in% names(paper))) {
    logger(
      paste0("Page size of: ", size, " is not accepted. Options: A3p A3l A4p A4l. See help for further info"),
      "critical"
    )
    quit()
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

calc_text_size = function(height, phylo, buffer = 5) {
  if (class(phylo) == "treedata") {
    taxa = phylo@phylo$Nnode
  }
  if (class(phylo) == "phylo") {
    taxa = phylo$Nnode
  } else {
    logger("Unable to get num of tips in input tree. Check tree type is nwk or nhx", "critical")
    quit()
  }

  max_line_height = (height + buffer) / taxa # mm
  return(max_line_height)
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
    quit()
  }

  mismatch = taxa_labels[!taxa_labels %in% meta_desig]

  if (!length(mismatch)) {
    logger(text = "All tip labels present in metafile", level = "info")
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
  # returns FALSE when empty
  if (is.null(var)) {
    return(FALSE)
  }
  if (var == '') {
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
  logger("CLI Inputs:", "info", TRUE)
  arg_names = names(arguments)
  for (n in seq_along(arg_names)) {
    if (check_empty(arguments[n]) & !(arg_names[n] == "help")) {
      message = paste0("  ", arg_names[n], ": ", arguments[n])
      logger(message, "info", TRUE)
    }
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
    c("-m", "--meta"),
    help = "Required: tsv file containing metafile.",
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
    c("--clades-file"),
    help = "TSV - 'clade\tnode number' - used for vertical lines on a tree.",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("--colors-file"),
    help = "TSV file - 'category\tcolor' for coloring taxa names.",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("--color-var"),
    help = "Variable in metafile to control the color of taxa labels. Ensure all categories are in taxa/metafile.",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("--shapes-file"),
    help = "TSV file - See below for information",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("--shape-var"),
    help = "Variable in file to control the tip points shape and color. Ensure all categories are in taxa/metafile.",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("--title"),
    help = "Optional: The title of the final output",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("--paper-size"),
    help = "Optional: Output size - options: A3p, A3l, A4l, [A4p] - p/l is the orientation. [default %default]",
    action = "store",
    type = "character",
    default = "A4p"
  )
)

parser = OptionParser(
  epilogue = c(
    "Metafile must be a tab separated file, with the first column containing the sample id matching taxa labels on tree.",
    "All subsequent columns are optional, any extra columns required must match --colors-file and --tip-point inputs.",
    "--shapes-file enables control of the tippoint shapes and colors.",
    "\n",
    "Common Issues and solutions:",
    "\t- Text is squished together; Increase output side from A4 to A3."
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
      quit()
    }
  }
)

#########################################
############# Input checks ##############
#########################################
print_inputs((arguments))
tree = reader(arguments$tree, "tree")
metadata = reader(arguments$meta, "tsv")

if (check_empty(arguments$clades)) {
  clades = reader(arguments$clades, "tsv")
}
if (check_empty(arguments$`colors-file`)) {
  ucolors = reader(arguments$`colors-file`)
}
if (check_empty(arguments$`shapes-file`)) {
  ushapes = reader(arguments$`shapes-file`)
  names(ushapes$shapes_type) = ushapes$`shape-cats`
}
if (!check_empty(arguments$colorTaxa)) {
  arguments$colorTaxa = "black"
}
if (check_empty(arguments$`shape-var`)) {
  check_var_in_meta(metadata, arguments$tipPoint)
}

check_taxa_names(tree, metadata[, 1])
output_size = get_paper_size(arguments$`paper-size`)
tipLabSize = calc_text_size(height = output_size[[1]], phylo = tree)

logger("----- All Checks Okay - Plotting tree -----", "info")

##############################################
################ Tree plotting ###############
##############################################

plot = ggtree(tree) %<+%
  metadata +
  geom_tiplab(
    aes(color = !!sym(arguments$colorTaxa)),
    geom = "text",
    size = tipLabSize,
    key_glyph = rectangle_key_glyph(
      fill = color,
      padding = margin(0, 0, 0, 0),
      color = "black",
      linetype = 3
    ),
    offset = 0.00009,
    family = "Arial"
  ) +

  geom_tippoint(
    size = tipLabSize,
    aes(
      fill = !!sym(arguments$`shape-var`),
      shape = !!sym(arguments$`shape-var`)
    )
  ) +

  scale_color_manual(
    arguments$colorTaxa,
    limits = ucolors$category,
    values = ucolors$color,
    na.value = "#000000"
  ) +

  scale_fill_manual(
    arguments$tipPoint,
    values = ushapes$colors,
    limits = ushapes$category,
    na.value = "#000000"
  ) +

  scale_shape_manual(
    "Legend",
    values = shapes$shape,
    breaks = shapes$shape
  ) +

  guides(
    fill = guide_legend(
      override.aes = list(
        size = tipLabSize * 1.5,
        label = "",
        shape = shapes_file$shapes_type
      )
    )
  ) +

  guides(
    color = guide_legend(
      override.aes = list(
        size = tipLabSize * 5,
        label = "\u25A0",
        linetype = 3
      )
    )
  ) +

  ggtitle(tree_title(arguments$title)) +

  theme(
    legend.position = c(0.1, 0.65),
    legend.key.size = unit(tipLabSize * 2, "mm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    plot.title = element_text(hjust = 0.06, vjust = -15, size = 20),
    plot.subtitle = element_text(hjust = 0.02, vjust = -12, size = 20)
  )

# mutations on branches
nudge = ggplot_build(treeplot)$layout$panel_scales_y[[1]]$range$range[1] / 3

treeplot = treeplot +
  geom_text(
    aes(x = branch, label = aa_muts),
    size = tipLabSize - 0.5,
    nudge_y = nudge
  )

# Fix tip label clipping
plot_dim_x = ggplot_build(treeplot)$layout$panel_scales_x[[1]]$range$range[2]

treeplot = treeplot +
  coord_cartesian(clip = "off", expand = FALSE) +
  xlim(NA, ((0.40 * plot_dim_x) + plot_dim_x))

# add clades
if (arguments$clades) {
  treeplot = treeplot + add_clades(cladesFile, tree@data, plot_dim_x)
}

ggsave(
  filename = arguments$output,
  plot = treeplot,
  device = set_device(arguments$outfile),
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
