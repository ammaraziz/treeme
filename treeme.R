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

option_list = list(
  make_option(
    c("-t", "--tree"),
    help = "input nhx tree",
    action = "store",
    type = "character",
    default = NA
  ),
  make_option(
    c("-o", "--output"),
    help = "output name - must end in pdf or svg",
    action = "store",
    type = "character",
    default = NA
  ),
  make_option(
    c("-c", "--clades"),
    help = "clade(str) and node number in tsv format",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("-m", "--meta"),
    help = "tsv file containing meta data",
    action = "store",
    type = "character",
    default = NA
  ),
  make_option(
    c("--colors"),
    help = "TSV file - 'category\tcolor' for coloring ",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("--shapes"),
    help = "TSV file - 'category\tcolor' for tip point shapes ",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("-l", "--colorTaxa"),
    help = "which variable (in meta file) to color taxa labels",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("-p", "--tipPoint"),
    help = "which variable (in meta file) to color/shape tip points",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("-g", "--title"),
    help = "The title of the final output",
    action = "store",
    type = "character",
    default = ""
  ),
  make_option(
    c("-s", "--paperSize"),
    help = "output size - options: A3p, A3l, A4l, [A4p] - p/l is the orientation",
    action = "store",
    type = "character",
    default = "A4p"
  )
)

parser = OptionParser(
  # epilogue = "Epilogue goes here",
  option_list = option_list,
  usage = "treeme.r -t {tree} -o {pdf} -c {file} -m {meta} -l {variable} -p {variable} -g {title} -s {size}"
)

###############################
########## Functions ##########
###############################

logger = function(text, level) {
  time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  off = "\033[0m"
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
  cat(levels[level], time, " | ", toupper(level), " | ", text, off)
}

reader = function(infile, type) {
  # safely read in files
  tryCatch(
    expr = {
      if (type == "tree") {
        t = read.beast(infile)
        logger("Tree file successfully read in.", "info")
        return(t)
      }
      if (type == "tsv") {
        m = read.delim(
          infile,
          sep = "\t",
          na.strings = "",
          stringsAsFactors = F
        )
        logger(paste0(infile, " successfully read in."), "info")
        return(m)
      }
    },
    error = function(e) {
      logger(
        text = paste0("Unable to read infile: ", infile),
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
      warning = "warning"
    )
    quit()
  }
}

set_device = function() {
  ext = strsplit(arguments$outfile, ".", fixed = T)[[1]][-1]
  if (ext == "pdf") {
    device = "cairo_pdf"
  }
  if (ext == "svg") {
    device = "svg"
  }
  return(device)
}

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
      cat("Missing arguments:", paste(names(which(is.na(arguments)))))
      quit()
    }
  }
)

###############################################
########### Functions for Tree ################
###############################################

tree_title = function(title) {
  today = format(Sys.time(), format = "%d %b %Y")
  out = paste0(as.character(title), "\n", today)
}

calc_text_size = function(height, lines, buffer) {
  max_line_height = (height) / lines # mm
  return(max_line_height)
}

check_taxa_names = function(tree_labs, meta_desig) {
  # tree_labs  :  tree@phylo$tip.label
  # meta_desig :  meta[, 'designation']
  mismatch = tree_labs[!tree_labs %in% meta_desig]
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

get_paper_size = function(size_argument) {
  paperSizes = list(
    "A2p" = c(297 * 2, 420 * 2),
    "A2l" = c(420 * 2, 297 * 2),
    "A3p" = c(297, 420),
    "A3l" = c(420, 297),
    "A4p" = c(210, 297),
    "A4l" = c(297, 210)
  )
  if (!any(size_argument %in% names(paperSizes))) {
    stop(red("Paper size unknown. Options:  A3p A3l A4p A4l"))
  }
  return(paperSizes[[size_argument]])
}

#########################################
############# Input checks ##############
#########################################

tree = reader(arguments$tree, "tree")
meta = reader(arguments$meta, "tsv")
clades = reader(arguments$clades, "tsv")
cols = reader(arguments$colors)
shapes = reader(arguments$shapes)
names(shapes$shapes_type) = shapes$shape_cats

check_var_in_meta(meta, arguments$colorTaxa)
check_var_in_meta(meta, arguments$tipPoint)
check_taxa_names(tree@phylo$tip.label, meta[, "strain"])
logger("~~~All Checks Okay - Plotting tree~~~", "info")

##############################################
################ Tree plotting ###############
##############################################

output_size = get_paper_size(arguments$paperSize)

tipLabSize = calc_text_size(
  height = output_size[[1]],
  lines = tree@phylo$Nnode,
  buffer = 15
)

treeplot = ggtree(tree) %<+%
  meta +
  geom_tiplab(
    geom = "text",
    size = tipLabSize,
    aes(color = !!sym(arguments$colorTaxa)),
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
      fill = !!sym(arguments$tipPoint),
      shape = !!sym(arguments$tipPoint)
    )
  ) +

  scale_color_manual(
    arguments$colorTaxa,
    limits = cols_file$category,
    values = cols_file$color,
    na.value = "#000000"
  ) +

  scale_fill_manual(
    arguments$tipPoint,
    values = shapes_file$shape_colors,
    limits = shapes_file$shape_cats,
    na.value = "#000000"
  ) +

  scale_shape_manual(
    "Legend",
    values = shapes_file$shapes_type,
    breaks = shapes_file$shapes_type
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
  device = set_device(),
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
