# PLOT THEMES -------------------------------------------------------------


#' Generic theme for ggplot2 output.
#'
#' Uses 'theme_bw()' and puts legend at the bottom.
#'
#' @keywords theme_bw
#' @return A generic theme with legend at the bottom.
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalised counts, according to the covariate information in meta
#' \dontrun{plot <- ggplot2::boxplot(data) + theme_new}
theme_new <- ggplot2::theme_bw() +
  # adjustments for the legend
  ggplot2::theme(legend.position="bottom",
                 legend.text = ggplot2::element_text(color = "black", size=12),
                 legend.key = ggplot2::element_rect(size=12),
                 axis.text.x=ggplot2::element_text(angle=35, hjust=1))


#' Poster theme for ggplot2 output.
#'
#' Adjusts plot to use University of Luebeck (UzL) colors.
#'
#' @keywords theme poster
#' @param base_size set base font size - default is 18
#' @param base_family set base font family - default is Myriad Pro
#' @return A poster theme with legend at the bottom.
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalised counts, according to the covariate information in meta
#' \dontrun{plot <- ggplot2::boxplot(data) + theme_poster(base_size = 18,
#' base_family = "Myriad Pro")}
theme_poster <- function (base_size = 18, base_family = "Myriad Pro") {
  ## dummy ##
  main_color <- "#014B5A"

  ggplot2::theme(
    line = ggplot2::element_line(colour = "#7F7F7F", size = 0.5,
                                 linetype = 1, lineend = "butt"),
    rect = ggplot2::element_rect(fill = "transparent", colour = NA,
                                 size = 0.5, linetype = 1),
    text = ggplot2::element_text(family = base_family, face = "plain",
                                 colour = eval(main_color), size = base_size,
                                 lineheight = 0.9,  hjust = 0.5,
                                 vjust = 0.5, angle = 0,
                                 margin = margin(), debug = FALSE),

    axis.text.x=ggplot2::element_text(angle=35, hjust=1, color = eval(main_color), size=18),
    axis.text.y=ggplot2::element_text(color = eval(main_color), size=12),
    axis.ticks = ggplot2::element_line(color = "#7F7F7F"),
    axis.line = ggplot2::element_line(color = "#7F7F7F"),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = rel(1.8), angle = 90),

    ### FACET-GRID/WRAP THEME ###
    strip.text.x = ggplot2::element_text(size = eval(base_size), color = eval(main_color), face = "bold"),
    strip.text.y = ggplot2::element_text(size = base_size, color = eval(main_color), face = "bold"),
    strip.background = ggplot2::element_rect(color=eval(main_color), fill="transparent", size=1.5, linetype="solid"),

    legend.position="bottom",
    legend.text = ggplot2::element_text(color = eval(main_color), size=18),
    legend.key = ggplot2::element_rect(size=20),
    legend.justification = "center",
    legend.background = ggplot2::element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = ggplot2::element_blank(), # get rid of legend panel bg

    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA), # bg of the panel
    panel.border = ggplot2::element_blank(),
    #panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = ggplot2::element_blank(), # get rid of minor grid

    ### ALTERNATIV ###
    panel.grid.major = ggplot2::element_line(colour = "#7F7F7F", size = 0.25),
    #panel.grid.minor = element_line(colour = eval(main_color), size = 0.25),
    ### ALTERNATIV ###

    plot.background = ggplot2::element_rect(fill = "transparent", color = NA) # bg of the plot
  )
}



#' RLE theme for ggplot2 output.
#'
#' Default theme for plotting the packages RLE-plots.
#'
#' @keywords RLE theme
#' @param base_size set base font size - default is 18
#' @param base_family set base font family - default is Myriad Pro
#' @param legend_position on of "top", "bottom", "left", "right" determines legend position
#' @return A poster theme with legend at the bottom.
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalised counts, according to the covariate information in meta
#' \dontrun{plot <- mbecRLE(input.obj, model.vars=c("group","batch"),
#' return.data=FALSE) + theme_rle(base_size = 18, base_family = "Myriad Pro",
#' legend_position = 'bottom')}
theme_rle <- function(base_size = 18, base_family = "Myriad Pro",legend_position = 'bottom') {
  ggplot2::theme_bw() +
    # adjustments for the legend
    ggplot2::theme(legend.position=eval(legend_position),
          legend.text = ggplot2::element_text(color = "black", size=12),
          legend.key = ggplot2::element_rect(size=12),
          axis.text.x=ggplot2::element_text(angle=35, hjust=1),
          axis.title.x = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          #axis.text.x=element_blank(),
          axis.ticks.x=ggplot2::element_blank())

}


#' PCA theme for ggplot2 output.
#'
#' Default theme for plotting the packages PCA-plots.
#'
#' @keywords PCA theme
#' @param base_size set base font size - default is 18
#' @param base_family set base font family - default is Myriad Pro
#' @param legend_position on of "top", "bottom", "left", "right" determines legend position
#' @return A poster theme with legend at the bottom.
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalised counts, according to the covariate information in meta
#' \dontrun{plot <- mbecPCA(input.obj, model.vars=c("group","batch"), return.data=FALSE) +
#' theme_pca(base_size = 18, base_family = "Myriad Pro",legend_position = 'bottom')}
theme_pca <- function(base_size = 18, base_family = "Myriad Pro",legend_position = 'top') {
  # ToDo: get rid of this
  x.angle = 0
  x.hjust = 0.5
  # density.lwd = 0.2
  # title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex =0.75

  ggplot2::theme_bw() +
    ggplot2::theme(
      #legend.position=eval(legend_position), I am not sure what I was trying to do here..
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = 'right', legend.box = 'horizontal',
      legend.direction = 'vertical',
      legend.key.height = ggplot2::unit(0.2, 'cm'),
      legend.key.width = ggplot2::unit(0.1, 'cm'),
      legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
      legend.spacing.x = ggplot2::unit(0.1, 'cm'),
      legend.spacing.y = ggplot2::unit(0.1, 'cm'),
      legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))
    )
}


#' BOX theme for ggplot2 output.
#'
#' Default theme for plotting the packages BOX-plots.
#'
#' @keywords BOX theme
#' @param base_size set base font size - default is 18
#' @param base_family set base font family - default is Myriad Pro
#' @param legend_position on of "top", "bottom", "left", "right" determines legend position
#' @return A poster theme with legend at the bottom.
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalized counts, according to the covariate information in meta
#' \dontrun{plot <- mbecPCA(input.obj, model.vars=c("group","batch"), return.data=FALSE) +
#' theme_box(base_size = 18, base_family = "Myriad Pro",legend_position = 'bottom')}
theme_box <- function(base_size = 18, base_family = "Myriad Pro",legend_position = 'bottom') {
  # ToDo: get rid of this
  x.angle = 0
  x.hjust = 0.5
  density.lwd = 0.2
  title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex =0.75

  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = x.angle, hjust = x.hjust), panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5),

      legend.position = eval(legend_position), legend.box = 'horizontal',
      legend.direction = 'horizontal',
      legend.key.height = ggplot2::unit(0.8, 'cm'),
      legend.key.width = ggplot2::unit(0.4, 'cm'),
      legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
      legend.spacing.x = ggplot2::unit(0.4, 'cm'),
      legend.spacing.y = ggplot2::unit(0.4, 'cm'),
      legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))
    )
}


#' MOSAIC theme for ggplot2 output.
#'
#' Default theme for plotting the packages BOX-plots.
#'
#' @keywords MOSAIC theme
#' @param base_size set base font size - default is 18
#' @param base_family set base font family - default is Myriad Pro
#' @param legend_position on of "top", "bottom", "left", "right" determines legend position
#' @return A poster theme with legend at the bottom.
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalised counts, according to the covariate information in meta
#' \dontrun{plot <- mbecPCA(input.obj, model.vars=c("group","batch"), return.data=FALSE) +
#' theme_box(base_size = 18, base_family = "Myriad Pro",legend_position = 'bottom')}
theme_mosaic <- function (base_size = 18, base_family = "Myriad Pro",legend_position = 'top') {

  ## ToDo: figure out what to do with these options
  main_color = "#004B5A"
  x.angle = 0
  x.hjust = 0.5
  density.lwd = 0.2
  title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex =0.75

  ggplot2::theme(
    axis.text.x=ggplot2::element_blank(),
    axis.text.y=ggplot2::element_text(color = eval(main_color), size=12),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "#7F7F7F"),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), angle = 90),
    legend.position = eval(legend_position), legend.box = 'horizontal',
    legend.direction = 'horizontal',
    legend.key.height = ggplot2::unit(0.2, 'cm'),
    legend.key.width = ggplot2::unit(0.1, 'cm'),
    legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
    legend.spacing.x = ggplot2::unit(0.1, 'cm'),
    legend.spacing.y = ggplot2::unit(0.1, 'cm'),
    legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))
  )
}
