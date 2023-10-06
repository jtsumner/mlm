nature_theme <- function(x_axis_labels, y_label) {
    # set default text format based on categorical and length
    # code modified from Maaslin2 source code (Mallick et al., (2021). PLoS Computational Biology)
    angle = NULL
    hjust = NULL
    size = 8
    if (max(nchar(x_axis_labels), na.rm=TRUE) > 5) {
        angle = 45
        hjust = 1
        size = 6
    }
    axis_title_size = 10
    if (nchar(y_label) > 15) {
        axis_title_size = 8
    }
    if (nchar(y_label) > 25) {
        axis_title_size = 6
    }
    return ( ggplot2::theme_bw() + ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = size, vjust = 1, hjust = hjust, angle = angle, colour = "black"),
        axis.text.y = ggplot2::element_text(size = 8, hjust = 1, colour = "black"),
        axis.title = ggplot2::element_text(size = axis_title_size),
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
        legend.title = ggplot2::element_text(size = 6, face = 'bold'),
        legend.text = ggplot2::element_text(size = 6),
        axis.line = ggplot2::element_line(colour = 'black', size = .5),
        axis.line.x = ggplot2::element_line(colour = 'black', size = .5),
        axis.line.y = ggplot2::element_line(colour = 'black', size = .5),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.ticks = element_line(color="black"),
        text = element_text(family = "Arial"))
   )
}