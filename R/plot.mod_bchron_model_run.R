#' This function generates  age model plots given the output of an ageModel run.
#'
#' @param model Output of the \code{ageModel} function
#' @param type c('pdf', rectangle', 'posterior')
#'
#' @import "ggplot2"
#' @import "dplyr"
#' @import "tidyr"
#' @import "purrr"
#' @import "ggridges"

#' @export

plot.mod_bchron_model_run <- function(model, type = c('pdf',
                                                          'rectangle',
                                                          'posterior')) {

  # first format the model ----------------------------------------------------
  model$HDI <- model$HDI |>
    t() |>
    as.data.frame() |>
    add_column(position = model$predictPositions)

  model$likelihoods <- model$likelihoods |>
    as.data.frame() |>
    set_names(nm = model$ids) |>
    add_column(age = model$ageGrid) |>
    pivot_longer(cols = model$ids,
                 names_to = 'id',
                 values_to = 'probability')
  # add the positions
  positions <- data.frame(id = model$ids,
                          position = model$masterPositions)

  model$likelihoods <- model$likelihoods |>
    full_join(positions, by = 'id') |>
    drop_na() |>
    mutate(id = factor(id, levels = model$ids[order(model$masterPositions)]))

  # plot based on type --------------------------------------------------------
  p <- switch(type,
              pdf = plot_pdf_model(model),
              rectangle = plot_rectangle_model(model),
              posterior = 'posterior')
  return(p)
}


###############################################################################
# first format the model output to be more useful for ggplot

###############################################################################
# basic PDF model plot
plot_pdf_model <- function(model) {
  # format the HDI for plotting -----------------------------
  p <- model$HDI |>
    ggplot(mapping = aes(y = position)) +
    geom_ribbon(mapping = aes(xmin = `2.5%`,
                              xmax = `97.5%`),
                alpha = 0.5) +
    geom_line(mapping = aes(x = `50%`)) +
    xlab('Age') +
    ylab('position') +
    geom_density_ridges2(data = model$likelihoods,
                         mapping = aes(x = age,
                                       y = position,
                                       height = probability,
                                       group = id,
                                       fill = id),
                         stat = 'identity',
                         color = NA,
                         scale = 1) +
    theme(legend.position = 'top')
  return(p)
}


###############################################################################
# time scale style rectangle plots

plot_rectangle_model <- function(model) {

  # get the 95% range of the likelihoods
  ranges <- model$likelihoods |>
    group_by(id) |>
    arrange(age) |>
    mutate(cdf = cumsum(probability)) |>
    filter(between(cdf, 0.025, 0.975)) |>
    ungroup() |>
    group_by(id) |>
    summarise(min = min(age),
              max = max(age)) |>
    ungroup()
  # get thicknesses and positions  --------------------------
  positions <- data.frame(id = model$ids,
                          position = model$masterPositions)

  rectangles <- data.frame(id = rep(model$ids,
                                    times = model$nAges),
                           thickness = model$positionThicknesses) |>
    distinct() |>
    full_join(positions, by = 'id') |>
    mutate(top = position + thickness,
           bottom = position - thickness) |>
    full_join(ranges, by = 'id') |>
    mutate(id = factor(id, levels = model$ids[order(model$masterPositions)]))

  p <- model$HDI |>
    ggplot(mapping = aes(y = position)) +
    geom_ribbon(mapping = aes(xmin = `2.5%`,
                              xmax = `97.5%`),
                alpha = 0.5) +
    geom_line(mapping = aes(x = `50%`)) +
    xlab('Age') +
    ylab('position') +
    geom_rect(data = rectangles,
              mapping = aes(xmin = min,
                            xmax = max,
                            ymin = bottom,
                            ymax = top,
                            group = id,
                            color = id),
              fill = NA) +
    theme(legend.position = 'top')
  p <- p +
    geom_text(data = rectangles,
              mapping = aes(x = min,
                            y = position,
                            label = id,
                            color = id),
              hjust = 1.1) +
    theme(legend.position = 'none')
  return(p)
}

