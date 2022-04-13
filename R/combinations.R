#' Generate levels of all combinations
#'
#' @param nms Character vector to generate combinations
#' @param collapse Character string to separate the combinations
#' @importFrom utils head
#' @importFrom grDevices rgb
#' @importFrom dplyr select filter mutate_at mutate
#' @importFrom purrr set_names
#' @import ggplot2
gen_combinations <- function(nms, collapse = "/") {

  # first lets create all the possible combinations
  levels <- parse(text = paste0("expand.grid(",
                                paste0(rep(
                                  "c(TRUE,FALSE)", length(nms)
                                ), collapse = ","),
                                ")"))
  levels <- head(eval(levels), 2 ^ length(nms) - 1)

  # tidy the names of the levels
  levels <-
    (apply(levels, 1, function(x) {
      paste0(nms[x], collapse = collapse)
    }))
  levels <- levels[order(nchar(levels), levels)]

  # return this as an empty numeric vector with these names
  return(as.character(levels))


}

gen_color_combinations <- function(nms, colors, collapse = "/") {

  # first lets create all the possible combinations
  levels <- parse(text = paste0("expand.grid(",
                                paste0(rep(
                                  "c(TRUE,FALSE)", length(nms)
                                ), collapse = ","),
                                ")"))

  levels <- head(eval(levels), 2 ^ length(nms) - 1)

  # from these generate the color mixes
  color_levels <- apply(levels, MARGIN = 1, FUN = function(x) {
    hex_mix(colors[x])
  })

  # now bring back to the levels of the nms
  levels <-
    (apply(levels, 1, function(x) {
      paste0(nms[x], collapse = collapse)
    }))

  color_levels <- color_levels[order(nchar(levels), levels)]
  levels <- levels[order(nchar(levels), levels)]

  return(set_names(color_levels, levels))

}

rgb2hex <- function(r,g,b, maxColorValue = 255) {
  rgb(r, g, b, maxColorValue = maxColorValue)
}

mix_two_hex <- function(col1, col2, alpha = 0.5) {

  rgb1 <- colorspace::hex2RGB(col1)
  rgb2 <- colorspace::hex2RGB(col2)
  rgmix <- colorspace::mixcolor(alpha, rgb1, rgb2)

  rgb2hex(rgmix@coords[1], rgmix@coords[2], rgmix@coords[3], 1)

}

hex_mix <- function(colors) {

  if(length(colors) < 2) {

    mix <- colors

  } else {

    mix <- colors[1]
    i <- 2

    while(i <= length(colors)) {

      mix <- mix_two_hex(mix, colors[i], 1/i)
      i <- i + 1

    }
  }

  return(mix)

}
