#' Odds Ratio Plot Helper
#'
#' @param mod Linear model output
#' @param x_pos X axis position of ORs column
#' @param x_pos_sh X axis position start and end of underline in OR columns
#' @param breaks X axis breaks
#' @param nms Named vector for term labels in model
#' @param order Vector for order of terms
or_plot <- function(mod,
                    x_pos = 2,
                    x_pos_sh = rep(0.5,2),
                    breaks = c(0.5, 1, 1.5),
                    nms = c(),
                    order = NULL){

  # Tidy up mixed model data to be ORs from binom model
  td <- broom.mixed::tidy(mod,conf.int=TRUE,) %>%
    mutate_at(vars(dplyr::matches("est|conf")),.funs = exp)

  # label our covariates more clearly
  td <- td %>% filter(!grepl("intercept", .data$term, ignore.case = TRUE))
  td <- td %>% mutate(term = replace(.data$term, which(.data$term %in% names(nms)), nms[match(.data$term, names(nms))]))

  # order covariates by length of label for easy plotting
  if(!is.null(order)){
    td <- td[order,]
    td <- td[rev(seq_len(nrow(td))),]
  }

  # format the labels and title
  td$label <- paste0(round(td$.data$estimate,2), " (",
                     round(td$conf.low,digits = 2), "-",
                     round(td$conf.high,digits = 2), ")")
  td <- rbind(td, td[1,])
  td[nrow(td),2:ncol(td)] <- NA
  td$term[nrow(td)] <- ""
  td$label[nrow(td)] <- "aOR (95% CI)"

  td$term <- factor(td$term, levels = td$term)
  if("effect" %in% names(td)) {
    td <- td[td$effect == "fixed",]
  }

  # create OR plot
  gg <- ggplot(td, aes(.data$estimate, .data$term, color = as.logical(.data$p.value<0.05))) +
    geom_point() +
    geom_errorbarh(aes(xmin = .data$conf.low, xmax = .data$conf.high, height=0.25)) +
    geom_vline(xintercept = 1) + theme_bw() +
    scale_color_manual(name="p-value < 0.05", values = c("FALSE"=rb[2],"TRUE"=rb[1])) +
    theme(text = element_text(size = 14), panel.spacing = unit(10,units = "pt"),) +
    xlab("Odds Ratio") + ylab("") +
    geom_text(inherit.aes = FALSE, data = td,
              mapping = aes(y = .data$term, x = .data$x_pos, label = .data$label)) +
    theme(legend.position = "none",
          panel.grid.major.y =  element_blank(),
          panel.grid.minor.y =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(hjust = 0.4)) +
    scale_x_continuous(breaks = breaks, limits = c(0,x_pos + 0.5))

  # Add OR desc to right
  gg + geom_curve(data = data.frame(x = x_pos - x_pos_sh[1], y = nrow(td)-0.25,
                                    xend = x_pos + x_pos_sh[2], yend = nrow(td)-0.25),
                  mapping = aes(x = .data$x,y = .data$y, xend = .data$xend, yend = .data$yend),
                  angle = 0L, arrow = arrow(30L, unit(0L, "inches"), "last", "closed"),
                  alpha = 1, inherit.aes = FALSE)

}


rb <- c("#a63924", "#2a4e91")
