#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


#' Bind list of data frames together
#'
#' @param x List of data.frames to bind together
#' @examples {
#'
#' df1 <- data.frame("a" = 1, "b" = 2)
#' df2 <- data.frame("a" = 2, "b" = 4)
#' l <- list(df1, df2)
#' df3 <- rbind_list_base(l)
#' df3
#'
#' }
rbind_list_base <- function(x) {
  x2 <- do.call(
    rbind.data.frame,
    c(x, stringsAsFactors = FALSE, make.row.names = FALSE)
  )
  rownames(x2) <- seq_len(dim(x2)[1])
  x2
}
