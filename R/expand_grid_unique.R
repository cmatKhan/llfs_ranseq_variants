#' expand_grid_unique
#'
#' @description Expand grid unique
#'
#' @param x vector x.
#' @param y vector y
#' @param equals If TRUE also same combinations are included
#'
#' @details
#' The function returns all possible combinations of two vectors.
#'
#' @return vector
#'
#' @examples
#' x <- 1:5
#' y <- 1:5
#' expand_grid_unique(x, y)
#'
#' @aliases expand_grid_unique
#' @rdname expand_grid_unique
#' @author https://github.com/mhesselbarth/suppoRt/blob/HEAD/R/expand_grid_unique.R
expand_grid_unique <- function(x, y, equals = FALSE)
{
  # get unique vectors for x and y
  x <- unique(x)
  y <- unique(y)

  # For an element in x, this function creates a matrix of comparisons against
  # y. If equals is true then x1 y2 and y2 x1 are both returned. If equals is
  # false, only one of those two pairs is returned.
  g <- function(i) {
    # if the length of x is 3 and equals is FALSE, the set diffs are:
    # setdiff(y,x[1])
    # setdiff(y,x[1:2])
    # setdiff(y,x[1:3])
    # if equals is TRUE, the setdiffs are
    # setdiff(y,x[0])
    # setdiff(y,x[1])
    # setdiff(y,x[1:2])
    # which results in this:
    # > setdiff(y,x[seq_len(0)])
    # [1] "PB.6.2"    "PB.7.3"    "PB.7566.3"
    # > setdiff(y,x[seq_len(1)])
    # [1] "PB.7.3"    "PB.7566.3"
    # > setdiff(y,x[seq_len(2)])
    # [1] "PB.7566.3"
    z <- setdiff(y, x[seq_len(i - equals)])
    # if z exists, then a matrix is made where the left column is an element
    # of x and the right elements are the result of the setdiff function above
    if(length(z)){
      cbind(x[i], z, deparse.level = 0)
    }
  }
  # apply function g over every element of x
  do.call(rbind, lapply(seq_along(x), g)) %>%
    as.data.frame()
}
