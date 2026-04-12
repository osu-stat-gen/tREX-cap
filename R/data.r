#' Simulated Hi-C contact matrix
#'
#' @format A \eqn{n \times n} matrix where the (i, j) element is the interaction count between bin i and bin j.
#' The matrix is symmetric with 0 on the diagonal.
"sim_hic"



#' Simulated bias matrix
#'
#' @format A matrix with n rows and 3 columns: 
#' \describe{
#'   \item{column 1}{effective fragment information.}
#'   \item{column 2}{GC content}
#'   \item{column 3}{mappability}
#' }
"sim_bias"
