
#' Asymptotic Score Continuity-corrected Confidence Interval for Independent RR
#'
#' @param x A four-element vector `c(x11, x12, x21, x22)` represents the numbers 
#'          in a 2-by-2 contingency table. See details.
#' @param m An number indicates the strength of the continuity correction. The 
#'          default is `NULL`, which means no continuity correction 
#'          \insertCite{miettinen1985comparative}{ASCCIRR}. We recommend
#'           `2` if the continuity correction is needed.
#' @param level Confidence level. Default is `0.95`.
#' @param tol The desired accuracy (convergence tolerance) for the iterative root 
#'            finding procedure when finding intevals. Default is `1e-10`.
#'
#' @return A list of two objects: `result` and `method`.
#' 
#' @details
#' The contingency table is considered as:
#' |           | Event     | Non-Event  |
#' | ----------|:---------:| ----------:|
#' | Group 1   | x11       | x12        |
#' | Group 2   | x21       | x22        |
#' 
#' @references{
#'   \insertRef{miettinen1985comparative}{ASCCIRR}
#' }
#' 
#' @importFrom Rdpack reprompt
#' @importFrom stats qnorm uniroot
#'
#' @examples
#' x <- c(7, 27, 1, 33)
#' 
#' ASCIndep(x)
#' ASCIndep(x, m=2)   # ASC
#' 
#' @export
ASCIndep <- function(x, m=NULL, level=0.95, tol=1e-10) {
    
    x11 <- x[1]; x12 <- x[2]; x21 <- x[3]; x22 <- x[4]
    x1p <- x11 + x12
    x2p <- x21 + x22
    N <- sum(x)
    
    r1 <- x11 / x1p
    r2 <- x21 / x2p
    
    estimate <- r1 / r2
    
    alpha <- 1 - level
    z <- qnorm(1 - alpha/2)
    
    if (is.null(m)) {
        
        scoreform <- function(RR, lo=TRUE) {
            A <- N * RR
            B <- - (x1p * RR + x11 + x2p + x21 * RR)
            C <- x11 + x21
            R2 <- (- B - sqrt(B^2 - 4*A*C)) / (2 * A)
            R1 <- R2 * RR
            V <- (R1*(1-R1)/x1p + RR^2*R2*(1-R2)/x2p)*N/(N-1)*10000
            if (lo) {
                (r1 - r2*RR)*100 / sqrt(V) - z
            } else {
                (r2*RR - r1)*100 / sqrt(V) - z
            }
        }
        
        method <- data.frame(Methods = "Score", level = level)
        
    } else {
        
        scoreform <- function(RR, lo=TRUE) {
            A <- N * RR
            B <- - (x1p * RR + x11 + x2p + x21 * RR)
            C <- x11 + x21
            R2 <- (- B - sqrt(B^2 - 4*A*C)) / (2 * A)
            R1 <- R2 * RR
            V <- (R1*(1-R1)/x1p + RR^2*R2*(1-R2)/x2p)*N/(N-1)
            if (lo) {
                (r1 - r2*RR - r2/(m*N)) / sqrt(V) - z
            } else {
                (r2*RR - r1 - r2/(m*N)) / sqrt(V) - z
            }
        }
        
        method <- data.frame(Methods = "Score with CC", CC = m, level = level)
        
    }
    
    
    if (x11 == 0) {
        lower <- 0
    } else {
        lower <- uniroot(scoreform, c(1e-10, 1e+10), tol=tol)$root
    }
    
    if (x21 == 0) {
        upper <- Inf
    } else {
        upper <- uniroot(scoreform, c(1e-10, 1e+10), lo=FALSE, tol=tol)$root
    }
    
    result <- data.frame(Lower = lower, Estimate = estimate, Upper = upper)
    
    return(list(result=result, method=method))
    
}