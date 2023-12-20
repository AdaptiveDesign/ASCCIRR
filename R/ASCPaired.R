
#' Asymptotic Score Continuity-corrected Confidence Interval for Paired RR
#'
#' @param x A four-element vector `c(x11, x12, x21, x22)` represents the numbers 
#'          in a 2-by-2 contingency table. See details.
#' @param m An number indicates the strength of the continuity correction. The 
#'          default is `NULL`, which means no continuity correction 
#'          \insertCite{tang2012confidence}{ASCCIRR}. And the values 
#'          of `2`, `4`, and `8` represent the continuity corrections of high, 
#'          medium, and low strength, respectively \insertCite{delrocco2023new}{ASCCIRR}.
#' @param level Confidence level. Default is `0.95`.
#' @param ctrl A very small number will be added to the contingency table if `0` exists.
#'             Default is `1e-10`.
#'
#' @return A list of two objects: `result` and `method`.
#' 
#' @details
#' The contingency table is considered as:
#' |                   | Test2: Event   | Test2: Nonevent  |
#' | ------------------|:--------------:| ----------------:|
#' | Test1: Event      | x11            | x12              |
#' | Test1: Nonevent   | x21            | x22              |
#' 
#' @references{
#'   \insertRef{tang2012confidence}{ASCCIRR}
#'   \insertRef{delrocco2023new}{ASCCIRR}
#' }
#' 
#' @importFrom Rdpack reprompt
#' @importFrom utils install.packages installed.packages
#' @importFrom polynom polynomial
#'
#' @examples
#' x <- c(9, 36, 14, 1371)
#' 
#' ASCPaired(x)
#' ASCPaired(x, m=2)  # ASCC-H
#' ASCPaired(x, m=4)  # ASCC-M
#' ASCPaired(x, m=8)  # ASCC-L
#' 
#' @export
ASCPaired <- function(x, m=NULL, level=0.95, ctrl=1e-10) {
    
    if (any(x == 0)) {
        x <- x + ctrl
    }
    x11 <- x[1]
    x12 <- x[2]
    x21 <- x[3]
    x22 <- x[4]
    
    xp1 <- x11 + x21
    x1p <- x11 + x12
    estimate <- x1p / xp1
    N <- sum(x)
    
    alpha <- 1 - level
    z <- qnorm(1 - alpha/2)
    
    if (is.null(m)) {
        
        a <- xp1^4 + z^2 * xp1^3
        b <- -(2*xp1^2 + z^2*xp1) * (2*xp1*x1p + z^2*(x11+x12+x21))
        c <- 6*xp1^2*x1p^2 + z^4*(x1p+xp1)*(x11+x12+x21) + z^2*xp1*x1p*(6*x11+5*x12+5*x21)
        d <- - (2*x1p^2 + z^2*x1p) * (2*xp1*x1p + z^2*(x11+x12+x21))
        e <- x1p^4 + z^2 * x1p^3
        
        b <- b / a
        c <- c / a
        d <- d / a
        e <- e / a
        a <- a / a
        
        theta <- solve(polynom::polynomial(c(e, d, c, b, a)))
        theta <- as.numeric(ifelse(Im(theta) == 0, Re(theta), NA))
        
        flag <- calculate_score(theta, x, level=level)
        
        if (sum(flag) == 2) {
            root <- theta[which(flag == 1)]
        } else {
            print("More than or less than 2 roots have the score close to z_0.975")
            root <- rep(NA, 2)
        }
        
        lower <- min(root)
        upper <- max(root)
        
        method <- data.frame(Methods = "Score Asymptotic CI", level = level)
        
    } else {
        
        a_lower <- a_upper <- xp1^4 + z^2 * xp1^3
        
        b_lower <- -(2*xp1^2 + z^2*xp1) * (2*xp1*(x1p-xp1/(m*N)) + z^2*(x11+x12+x21))
        b_upper <- -(2*xp1^2 + z^2*xp1) * (2*xp1*(x1p+xp1/(m*N)) + z^2*(x11+x12+x21))
        
        c_lower <- 6*xp1^2*(x1p-xp1/(m*N))^2 + z^4*(x1p+xp1)*(x11+x12+x21) + z^2*xp1*((x1p-xp1/(m*N))^2 + 4*(x11+x12+x21)*(x1p-xp1/(m*N)) + x1p*xp1)
        c_upper <- 6*xp1^2*(x1p+xp1/(m*N))^2 + z^4*(x1p+xp1)*(x11+x12+x21) + z^2*xp1*((x1p+xp1/(m*N))^2 + 4*(x11+x12+x21)*(x1p+xp1/(m*N)) + x1p*xp1)
        
        d_lower <- - (2*(x1p-xp1/(m*N))^2 + z^2*x1p) * (2*xp1*(x1p-xp1/(m*N)) + z^2*(x11+x12+x21))
        d_upper <- - (2*(x1p+xp1/(m*N))^2 + z^2*x1p) * (2*xp1*(x1p+xp1/(m*N)) + z^2*(x11+x12+x21))
        
        e_lower <- (x1p-xp1/(m*N))^4 + z^2 * x1p * (x1p-xp1/(m*N))^2
        e_upper <- (x1p+xp1/(m*N))^4 + z^2 * x1p * (x1p+xp1/(m*N))^2
        
        b_lower <- b_lower / a_lower
        b_upper <- b_upper / a_upper
        c_lower <- c_lower / a_lower
        c_upper <- c_upper / a_upper
        d_lower <- d_lower / a_lower
        d_upper <- d_upper / a_upper
        e_lower <- e_lower / a_lower
        e_upper <- e_upper / a_upper
        a_lower <- a_lower / a_lower
        a_upper <- a_upper / a_upper
        
        lowertheta <- solve(polynom::polynomial(c(e_lower, d_lower, c_lower, b_lower, a_lower)))
        lowertheta <- as.numeric(ifelse(Im(lowertheta) == 0, Re(lowertheta), NA))
        uppertheta <- solve(polynom::polynomial(c(e_upper, d_upper, c_upper, b_upper, a_upper)))
        uppertheta <- as.numeric(ifelse(Im(uppertheta) == 0, Re(uppertheta), NA))
        
        lowerflag <- calculate_score(lowertheta, x, m=m, lower=TRUE, level=level)
        upperflag <- calculate_score(uppertheta, x, m=m, lower=FALSE, level=level)
        
        if (sum(lowerflag) == 1) {
            lower <- lowertheta[which(lowerflag == 1)]
        } else {
            print("No root for lower bound has the score close to z_0.975")
            lower <- NA
        }
        
        if (sum(upperflag) == 1) {
            upper <- uppertheta[which(upperflag == 1)]
        } else {
            print("No root for upper bound has the score close to z_0.975")
            upper <- NA
        }
        
        method <- data.frame(Methods = "Continuity-corrected Score Asymptotic CI", CC = m, level = level)
        
    }
    
    result <- data.frame(Lower = lower, Estimate = estimate, Upper = upper)
    
    return(list(result=result, method=method))
    
}





calculate_score <- function(theta, x, m=NULL, lower=TRUE, level) {
    
    x11 <- x[1]
    x12 <- x[2]
    x21 <- x[3]
    x22 <- x[4]
    xp1 <- x11 + x21
    x1p <- x11 + x12
    N <- sum(x)
    
    alpha <- 1 - level
    z <- qnorm(1 - alpha/2)
    
    if (is.null(m)) {
        Stheta <- x1p - xp1 * theta
    } else if (lower) {
        Stheta <- (x1p - xp1 * theta) - xp1 / (m * N)
    } else {
        Stheta <- (x1p - xp1 * theta) + xp1 / (m * N)
    }
    
    A <- N * (1 + theta)
    B <- xp1 * theta^2 - (x1p + 2 * x21)
    C_ <- x21 * (1 - theta) * (x1p + x21)/N
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    q21 <- ifelse(num == 0, 0, num/(2 * A))
    V <- pmax(0, N * (1 + theta) * q21 + (x1p + x21) * (theta - 1))
    score <- ifelse(Stheta == 0, 0, Stheta/sqrt(V))
    
    if (is.null(m)) {
        flag <- ifelse((!is.na(score)) & (abs(abs(score) - z) <= 0.00005), 1, 0)
    } else if (lower) {
        flag <- rep(0, 4)
        if (sum((!is.na(score)) & (abs(score - z) <= 0.00005)) >= 1) {
            flag[which.min(abs(score - z))] <- 1
        }
    } else {
        flag <- rep(0, 4)
        if (sum((!is.na(score)) & (abs(score + z) <= 0.00005)) >= 1) {
            flag[which.min(abs(score + z))] <- 1
        }
    }
    
    return(flag)
}
