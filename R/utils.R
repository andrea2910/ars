#' @title check_positive
#' @description Check whether a given function is only positive
#' @param f a function to check whether it is positive on the y axis
#' @param var_min the lower bound domain to check with f
#' @param var_max the upper bound domain to check with f
#' @return a boolean TRUE or FALSE. A TRUE means the function is only positive
#' @examples
#' check_positive(dnorm, -5, 5)
#' @export
check_positive = function(f, var_min, var_max) {
  # quick checks
  if (class(f) != "function")
    stop("f is not a function")

  if (class(var_min) != "numeric" |
      length(var_min) != 1)
    stop("var_min is not numeric or not length 1")

  if (class(var_max) != "numeric" |
      length(var_max) != 1)
    stop("var_max is not numeric or not length 1")

  # choose any point in interval (var_min, var_max)
  # and check whether it is negative
  test_pts = max(var_max - 1, var_min)
  if (f(test_pts) < 0)
    return(FALSE)

  # find f function values at boundaries
  f_var_min <- f(var_min)
  f_var_max <- f(var_max)

  # if the bound values have different signs, then return false
  # because it must pass through 0
  if (sign(f_var_min) != sign(f_var_max))
    return(FALSE)

  # update var_min and var_max if var_min or var_max is infinite
  if (is.infinite(var_min))
    var_min = -10 ^ 8
  if (is.infinite(var_max))
    var_max = 10 ^ 8

  # try to find 100 root values between lower and upper
  results <-
    try(uniroot.all(Vectorize(f), lower = var_min, upper = var_max))

  if (class(results) == "try-error")
    stop("Error using uniroot.all. Try different values")

  if (length(results) == 0) {
    # no roots found
    if (f_var_min > 0)
      return(TRUE) # if f_var_min is positive, then we return TRUE
    if (f_var_min <= 0)
      return(FALSE) # if it isn't, then we return FASLE
  } else{
    if (sum(results == 0) == 0)
      return(TRUE) # sum up booleans
    else
      return(FALSE)
  }
}

#' @title calc_deriv
#' @description Calculate the derivative of function f: f'(x)
#' @param x a specific x point to find the derivative at
#' @param f a function to calculate the derivative of
#' @param lower the lower bound domain to check with f
#' @param upper the upper bound domain to check with f
#' @param ... additional variables needed in the f function
#' @return a numeric value that represents the derivative at a certain point
#' @examples
#' calc_deriv(5, dnorm, -5, 5, mean=5, sd=1)
#' @export
calc_deriv = function(x, f, lower, upper, ...) {
  if(class(f) != "function") stop("f must be a function")
  if(class(x) != "numeric") stop("x must be a number")
  if(length(x) != 1) stop("x must be only 1 number")
  if(class(lower) != "numeric") stop("lower must be a number")
  if(class(upper) != "numeric") stop("upper must be a number")

  eps <- (.Machine$double.eps) ^ (1 / 4) # find a small difference value
  d <- NA
  if (x == lower)
    d <- (f(x + eps, ...) - f(x, ...)) / eps
  else if (x == upper)
    d <-  (f(x, ...) - f (x - eps, ...)) / eps
  else if (lower <= x &&
           x <= upper)
    d <- (f(x + eps, ...) - f(x - eps, ...)) / (2 * eps)
  else
    stop("Out of bounds")

  # if limit doesn't exist then we need to stop
  if (is.na(d) |
      is.infinite(d) | is.nan(d))
    stop("The derivative does not exist.
         Please make sure your function has a derivative.")
  return(d)
}

#' @title check_interpolconcave
#' @description given a (sorted) set of points,
#' check if the interpolating function is concave
#' @param f represent f(x)
#' @param x x points
#' @return a boolean TRUE or FALSE. A TRUE means the function is concave
#' @export
#' @examples
#' log_f <- function(x) log(dnorm(x))
#' x <- seq(-3, 3, by=1)
#' y <- log_f(x)
#' check_interpolconcave(x, y)
check_interpolconcave = function(x, f) {
  if (length(x) != length(f)) {
    stop('Length of x and f should be equal.')
  }
  sig = TRUE
  for (i in 1:(length(x) - 2)) {
    inter_f = f[i] + (x[i + 1] - x[i]) * (f[i + 2] - f[i]) / (x[i + 2] - x[i])
    if (inter_f > f[i + 1]) {
      sig = FALSE
      break
    }
  }
  sig
}


#' @title create_upphull
#' @description given a (sorted) set of points, create the upper hull function
#' @param f represent f(x)
#' @param x x points
#' @return a function that calculates the upperhull for a function
#' @examples
#' log_f <- function(x) log(dnorm(x))
#' x <- seq(-3, 3, by=1)
#' y <- log_f(x)
#' u <- create_upphull(x, y)
#' u(0)
#' @export
create_upphull = function(x, f) {
  if (length(x) != length(f)) {
    stop('Length of x and f should be equal.')
  }
  k = length(x)
  u = function(p) {
    val = NA
    j = findInterval(p, x)

    if (j == 0) {
      val = f[1] + (p - x[1]) * (f[2] - f[1]) / (x[2] - x[1])
    } else if (j == k) {
      val = f[k] + (p - x[k - 1]) * (f[k] - f[k - 1]) / (x[k] - x[k - 1])
    } else {
      val = f[j] + (p - x[j]) * (f[j + 1] - f[j]) / (x[j + 1] - x[j])
    }
    val
  }
  u
}

#' @title create_lowhull
#' @description given a (sorted) set of points, create the lower hull function
#' @param f represent f(x)
#' @param x x points
#' @return a function that calculates the lowerhull for a function
#' @examples
#' log_f <- function(x) log(dnorm(x))
#' x <- seq(-3, 3, by=1)
#' y <- log_f(x)
#' l <- create_lowhull(x, y)
#' l(0)
#' @export
create_lowhull = function(x, f) {
  if (length(x) != length(f)) {
    stop('Length of x and f should be equal.')
  }
  k = length(x)
  l = function(p) {
    val = NA
    j = findInterval(p, x)
    if (j == 0) {
      val = -Inf
    } else if (j == k) {
      val = Inf
    } else {
      val = f[j] + (p - x[j]) * (f[j + 1] - f[j]) / (x[j + 1] - x[j])
    }
    val
  }
  l
}


#' @title calc_uz
#' @description calculate the z points given a (sorted) set of points and their h values and h_prim values
#' @param y_l y points to the left
#' @param y_r y points to the right
#' @param hy_l h(y_l) values
#' @param hy_r h(y_r) values
#' @param hyprim_l h'(y_l)
#' @param hyprim_r h'(y_r)
#' @return a namedVector object with names=z and values=uz
#' @export
calc_uz = function(y_l, y_r, hy_l, hy_r, hyprim_l, hyprim_r) {
  z = y_l + (hy_l - hy_r + (y_r - y_l) * hyprim_r) / (hyprim_r - hyprim_l)
  uz = hy_l + (z - y_l) * hyprim_l

  namedVector$new(z, uz)
}


#' @title integ_expinterpol
#' @description integrate the exponential of interpolating function for each interval given a (sorted) set of points
#' @param f represent f(x)
#' @param x x points
#' @param y may represent different y values from f(x)
#' @return a function that calculates the lowerhull for a function
#' @export
integ_expinterpol = function(x, f, y = x) {
  # detail: given (x1,f1) and (x2,f2), ax+b is the line that goes through both of the two points,
  ## this function calculates I_y1y2 = Int_x1^x2 exp(ax+b) for each interval (y1,y2), and stores them in
  ## namedVector as (names: y2, values: I_y1y2). (names[1] = y[1], values[1] = 0).
  ## Note that the points used to calculate a and b might not be the endpoints of interval.
  if (length(x) != length(y)) {
    stop('Length of inputs should be the same.')
  }
  if (length(f) != length(y)) {
    stop('Length of inputs should be the same.')
  }

  y_l = head(y, n = -1)
  y_r = tail(y, n = -1)
  x_l = head(x, n = -1)
  x_r = tail(x, n = -1)
  f_l = head(f, n = -1)
  f_r = tail(f, n = -1)

  a = (f_r - f_l) / (x_r - x_l)
  b = f_l - a * x_l
  vals = (exp(a * y_r + b) - exp(a * y_l + b)) / a

  namedVector$new(names = y, values = c(0, vals))
}
