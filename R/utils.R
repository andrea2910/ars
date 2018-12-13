#' Class providing object with methods for communication with R6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @return Object of \code{\link{R6Class}} with methods to manipulate data points
#' @format \code{\link{R6Class}} object.
#' @examples
#' ex <- namedVector$new(names = c(1,2,3), values = c(4,5,6))
#' @field names Stores points
#' @field values Stores y values associated with points
#' @section Methods:
#' \describe{
#'   \item{\code{add(nams, vals)}}{This method adds names and values to the \code{namedVector}}
#'
#'   \item{\code{delete_byname(nams)}}{This method deletes both the name and value of inputted nams from \code{namedVector}}
#'   \item{\code{sort_byname(decreasing=FALSE)}}{This method sorts the \code{namedVector} by name values.}}
namedVector = R6Class(
  "namedVector",
  lock_objects = FALSE,
  public = list(
    ### functions ###
    # initialization
    initialize = function(names, values) {
      if (length(names) != length(values)) {
        stop('Length of names and values should be equal.')
      }
      self$names = names
      self$values = values
    },

    add = function(nams, vals) {
      # delete repeated names in nams
      duplicated_ind = duplicated(nams)
      nams = nams[!duplicated_ind]
      vals = vals[!duplicated_ind]
      # delete names in nams that are already in self$names
      index_names = !(nams %in% self$names)
      nams = nams[index_names]
      vals = vals[index_names]

      # delete values that are NaN
      #no_nan_ind = !is.nan(vals)
      #nams = nams[no_nan_ind]
      #vals = vals[no_nan_ind]

      # add
      self$names = c(self$names, nams)
      self$values = c(self$values, vals)
    },

    # delete an element by name
    delete_byname = function(nams) {
      if (length(levels(factor(self$names))) != length(self$names) ||
          length(levels(factor(self$nams))) != length(self$nams)) {
        stop('For now this method only functions when names are unique.')
      }
      index = match(nams, self$names)
      index = index[!is.na(index)] # drop NAs
      self$names = self$names[-index]
      self$values = self$values[-index]
    },

    # sort namedVector by name
    sort_byname = function(decreasing = FALSE) {
      if (length(levels(factor(self$names))) != length(self$names)) {
        stop('For now this method only functions when names are unique.')
      }
      sorted_names = sort(self$names, decreasing = decreasing)
      sorted_values = self$values
      for (i in 1:length(sorted_names)) {
        index = match(sorted_names[i], self$names)
        sorted_values[i] = self$values[index]
      }
      self$names = sorted_names
      self$values = sorted_values
    },

    ### variables ###
    names = c(),
    values = c()
  )
)

#' @title check_positive
#' @description Check whether a function has negative values
#' @importFrom rootSolve uniroot.all
#' @param f a function to search across
#' @param var_min the lower bound domain to check with f
#' @param var_max the upper bound domain to check with f
#' @return a boolean value that represents whether a function is only positive
#' @examples
#' check_positive(dnorm, -Inf, Inf)
#' @export
check_positive = function(f, var_min, var_max) {
  # choose any point in interval (var_min, var_max)
  new_f = function(x) f(x) + .Machine$double.eps

  # find root
  f_var_min <- new_f(var_min)
  f_var_max <- new_f(var_max)

  if (sign(f_var_min) != sign(f_var_max) || sign(f_var_min) < 1)
    return(FALSE) # if the bound valeus differ

  # update var_min and var_max because we need bounds
  if (is.infinite(var_min))
    var_min = - .Machine$double.xmax
  if (is.infinite(var_max))
    var_max = .Machine$double.xmax

  results <-
    try(uniroot.all(Vectorize(new_f), lower = var_min, upper = var_max))
  if (class(results) == "try-error")
    stop("Error using uniroot.all. Try different values")
  if (length(results) == 0) {
    # no roots found
    return(TRUE)
  } else{
    return(FALSE)
  }
}

#' @title calc_deriv
#' @description Calculate the derivative of function f: f'(x) with an increment of machine epsilon
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

  # find a small difference value
  # note: when x is significantly large, we need to have larger eps to avoid underflow
  eps <- (.Machine$double.eps) ^ (1/4) * max(abs(x)/10e8, 1)
  if (abs(x) > 10e8)
    warning("Very hard to get good approximation of the derivatives at this scale.")

  d <- NA
  if (x >= lower && x <= lower + eps)
    d <- (f(x + eps, ...) - f(x, ...)) / eps
  else if (x <= upper && x >= upper - eps)
    d <-  (f(x, ...) - f (x - eps, ...)) / eps
  else if (x > lower + eps && x < upper - eps)
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
#' check if the interpolating function is concave he function will loop through all points in x and calculate the interpolated value
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
  if(class(x)!="numeric") stop("Input must be numeric")
  if(class(f)!="numeric") stop("Input must be numeric")

  sig <- TRUE
  for (i in 1:(length(x) - 2)) {
    inter_f = f[i] + (x[i + 1] - x[i]) * (f[i + 2] - f[i]) / (x[i + 2] - x[i])
    if (inter_f > f[i + 1] + .Machine$double.eps*abs(f[i+1])) {
      sig <- FALSE
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
      val = -Inf
    } else {
      val = f[j] + (p - x[j]) * (f[j + 1] - f[j]) / (x[j + 1] - x[j])
    }
    val
  }
  l
}


#' @title calc_uz
#' @description calculate the z points (or intersection between the tangent lines at given points) given a (sorted) set of points and their h values and h_prim values
#' @param y_l y points to the left
#' @param y_r y points to the right
#' @param hy_l h(y_l) values
#' @param hy_r h(y_r) values
#' @param hyprim_l h'(y_l)
#' @param hyprim_r h'(y_r)
#' @return a namedVector object with names=z and values=uz
#' @export
calc_uz = function(y_l, y_r, hy_l, hy_r, hyprim_l, hyprim_r) {
  z = y_l + (hy_l - hy_r + (y_r - y_l) * hyprim_r) / round(hyprim_r - hyprim_l, 8)
  uz = hy_l + (z - y_l) * hyprim_l

  # additional check to catch certain cases
  replace_w_nan = round(hyprim_r - hyprim_l, 8) == 0
  z[replace_w_nan] = NaN
  uz[replace_w_nan] = NaN

  namedVector$new(z, uz)
}


#' @title integ_expinterpol
#' @importFrom utils head tail
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

  ind = is.infinite(vals) | is.nan(vals)
  vals[ind] = exp(b[ind]) * (y_r[ind] - y_l[ind])

  namedVector$new(names = y, values = c(0, vals))
}
