library(R6)

#' Class providing object with methods for communication with R6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @return Object of \code{\link{R6Class}} with methods to manipulate object
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
      # delete names in nams taht are already in self$names
      index_names = !(nams %in% self$names)
      nams = nams[index_names]
      vals = vals[index_names]
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
      if (NA %in% index) {
        stop("Couldn\'t find some of the passed names.")
      }
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

#' Class providing object with methods for communication with R6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @return Object of \code{\link{R6Class}} with methods to generate adaptive rejection sampling
#' @param funx a function we want to sample from
#' @param D a vector of length 2 of numbers for bounds of function
#' @format \code{\link{R6Class}} object.
#' @examples
#' ex <- ARS$new(dnorm, c(-Inf,Inf),mean=3, sd=1)
#' samples <- ex$sample(n=2000)
#' ex$plot_samples()
#' @field self$f Function we want to sample from
#' @field self$var_max Upper bound on function we want to sample from
#' @field self$var_min Lower bound on function we want to sample from
#' @field private$h_prims A R6 namedVector that contains h'(x)
#' @field private$h_vals A R6 namedVector that contains h(x)
#' @field private$s_cdfs A R6 namedVector that contains each interval's CDF
#' @field private$u_vals A R6 namedVector that contains u(x)
#' @field private$x a vector of x values that are our sample values. Empty unless sample() run
#' @field private$y a vector of y values that begins to increase once sample() is run
#' @field private$z a vector of z values, or intersection points, that begins to increase once sample() is run
#' @section Public Methods:
#' \describe{
#'   \item{\code{f(x)}}{This method calculates our funx observation at a given x }
#'   \item{\code{plot_sampdist()}}{This method plots our envelope function }
#'   \item{\code{calc_sampdist()}}{This method calculates the sampling distribution from a uniform distribution}
#'   \item{\code{plot_samples()}}{This method plots our samples in a histogram}
#'   \item{\code{s()}}{This method calculates the envelope or sampling function }
#'   \item{\code{sample(n)}}{This method samples n points from using adaptive rejection sampling}}
#' @section Private Methods:
#' \describe{
#'   \item{\code{h(x)}}{This method calculates the log of our funcx at a given point x}
#'   \item{\code{init_hprim()}}{This method intializes our hprim values}
#'   \item{\code{init_hval()}}{This method initalizes our h values}
#'   \item{\code{init_scdf(}}{This method initalizes our cdf under the upper hull}
#'   \item{\code{init_u()}}{This method initializes our upper hull values}
#'   \item{\code{init_y()}}{This method initalizes our y values }
#'   \item{\code{l(p)}}{This method is a function that returns the slope of a lowerhull line given a point p}
#'   \item{\code{u(p)}}{This method is a function that returns the slope of a upperhull line given a point p}
#'   \item{\code{sampl_exph(}}{This method samples from our exp(h values)}
#'   \item{\code{update(y, hy, hy_prim)}}{This method updates our private variables if we reject a sample x_star}}
ARS = R6Class(
  "ARS",
  lock_objects = FALSE,
  public = list(
    ### functions ###
    # initialization
    initialize = function(funx, D = c(-Inf, Inf), ...) {
      # input: function and its support, initial points, minimum accepted number of points
      # note: should we give default value of range?

      if (class(funx) != "function")
        stop("funx must be a function")

      f = function(x)
        funx(x, ...)

      if (length(D) != 2)
        stop("Input dimension is not a length of 2")

      var_min = min(D)
      var_max = max(D)

      if (var_min == var_max)
        stop("Input dimension must differ")

      # check if the input function is valid
      # note: anything else to check here?
      if (!check_positive(f, var_min, var_max)) {
        stop('Input function is not positive in the given support.')
      }

      Int_f = integrate(Vectorize(f), lower = var_min, upper =
                          var_max)
      if (abs(Int_f$value - 1.) > 10 * Int_f$abs.error) {
        # note: change tolerance level?
        warning('Input function is not normalized. The pdf does not integrate to 1.')
      }

      # set values
      self$f = f

      private$h = function(x) {
        val = f(x)
        if (val <= 0) {
          stop('Input function is not positive in the given support.')
        }
        log(f(x))
      }

      self$var_min = var_min
      self$var_max = var_max
    },
    # initalize ends here

    calc_sampdist = function() {
      exp_u = function(x)
        exp(private$u(x))
      self$s = function(x)
        exp_u(x) / integrate(Vectorize(exp_u),
                             lower = self$var_min,
                             upper = self$var_max)$value
    },

    # main function. sample points from distribution
    sample = function(n = 1000) {
      start_npts = as.integer(n ^ (1 / 10)) + 1
      # note: this is just an example of choosing the number of starting points

      # sample starting points
      private$init_y()
      if (self$var_min == -Inf && self$var_max == Inf) {
        private$y = c(private$y, rnorm(
          start_npts - length(private$y),
          mean = mean(private$y),
          sd = sd(private$y)
        ))
      } else if (self$var_min == -Inf) {
        private$y = c(private$y,-rexp(start_npts - length(private$y)) + self$var_max)
      } else if (self$var_max == Inf) {
        private$y = c(private$y, rexp(start_npts - length(private$y)) + self$var_min)
      } else {
        private$y = c(private$y,
                      runif(
                        start_npts - length(private$y),
                        min = self$var_min,
                        max = self$var_max
                      ))
      }
      # note: this is just an example of sampling starting points

      # construct the unnormalized cdf of s and everything needed for the construction,
      ## which is everything that is related to the upper hull
      private$y <- sort(private$y)
      private$init_scdf()
      private$construct()

      # sample and update
      private$x = vector(mode = 'numeric', length = n)
      samp = 0
      loop = 0
      while (samp < n && loop < 10 * n) {
        # note: l is the number of loops conducted.
        ## This is just an example of the criterion to stop the loop.

        ## sample ##
        x_star = private$samp_exph()
        u_star = runif(1)


        if (u_star <= exp(private$l(x_star) - private$u(x_star))) {
          samp = samp + 1
          private$x[samp] = x_star
        } else {
          hxstar = private$h(x_star)
          hxstar_prim = calc_deriv(x_star, private$h, self$var_min, self$var_max)
          if (u_star <= exp(hxstar - private$u(x_star))) {
            samp = samp + 1
            private$x[samp] = x_star
          }
          ## update ##
          private$update(x_star, hxstar, hxstar_prim)
          private$construct()
        }
        loop = loop + 1
      }
      self$calc_sampdist()
      private$x
    },

    plot_samples = function() {
      library(ggplot2)

      # plot median to ignore outliers
      ggplot() + geom_histogram(aes(x = private$x)) +
        geom_vline(xintercept = median(private$x), color = "blue") +
        theme_bw() +
        ggtitle("Sample distribution") + xlab("Samples")

    },

    plot_sampdist = function() {
      library(ggplot2)

      self$calc_sampdist()
      xs <-
        seq(private$y[1], private$y[length(private$y)],
            length.out = 100)
      ys <- sapply(xs, self$s)

      ggplot() + geom_line(aes(x = xs, y = ys)) + theme_bw() +
        ggtitle("S distribution") + xlab("x") +
        ylab("P(S(X)=x)")

    },
    # print = function(){
    #   #invisible(self)
    # },
    ### variables ###
    # density function
    f = NULL,
    # min and max value of the support of the function
    var_min = NA,
    var_max = NA,
    # sampling function
    s = NULL
  ),

  private = list(
    ### functions ###
    # choose initial points
    init_y = function() {
      # (re)initialize sample points vector
      private$y = vector(mode = 'numeric')
      # if unbounded on the left
      if (self$var_min == -Inf) {
        y_min = -1
        iter = 0
        while (calc_deriv(y_min,
                          private$h,
                          lower = self$var_min,
                          upper = self$var_max) <= 0 && iter < 1000) {
          y_min = y_min - runif(1) * abs(private$h(y_min))
          # note: explain the reason to choose a random number
          ## and the reason to multiply by the function value.
          iter = iter + 1
        }
        if (iter >= 1000) {
          stop('Unlucky. Have searched over 1,000 times to find a h prime less than or equal to 0')
        }
        private$y = c(private$y, y_min)
      }
      # if unbounded on the right
      if (self$var_max == Inf) {
        y_max = 1
        iter = 0
        while (calc_deriv(y_max,
                          private$h,
                          lower = self$var_min,
                          upper = self$var_max) >= 0 && iter < 1000) {
          y_max = y_max + runif(1) * abs(private$h(y_max))
          # note: explain the reason to choose a random number
          ## and the reason to multiply by the function value.
          iter = iter + 1
        }
        if (iter >= 1000) {
          stop('Unlucky. Have searched over 1,000 tims to find a h prime larger than or equal to 0')
        }
        private$y = c(private$y, y_max)
      }
    },

    # initialize values of sampled points, called by init_u
    init_hval = function() {
      # calculate derivatives
      values = sapply(private$y, private$h)
      # (re)initialize derivatives of sample points
      private$h_vals = namedVector$new(names = private$y, values =
                                         values)
    },

    # initialize derivatives of sampled points, called by init_u
    init_hprim = function() {
      # calculate derivatives
      derivs = sapply(
        private$y,
        FUN = calc_deriv,
        f = private$h,
        lower = self$var_min,
        upper = self$var_max
      )
      # (re)initialize derivatives of sample points
      private$h_prims = namedVector$new(names = private$y, values =
                                          derivs)
    },

    # initialize upper hull
    init_u = function() {
      # set points
      m = length(private$y)
      y = private$y
      # set h values
      private$init_hval()
      hy = private$h_vals$values
      # set h derivatives
      private$init_hprim()
      hy_prim = private$h_prims$values
      # calculate Z
      uz = calc_uz(y[1:(m - 1)], y[2:m], hy[1:(m - 1)], hy[2:m], hy_prim[1:(m -
                                                                              1)], hy_prim[2:m])
      # check concavity
      u_pts = c(y[1], uz$names, y[m])
      u_vals = c(hy[1], uz$values, hy[m])
      if (!check_interpolconcave(u_pts, u_vals)) {
        stop('Input function is not log-concave.')
      }
      # (re)initialize z
      private$z = uz$names
      # (re)initialize upper hull values
      private$u_vals = namedVector$new(names = u_pts, values =
                                         u_vals)
      # (re)create upper hull function
      private$u = create_upphull(u_pts, u_vals)
    },

    # initialize cdfs of u
    init_scdf = function() {
      # set u
      private$init_u()
      u_pts = private$u_vals$names
      u_vals = private$u_vals$values
      # (re)initialize cdf
      cdf_pts = c(self$var_min, private$z, self$var_max)
      private$s_cdfs = integ_expinterpol(u_pts, u_vals, cdf_pts)
      # should it have 3?
    },

    # construct function u and l
    construct = function() {
      private$u = create_upphull(private$u_vals$names, private$u_vals$values)
      private$l = create_lowhull(private$h_vals$names, private$h_vals$values)
    },
    # essential function 1 of main function: sample a point from s
    samp_exph = function() {
      # generate u
      rand = runif(1, min = 0, max = sum(private$s_cdfs$values))
      # find the interval of x where F(x) = u

      cumm_cdfs = cumsum(private$s_cdfs$values)

      interv = findInterval(rand, cumm_cdfs)

      cdf = cumm_cdfs[interv]


      # calculate x and return the value
      u_pts = private$u_vals$names
      u_vals = private$u_vals$values
      cdf_pts = private$s_cdfs$names

      a = (u_vals[interv + 1] - u_vals[interv]) / (u_pts[interv +
                                                           1] - u_pts[interv])
      b = u_vals[interv] - a * u_pts[interv]

      n_test = length(cdf_pts)
      a_test = (u_vals[2:n_test] - u_vals[1:(n_test - 1)]) / (u_pts[2:n_test] - u_pts[1:(n_test -
                                                                                           1)])
      b_test = u_vals[1:(n_test - 1)] - a_test * u_pts[1:(n_test -
                                                            1)]

      (log(a * (rand - cdf) + exp(a * cdf_pts[interv] + b)) - b) /
        a
    },

    # essential function 2 of main function: update
    ## update everything that needs to be updated given a new y and its h value and h_prim value
    update = function(y, hy, hy_prim) {
      # note: updating needs identifying the interval the new value is in and a certain sequence,
      ## it is probably better to write it in a single function.

      # 1. update z, u_vals, s_cdfs backwards
      l_y = length(private$y)
      l_z = length(private$z)

      int = findInterval(y, private$y)

      if (int == 0) {
        # calculate new z
        u_z = calc_uz(y,
                      private$y[1],
                      hy,
                      private$h_vals$values[1],
                      hy_prim,
                      private$h_prims$values[1])
        z_new = u_z$names
        uz_new = u_z$values
        # check concavity
        uvals_new = private$u_vals$clone()
        uvals_new$delete_byname(private$y[1])
        uvals_new$add(c(y, z_new), c(hy, uz_new))
        uvals_new$sort_byname()
        if (!check_interpolconcave(uvals_new$names, uvals_new$values)) {
          stop('Input function is not log-concave.')
        }
        # calculate new cdf
        calc_pts = c(y, z_new, private$z[1])
        calc_vals = c(hy, uz_new, private$u_vals$values[2])
        interv_pts = c(self$var_min, z_new, private$z[1])

        scdf_new = integ_expinterpol(calc_pts, calc_vals, interv_pts)
        # update s_cdfs
        private$s_cdfs$delete_byname(private$z[1])
        private$s_cdfs$add(scdf_new$names, scdf_new$values)
        private$s_cdfs$sort_byname()
        # update u_vals
        private$u_vals$delete_byname(private$y[1])
        private$u_vals$add(c(y, z_new), c(hy, uz_new))
        private$u_vals$sort_byname()
        # update z
        private$z = sort(c(private$z, z_new))
      } else if (int == l_y) {
        # calculate new z
        u_z = calc_uz(private$y[l_y],
                      y,
                      private$h_vals[l_y],
                      hy,
                      private$h_prims[l_y],
                      hy_prim)
        z_new = u_z$names
        uz_new = u_z$values
        # check concavity
        uvals_new = private$u_vals$clone()
        uvals_new$delete_byname(private$y[l_y])
        uvals_new$add(c(z_new, y), c(uz_new, hy))
        uvals_new$sort_byname()
        if (!check_interpolconcave(uvals_new$names, uvals_new$values)) {
          stop('Input function is not log-concave.')
        }
        # calculate new cdf
        calc_pts = c(private$z[l_z], z_new, y)
        calc_vals = c(private$u_vals$values[l_z + 1], uz_new, hy)
        interv_pts = c(private$z[l_z], z_new, self$var_max)

        scdf_new = integ_expinterpol(calc_pts, calc_vals, interv_pts)
        # update s_cdf
        private$s_cdfs$delete_byname(self$var_max)
        private$s_cdfs$add(scdf_new$names, scdf_new$values)
        private$s_cdfs$sort_byname()
        # update u_vals
        private$u_vals$delete_byname(private$y[l_y])
        private$u_vals$add(c(z_new, y), c(uz_new, hy))
        private$u_vals$sort_byname()
        # update z
        private$z = sort(c(private$z, z_new))
        # note: this could potentially merge with the previous situation
      } else {
        # calculate z

        y_l = c(private$y[int], y)
        y_r = c(y, private$y[int + 1])

        hy_l = c(private$h_vals$values[int], hy)
        hy_r = c(hy, private$h_vals$values[int + 1])
        hyprim_l = c(private$h_prims$values[int], hy_prim)
        hyprim_r = c(hy_prim, private$h_prims$values[int + 1])

        u_z = calc_uz(y_l, y_r, hy_l, hy_r, hyprim_l, hyprim_r)
        z_new = u_z$names
        uz_new = u_z$values

        # check concavity
        uvals_new = private$u_vals$clone()

        uvals_new$delete_byname(private$z[int])
        uvals_new$add(z_new, uz_new)
        uvals_new$sort_byname()
        if (!check_interpolconcave(uvals_new$names, uvals_new$values)) {
          stop('Input function is not log-concave.')
        }
        # calculate new cdf
        calc_pts = c(private$u_vals$names[int], z_new, private$u_vals$names[int +
                                                                              2])
        calc_vals = c(private$u_vals$values[int], uz_new, private$u_vals$values[int +
                                                                                  2])


        interv_pts = c(private$s_cdfs$names[int], z_new, private$s_cdfs$names[int +
                                                                                2])

        scdf_new = integ_expinterpol(calc_pts, calc_vals, interv_pts)
        # update s_cdf
        private$s_cdfs$delete_byname(private$s_cdfs$names[(int +
                                                             1):(int + 2)])
        private$s_cdfs$add(scdf_new$names, scdf_new$values)
        private$s_cdfs$sort_byname()
        # update u_vals
        private$u_vals$delete_byname(private$z[int])
        private$u_vals$add(z_new, uz_new)
        private$u_vals$sort_byname()
        # update z
        private$z = private$z[-int]
        private$z = sort(c(private$z, z_new))
      }
      # 2. update y, h_vals, h_prims
      # update y
      private$y = c(private$y, y)
      private$y = sort(private$y)
      # update h_vals
      private$h_vals$add(y, hy)
      private$h_vals$sort_byname()
      # update h_prims
      private$h_prims$add(y, hy_prim)
      private$h_prims$sort_byname()
      # recreate u
      private$u = create_upphull(private$u_vals$names, private$u_vals$values)
    },

    ### variables ###
    # log of the input function
    h = NULL,
    # construction points
    y = c(),
    # values of sampled points
    h_vals = NULL,
    # derivatives of sampled points
    h_prims = NULL,
    # points of upper hull
    z = c(),
    # upper hull function
    u = NULL,
    # upper hull, stored pointwise(include the value of all z points and the first and last y points)
    u_vals = NULL,
    # cdf of each interval of exponential upper hull, unnormalized and stored pointwise
    s_cdfs = NULL,
    # lower hull function
    l = NULL,
    # sampled points
    x = c()
  )
)
