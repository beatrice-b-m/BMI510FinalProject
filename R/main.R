# This project has a total of 20 points.
#   -   The code for each function below is worth 1 point, with the exceptions of more complex
#       standardizeNames (2 points) and downloadRedcapReport (3 points).
#
#   -   The NAMESPACE , README.md , and .gitignore files in your package are each worth 1 point.
#
#   -   The documentation for each function below is worth 1 point.

# ----------------------------------------------------------------------------------------------
# logLikBernoulli = function(data)
# Write a function that takes a vector like
# data = c(1,0,0,0,1,1,1,) and calculates the parameter p that maximizes the log-likelihood.
# log(P(p|data))
# Use a grid-based search with p in steps of 0.001.

logLikBernoulli = function(data) {
  # create a vector of values between 0.001 and 0.999 to test for p
  p_vec = seq(0.001, 0.999, by = 0.001)

  # use sapply to extract all log likelihoods for the vec of tested p vals
  log_lik_vec = sapply(p_vec, function(p) .BernoulliLogLik(data, p))

  # find the idx corresponding to the max log likelihood
  max_idx = which.max(log_lik_vec)

  # return the value of p that maximizes the log likelihood for the sample
  return (p_vec[max_idx])
}

# using my BernoulliLogLik function from MT-02
.BernoulliLogLik = function(x, p) {
  x = as.integer(x)
  log_likelihood = sum(x * log(p) + (1 - x) * log(1 - p))
  return(log_likelihood)
}

# ----------------------------------------------------------------------------------------------
# survCurv = function(status,time)
# Write a function that takes a numerical vector status and a
# numerical vector time , and calculates and plots a survival
# curve S(t) . Test your function on the dataset here:
# (https://jlucasmckay.bmi.emory.edu/global/bmi510/Labs-Materials/survival.csv).

survCurv = function(status, time) {
  # format the data into a tibble
  data = dplyr::tibble(status=status, time=time)

  # count the number of initial individuals
  # assuming each person is accounted for
  init_n = nrow(data)

  # based on the s_t formula given here:
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059453/
  surv_data = data |>
    arrange(time) |> # sort by time
    mutate(cum_event = cumsum(status == 1)) |> # get the cum. sum of events
    mutate(cum_censor = cumsum(status == 0)) |> # get the cum. sum of censors
    mutate(time_unit = ceiling(time)) |> # round each time up to the nearest integer
    group_by(time_unit) |> # group by rounded time units
    summarise(
      n_event = max(cum_event), # get the max cum. events per time unit
      n_censor = max(cum_censor) # get the max cum. censors per time unit
    ) |>
    mutate(p_surv = (init_n - n_event) / init_n) # calculate the prob. of survival at each time step

  # plot the survival curve
  ggplot2::ggplot(surv_data, ggplot2::aes(x=time_unit,y=p_surv)) +
    ggplot2::geom_step() +
    ggplot2::labs(
      x = "Time",
      y = "Estimated Survival Probability",
      title = "Estimated Survival Over Time"
    ) +
    ggplot2::ylim(0.0, 1.0) +
    ggplot2::geom_point()  # Add points at each time step
}

# ----------------------------------------------------------------------------------------------
# unscale = function(x)
# Write a function that takes a vector that has been put through
# scale and reverses the centering/scaling, if any.

#' Undo Scaling and Centering of Data
#'
#' This function reverts the scaling and centering transformations applied to data.
#' It checks for attributes indicating scaling and centering transformations ('scaled:center'
#' and 'scaled:scale') and applies the inverse operations to restore the original data values.
#'
#' @param x A numeric vector or matrix that has been scaled and/or centered.
#'
#' @return The original data with scaling and centering transformations reversed.
#'         If no scaling or centering has been done, the input data `x` is returned unchanged.
#'
#' @examples
#' # Scale and then unscale data
#' data = 1:10
#' scaled_data = scale(data)
#' original_data = unscale(scaled_data)
#'
#' @export
unscale = function(x) {
  # extract the center/scaling factors from the scaled data
  center_factor = attr(x, "scaled:center")
  scale_factor = attr(x, "scaled:scale")

  # undo scaling if it's been done
  if (!is.null(scale_factor)) {
    x = x * scale_factor
  }
  # undo centering if it's been done
  if (!is.null(center_factor)) {
    x = x + center_factor
  }

  return (x)
}

# ----------------------------------------------------------------------------------------------
# pcApprox = function(x, npc)
# Write a function that returns an approximation to the data x
# based on npc PCs. (Note that the approximation should be
#                    rescaled and centered to match the original data).

#' Perform PCA and Approximate Original Data
#'
#' This function scales and centers the data, performs PCA, and then reconstructs
#' an approximation of the original data using a specified number of principal
#' components. The data is then rescaled and recentered to match the original scale.
#'
#' @param x A numeric matrix of data to be analyzed.
#' @param npc The number of principal components to use for the approximation.
#'
#' @return A matrix that represents an approximation of the original data,
#'         reconstructed using the specified number of principal components and
#'         rescaled to match the scale and center of the original data.
#'
#' @examples
#' test_data = matrix(rnorm(100), nrow = 10, ncol = 10)
#' approx_data = pcApprox(test_data, npc = 2)
#' print(approx_data)
#'
#' @export
pcApprox = function(x, npc) {
  # consulted the R pca resource here:
  # https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf
  # scale the data and extract the scale/center factors for later
  scaled_x = scale(x)
  center_factor = attr(scaled_x, "scaled:center")
  scale_factor = attr(scaled_x, "scaled:scale")

  # apply pca to the data
  pca = prcomp(scaled_x)

  # extract the x and rotation matrices from the pca object up to the npc-th pc
  pca_x = pca$x[, 1:npc]
  pca_rot = pca$rotation[, 1:npc]

  # reconstruct an approximation of the data
  approx_x = pca_x %*% t(pca_rot)

  # unscale/uncenter the approximation
  approx_x = (approx_x + center_factor) + scale_factor

  return (approx_x)
}

# ----------------------------------------------------------------------------------------------
# standardizeNames = function(data)
# Write a wrapper around dplyr::rename_with and
# janitor::make_clean_names that converts the variables in a tibble
# data to "small_camel" case (or another case if you like). The idea
# here is to have a reliable function that standardizes the variable
# names in data you’re dealing with. This function should import
# elements of the janitor and snakecase packages; Roxygen2 will handle
# these dependencies for you

#' Standardize Column Names to a Specific Case
#'
#' This function standardizes the column names of a dataframe to a specific naming convention.
#' It uses `janitor::make_clean_names` and `snakecase::to_any_case` to clean columns names and convert them
#' to the specified case format. The case format can be adjusted via the `target_case` parameter.
#'
#' @param data A dataframe or tibble whose column names are to be standardized.
#' @param target_case A string specifying the target case format for the column names.
#'        The default is "small_camel". Supports all case types supported by
#'        `snakecase::to_any_case` including "snake", "small_camel", "big_camel",
#'        and "screaming_snake".
#'
#' @return A dataframe or tibble with standardized column names.
#'
#' @examples
#' data = data.frame(`First name` = c("Alice", "Bob"), `Last.Name` = c("Smith", "Jones"))
#' new_data = standardizeNames(data, target_case = "snake")
#'
#' @export
standardizeNames = function(data, target_case="small_camel") {
  out_data = dplyr::rename_with(
    data, ~ janitor::make_clean_names(.) |> snakecase::to_any_case(case=target_case),
    .cols=everything()
  )
  return (out_data)
}

# ----------------------------------------------------------------------------------------------
# minimumN = function(x1,x2)
# Write a wrapper around pwr::pwr.t2n.test that takes either one ( x1 )
# or two ( x2 ) samples of preliminary data and returns the minimum sample
# size needed for a t-test of the null hypotheses that either mu_x1 == 0
# or mu_x1 == mu_x2 with 80% power at alpha=0.05.

#' Calculate Minimum Sample Size for Power Analysis
#'
#' This function calculates the minimum sample size required to achieve a given power for a one-sample
#' or two-sample t-test. For a one-sample test, the null hypothesis is that the mean of `x1` is zero.
#' For a two-sample test, the null hypothesis is that the means of `x1` and `x2` are equal.
#'
#' @param x1 Numeric vector for the first group or the single sample.
#' @param x2 Optional numeric vector of data for the second sample (default is `NULL` for one-sample tests).
#' @param alpha Desired significance level (default is 0.05).
#' @param power Desired power of the test (default is 0.80).
#'
#' @return The ceiling of the calculated minimum sample size necessary to achieve the specified power
#'         for the given significance level and effect size.
#'
#' @importFrom pwr pwr.t.test
#' @examples
#' # For a one-sample test
#' x1 <- rnorm(100, 1, 1)
#' minimumN(x1)
#'
#' # For a two-sample test
#' x2 <- rnorm(100, 1.2, 1)
#' minimumN(x1, x2)
#'
#' @export
# Define the function
minimumN = function(x1, x2=NULL, alpha=0.05, power=0.80) {
  if (is.null(x2)) {
    # if x2 is null do a one sample test
    # with null hypothesis mu_x1 == 0

    # calculate stats for x1 and x2
    stats1 = getStats(x1)

    # calculate effect size
    effect_size = stats1$Mean / stats1$SD

    test_type = "one.sample"

  } else {
    # otherwise do a two sample test
    # with null hypothesis mu_x1 == mu_x2

    # calculate stats for x1 and x2
    stats1 = getStats(x1)
    stats2 = getStats(x2)

    # based on the implementation here:
    # https://www.statology.org/pooled-standard-deviation-in-r/
    # calculate pooled standard deviation and effect size
    pooled_sd = sqrt(((stats1$N-1)*stats1$SD^2 + (stats2$N-1)*stats2$SD^2) / (stats1$N+stats2$N-2))
    effect_size = abs(stats1$Mean - stats2$Mean) / pooled_sd
    test_type = "two.sample"
  }

  # perform power analysis
  result = pwr::pwr.t.test(
    d=effect_size,
    sig.level=alpha,
    power=power,
    type=test_type,
    alternative="two.sided"
  )

  # extract the minimum sample size from the results and return it
  return (ceiling(result$n))
}

getStats = function(x) {
  return (list(Mean=mean(x), SD=sd(x), N=length(x)))
}

# ----------------------------------------------------------------------------------------------
# downloadRedcapReport = function(redcapTokenName,redcapUrl,redcapReportId)
# Using the block of RedCap template code below, write a function that:
#   -   uses Sys.getenv() to read an API token called redcapTokenName
#       from the users’ .REnviron file.
#
#   -   queries redcapUrl to return the Redcap Report redcapReportId .
#       (Notice these are the data from our simulated stroke study,
#       now nicely and securely hosted on RedCap.)
#
#   -   returns the contents as a tibble.

#' Download a Report from REDCap
#'
#' This function downloads a specified report from a REDCap project using the REDCap API.
#' The report is requested via POST request, retrieved in CSV format, and returned as a tibble.
#'
#' @param redcapTokenName A string specifying the name of the environment variable where
#'   the REDCap API token is stored. This token is used to authenticate the request.
#' @param redcapUrl A string containing the base URL of the REDCap instance.
#' @param redcapReportId An integer or string specifying the unique identifier of the report to be downloaded.
#'
#' @return Returns a tibble of the requested report.
#' @examples
#' # Assuming you have a REDCap token stored in an
#' # environment variable named "REDCAP_TOKEN" and
#' # the REDCap instance URL is: "https://redcap.example.com/api/"
#' # for a report with ID: 12345
#' report_tibble = downloadRedcapReport(
#'   "REDCAP_TOKEN",
#'   "https://redcap.example.com/api/",
#'   "12345"
#' )
#'
#' @export
downloadRedcapReport = function(redcapTokenName, redcapUrl, redcapReportId) {
  # retrieve the redcapTokenName from environment variables
  redcap_token = Sys.getenv(redcapTokenName)

  # format the formData
  formData = list(
    "token"=redcap_token,
    content='report',
    format='csv',
    report_id=redcapReportId,
    csvDelimiter='',
    rawOrLabel='raw',
    rawOrLabelHeaders='raw',
    exportCheckboxLabel='false',
    returnFormat='csv'
  )
  # submit the request and parse the content
  response = httr::POST(redcapUrl, body=formData, encode="form")
  result = httr::content(response)

  return (result)
}
