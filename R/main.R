# ----------------------------------------------------------------------------------------------
#' Calculate the Maximum Likelihood Estimate of Probability for a Bernoulli Distribution
#'
#' This function calculates the maximum likelihood estimate (MLE) of the probability parameter `p`
#' for a Bernoulli distribution given a binary data vector. It searches across a grid of potential
#' `p` values ranging from 0.001 to 0.999 in steps of 0.001.
#'
#' @param data A binary vector (0s and 1s) representing the trial outcomes.
#'
#' @return The value of `p` that maximizes the Bernoulli log-likelihood for the provided data.
#'
#' @examples
#' sample_data = c(1, 0, 0, 0, 1, 1, 1)
#' mle_p = logLikBernoulli(sample_data)
#'
#' @export
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

#' Estimate and Plot a Survival Curve
#'
#' This function creates a survival curve based on status and time data provided. It first formats the data into a
#' tibble, calculates cumulative events and censorships, and then computes the survival probability at each unique
#' time unit. Finally, it plots the resulting survival curve using ggplot2.
#'
#' @param status A numeric vector where '1' indicates an event (e.g., death) and '0' indicates censoring.
#' @param time A numeric vector representing the time at which each event or censoring occurred.
#'
#' @return A ggplot object displaying the estimated survival curve over time.
#'
#' @examples
#' status_vector <- c(1, 0, 1, 0, 0, 1)
#' time_vector <- c(5, 12, 15, 20, 22, 25)
#' survival_plot = survCurv(status_vector, time_vector)
#' print(survival_plot)
#'
#' @export
survCurv = function(status, time) {
  # format the data into a tibble
  data = dplyr::tibble(status=status, time=time)

  # count the number of initial individuals
  # assuming each person has one row
  init_n = nrow(data)

  # based on the s_t formula given here:
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059453/
  surv_data = data |>
    dplyr::arrange(time) |> # sort by time
    dplyr::mutate(cum_event = cumsum(status == 1)) |> # get the cum. sum of events
    dplyr::mutate(cum_censor = cumsum(status == 0)) |> # get the cum. sum of censors
    dplyr::mutate(time_unit = ceiling(time)) |> # round each time up to the nearest integer
    dplyr::group_by(time_unit) |> # group by rounded time units
    dplyr::summarise(
      n_event = max(cum_event), # get the max cum. events per time unit
      n_censor = max(cum_censor) # get the max cum. censors per time unit
    ) |>
    dplyr::mutate(p_surv = (init_n - n_event) / init_n) # calculate the prob. of survival at each time step

  # plot the survival curve
  survival_curve = ggplot2::ggplot(surv_data, ggplot2::aes(x=time_unit,y=p_surv)) +
    ggplot2::geom_step() +
    ggplot2::labs(
      x = "Time",
      y = "Survival Probability",
      title = "Estimated Survival Over Time"
    ) +
    ggplot2::ylim(0.0, 1.0) +
    ggplot2::geom_point()  # Add points at each time step

  return(survival_curve)
}

# ----------------------------------------------------------------------------------------------

#' Undo Scaling and Centering of Data
#'
#' This function reverts any scaling and centering transformations applied to data.
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
  approx_x = (approx_x * scale_factor) + center_factor

  return (approx_x)
}

# ----------------------------------------------------------------------------------------------

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
#' data = data.frame(`first name` = c("Alice", "Bob"), `Last.Name` = c("Smith", "Jones"))
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

  # format the form data
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
