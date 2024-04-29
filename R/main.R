devtools::create("/Users/beatrice/Library/Mobile Documents/com~apple~CloudDocs/Files/College/BMI 510; Biostatistics for ML/BMI510FinalProject/")

# This project has a total of 20 points.
#   -   The code for each function below is worth 1 point, with the exceptions of more complex
#       standardizeNames (2 points) and downloadRedcapReport (3 points).
#
#   -   The NAMESPACE , README.md , and .gitignore files in your package are each worth 1 point.
#
#   -   The documentation for each function below is worth 1 point.

# ---------------------------------------------------------------------
# logLikBernoulli = function(data)
# Write a function that takes a vector like
# data = c(1,0,0,0,1,1,1,) and calculates the parameter p that maximizes the log-likelihood.
# log(P(p|data))

# Use a grid-based search with p in steps of 0.001.

# existing func:
# BernoulliLogLik = function(x, p) {
#   x = as.integer(x)
#   log_likelihood = sum(x * log(p) + (1 - x) * log(1 - p))
#   return(log_likelihood)
# }

# ---------------------------------------------------------------------
# survCurv = function(status,time)
# Write a function that takes a numerical vector status and a
# numerical vector time , and calculates and plots a survival 
# curve S(t) . Test your function on the dataset here: 
# (https://jlucasmckay.bmi.emory.edu/global/bmi510/Labs-Materials/survival.csv).

# ---------------------------------------------------------------------
# unscale = function(x)
# Write a function that takes a vector that has been put through 
# scale and reverses the centering/scaling, if any.

# ---------------------------------------------------------------------
# pcApprox = function(x, npc)
# Write a function that returns an approximation to the data x 
# based on npc PCs. (Note that the approximation should be 
#                    rescaled and centered to match the original data).

# ---------------------------------------------------------------------
# standardizeNames = function(data)
# Write a wrapper around dplyr::rename_with and
# janitor::make_clean_names that converts the variables in a tibble 
# data to "small_camel" case (or another case if you like). The idea 
# here is to have a reliable function that standardizes the variable 
# names in data you’re dealing with. This function should import 
# elements of the janitor and snakecase packages; Roxygen2 will handle 
# these dependencies for you

# ---------------------------------------------------------------------
# minimumN = function(x1,x2)
# Write a wrapper around pwr::pwr.t2n.test that takes either one ( x1 )
# or two ( x2 ) samples of preliminary data and returns the minimum sample 
# size needed for a t-test of the null hypotheses that either mu_x1 == 0 
# or mu_x1 == mu_x2 with 80% power at alpha=0.05.

# ---------------------------------------------------------------------
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

print(Sys.getenv('redcapTokenName'))

downloadRedcapReport = function(redcapTokenName, redcapUrl, redcapReportId) {
  formData = list("token"=token,
                   content='report',
                   format='csv',
                   report_id='46524',
                   csvDelimiter='',
                   rawOrLabel='raw',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='false',
                   returnFormat='csv'
  )
  response = httr::POST(url, body=formData, encode="form")
  result = httr::content(response)
}

redcap_url = "https://redcap.emory.edu/api/"
print(result)