# BMI 510 Final Project
### Beatrice Brown-Mulry

This R package was completed for BMI 510 at Emory University. It includes a collection of functions ranging from power analysis to data analysis and visualization. Below is a description of each function included in the package.

## Functions

### logLikBernoulli

Calculates the maximum likelihood estimate (MLE) of the probability parameter `p` for a Bernoulli distribution given a binary data vector.

**Usage Example:**

```R
sample_data = c(1, 0, 0, 0, 1, 1, 1)
mle_p = logLikBernoulli(sample_data)
```
---

### survCurv

Creates and plots a survival curve based on provided status and time data using ggplot2.

**Usage Example:**

```R
status_vector = c(1, 0, 1, 0, 0, 1)
time_vector = c(5, 12, 15, 20, 22, 25)
survival_plot = survCurv(status_vector, time_vector)
print(survival_plot)
```

---

### unscale

Reverts scaling and centering transformations applied to data.

**Usage Example:**

```R
data = 1:10
scaled_data = scale(data)
original_data = unscale(scaled_data)
```

---

### pcApprox

Performs Principal Component Analysis (PCA) and reconstructs an approximation of the original data using a specified number of principal components.

**Usage Example:**

```R
test_data = matrix(rnorm(100), nrow = 10, ncol = 10)
approx_data = pcApprox(test_data, npc = 2)
print(approx_data)
```

---

### standardizeNames

Standardizes the column names of a dataframe or tibble to a specified case format.

**Usage Example:**

```R
data = data.frame(`first name` = c("Alice", "Bob"), `Last.Name` = c("Smith", "Jones"))
new_data = standardizeNames(data)
```
---

### minimumN

Calculates the minimum sample size required for a given power in a one-sample or two-sample t-test.

**Usage Example:**

```R
x1 = rnorm(100, 1, 1)
x2 = rnorm(100, 1.2, 1)
minN = minimumN(x1, x2)
```
---

### downloadRedcapReport

Downloads a specified report from a REDCap project using the REDCap API.

**Usage Example:**

```R
report_tibble = downloadRedcapReport(
  "REDCAP_TOKEN",
  "https://redcap.example.com/api/",
  "12345"
)
```

---

## Installation

This package can be installed from GitHub via the `devtools` R package.

```R
devtools::install_github("beatrice-b-m/BMI510FinalProject")
```
