###------------------------------------------------------------------###
###---------------Installing the required packages-------------------###
###------------------------------------------------------------------###


install.packages('MSwM',dependencies = TRUE)
install.packages(c('dplyr','moments','janitor','gtsummary',
                   'rstatix','fastDummies','flextable','tidyverse',
                   'officer', 'gridExtra', 'cowplot','modelsummary',
                   'brant'),dependencies = TRUE)
library(zoo)
library(vars)
library(factoextra)
library(corrplot)
library(dynlm)
library(forecast)
library(moments) # to get summary stats
library(janitor) # to get summary stats
library(gtsummary) # to get summary stats
library(rstatix) # to get summary stats
library(fastDummies)
library(flextable) # to create aesthetic table
library(tidyverse)
library(officer)
library(gridExtra) # to plot multiple graphs on same page
library(cowplot) # to plot multiple graphs on same page
library(broom) # convert model results to table
library(modelsummary) # convert model results to table, even compare models
library(dplyr) # Need to install first
library(MSwM) # Need to install first
library(lmtest)
library(stats)
options(scipen = 9) # Avoid scientific notation


###------------------------------------------------------------------###
###---------------Loading & transforming data------------------------###
###------------------------------------------------------------------###

# Loading the data
if('GVAR Data Malaysia.csv' %in% list.files()){
  dframe <- read.table("GVAR Data Malaysia.csv", sep =',', header = T)
} else { # file.choose allows pop-up to choose file
  dframe <- read.table(file.choose(), sep=',', header = T)
}

# Rename column "y" to "gdp"
colnames(dframe)[2] <- "gdp"

# Transform into time-series object
mymacro <- ts(dframe[,"gdp"],
              start = c(1979,2),
              end = c(2023,3),
              frequency = 4)

# Plot the gdp data
ts.plot(mymacro, 
        main = "Log of GDP",      
        col  = "black",        # Define colours
        lty  = "solid",        # Define line types
        lwd  = 2,              # Define line width
        ylab = "Log GDP", # Label of y-axis
        xlab = "Year"       # Label of x-axis
)

gdp <- mymacro[]  #select gdp series
gdp.lag4 <- xts::lag.xts(gdp, k = 4)  #define 4th lag
gdp.yoy <- (gdp - gdp.lag4)/gdp.lag4*100 #compute growth rate
gdp.yoy <- window(gdp.yoy, start = c(1980,2))

# Plot the gdp.yoy data
ts.plot(gdp.yoy, 
        main = "Y-o-Y GDP Growth",      
        col  = "black",        # Define colours
        lty  = "solid",        # Define line types
        lwd  = 2,              # Define line width
        ylab = "Growth Rate (%)", # Label of y-axis
        xlab = "Year"       # Label of x-axis
)

###------------------------------------------------------------------###
###---------------Creating train and test samples--------------------###
###------------------------------------------------------------------###

# Keep last 4 periods of observation and test sample
test <- window(gdp.yoy, start = c(2022,4), end = c(2023,3))

# Take out last observation from training sample
train <- window(gdp.yoy, end = c(2022,3))

###------------------------------------------------------------------###
###---------------Estimate & forecast AR(1) model--------------------###
###------------------------------------------------------------------###


# Set seed
set.seed(123)

# Estimating the AR(1) model with the arima function
est_ar1 <- arima(train, order = c(1,0,0), method = 'CSS')

summary(est_ar1) # Getting the coefficients


# 4-steps ahead point forecast
fc_ar1 <- forecast(est_ar1,
                   h = 4,
                   level = seq(5, 95, 10))

# Plot the forecast fan chart of AR(1) model
plot(fc_ar1,
     main = "Forecast fan chart of the AR(1) model: Y-o-Y GDP Growth",
     lwd = 2,
     ylab = 'Y-o-Y Growth (in %)',
     xlab = 'Periods',
     showgap = F,
     fcol = "red",
     flty = 2)

# Evaluation stats for the AR(1) model

loglik_ar1 <- est_ar1[["loglik"]] # Getting the loglikelihood value

# Improving the forecast error calculation section
# Extract all predictions and confidence intervals
pred_ar1 <- predict(est_ar1, n.ahead = 4)

# Calculate forecast errors (more clearly defined)
fe_ar1 <- test - pred_ar1$pred

# Deriving the squared forecast errors 
sfe_ar1 <- fe_ar1^2

# Calculate Root Mean Squared Forecast Error (RMSFE)
rmsfe_ar1 <- sqrt(mean(sfe_ar1))

# Calculate Mean Absolute Forecast Error (MAFE)
mafe_ar1 <- mean(abs(fe_ar1))

# Create a data frame to compare actual vs. predicted values
forecast_comparison <- data.frame(
  Period = time(test),
  Actual = as.numeric(test),
  Predicted = as.numeric(pred_ar1$pred),
  Error = as.numeric(fe_ar1),
  Squared_Error = as.numeric(sfe_ar1)
)

# Display forecast comparison
print(forecast_comparison)

# Ljung-Box test for autocorrelation
resid.ar1 <- residuals(est_ar1)
lb_result <- Box.test(resid.ar1, lag = 10, type = "Ljung-Box")
lb_result

# Shapiro-Wilk test for normality
shapiro_result <- shapiro.test(resid.ar1)
shapiro_result

# Breush-Pagan test for heteroscedasticity
# Create df of resid and fitted values
bp.test.data <- data.frame(
  residuals = resid.ar1,
  fitted = fitted(est_ar1),
  lag_residuals = c(NA, resid.ar1[-length(resid.ar1)])
)
# Remove NA values from lagged resids
bp.test.data <- na.omit(bp.test.data)
# Run test
bp_result <- bptest(residuals ~ fitted + lag_residuals,
                    data = bp.test.data)
bp_result

# Calculate the AIC
aic_ar1 <- (2*2)-(2*loglik_ar1)
aic_ar1

#Calculate the BIC
bic_ar1 <- (-2 * loglik_ar1) + (2*log(length(train)))
bic_ar1

# Enhanced evaluation metrics including RMSFE
eval_ar1 <- c("Log-likelihood" = loglik_ar1,
              "AIC" = aic_ar1,
              "BIC" = bic_ar1,
              "RMSFE" = rmsfe_ar1,
              "MAFE" = mafe_ar1,
              "Ljung-Box Test" = lb_result$statistic,
              'Shapiro-Wilk Test' = shapiro_result$statistic,
              'Breusch-Pagan Test' = bp_result$statistic)
eval_ar1

# Visualize forecast errors
plot(forecast_comparison$Period, forecast_comparison$Error, 
     type = "o", col = "blue", 
     main = "Forecast Errors for AR(1) Model",
     xlab = "Period", ylab = "Forecast Error")
abline(h = 0, lty = 2)

###------------------------------------------------------------------###
###---------------Estimating 2-state MS-AR(1) model------------------###
###------------------------------------------------------------------###

# Data set = gdp.yoy
# Set seed
set.seed(123)

# Create train/test split of dataset
# Keep last 4 periods of observation and test sample
test.ms <- window(gdp.yoy, start = c(2022,4), end = c(2023,3))

# Take out last observation from training sample
train.ms <- window(gdp.yoy, end = c(2022,3))

# define state and number of AR lags
k = 2 # number of states
p = 1 # number of AR lags

# create MS-AR(1) model

two_ms <- msmFit(lm(train.ms~1),
                 k = k,
                 sw = c(TRUE, TRUE, TRUE),
                 p = p)
summary(two_ms)

# Obtaining transition probabilities
transition_matrix <- two_ms@transMat
print("Transition Probability Matrix")
print(transition_matrix)

# create smoothed probability chart for both regimes

probr1 <- ts(two_ms@Fit@smoProb, end = c(2022,2), frequency = 4)

plot.ts(probr1[,1],
        main = 'Smoothed probabilities for both regimes',
        col = 'blue',
        xlab = '',
        ylab = '',
        lwd = 1.5)
lines(probr1[,2], col = 'red')
legend('left', legend = c('Regime 1','Regime 2'),
       lwd = c(2,1), col = c('blue','red'), bty = 'n', y.intersp = 1.5)


# Calculate expected duration of each regime
duration_regime1 <- 1 / (1 - transition_matrix[1,1])
duration_regime2 <- 1 / (1 - transition_matrix[2,2])
print(paste("Expected duration of regime 1 (expansion):", 
            round(duration_regime1, 1), "quarters"))
print(paste("Expected duration of regime 2 (recession):", 
            round(duration_regime2, 1), "quarters"))

###------------------------------------------------------------------###
###--------------Forecasting 2-state MS-AR(1) model------------------###
###------------------------------------------------------------------###

# Extract regime-specific parameters
cat("\nRegime-specific parameters:\n")
print(two_ms@Coef)

# Extract transition probabilities
cat("\nTransition probabilities:\n")
print(two_ms@transMat)

# -------------------------------------------------------------------------
# MS(2)-AR(1) Forecasting
# -------------------------------------------------------------------------

# We need to implement a forecasting function for MS-AR model
# Since MSwM doesn't have a built-in forecast function, we'll create our own

# Function to generate forecasts from MS-AR model
forecast_msar <- function(model, h=4) {
  # Debug: Print model structure
  cat("Model class:", class(model), "\n")
  cat("Available slots:", paste(slotNames(model), collapse=", "), "\n")
  
  # For MSM.lm models, data is in model@model
  model_data <- model@model
  cat("Model$model structure:", class(model_data), "\n")
  cat("Model$model components:", paste(names(model_data), collapse=", "), "\n")
  
  # Most likely, the data is in model@model$model
  if (is.null(model_data$model)) {
    stop("Cannot find data in model object. Please inspect the model structure.")
  }
  
  data_frame <- model_data$model
  cat("Data columns:", paste(names(data_frame), collapse=", "), "\n")
  
  # Get the response variable (likely the first column)
  response_var <- names(data_frame)[1]
  last_obs <- tail(data_frame[[response_var]], 1)
  cat("Last observation:", last_obs, "\n")
  
  # Get coefficient structure
  cat("Coefficient structure:\n")
  print(model@Coef)
  
  # Get transition matrix
  cat("Transition matrix:\n")
  print(model@transMat)
  
  # Get filtered probabilities
  cat("Smooth probabilities structure:\n")
  print(head(model@Fit@filtProb))
  cat("Last regime probabilities:", tail(model@Fit@filtProb, 1), "\n")
  
  
  # Get the original data from the model
  # For MSM.lm objects, the data is in model@model$model
  original_data <- model@model$model
  response_var <- names(original_data)[1]  # First column is response variable
  last_obs <- tail(original_data[[response_var]], 1)
  
  # Get parameters for each regime
  coef1 <- model@Coef[1,]  # Regime 1 coefficients
  coef2 <- model@Coef[2,]  # Regime 2 coefficients
  
  # Extract transition probabilities
  trans_mat <- model@transMat
  
  # Get filtered probabilities to determine current regime
  # In MSM.lm objects, it's stored in Fit@filtProb
  smooth_probs <- model@Fit@filtProb
  current_regime_prob <- tail(smooth_probs, 1)
  
  # Initialize forecast vectors
  forecasts <- numeric(h)
  
  # Generate forecasts
  prev_value <- last_obs
  
  for (i in 1:h) {
    # For AR(1) model with intercept
    # Check the actual names of coefficients
    int_name <- "(Intercept)"
    if (!(int_name %in% names(coef1))) {
      int_name <- names(coef1)[1]  # Use first coefficient as intercept
    }
    
    slope_name <- response_var
    if (!(slope_name %in% names(coef1))) {
      slope_name <- names(coef1)[2]  # Use second coefficient as slope
    }
    
    # Calculate forecasts for each regime
    forecast_r1 <- coef1[int_name]
    if (length(coef1) > 1) forecast_r1 <- forecast_r1 + coef1[slope_name] * prev_value
    
    forecast_r2 <- coef2[int_name]
    if (length(coef2) > 1) forecast_r2 <- forecast_r2 + coef2[slope_name] * prev_value
    
    # Weight forecasts by probability of being in each regime
    p_regime1 <- current_regime_prob[1] * trans_mat[1,1] + 
      current_regime_prob[2] * trans_mat[2,1]
    p_regime2 <- current_regime_prob[1] * trans_mat[1,2] + 
      current_regime_prob[2] * trans_mat[2,2]
    
    # Weighted forecast
    forecasts[i] <- p_regime1 * forecast_r1 + p_regime2 * forecast_r2
    
    # Update for next iteration
    prev_value <- forecasts[i]
    current_regime_prob <- c(p_regime1, p_regime2)
  }
  
  return(forecasts)
}

# Generate forecasts for the hold-out period (4 quarters of 2024)
fc_msar <- forecast_msar(two_ms, h=4)

# Create a data frame with forecast values
fc_msar_df <- data.frame(
  Point_Forecast = as.numeric(fc_msar)
)

# Print MS-AR forecasts
cat("\nMS-AR Model Forecasts for 2024:\n")
print(fc_msar)

# -------------------------------------------------------------------------
# MS(2)-AR(1) Forecast Evaluation
# -------------------------------------------------------------------------
# Check the structure of both objects
cat("Structure of test.ms:\n")
print(class(test.ms))
print(str(test.ms))

cat("\nStructure of fc_msar:\n")
print(class(fc_msar))
print(str(fc_msar))

# Convert forecasts to a time series object with same time attributes as test.ms
fc_msar_ts <- ts(fc_msar, 
                 start = time(test.ms)[1],
                 frequency = frequency(test.ms))

# Or alternatively, convert both to numeric vectors for comparison
test.ms_numeric <- as.numeric(test.ms)
fc_msar_numeric <- as.numeric(fc_msar)

# Calculate forecast errors using the compatible objects
msar_errors <- test.ms_numeric - fc_msar_numeric
# OR
msar_errors <- test.ms - fc_msar_ts

# Calculate RMSFE
msar_rmsfe <- sqrt(mean(msar_errors^2, na.rm = TRUE))

# Calculate MAFE
msar_mafe <- mean(abs(msar_errors), na.rm = TRUE)

##################################################################
# Calculate forecast errors
msar_errors <- test.ms - fc_msar

# Calculate RMSFE
msar_rmsfe <- sqrt(mean(msar_errors^2, na.rm = TRUE))

# Calculate MAFE
msar_mafe <- mean(abs(msar_errors), na.rm = TRUE)

# Add MS-AR results to the forecast accuracy data frame
forecast_accuracy <-data.frame(Model = "MS(2)-AR(1) Model",
                               RMSFE = msar_rmsfe,
                               MAFE = msar_mafe)
print(forecast_accuracy)


###------------------------------------------------------------------###
###---------------Estimating 3-state MS-AR(1) model------------------###
###------------------------------------------------------------------###

# Data set = gdp.yoy
# Set seed
set.seed(123)

# Create train/test split of dataset
# Keep last 4 periods of observation and test sample
test.ms <- window(gdp.yoy, start = c(2022,4), end = c(2023,3))

# Take out last observation from training sample
train.ms <- window(gdp.yoy, end = c(2022,3))

# define state and number of AR lags
k = 3 # number of states
p = 1 # number of AR lags

# create MS-AR(1) model

three_ms <- msmFit(lm(train.ms~1),
                 k = k,
                 sw = c(TRUE, TRUE, TRUE),
                 p = p)
summary(three_ms)

# Obtaining transition probabilities
transition_matrix <- three_ms@transMat
print("Transition Probability Matrix")
print(transition_matrix)

# create smoothed probability chart for both regimes

probr1 <- ts(three_ms@Fit@smoProb, end = c(2022,2), frequency = 4)
par(xpd=TRUE)
plot.ts(probr1[,1],
        main = 'Smoothed probabilities for all regimes',
        col = 'blue',
        xlab = '',
        ylab = '',
        lwd = 2)
lines(probr1[,2], col = 'red', lwd = 2)
lines(probr1[,3], col = 'black', lwd = 2)
legend('topright', legend = c('Regime 1','Regime 2', 'Regime 3'),
       lwd = c(2,2,2), col = c('blue','red','black'), bty = 'n', y.intersp = 1.5)

# Calculate expected duration of each regime
duration_regime1 <- 1 / (1 - transition_matrix[1,1])
duration_regime2 <- 1 / (1 - transition_matrix[2,2])
duration_regime3 <- 1 / (1 - transition_matrix[3,3])
print(paste("Expected duration of regime 1 (stable):", 
            round(duration_regime1, 1), "quarters"))
print(paste("Expected duration of regime 2 (recession):", 
            round(duration_regime2, 1), "quarters"))
print(paste("Expected duration of regime 3 (growth):", 
            round(duration_regime3, 1), "quarters"))

# -------------------------------------------------------------------------
# MS(3)-AR(1) Forecasting
# -------------------------------------------------------------------------

forecast_msar <- function(model, h=4) {
  # Debug: Print model structure
  cat("Model class:", class(model), "\n")
  cat("Available slots:", paste(slotNames(model), collapse=", "), "\n")
  
  # For MSM.lm models, data is in model@model
  model_data <- model@model
  cat("Model$model structure:", class(model_data), "\n")
  cat("Model$model components:", paste(names(model_data), collapse=", "), "\n")
  
  # Get the data from the model
  if (is.null(model_data$model)) {
    stop("Cannot find data in model object. Please inspect the model structure.")
  }
  
  data_frame <- model_data$model
  cat("Data columns:", paste(names(data_frame), collapse=", "), "\n")
  
  # Get the response variable (likely the first column)
  response_var <- names(data_frame)[1]
  last_obs <- tail(data_frame[[response_var]], 1)
  cat("Last observation:", last_obs, "\n")
  
  # Get coefficient structure
  cat("Coefficient structure:\n")
  print(model@Coef)
  
  # Get number of regimes (for MS(3)-AR(1), k=3)
  k <- nrow(model@Coef)
  cat("Number of regimes:", k, "\n")
  
  # Get transition matrix
  cat("Transition matrix:\n")
  print(model@transMat)
  
  # Get filtered probabilities
  cat("Last regime probabilities:", tail(model@Fit@filtProb, 1), "\n")
  
  # Get parameters for each regime - now handling any number of regimes
  coef_list <- list()
  for (i in 1:k) {
    coef_list[[i]] <- model@Coef[i,]
  }
  
  # Extract transition probabilities
  trans_mat <- model@transMat
  
  # Get filtered probabilities to determine current regime
  smooth_probs <- model@Fit@filtProb
  current_regime_prob <- as.numeric(tail(smooth_probs, 1))  # Ensure numeric
  cat("Current regime probabilities (numeric):", current_regime_prob, "\n")
  
  # Initialize forecast vectors
  forecasts <- numeric(h)
  
  # Generate forecasts
  prev_value <- as.numeric(last_obs)  # Ensure numeric
  cat("Initial previous value (numeric):", prev_value, "\n")
  
  for (i in 1:h) {
    # For AR(1) model with intercept
    # Check the actual names of coefficients
    int_name <- "(Intercept)"
    if (!(int_name %in% names(coef_list[[1]]))) {
      int_name <- names(coef_list[[1]])[1]  # Use first coefficient as intercept
    }
    
    slope_name <- response_var
    if (length(coef_list[[1]]) > 1 && !(slope_name %in% names(coef_list[[1]]))) {
      slope_name <- names(coef_list[[1]])[2]  # Use second coefficient as slope
    }
    
    cat("\nIteration", i, "- Using coefficient names:", int_name, "and", slope_name, "\n")
    
    # Calculate forecasts for each regime
    regime_forecasts <- numeric(k)
    for (j in 1:k) {
      regime_forecasts[j] <- as.numeric(coef_list[[j]][int_name])
      if (length(coef_list[[j]]) > 1) {
        regime_forecasts[j] <- regime_forecasts[j] + as.numeric(coef_list[[j]][slope_name]) * prev_value
      }
      cat("Regime", j, "forecast:", regime_forecasts[j], "\n")
    }
    
    # Calculate transition probabilities for each regime
    new_regime_probs <- numeric(k)
    for (j in 1:k) {
      new_regime_probs[j] <- sum(current_regime_prob * trans_mat[, j])
      cat("New probability for regime", j, ":", new_regime_probs[j], "\n")
    }
    
    # Debug: Check if vectors are numeric before multiplication
    cat("regime_forecasts is numeric:", is.numeric(regime_forecasts), "\n")
    cat("new_regime_probs is numeric:", is.numeric(new_regime_probs), "\n")
    
    # Weighted forecast across all regimes - element-wise multiplication then sum
    forecasts[i] <- sum(new_regime_probs * regime_forecasts)
    cat("Weighted forecast for step", i, ":", forecasts[i], "\n")
    
    # Update for next iteration
    prev_value <- forecasts[i]
    current_regime_prob <- new_regime_probs
  }
  
  return(forecasts)
}

# Generate forecasts for the hold-out period (4 quarters of 2024)
fc_msar <- forecast_msar(three_ms, h=4)

# Create a data frame with forecast values
fc_msar_df <- data.frame(
  Point_Forecast = as.numeric(fc_msar)
)

# Print MS-AR forecasts
cat("\nMS-AR Model Forecasts for 2024:\n")
print(fc_msar)

# -------------------------------------------------------------------------
# MS(3)-AR(1) Forecast Evaluation
# -------------------------------------------------------------------------
# Check the structure of both objects
cat("Structure of test.ms:\n")
print(class(test.ms))
print(str(test.ms))

cat("\nStructure of fc_msar:\n")
print(class(fc_msar))
print(str(fc_msar))

# Convert forecasts to a time series object with same time attributes as test.ms
fc_msar_ts <- ts(fc_msar, 
                 start = time(test.ms)[1],
                 frequency = frequency(test.ms))

# Or alternatively, convert both to numeric vectors for comparison
test.ms_numeric <- as.numeric(test.ms)
fc_msar_numeric <- as.numeric(fc_msar)

# Calculate forecast errors using the compatible objects
msar_errors <- test.ms_numeric - fc_msar_numeric
# OR
msar_errors <- test.ms - fc_msar_ts

# Calculate RMSFE
msar_rmsfe <- sqrt(mean(msar_errors^2, na.rm = TRUE))

# Calculate MAFE
msar_mafe <- mean(abs(msar_errors), na.rm = TRUE)

##################################################################
# Calculate forecast errors
msar_errors <- test.ms - fc_msar

# Calculate RMSFE
msar_rmsfe <- sqrt(mean(msar_errors^2, na.rm = TRUE))

# Calculate MAFE
msar_mafe <- mean(abs(msar_errors), na.rm = TRUE)

# Add MS-AR results to the forecast accuracy data frame
forecast_accuracy <-data.frame(Model = "MS(3)-AR(1) Model",
                               RMSFE = msar_rmsfe,
                               MAFE = msar_mafe)
print(forecast_accuracy)
