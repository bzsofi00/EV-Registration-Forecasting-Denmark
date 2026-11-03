# Load necessary libraries
library(fpp3)
library(readr)
library(urca)
library(tseries)
library(forecast)
library(fable)
#install.packages("strucchange")
library(strucchange)
library(purrr)

# Load the cleaned data
cars_data <- read_csv("cleaned_data.csv")

# Filter for electric vehicle sales, convert time column and store as tsibble
ev_data <- cars_data %>%
  filter(Propellant == "Electricity") %>%
  mutate(Date = yearmonth(ymd(Date))) %>%
  as_tsibble(index = Date)

# Plot the EV sales over time
plot_ev_sales <- ev_data %>%
  autoplot(Registrations) +
  labs(title = "Monthly New Electric Vehicle Registrations in Denmark (2011–2025)",
       x = "Date", y = "Number of Registrations") +
  theme_minimal()
plot_ev_sales

# seasonal plot
ev_data %>%
  gg_season(Registrations, period = "year") +
  labs(title = "Seasonal Plot of Electric Car Registrations",
       x = "Month", y = "Registrations") +
  theme_minimal()

# Plot the autocorrelation function
plot_acf_ev <- ev_data %>%
  ACF(Registrations) %>%
  autoplot() +
  labs(title = "ACF of Electric Car Registrations",
       x = "Lag", y = "ACF") +
  theme_minimal()
plot_acf_ev

# STL decomposition
ev_data %>%
  model(STL(Registrations ~ season(window = "periodic"))) %>%
  components() %>%
  autoplot() +
  labs(title = "STL Decomposition of EV Registrations") +
  theme_minimal()

# Log transformation
log_ev_data <- ev_data %>%
  mutate(Registrations_log = log(Registrations)) %>%
  select(-Registrations)

# Plot log-transformed data
plot_log_ev <- log_ev_data %>%
  autoplot(Registrations_log) +
  labs(title = "Log-Transformed Electric Car Registrations",
       x = "Date", y = "Log of Registrations") +
  theme_minimal()
plot_log_ev

# Box-Cox transformation: determine lambda
lambda <- ev_data %>%
  features(Registrations, features = guerrero) %>%
  pull(lambda_guerrero)
lambda

# Use lambda found above
bc_ev_data <- ev_data %>%
  mutate(Registrations_bc = box_cox(Registrations, lambda)) %>%
  select(-Registrations)

bc_ev_data %>%
  autoplot(Registrations_bc) +
  labs(title = "Box-Cox Transformed Electric Car Sales", x = "Date", y = "Box-Cox of Registrations") +
  theme_minimal()

# Plot the autocorrelation function of Box-Cox transformed data
plot_acf_bc_ev <- bc_ev_data %>%
  ACF(Registrations_bc) %>%
  autoplot() +
  labs(title = "ACF of Box-Cox Transformed EV Registrations",
       x = "Lag", y = "ACF") +
  theme_minimal()
plot_acf_bc_ev

# STL decomposition of Box-Cox transformed series
bc_ev_data %>%
  model(STL(Registrations_bc ~ season(window = "periodic"))) %>%
  components() %>%
  autoplot() +
  labs(title = "STL Decomposition of Box-Cox Transformed EV Sales") +
  theme_minimal()

# ADF test
summary(ur.df(bc_ev_data$Registrations_bc, type = "trend"))
# test statistic -4.3202 < critical value -3.99, so we reject the null hypothesis of a unit root at the 1% level
# The box-cox transformed series is stationary

# KPSS test
summary(ur.kpss(bc_ev_data$Registrations_bc, type = "tau"))
# test statistic 0.5606 > critical values at all levels , so we reject the null hypothesis of stationarity, meaning the series is not stationary

-----------------------------------
# STRUCTURAL BREAK DETECTION: QLR (supF) Test

# Convert to ts object
start_date <- yearmonth(min(bc_ev_data$Date))
ev_ts <- ts(bc_ev_data$Registrations_bc, 
            start = c(year(start_date), month(start_date)), 
            frequency = 12)

# Compute rolling F-statistics across all eligible break points (excluding first and last 15% of observations)
f_stats <- Fstats(ev_ts ~ 1, from = 0.15)
# Compute rolling F-statistics allowing for a trend
f_stats_trend <- Fstats(ev_ts ~ time(ev_ts), from = 0.15)

# supF test for structural break at unknown date
supf_result <- sctest(f_stats, type = "supF")
print(supf_result)
# The p-value is extremely small (essentially 0), so we reject the null hypothesis of no structural break
# supF test for structural break in the trend
supf_result_trend <- sctest(f_stats_trend, type = "supF")
print(supf_result_trend)

# Estimate most likely breakpoints
breakpoints_estimates <- breakpoints(f_stats, breaks = 3, alpha = 0.05)
summary(breakpoints_estimates)
# Optimal 2-segment partition, Breakpoints at observation number 115, which corresponds to July 2020 (2020(7))
breakpoints_estimates_trend <- breakpoints(f_stats_trend, breaks = 3, alpha = 0.05)

# Plot the F-statistics and highlight the breakpoints
plot(f_stats, alpha = 0.1, main = "F Statistics")
lines(breakpoints_estimates)

# Access the breakpoints from the breakpoints_estimates object
break_index <- breakpoints_estimates$breakpoints[1]
break_date <- bc_ev_data$Date[break_index]
print(break_date)

# Split the data into two segments based on the breakpoint
bc_ev_post_break <- bc_ev_data %>%
  filter(Date > break_date)

bc_ev_pre_break <- bc_ev_data %>%
  filter(Date <= break_date)
-----------------------------
# Differencing the whole series

# Check how many times we need to difference the whole series
bc_ev_data %>%
  features(Registrations_bc, unitroot_ndiffs)

# Check whether seasonal differencing is required
bc_ev_data %>%
  features(Registrations_bc, unitroot_nsdiffs)

# Difference the whole series
bc_ev_data_diff <- bc_ev_data %>%
  mutate(Registrations_bc_diff = difference(Registrations_bc)) %>%
  filter(!is.na(Registrations_bc_diff)) %>%
  select(-Registrations_bc)

# Plot the differenced whole series
plot_diff_ev <- bc_ev_data_diff %>%
  autoplot(Registrations_bc_diff) +
  labs(title = "Differenced Box-Cox Transformed Electric Car Registrations",
       x = "Date", y = "Differenced Box-Cox of Registrations") +
  theme_minimal()
plot_diff_ev

# Check for stationarity of the differenced whole series
# ADF test  
summary(ur.df(bc_ev_data_diff$Registrations_bc_diff, type = "trend"))
# test statistic -17.056 < critical value -3.99, so we reject the null hypothesis of a unit root at the 1% level
# The differenced series is stationary

# KPSS test
summary(ur.kpss(bc_ev_data_diff$Registrations_bc_diff, type = "tau"))
# test statistic 0.022 < critical values at all levels , so we fail to reject the null hypothesis of no unit root
# the differenced series is stationary

# ACF and PACF of the differenced series
plot_acf_diff_ev <- bc_ev_data_diff %>%
  ACF(Registrations_bc_diff) %>%
  autoplot() +
  labs(title = "ACF of Differenced Box-Cox Transformed EV Registrations",
       x = "Lag", y = "ACF") +
  theme_minimal()
plot_acf_diff_ev

plot_pacf_diff_ev <- bc_ev_data_diff %>%
  PACF(Registrations_bc_diff) %>%
  autoplot() +
  labs(title = "PACF of Differenced Box-Cox Transformed EV Registrations",
       x = "Lag", y = "PACF") +
  theme_minimal()
plot_pacf_diff_ev

# Identify an ARIMA model
# ARIMA(2,1,3)(1,0,1)[12]

# Check how many times we need to difference the post-break series
bc_ev_post_break %>%
  features(Registrations_bc, unitroot_ndiffs)
# 1 difference is needed
# Check whether seasonal differencing is required
bc_ev_post_break %>%
  features(Registrations_bc, unitroot_nsdiffs)
# 1 seasonal difference is needed

# Difference the post-break series
bc_ev_post_break_diff <- bc_ev_post_break %>%
  mutate(
    # seasonal differencing
    Registrations_bc_seasonal = difference(Registrations_bc, lag = 12),  
    
    # first-order differencing
    Registrations_bc_diff = difference(Registrations_bc_seasonal)  
  ) %>%
  filter(!is.na(Registrations_bc_diff)) %>%
  select(-Registrations_bc)  # Remove original and seasonally differenced columns

# Plot the differenced post-break series
plot_bc_ev_post_break_diff <- bc_ev_post_break_diff %>%
  autoplot(Registrations_bc_diff) +
  labs(title = "Differenced Box-Cox Transformed Electric Car Registrations (Post-Break)",
       x = "Date", y = "Differenced Box-Cox of Registrations") +
  theme_minimal()
plot_bc_ev_post_break_diff

# Check for stationarity of the differenced post-break series
# ADF test
summary(ur.df(bc_ev_post_break_diff$Registrations_bc_diff, type = "trend"))
# test statistic -7.8769 < critical value at all significance levels, so we reject the null hypothesis of a unit root
# The differenced post-break series is stationary

# KPSS test
summary(ur.kpss(bc_ev_post_break_diff$Registrations_bc_diff, type = "tau"))
# test statistic 0.0622 < critical values at all significance levels, so we fail to reject the null hypothesis of no unit root
# The differenced post-break series is stationary

# ACF and PACF of the differenced post-break series
plot_acf_post_break_diff <- bc_ev_post_break_diff %>%
  ACF(Registrations_bc_diff) %>%
  autoplot() +
  labs(title = "ACF of Differenced Box-Cox Transformed EV Registrations (Post-Break)",
       x = "Lag", y = "ACF") +
  theme_minimal()
plot_acf_post_break_diff

plot_pacf_post_break_diff <- bc_ev_post_break_diff %>%
  PACF(Registrations_bc_diff) %>%
  autoplot() +
  labs(title = "PACF of Differenced Box-Cox Transformed EV Registrations (Post-Break)",
       x = "Lag", y = "PACF") +
  theme_minimal()
plot_pacf_post_break_diff

# Identify an ARIMA model for the post-break series
# ARIMA(2, 1, 1)(0, 1, 0)[12]

# reapply the QLR test on the differenced data
# Convert to ts object 
start_date_diff <- yearmonth(min(bc_ev_data_diff$Date))
ev_ts_diff <- ts(bc_ev_data_diff$Registrations_bc_diff, 
                 start = c(year(start_date_diff), month(start_date_diff)), 
                 frequency = 12)

# Compute rolling F-statistics across all eligible break points (excluding first and last 15% of observations)
f_stats_diff <- Fstats(ev_ts_diff ~ 1, from = 0.15)
f_stats_diff2 <- Fstats(ev_ts_diff ~ stats::lag(ev_ts_diff, 1), from = 0.15)

# supF test for structural break at unknown date
supf_result_diff <- sctest(f_stats_diff, type = "supF")
print(supf_result_diff)
# p-value is 1, so we fail to reject the null hypothesis of no structural break
# After differencing the data, the series does not exhibit any significant structural breaks in its stationary component

supf_result_diff2 <- sctest(f_stats_diff2, type = "supF")
print(supf_result_diff2)
----------------------------------------------------------
  # create a dummy variable for the break
  bc_ev_data <- bc_ev_data %>%
  mutate(
    # Create a dummy variable for the break
    Break = if_else(Date > break_date, 1, 0), # level shift
    # Create a dummy variable for the trend change
    Trend = if_else(Date > break_date, Date - break_date, 0) # trend shift
  ) %>%
  as_tsibble(index = Date)

# fit the model with only level shift
model_level_shift <- bc_ev_data %>%
  model(ARIMA(Registrations_bc ~ Break)) %>%
  report()
# fit the model with level and trend shift
model_level_trend_shift <- bc_ev_data %>%
  model(ARIMA(Registrations_bc ~ Break + Trend)) %>%
  report()
# The trend dummy is not significant, Model fit statistics (AIC/BIC) are actually worse when it’s added
# This aligns with earlier supF test, which also found no significant structural break after differencing

# compare models
model_comparison <- bc_ev_data %>%
  model(
    level_shift = ARIMA(Registrations_bc ~ Break),
    level_trend_shift = ARIMA(Registrations_bc ~ Break + Trend)
  ) %>%
  glance() %>%
  select(.model, AICc, BIC) %>%
  arrange(AICc)

# Print the model comparison
print(model_comparison)

# fit multiple ARIMA models with the dummy variable
# Fit ARIMA models with the dummy variable
fit_with_break <- bc_ev_data %>%
  model(
    arima_b1 = ARIMA(Registrations_bc ~ Break + 1 + pdq(1, 1, 2) + PDQ(1, 0, 1, 12)),
    arima_b2 = ARIMA(Registrations_bc ~ Break + 1 + pdq(1, 1, 3) + PDQ(1, 0, 1, 12)),
    arima_b3 = ARIMA(Registrations_bc ~ Break + 1 + pdq(2, 1, 2) + PDQ(1, 0, 1, 12)),
    arima_b4 = ARIMA(Registrations_bc ~ Break + 1 + pdq(2, 1, 3) + PDQ(1, 0, 1, 12)),
    
    arima_b5 = ARIMA(Registrations_bc ~ Break + 0 + pdq(1, 1, 2) + PDQ(1, 0, 1, 12)),
    arima_b6 = ARIMA(Registrations_bc ~ Break + 0 + pdq(1, 1, 3) + PDQ(1, 0, 1, 12)),
    arima_b7 = ARIMA(Registrations_bc ~ Break + 0 + pdq(2, 1, 2) + PDQ(1, 0, 1, 12)),
    arima_b8 = ARIMA(Registrations_bc ~ Break + 0 + pdq(2, 1, 3) + PDQ(1, 0, 1, 12))
  )

# 2. Then get a summary separately
fit_with_break_summary <- fit_with_break %>%
  glance() %>%
  arrange(AICc) %>%
  select(.model, AICc, BIC)
# Print the summary
print(fit_with_break_summary)

dof_map_break <- c(
  arima_b1 = 1 + 2 + 1 + 1 + 1 + 1,  # p + q + P + Q + Break + intercept
  arima_b2 = 1 + 3 + 1 + 1 + 1 + 1,
  arima_b3 = 2 + 2 + 1 + 1 + 1 + 1,
  arima_b4 = 2 + 3 + 1 + 1 + 1 + 1,
  
  arima_b5 = 1 + 2 + 1 + 1 + 1,      # no intercept
  arima_b6 = 1 + 3 + 1 + 1 + 1,
  arima_b7 = 2 + 2 + 1 + 1 + 1,
  arima_b8 = 2 + 3 + 1 + 1 + 1
)


ljung_break_results <- purrr::imap_dfr(dof_map_break, function(dof, model_name) {
  fit_with_break %>%
    select(!!model_name) %>%
    augment() %>%
    features(.innov, ljung_box, lag = 24, dof = dof) %>%
    mutate(model = model_name)
})

ljung_break_results %>%
  select(model, lb_stat, lb_pvalue) %>%
  arrange(desc(lb_pvalue))

best_arima_break <- fit_with_break %>%
  select(arima_b5)
report(best_arima_break)

# check the residuals
best_arima_break %>% gg_tsresiduals()

# Check the residuals for autocorrelation
augment(best_arima_break) %>%
  features(.innov, ljung_box, lag = 24, dof = 6)
# p-value 0.0661 > 0.05, so we can accept the null hypothesis of no autocorrelation

final_model <- bc_ev_data %>%
  model(arima_b5 = ARIMA(Registrations_bc ~ Break + 0 + pdq(1, 1, 2) + PDQ(1, 0, 1, 12)))
report(final_model)

-------------------
  # Fit ARIMA on post-break data based on ACF and PACF plots
  # ARIMA(2, 1, 1)(0, 1, 0)[12]
  manual_post_break_arima <- bc_ev_post_break %>%
  model(ARIMA(Registrations_bc ~ 0 + pdq(2, 1, 1) + PDQ(0, 1, 0, period=12)))
report(manual_post_break_arima)

# Check the residuals
manual_post_break_arima %>%
  gg_tsresiduals() +
  labs(title = "Residuals of Manual ARIMA Model on Post-Break Data")

# Check the residuals for autocorrelation
augment(manual_post_break_arima) %>%
  features(.innov, ljung_box, lag = 24, dof = 3)
# p-value 0.973 > 0.05, so we fail to reject the null hypothesis of no autocorrelation
# there is no autocorrelation in the residuals

# Fit ARIMA on post-break data automatically
auto_post_break_arima <- bc_ev_post_break %>%
  model(ARIMA(Registrations_bc))
report(auto_post_break_arima)
# Model: ARIMA(0,0,1)(1,1,0)[12] w/ drift 

# Check the residuals
auto_post_break_arima %>%
  gg_tsresiduals() +
  labs(title = "Residuals of Automatic ARIMA Model on Post-Break Data")

# Check the residuals for autocorrelation
augment(auto_post_break_arima) %>%
  features(.innov, ljung_box, lag = 24, dof = 3)
# p-value 0.395 > 0.05, so we fail to reject the null hypothesis of no autocorrelation
# there is no autocorrelation in the residuals

# Explore different ARIMA models on post-break data
post_fit <- bc_ev_post_break %>%
  model(
    #arima0 = ARIMA(Registrations_bc ~ 0 + pdq(0, 0, 1) + PDQ(1, 1, 0, period=12)),
    arima1 = ARIMA(Registrations_bc ~ 1 + pdq(0, 0, 1) + PDQ(1, 1, 0, period=12)),
    #arima2 = ARIMA(Registrations_bc ~ 0 + pdq(0, 0, 1) + PDQ(0, 1, 0, period=12)),
    arima2 = ARIMA(Registrations_bc ~ 0 + pdq(0, 1, 1) + PDQ(0, 1, 0, period=12)),
    arima3 = ARIMA(Registrations_bc ~ 0 + pdq(1, 1, 1) + PDQ(0, 1, 0, period=12)),
    arima4 = ARIMA(Registrations_bc ~ 0 + pdq(1, 1, 1) + PDQ(1, 1, 0, period=12)),
    arima5 = ARIMA(Registrations_bc ~ 0 + pdq(0, 1, 2) + PDQ(0, 1, 0, period=12)),
    arima6 = ARIMA(Registrations_bc ~ 0 + pdq(1, 1, 2) + PDQ(0, 1, 0, period=12)),
    arima7 = ARIMA(Registrations_bc ~ 0 + pdq(2, 1, 0) + PDQ(0, 1, 0, period=12)),
    arima8 = ARIMA(Registrations_bc ~ 0 + pdq(2, 1, 0) + PDQ(1, 1, 0, period=12)),
    arima9 = ARIMA(Registrations_bc ~ 0 + pdq(2, 1, 1) + PDQ(0, 1, 0, period=12)),
    arima10 = ARIMA(Registrations_bc ~ 0 + pdq(2, 1, 2) + PDQ(0, 1, 0, period=12))
  )

post_fit %>%
  glance() %>%
  arrange(AICc) %>%
  select(.model, AICc, BIC)

# Check the autocorrelation of the residuals
# Define degrees of freedom for each model
dof_map_post <- c(
  #arima0 = 0 + 1 + 1 + 0,
  arima1 = 0 + 1 + 1 + 0 + 1,
  #arima3 = 0 + 1 + 0 + 0,
  arima2 = 0 + 1 + 0 + 0,
  arima3 = 1 + 1 + 0 + 0,
  arima4 = 1 + 1 + 1 + 0,
  arima5 = 0 + 2 + 0 + 0,
  arima6 = 1 + 2 + 0 + 0,
  arima7 = 2 + 0 + 0 + 0,
  arima8 = 2 + 0 + 1 + 0,
  arima9 = 2 + 1 + 0 + 0,
  arima10 = 2 + 2 + 0 + 0
)

ljung_post_results <- purrr::imap_dfr(dof_map_post, function(dof, model_name) {
  post_fit %>%
    select(!!model_name) %>%
    augment() %>%
    features(.innov, ljung_box, lag = 24, dof = dof) %>%
    mutate(model = model_name)
})

ljung_post_results %>%
  select(model, lb_stat, lb_pvalue) %>%
  arrange(desc(lb_pvalue))

# check best arima
best_arima_post <- post_fit %>%
  select(arima1) %>%
  report(best_arima_post)

# Check the accuracy of the models
accuracy(post_fit)

# Choose the best model
# The best model is arima1
# ARIMA(0, 0, 1)(1, 1, 0)[12] with drift
best_fit_post <- bc_ev_post_break %>%
  model(
    arima2 = ARIMA(Registrations_bc ~ 1 + pdq(0, 0, 1) + PDQ(1, 1, 0, period=12)))
best_fit_post %>% report()

# Check the residuals of the best model
best_fit_post %>%
  gg_tsresiduals() +
  labs(title = "Residuals of Best ARIMA Model (ARIMA(0, 0, 1)(1, 1, 0)[12]) on Post-Break Data")

# Check the residuals for autocorrelation
augment(best_fit_post) %>%
  features(.innov, ljung_box, lag = 24, dof = 3)
# p-value 0.395 > 0.05, so we fail to reject the null hypothesis of no autocorrelation
# there is no autocorrelation in the residuals

-------------------------
# Forecast on training data with full model + break dummy
# Split the data into training and testing sets
# Train: Jan 2011 – Dec 2022
# Test: Jan 2023 – Feb 2025
# 12-year training window and a 2-year test window

train_data <- bc_ev_data %>%
  filter(Date <= yearmonth("2022 Dec"))

test_data <- bc_ev_data %>%
  filter(Date >= yearmonth("2023 Jan"))

nrow(train_data)
nrow(test_data)

test_data %>% select(Date, Break)

# Fit the model on the training set
fit_train <- train_data %>%
  model(arima_b5 = ARIMA(Registrations_bc ~ Break + 0 + pdq(1,1,2) + PDQ(1,0,1,12)))
report(fit_train)

# give prediction intervals for your forecasts
fc <- fit_train %>% forecast(new_data = test_data, level = 95)

# Plot actual vs forecast
fc %>%
  autoplot(bc_ev_data) +
  labs(title = "Forecasts using ARIMA5(1,1,2)(1,0,1)[12] with Break dummy", y = "Registrations (Box-Cox)",
       x = "Date") +
  autolayer(test_data, Registrations_bc) +
  theme_minimal()

# compare accuracy
accuracy(fc, bc_ev_data)

# ME    = 1.77  Positive bias—forecasts are, on average, overestimating slightly
# RMSE  = 2.05  Penalizes larger errors more. Lower is better. This is in the Box-Cox scale, so its absolute meaning is less interpretable unless back-transformed
# MAE   = 1.88  On average, forecasts are off by about 1.88 units (again in Box-Cox space)
# MPE   = 10.9%  
# MAPE  = 11.7%  This is a relatively acceptable error margin for time series forecasting. Anything below 10–15% is often considered good in applied settings
# MASE  = 1.13  This is just above 1, meaning it’s slightly worse than a naive seasonal benchmark. You ideally want MASE < 1 for your model to outperform the naive baseline

# 2. Extract prediction intervals (convert hilo to columns)
fc_extracted <- fc %>%
  mutate(
    intervals = hilo(Registrations_bc, 95),
    lower = map_dbl(intervals, ~ .x$lower),
    upper = map_dbl(intervals, ~ .x$upper)
  )

# 3. Invert Box-Cox transformation manually
fc_inverted <- fc_extracted %>%
  mutate(
    mean_orig = inv_box_cox(.mean, lambda = lambda),
    lower_orig = inv_box_cox(lower, lambda = lambda),
    upper_orig = inv_box_cox(upper, lambda = lambda)
  )

# Historical original series
original_data <- bc_ev_data %>%
  mutate(Registrations_orig = inv_box_cox(Registrations_bc, lambda))

# Plot forecasts with historical data
fc_inverted %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data, aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Forecasts using ARIMA5(1,1,2)(1,0,1)[12] with Break dummy",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

# plot the forecast from 2020
fc_inverted %>%
  filter(Date >= yearmonth("2020 Jul")) %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data %>% filter(Date >= yearmonth("2020 Jul")), aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Forecasts using ARIMA5(1,1,2)(1,0,1)[12] with Break dummy",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

--------------------
# forecast on training data with post-break model
  # Split the post-break data into training and testing sets
 # split at end of 2023
bc_ev_post_break_train <- bc_ev_post_break %>%
  filter(Date <= yearmonth("2023 Dec"))
bc_ev_post_break_test <- bc_ev_post_break %>%
  filter(Date > yearmonth("2023 Dec"))

# 2. Refit best post-break ARIMA model on the training data only
# ARIMA(0, 0, 1)(1, 1, 0)[12] with drift
best_fit_post_train <- bc_ev_post_break_train %>%
  model(
    arima2 = ARIMA(Registrations_bc ~ 1 + pdq(0, 0, 1) + PDQ(1, 1, 0, period = 12))
  )
best_fit_post_train %>% report()

# 3. Forecast into the test period
fc_post_break <- forecast(best_fit_post_train, h = nrow(bc_ev_post_break_test))

# 4. Plot actual vs forecast (post-break data)
plot_forecast_arima <- fc_post_break %>%
  autoplot(bc_ev_data) +
  labs(title = "Forecast using ARIMA1(0,0,1)(1,1,0)[12] with drift on post-break series",
       x = "Date", y = "Box-Cox of Registrations") +
  theme_minimal()
plot_forecast_arima

# Plot actual vs forecast (post-break data) with test data
plot_postbreak_forecast_arima <- fc_post_break %>%
  autoplot(bc_ev_data %>% filter(Date >= yearmonth("2020 Jul"))) +
  labs(title = "Post-Break Forecast vs Actual",
       x = "Date", y = "Box-Cox of Registrations") +
  theme_minimal()
plot_postbreak_forecast_arima

# 5. Accuracy metrics: Compare forecast vs actual
accuracy(fc_post_break, bc_ev_data)
# ME: 0.089
# RMSE: 0.502
# MAE: 0.404
# MPE: 0.535
# MAPE: 2.49
# MASE: 0.240

  # extract prediction intervals
  fc_extracted <- fc_post_break %>%
  mutate(
    intervals = hilo(Registrations_bc, level = 0.95),
    lower = map_dbl(intervals, ~ .x$lower),
    upper = map_dbl(intervals, ~ .x$upper)
  )

# invert box cox manually
fc_inverted <- fc_extracted %>%
  mutate(
    mean_orig = inv_box_cox(.mean, lambda = lambda),
    lower_orig = inv_box_cox(lower, lambda = lambda),
    upper_orig = inv_box_cox(upper, lambda = lambda)
  )

# historical original series
original_data <- bc_ev_data %>%
  mutate(Registrations_orig = inv_box_cox(Registrations_bc, lambda))

# plot the forecast with original data
fc_inverted %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data, aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Forecasts using ARIMA1(0, 0, 1)(1, 1, 0)[12] with drift on post-break series",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

# plot the forecast from 2020
fc_inverted %>%
  filter(Date >= yearmonth("2020 Jul")) %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data %>% filter(Date >= yearmonth("2020 Jul")), aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Forecasts using ARIMA1(0, 0, 1)(1, 1, 0)[12] with drift on post-break series (2020 onwards)",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

--------------------

# forecast into the future
# with final model (ARIMA5(1,1,2)(1,0,1)[12] with Break dummy)
# create a future new_data frame
future_dates <- yearmonth("2025 Mar") + 0:23

# assume break = 1 for all future periods 
future_data <- tibble(
  Date = future_dates,
  Break = 1
) %>% as_tsibble(index = Date)

# generate the forecast
future_fc <- forecast(final_model, new_data = future_data)

future_fc <- final_model %>%
  forecast(new_data = future_data, level = c(80, 95))  # Try both

# plot forecast with historical data
future_fc %>%
  autoplot(bc_ev_data) +
  labs(title = "Future Forecast using ARIMA5(1,1,2)(1,0,1)[12] with Break dummy",
       y = "Box-Cox Transformed Registrations",
       x = "Date") +
  theme_minimal()

names(future_fc)

# Extract prediction intervals (convert hilo to columns)
future_fc_extracted <- future_fc %>%
  mutate(
    intervals = hilo(Registrations_bc, 95),
    lower = map_dbl(intervals, ~ .x$lower),
    upper = map_dbl(intervals, ~ .x$upper)
  )

# invert box cox manually
future_fc_inverted <- future_fc_extracted %>%
  mutate(
    mean_orig = inv_box_cox(.mean, lambda = lambda),
    lower_orig = inv_box_cox(lower, lambda = lambda),
    upper_orig = inv_box_cox(upper, lambda = lambda)
  )

# historical original series
original_data <- bc_ev_data %>%
  mutate(Registrations_orig = inv_box_cox(Registrations_bc, lambda))

# plot forecasts with historical data
future_fc_inverted %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data, aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Future Forecasts using ARIMA5(1,1,2)(1,0,1)[12] with Break dummy",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

# With final model ARIMA1(0, 0, 1)(1, 1, 0)[12] with drift on post-break series
# forecast into 3 years
  # 1. Refit the model on the full post-break data
  final_post_model <- bc_ev_post_break %>%
  model(
    arima_final = ARIMA(Registrations_bc ~ 1 + pdq(0, 0, 1) + PDQ(1, 1, 0, 12))
  )

# 2. Forecast 36 months into the future starting after Feb 2025
future_forecast <- final_post_model %>%
  forecast(h = 36)

# 3. Plot the forecast
plot_future_forecast_final <- future_forecast %>%
  autoplot(bc_ev_data) +
  labs(title = "Future Forecast using ARIMA1(0, 0, 1)(1, 1, 0)[12] with drift on post-break series",
       x = "Date", y = "Box-Cox of Registrations") +
  theme_minimal()
plot_future_forecast_final

# Back-transform

fc_future_extracted <- future_forecast %>%
  mutate(
    intervals = hilo(Registrations_bc, level = 0.95),
    lower = map_dbl(intervals, ~ .x$lower),
    upper = map_dbl(intervals, ~ .x$upper)
  )

# invert box cox manually
fc_future_inverted <- fc_future_extracted %>%
  mutate(
    mean_orig = inv_box_cox(.mean, lambda = lambda),
    lower_orig = inv_box_cox(lower, lambda = lambda),
    upper_orig = inv_box_cox(upper, lambda = lambda)
  )

# historical original series
original_data <- bc_ev_data %>%
  mutate(Registrations_orig = inv_box_cox(Registrations_bc, lambda))

# plot the forecast with original data
fc_future_inverted %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data, aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Future Forecasts using ARIMA1(0, 0, 1)(1, 1, 0)[12] with drift on post-break series",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

# plot the forecast from 2020
fc_future_inverted %>%
  filter(Date >= yearmonth("2020 Jul")) %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data %>% filter(Date >= yearmonth("2020 Jul")), aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Future Forecasts with Final ARIMA Model (2020 onwards)",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

----------------------------------------------------------
  # ETS
  # Fit ETS model on training data (post-break data)
  fit_ets_post_train <- bc_ev_post_break_train %>%
  model(ETS(Registrations_bc))
# report
report(fit_ets_post_train)
# ETS(A,N,A)

# fit multiple ets models with different components
ets_models <- bc_ev_post_break_train %>%
  model(
    ETS_AAA = ETS(Registrations_bc ~ error("A") + trend("A") + season("A")),
    ETS_ANA = ETS(Registrations_bc ~ error("A") + trend("N") + season("A")),
    ETS_AAN = ETS(Registrations_bc ~ error("A") + trend("A") + season("N")),
    ETS_MAM = ETS(Registrations_bc ~ error("M") + trend("A") + season("M")),
    ETS_MNM = ETS(Registrations_bc ~ error("M") + trend("N") + season("M"))
  )

ets_models %>%
  glance() %>%
  arrange(AICc) %>%
  select(.model, AICc, BIC)
# The best model is ETS(A,N,A) with AICc = 163, BIC = 169

# Map dof for each model
dof_map_ets <- c(
  ETS_AAA = 3,
  ETS_ANA = 2,
  ETS_AAN = 2,
  ETS_MAM = 3,
  ETS_MNM = 2
)

# Apply Ljung-Box test using correct dof
ljung_ets_results <- purrr::imap_dfr(dof_map_ets, function(dof, model_name) {
  ets_models %>%
    select(!!model_name) %>%
    augment() %>%
    features(.innov, ljung_box, lag = 12, dof = dof) %>%
    mutate(model = model_name)
})

# View results
ljung_ets_results %>%
  select(model, lb_stat, lb_pvalue) %>%
  arrange(desc(lb_pvalue))

# report for ets_ANA
ets_models %>%
  select(ETS_ANA) %>%
  report()

ets_models %>%
  select(ETS_ANA) %>%
  gg_tsresiduals()

ets_models %>%
  select(ETS_ANA) %>%
  augment() %>%
  features(.innov, ljung_box, lag = 12, dof = 2)
# p-value 0.167 > 0.05, so we fail to reject the null hypothesis of no autocorrelation
# there is no autocorrelation in the residuals

# Forecast using ETS model
forecast_ets_post <- fit_ets_post_train %>%
  forecast(h = nrow(bc_ev_post_break_test))

forecast_ets_best <- ets_models %>%
  select(ETS_ANA) %>%
  forecast(h = nrow(bc_ev_post_break_test))

# Plot the forecast
plot_post_forecast_ets <- forecast_ets_best %>%
  autoplot(bc_ev_data) +
  autolayer(bc_ev_post_break_test, Registrations_bc) +
  labs(title = "Forecast using ETS(A,N,A) on post-break data",
       y = "Box-Cox Transformed Registrations",
       x = "Date") +
  theme_minimal()
plot_post_forecast_ets

# Compare accuracy
accuracy(forecast_ets_best, bc_ev_data)

# extract prediction intervals (convert hilo to columns)
forecast_ets_extracted <- forecast_ets_post %>%
  mutate(
    intervals = hilo(Registrations_bc, 95),
    lower = map_dbl(intervals, ~ .x$lower),
    upper = map_dbl(intervals, ~ .x$upper)
  )

# Invert Box-Cox transformation manually
forecast_ets_inverted <- forecast_ets_extracted %>%
  mutate(
    mean_orig = inv_box_cox(.mean, lambda = lambda),
    lower_orig = inv_box_cox(lower, lambda = lambda),
    upper_orig = inv_box_cox(upper, lambda = lambda)
  )

# Historical original series
original_data_ets <- bc_ev_data %>%
  mutate(Registrations_orig = inv_box_cox(Registrations_bc, lambda))
# Plot forecasts with historical data
forecast_ets_inverted %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data_ets, aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Forecasts using ETS(A,N,A) on post-break data",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

# Plot forecast from 2020
# plot the forecast from 2020
forecast_ets_inverted %>%
  filter(Date >= yearmonth("2020 Jul")) %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data_ets %>% filter(Date >= yearmonth("2020 Jul")), aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Forecasts using ETS(A,N,A) on post-break data",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

# Forecast into the future
# Combine training and test data (up to Feb 2025)
bc_ev_post_break_full <- bc_ev_data %>%
  filter(Date <= yearmonth("2025 Feb"))

# Fit ETS model on all available data up to Feb 2025
fit_ets_post_full <- bc_ev_post_break_full %>%
  model(ETS(Registrations_bc ~ error("A") + trend("N") + season("A")))


# 1. Set the Forecast Horizon (e.g., 12 months into the future)
forecast_horizon_ets <- 36

# 2. Make Forecast for the Future
#fc_future_ets <- forecast(fit_ets_post_train, h = forecast_horizon_ets)
fc_future_ets <- forecast(fit_ets_post_full, h = forecast_horizon_ets)

# 3. Plot the forecast
plot_future_forecast_ets <- fc_future_ets %>%
  autoplot(bc_ev_data) +
  labs(title = "Future Forecast using ETS(A,N,A) on post-break data",
       x = "Date", y = "Box-Cox of Registrations") +
  theme_minimal()
plot_future_forecast_ets

# Back-transform the forecast
back_transform_forecast <- function(forecast, lambda) {
  if (lambda == 0) {
    exp(forecast) # If lambda is 0, use the exponential transformation
  } else {
    (forecast * lambda + 1)^(1/lambda) # Box-Cox inverse transformation
  }
}
# Apply back transformation to the forecast
fc_future_ets_extracted <- fc_future_ets %>%
  mutate(
    intervals = hilo(Registrations_bc, 95),
    lower = map_dbl(intervals, ~ .x$lower),
    upper = map_dbl(intervals, ~ .x$upper)
  )
fc_future_ets_inverted <- fc_future_ets_extracted %>%
  mutate(
    mean_orig = back_transform_forecast(.mean, lambda),
    lower_orig = back_transform_forecast(lower, lambda),
    upper_orig = back_transform_forecast(upper, lambda)
  )
# Historical original series
original_data_future <- bc_ev_data %>%
  mutate(Registrations_orig = inv_box_cox(Registrations_bc, lambda))
# Plot forecasts with historical data
fc_future_ets_inverted %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data_future, aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Future Forecasts using ETS(A,N,A) on post-break data",
       y = "Registrations",
       x = "Date") +
  theme_minimal()

fc_future_ets_inverted %>%
  filter(Date >= yearmonth("2020 Jul")) %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = lower_orig, ymax = upper_orig), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_orig), color = "blue", size = 1) +
  geom_line(data = original_data_future %>% filter(Date >= yearmonth("2020 Jul")), aes(x = Date, y = Registrations_orig), color = "black") +
  labs(title = "Future Forecasts using ETS(A,N,A) on post-break data",
       y = "Registrations",
       x = "Date") +
  theme_minimal()
