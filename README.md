# Forecasting Electric Vehicle Registrations in Denmark: A Structural Time Series Analysis

## Overview

This project investigates the forecasting of monthly electric vehicle (EV) registrations in Denmark using statistical time series methods. The analysis compares ARIMA and Exponential Smoothing (ETS) models to identify the most accurate approach for forecasting EV registration trends.

The study places particular emphasis on identifying structural changes in the registration series and evaluating whether modelling the period after the identified structural break produces more accurate forecasts.

## Research Question

Which statistical time series model provides the most accurate forecasts for monthly electric vehicle registrations in Denmark?

## Data

The analysis uses monthly passenger car registration data from Statistics Denmark's BIL51 dataset, covering January 2011 to February 2025, with a total of 170 observations.

The dataset contains new passenger car registrations by ownership and fuel type. The analysis focuses on electric vehicle registrations.

## Methodology

The project follows a structured time series forecasting approach:

- **Exploratory Data Analysis:** Examined trends, seasonality and autocorrelation patterns using visualisations, seasonal plots and STL decomposition.
- **Data Transformation:** Evaluated logarithmic and Box-Cox transformations to address potential variance instability.
- **Stationarity Testing:** Applied Augmented Dickey-Fuller (ADF) and KPSS tests to assess stationarity.
- **Structural Break Detection:** Used the QLR (supF) test to identify significant changes in the underlying registration trend.
- **Model Development:** Estimated and compared ARIMA and ETS models, considering both the full time series and the period following the identified structural break.
- **Model Diagnostics:** Evaluated residuals using Ljung-Box tests and compared model performance using information criteria and forecast accuracy measures.
- **Forecasting:** Generated 36-month forecasts and transformed predictions back to the original registration scale.

## Key Findings and Results

The structural break analysis identified **July 2020** as a significant turning point in the EV registration series.

The analysis compared three main forecasting approaches:

| Model | MAPE | MAE | RMSE | MASE |
|---|---:|---:|---:|---:|
| Full-series ARIMA with break dummy | 11.70% | 1.88 | 2.05 | 1.13 |
| Post-break ARIMA with drift | **2.49%** | **0.404** | **0.502** | **0.24** |
| ETS(A,N,A) | 8.07% | – | – | 0.778 |

The post-break ARIMA model achieved the best forecasting accuracy among the evaluated models.

**Selected model:** ARIMA(0,0,1)(1,1,0)[12] with drift, estimated using data from the post-July 2020 period.

The Ljung-Box test returned a p-value of 0.395 for the selected model, indicating no statistically significant residual autocorrelation at conventional significance levels.

## Final Forecast

The selected post-break ARIMA model was used to generate a **36-month forecast**, covering March 2025 to February 2028.

The results demonstrate the importance of accounting for structural changes in time series data when developing forecasting models. In this case, focusing on the period following July 2020 substantially improved forecast accuracy compared with modelling the entire historical series.

## Tools and Technologies

- **R:** Statistical analysis and time series forecasting
- **Python:** Initial data cleaning and preparation
- **fpp3, forecast, fable:** Time series modelling and forecasting
- **strucchange:** Structural break detection
- **urca, tseries:** Stationarity testing
- **readr, purrr:** Data handling and workflow support

## Project Files

- **R Analysis Script:** `EV_Registration_Forecasting.R`
- **Project Report:** `EV_Registration_Forecasting_Report.pdf`

## Academic Context

- **Course:** Predictive Analytics (CDSCO1005E)
- **Programme:** MSc Business Administration and Data Science
- **Institution:** Copenhagen Business School
- **Project Type:** Home Assignment
