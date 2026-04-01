# Chapter 15 Lecture -- Nonlinear Regression

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 15

> **Prerequisites**: [Statistics] Linear regression (Ch 14).
>
> **Learning Objectives**:
> 1. Apply polynomial and multiple linear regression
> 2. Implement Gauss-Newton method for nonlinear regression
> 3. Compare model selection criteria

---

<br>

## Table of Contents

- [1. Review: Standard Error and Coefficient of Determination](#1-review-standard-error-and-coefficient-of-determination)
  - [1.1 Sample Standard Deviation (SD)](#11-sample-standard-deviation-sd)
  - [1.2 Standard Error of the Mean (SE)](#12-standard-error-of-the-mean-se)
  - [1.3 Standard Error of Regression](#13-standard-error-of-regression)
  - [1.4 Coefficient of Determination (R-squared)](#14-coefficient-of-determination-r-squared)
- [2. Linearization of Nonlinear Relationships](#2-linearization-of-nonlinear-relationships)
  - [2.1 Exponential Model](#21-exponential-model)
  - [2.2 Power Model](#22-power-model)
  - [2.3 Saturation-Growth-Rate Model](#23-saturation-growth-rate-model)
- [3. Polynomial Regression](#3-polynomial-regression)
  - [3.1 Problem Formulation](#31-problem-formulation)
  - [3.2 Deriving the Normal Equations](#32-deriving-the-normal-equations)
  - [3.3 Normal Equations in Matrix Form](#33-normal-equations-in-matrix-form)
  - [3.4 Standard Error of the Estimate](#34-standard-error-of-the-estimate)
- [4. Multiple Linear Regression](#4-multiple-linear-regression)
  - [4.1 Problem Formulation](#41-problem-formulation)
  - [4.2 Normal Equations](#42-normal-equations)
  - [4.3 Standard Error of the Estimate](#43-standard-error-of-the-estimate)
  - [4.4 Extension to Multivariable Power Models](#44-extension-to-multivariable-power-models)
- [5. General Linear Least Square Regression](#5-general-linear-least-square-regression)
  - [5.1 General Model Formulation](#51-general-model-formulation)
  - [5.2 Matrix Formulation](#52-matrix-formulation)
  - [5.3 Deriving the Normal Equations (Matrix Calculus)](#53-deriving-the-normal-equations-matrix-calculus)
  - [5.4 Computing the Solution](#54-computing-the-solution)
  - [5.5 Goodness-of-Fit Statistics](#55-goodness-of-fit-statistics)
  - [5.6 Covariance Matrix of Parameter Estimates](#56-covariance-matrix-of-parameter-estimates)
- [6. Nonlinear Regression (Gauss-Newton)](#6-nonlinear-regression-gauss-newton)
  - [6.1 Motivation: The Antoine Equation](#61-motivation-the-antoine-equation)
  - [6.2 Python Implementation with scipy.optimize](#62-python-implementation-with-scipyoptimize)
- [7. Summary Table](#7-summary-table)

---

<br>

## 1. Review: Standard Error and Coefficient of Determination

### 1.1 Sample Standard Deviation (SD)

**Sample Standard Deviation** (SD) measures how much individual data points vary (spread out) from the mean of a dataset.

$$
s_y = \sqrt{\frac{1}{n-1} \sum_{j=1}^{n} (y_j - \bar{y})^2}
$$

where:
- $n$ = number of data points
- $y_j$ = individual observation
- $\bar{y}$ = sample mean

**Example:** Heights = [170, 172, 168, 174, 176] cm

- Sample mean $\bar{y} = 172$ cm
- $s_y = 3.16$ cm

> **Interpretation:** On average, individual height measurements differ from the mean by about 3.16 cm.

### 1.2 Standard Error of the Mean (SE)

**Standard Error** (SE) measures how precisely the sample mean estimates the true population mean.

$$
SE = \frac{s_y}{\sqrt{n}}
$$

**Example (continued):**
- $SE = \frac{3.16}{\sqrt{5}} = 1.41$ cm

> **Interpretation:** The average height is estimated to be 172 cm, with a standard error of $\pm 1.41$ cm.

### 1.3 Standard Error of Regression

**Standard Error** (SE) of **Regression** measures how much observed values deviate from the linear regression line.

$$
s_e = \sqrt{\frac{1}{n-2} \sum_{j=1}^{n} (y_j - \hat{y}_j)^2} = \sqrt{\frac{SSE}{n-2}}
$$

where:
- $\hat{y}_j$ = predicted value from the regression line
- $n - 2$ = degrees of freedom (two parameters $\beta_0, \beta_1$ are estimated)

**Example:** Heights = [170, 172, 168, 174, 176] cm, $x = [1, 2, 3, 4, 5]$

- $SE = \sqrt{\frac{SSE}{n-2}} = 2.61$ cm

> **Interpretation:** The observed $y$-values deviate from the regression line by about 2.61 cm.

### 1.4 Coefficient of Determination (R-squared)

**Coefficient of Determination** $R^2$ measures how well the regression model explains the variability of the dependent variable.

The total sum of squares is decomposed as:

$$
SST = SSR + SSE
$$

Dividing both sides by $SST$:

$$
1 = \frac{SSR}{SST} + \frac{SSE}{SST}
$$

Therefore:

$$
R^2 = \frac{SSR}{SST} = 1 - \frac{SSE}{SST}
$$

where:
- $SST = \sum (y_i - \bar{y})^2$ -- Total Sum of Squares
- $SSR = \sum (\hat{y}_i - \bar{y})^2$ -- Regression Sum of Squares (explained variation)
- $SSE = \sum (y_i - \hat{y}_i)^2$ -- Error Sum of Squares (unexplained variation)

| $R^2$ Value | Meaning |
|---|---|
| $R^2 = 1$ | Perfect fit |
| $R^2 = 0$ | No fit |
| $0 < R^2 < 1$ | Partial fit |

**Example (ex 14.6):** $R^2 = 88\%$

> **Interpretation:** The relationship between the two variables explains 88% (approximately 90%) of the variation in the data.

---

<br>

## 2. Linearization of Nonlinear Relationships

The key idea is to **transform a nonlinear function into a straight line form** so that linear regression techniques can be applied.

### 2.1 Exponential Model

**Original model:**

$$
y = \alpha_1 e^{\beta_1 x}
$$

**Linearization:** Take the natural logarithm of both sides:

$$
\ln y = \ln(\alpha_1 e^{\beta_1 x}) = \ln \alpha_1 + \ln e^{\beta_1 x} = \ln \alpha_1 + \beta_1 x
$$

This gives a linear form: $\ln y = \underbrace{\ln \alpha_1}_{\text{intercept}} + \underbrace{\beta_1}_{\text{slope}} \cdot x$

> **Graphical Effect:** A plot of $y$ vs $x$ shows an exponential curve. A plot of $\ln y$ vs $x$ becomes a straight line with slope $\beta_1$ and intercept $\ln \alpha_1$.

### 2.2 Power Model

**Original model:**

$$
y = \alpha_2 x^{\beta_2}
$$

**Linearization:** Take the natural logarithm of both sides:

$$
\ln y = \ln(\alpha_2 x^{\beta_2}) = \ln \alpha_2 + \ln x^{\beta_2} = \ln \alpha_2 + \beta_2 \ln x
$$

This gives a linear form: $\ln y = \underbrace{\ln \alpha_2}_{\text{intercept}} + \underbrace{\beta_2}_{\text{slope}} \cdot \ln x$

> **Graphical Effect:** A plot of $y$ vs $x$ shows a power curve. A plot of $\ln y$ vs $\ln x$ becomes a straight line with slope $\beta_2$ and intercept $\ln \alpha_2$.

### 2.3 Saturation-Growth-Rate Model

**Original model:**

$$
y = \frac{\alpha_3 x}{\beta_3 + x}
$$

**Linearization:** Invert both sides:

$$
\frac{1}{y} = \frac{\beta_3 + x}{\alpha_3 x} = \frac{\beta_3}{\alpha_3} \cdot \frac{1}{x} + \frac{1}{\alpha_3}
$$

This gives a linear form: $\frac{1}{y} = \underbrace{\frac{\beta_3}{\alpha_3}}_{\text{slope}} \cdot \frac{1}{x} + \underbrace{\frac{1}{\alpha_3}}_{\text{intercept}}$

> **Graphical Effect:** A plot of $y$ vs $x$ shows a saturation curve. A plot of $\frac{1}{y}$ vs $\frac{1}{x}$ becomes a straight line with slope $\frac{\beta_3}{\alpha_3}$ and intercept $\frac{1}{\alpha_3}$.

| Model | Original Form | Transformed Variables | Slope | Intercept |
|---|---|---|---|---|
| Exponential | $y = \alpha_1 e^{\beta_1 x}$ | $\ln y$ vs $x$ | $\beta_1$ | $\ln \alpha_1$ |
| Power | $y = \alpha_2 x^{\beta_2}$ | $\ln y$ vs $\ln x$ | $\beta_2$ | $\ln \alpha_2$ |
| Saturation-Growth | $y = \frac{\alpha_3 x}{\beta_3 + x}$ | $\frac{1}{y}$ vs $\frac{1}{x}$ | $\frac{\beta_3}{\alpha_3}$ | $\frac{1}{\alpha_3}$ |

---

<br>

## 3. Polynomial Regression

### 3.1 Problem Formulation

When data shows a curvilinear pattern (e.g., parabolic), a linear line is not sufficient. We extend the linear regression idea to **polynomial regression**.

**Goal:**
- Extend the linear regression idea to polynomial regression
- Estimate confidence intervals

Given a set of data $\{x_i, y_i\}_{i=1}^{n}$, find $\beta_0, \beta_1, \beta_2$ such that:

$$
\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2
$$

The objective is to minimize the sum of squared errors:

$$
\min_{\beta_0, \beta_1, \beta_2} SSE
$$

where:

$$
SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)^2
$$

### 3.2 Deriving the Normal Equations

To determine $\beta_0, \beta_1, \beta_2$, take partial derivatives and set them to zero:

$$
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-1) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i^2) = 0
$$

After simplification, the **normal equations** become:

$$
\sum y_i = \beta_0 \cdot n + \beta_1 \sum x_i + \beta_2 \sum x_i^2
$$

$$
\sum y_i x_i = \beta_0 \sum x_i + \beta_1 \sum x_i^2 + \beta_2 \sum x_i^3
$$

$$
\sum y_i x_i^2 = \beta_0 \sum x_i^2 + \beta_1 \sum x_i^3 + \beta_2 \sum x_i^4
$$

### 3.3 Normal Equations in Matrix Form

The normal equations can be written compactly as:

$$
\begin{bmatrix} \sum y_i \\ \sum y_i x_i \\ \sum y_i x_i^2 \end{bmatrix}
=
\begin{bmatrix} n & \sum x_i & \sum x_i^2 \\ \sum x_i & \sum x_i^2 & \sum x_i^3 \\ \sum x_i^2 & \sum x_i^3 & \sum x_i^4 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
$$

> **Key Observation:** This is a $3 \times 3$ linear system that can be solved using methods from Part 3 (e.g., Gaussian elimination, LU decomposition).

### 3.4 Standard Error of the Estimate

For polynomial regression with 3 parameters ($\beta_0, \beta_1, \beta_2$):

$$
e_i = y_i - \hat{y}_i
$$

$$
SSE = \sum e_i^2
$$

$$
s_e = \sqrt{\frac{SSE}{n - 3}}
$$

> **Degrees of Freedom:** The denominator is $n - 3$ because three parameters ($\beta_0, \beta_1, \beta_2$) are used to construct $\hat{y}$. In general, for an $m$-th degree polynomial, the degrees of freedom are $n - (m+1)$.

---

<br>

## 4. Multiple Linear Regression

Multiple linear regression extends the concept from 1D to 2D (and beyond): fitting data with **two or more independent variables**.

### 4.1 Problem Formulation

The model fits a plane (or hyperplane) to the data:

$$
\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2
$$

The objective is:

$$
SSE = \sum e_i^2 = \sum (y_i - \hat{y}_i)^2 = \sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})^2
$$

### 4.2 Normal Equations

Seek $\beta_0, \beta_1, \beta_2$ such that $\min_{\beta_0, \beta_1, \beta_2} SSE$:

$$
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-1) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{1i}) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{2i}) = 0
$$

In matrix form:

$$
\begin{bmatrix} \sum y_i \\ \sum y_i x_{1i} \\ \sum y_i x_{2i} \end{bmatrix}
=
\begin{bmatrix} n & \sum x_{1i} & \sum x_{2i} \\ \sum x_{1i} & \sum x_{1i}^2 & \sum x_{1i} x_{2i} \\ \sum x_{2i} & \sum x_{1i} x_{2i} & \sum x_{2i}^2 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
$$

### 4.3 Standard Error of the Estimate

For multiple linear regression with two independent variables (3 parameters):

$$
s_e = \sqrt{\frac{SSE}{n - 3}}
$$

### 4.4 Extension to Multivariable Power Models

Multiple linear regression can be extended to **multivariable power models**:

$$
y = \beta_0 x_1^{\beta_1} x_2^{\beta_2} \cdots x_m^{\beta_m} + \varepsilon
$$

Taking the natural logarithm:

$$
\ln y = \ln \beta_0 + \beta_1 \ln x_1 + \beta_2 \ln x_2 + \cdots + \beta_m \ln x_m
$$

> **Key Insight:** After the log transformation, this becomes a standard multiple linear regression problem in the transformed variables $\ln y$, $\ln x_1$, $\ln x_2$, ..., $\ln x_m$.

---

<br>

## 5. General Linear Least Square Regression

### 5.1 General Model Formulation

The general linear least square model unifies all previous cases. The model is "linear in the parameters" $\beta_j$, but the basis functions $f_j$ can be arbitrary (e.g., trigonometric):

**Example -- Fourier-type model:**

$$
y = \beta_0 + \beta_1 \cos(\omega x) + \beta_2 \sin(\omega x)
$$

**General form:**

$$
y = \beta_0 f_0(x_1, x_2, \ldots, x_p) + \beta_1 f_1(x_1, x_2, \ldots, x_p) + \cdots + \beta_m f_m(x_1, x_2, \ldots, x_p)
$$

> **Important:** "Linear" here means linear in the **parameters** $\beta_j$, not necessarily in $x$. The functions $f_j$ can be nonlinear in $x$.

### 5.2 Matrix Formulation

Collect $n$ data points $\{y_i, x_{1i}, x_{2i}, \ldots, x_{pi}\}_{i=1}^{n}$.

For each data point $i = 1, 2, \ldots, n$:

$$
y_i = \beta_0 f_0(x_{1i}, x_{2i}, \ldots, x_{pi}) + \beta_1 f_1(x_{1i}, x_{2i}, \ldots, x_{pi}) + \cdots + \beta_m f_m(x_{1i}, x_{2i}, \ldots, x_{pi})
$$

In matrix form:

$$
\underset{(n \times 1)}{\mathbf{y}}
=
\underset{(n \times (m+1))}{\mathbf{X}}
\underset{((m+1) \times 1)}{\boldsymbol{\beta}}
+ \underset{(n \times 1)}{\mathbf{e}}
$$

where:

$$
\mathbf{y} = \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{bmatrix}, \quad
\mathbf{X} = \begin{bmatrix} f_{01} & f_{11} & f_{21} & \cdots & f_{m1} \\ f_{02} & f_{12} & f_{22} & \cdots & f_{m2} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ f_{0n} & f_{1n} & f_{2n} & \cdots & f_{mn} \end{bmatrix}, \quad
\boldsymbol{\beta} = \begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \\ \vdots \\ \beta_m \end{bmatrix}, \quad
\mathbf{e} = \begin{bmatrix} e_1 \\ e_2 \\ \vdots \\ e_n \end{bmatrix}
$$

Here $f_{ji}$ denotes the value of basis function $f_j$ evaluated at data point $i$, and $\mathbf{e}$ is the residual (error) vector.

> **Notation:** $\mathbf{X}$ has $n$ rows (data points) and $m+1$ columns (number of parameters including $\beta_0$).

### 5.3 Deriving the Normal Equations (Matrix Calculus)

The sum of squared errors in matrix form:

$$
SSE = \sum e_i^2 = \mathbf{e}^T \mathbf{e} = (\mathbf{y} - \hat{\mathbf{y}})^T (\mathbf{y} - \hat{\mathbf{y}})
$$

To minimize, take the derivative with respect to $\boldsymbol{\beta}$ and set to zero:

$$
\frac{\partial SSE}{\partial \boldsymbol{\beta}} = \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{e}^T \mathbf{e})
= \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})
= (-2\mathbf{X})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta}) = \mathbf{0}
$$

This yields the **normal equation**:

$$
\mathbf{X}^T \mathbf{X} \boldsymbol{\beta} = \mathbf{X}^T \mathbf{y}
$$

Solving for $\boldsymbol{\beta}$:

$$
\boxed{\boldsymbol{\beta} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}}
$$

> **Key Formula:** This is the most important result of the chapter. It provides a closed-form solution for the parameter vector $\boldsymbol{\beta}$ that minimizes the sum of squared errors, applicable to any general linear model.

### 5.4 Computing the Solution

Once $\boldsymbol{\beta}$ is computed:

| Quantity | Formula |
|---|---|
| Residual vector | $\mathbf{e} = \mathbf{y} - \mathbf{X}\boldsymbol{\beta}$ |
| Sum of squared errors | $SSE = \mathbf{e}^T \mathbf{e}$ |
| Predicted values | $\hat{\mathbf{y}} = \mathbf{X}\boldsymbol{\beta}$ |

### 5.5 Goodness-of-Fit Statistics

$$
SST = \sum (y_i - \bar{y})^2
$$

$$
SSR = \sum (\hat{y}_i - \bar{y})^2
$$

$$
R^2 = \frac{SSR}{SST}
$$

$$
SSE = SST - SSR
$$

**Standard error of the estimate** (general case with $m+1$ parameters):

$$
s_e = \sqrt{\frac{SSE}{n - (m+1)}}
$$

> **Degrees of Freedom:** The denominator $n - (m+1)$ accounts for the $m+1$ parameters used to fit the model.

### 5.6 Covariance Matrix of Parameter Estimates

The covariance matrix of the estimated parameters:

$$
\text{cov}(\hat{\boldsymbol{\beta}}) = \sigma_e^2 (\mathbf{X}^T \mathbf{X})^{-1}
$$

where $\sigma_e^2$ is estimated by $s_e^2 = \frac{SSE}{n-(m+1)}$.

> **Practical Use:** The diagonal elements of $\text{cov}(\hat{\boldsymbol{\beta}})$ give the variance of each individual parameter estimate, and their square roots provide the standard errors of the parameters. These can be used to construct confidence intervals.

---

<br>

## 6. Nonlinear Regression (Gauss-Newton)

### 6.1 Motivation: The Antoine Equation

Some models are **inherently nonlinear** in the parameters and cannot be linearized by transformation. For example, the **Antoine Equation** for vapor pressure:

$$
\log_{10} P_v = A - \frac{B}{C + T}
$$

where:
- $P_v$ = vapor pressure
- $T$ = temperature
- $A, B, C$ = adjustable parameters

Since $C$ appears inside a nonlinear expression ($C + T$ in the denominator), this equation **cannot** be expressed as a general linear model. We need nonlinear regression (optimization) to fit it.

### 6.2 Python Implementation with scipy.optimize

```python
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Sample vapor pressure data (e.g., for water) in mmHg at different temperatures (deg C)
T_data = np.array([40, 50, 60, 70, 80])    # Temperature in Celsius
P_data = np.array([55.3, 92.5, 149.4, 233.7, 355.1])  # Vapor pressure in mmHg

# Convert pressure to log base 10
logP_data = np.log10(P_data)

# Antoine equation model: log10(P) = A - B / (C + T)
def antoine_model(params, T):
    A, B, C = params
    return A - B / (C + T)

# Objective function: sum of squared errors between predicted and actual logP
def objective(params):
    predictions = antoine_model(params, T_data)
    return np.sum((logP_data - predictions) ** 2)

# Initial guess for A, B, C
initial_guess = [8, 1500, 200]

# Minimize the objective function
result = minimize(objective, initial_guess, method='L-BFGS-B')

# Extract optimized parameters
A_opt, B_opt, C_opt = result.x
```

> **Method:** The `scipy.optimize.minimize` function with `method='L-BFGS-B'` uses a quasi-Newton optimization algorithm. An initial guess is required, and the algorithm iteratively refines the parameters to minimize $SSE$.

> **Initial Guess Sensitivity:** Nonlinear regression methods depend on the initial guess. A poor initial guess may lead to convergence to a local minimum rather than the global minimum.

---

<br>

## 7. Summary Table

| Topic | Model | Parameters | Normal Equation / Method | $s_e$ |
|---|---|---|---|---|
| **Linear Regression** | $\hat{y} = \beta_0 + \beta_1 x$ | 2 | Closed-form (Ch.14) | $\sqrt{\frac{SSE}{n-2}}$ |
| **Polynomial Regression** | $\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2$ | 3 | $3 \times 3$ normal equations | $\sqrt{\frac{SSE}{n-3}}$ |
| **Multiple Linear Regression** | $\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2$ | 3 | $3 \times 3$ normal equations | $\sqrt{\frac{SSE}{n-3}}$ |
| **General Linear Least Squares** | $\hat{y} = \sum \beta_j f_j(\mathbf{x})$ | $m+1$ | $\boldsymbol{\beta} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\mathbf{y}$ | $\sqrt{\frac{SSE}{n-(m+1)}}$ |
| **Linearized Models** | Exponential / Power / Saturation | 2 | Transform then linear regression | depends on transform |
| **Nonlinear Regression** | e.g., Antoine eq. | varies | Iterative optimization (e.g., L-BFGS-B) | $\sqrt{\frac{SSE}{n-p}}$ |
