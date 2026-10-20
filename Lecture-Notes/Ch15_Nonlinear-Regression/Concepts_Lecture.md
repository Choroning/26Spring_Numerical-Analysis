# Chapter 15 — Nonlinear Regression

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 15

> **Learning Objectives:**
> 1. Linearize exponential, power, and saturation-growth-rate models and fit them by linear regression
> 2. Derive and apply the normal equations of polynomial and multiple linear regression
> 3. Derive the normal equations of the general linear least-squares model in matrix form and compute goodness-of-fit statistics
> 4. Fit models that are inherently nonlinear in their parameters by iterative optimization

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
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Review: Standard Error and Coefficient of Determination

### 1.1 Sample Standard Deviation (SD)

**Sample Standard Deviation** (SD) measures how much individual data points vary (spread out) from the mean of a dataset.

```math
s_y = \sqrt{\frac{1}{n-1} \sum_{j=1}^{n} (y_j - \bar{y})^2}
```

where:
- $`n`$ = number of data points
- $`y_j`$ = individual observation
- $`\bar{y}`$ = sample mean

**Example:** Heights = [170, 172, 168, 174, 176] cm

- Sample mean $`\bar{y} = 172`$ cm
- $`s_y = 3.16`$ cm

> **Note:** **Interpretation.** Individual height measurements differ from the mean by roughly 3.16 cm.

### 1.2 Standard Error of the Mean (SE)

**Standard Error** (SE) measures how precisely the sample mean estimates the true population mean.

```math
SE = \frac{s_y}{\sqrt{n}}
```

**Example (continued):**
- $`SE = \frac{3.16}{\sqrt{5}} = 1.41`$ cm

> **Note:** **Interpretation.** The average height is estimated to be 172 cm, with a standard error of $`\pm 1.41`$ cm.

### 1.3 Standard Error of Regression

The **standard error of regression** ($`s_e`$) measures how much observed values deviate from the linear regression line.

```math
s_e = \sqrt{\frac{1}{n-2} \sum_{j=1}^{n} (y_j - \hat{y}_j)^2} = \sqrt{\frac{SSE}{n-2}}
```

where:
- $`\hat{y}_j`$ = predicted value from the regression line
- $`n - 2`$ = degrees of freedom (two parameters $`\beta_0, \beta_1`$ are estimated)

**Example:** Heights = [170, 172, 168, 174, 176] cm, $`x = [1, 2, 3, 4, 5]`$

- $`s_e = \sqrt{\frac{SSE}{n-2}} = 2.61`$ cm

> **Note:** **Interpretation.** The observed $`y`$-values deviate from the regression line by about 2.61 cm.

### 1.4 Coefficient of Determination (R-squared)

**Coefficient of Determination** $`R^2`$ measures how well the regression model explains the variability of the dependent variable.

The total sum of squares is decomposed as:

```math
SST = SSR + SSE
```

Dividing both sides by $`SST`$:

```math
1 = \frac{SSR}{SST} + \frac{SSE}{SST}
```

Therefore:

```math
R^2 = \frac{SSR}{SST} = 1 - \frac{SSE}{SST}
```

where:
- $`SST = \sum (y_i - \bar{y})^2`$: Total Sum of Squares
- $`SSR = \sum (\hat{y}_i - \bar{y})^2`$: Regression Sum of Squares (explained variation)
- $`SSE = \sum (y_i - \hat{y}_i)^2`$: Error Sum of Squares (unexplained variation)

| $`R^2`$ Value | Meaning |
|---|---|
| $`R^2 = 1`$ | Perfect fit |
| $`R^2 = 0`$ | No fit |
| $`0 < R^2 < 1`$ | Partial fit |

*Table 1. Meaning of coefficient of determination values*

**Example (textbook Example 14.6):** $`R^2 = 88\%`$

> **Note:** **Interpretation.** The relationship between the two variables explains 88% (approximately 90%) of the variation in the data.

---

<br>

## 2. Linearization of Nonlinear Relationships

The key idea is to **transform a nonlinear function into a straight line form** so that linear regression techniques can be applied.

### 2.1 Exponential Model

**Original model:**

```math
y = \alpha_1 e^{\beta_1 x}
```

**Linearization:** Take the natural logarithm of both sides:

```math
\ln y = \ln(\alpha_1 e^{\beta_1 x}) = \ln \alpha_1 + \ln e^{\beta_1 x} = \ln \alpha_1 + \beta_1 x
```

This gives a linear form: $`\ln y = \underbrace{\ln \alpha_1}_{\text{intercept}} + \underbrace{\beta_1}_{\text{slope}} \cdot x`$

> **Note:** **Graphical Effect.** A plot of $`y`$ vs $`x`$ shows an exponential curve. A plot of $`\ln y`$ vs $`x`$ becomes a straight line with slope $`\beta_1`$ and intercept $`\ln \alpha_1`$.

### 2.2 Power Model

**Original model:**

```math
y = \alpha_2 x^{\beta_2}
```

**Linearization:** Take the natural logarithm of both sides:

```math
\ln y = \ln(\alpha_2 x^{\beta_2}) = \ln \alpha_2 + \ln x^{\beta_2} = \ln \alpha_2 + \beta_2 \ln x
```

This gives a linear form: $`\ln y = \underbrace{\ln \alpha_2}_{\text{intercept}} + \underbrace{\beta_2}_{\text{slope}} \cdot \ln x`$

> **Note:** **Graphical Effect.** A plot of $`y`$ vs $`x`$ shows a power curve. A plot of $`\ln y`$ vs $`\ln x`$ becomes a straight line with slope $`\beta_2`$ and intercept $`\ln \alpha_2`$.

### 2.3 Saturation-Growth-Rate Model

**Original model:**

```math
y = \frac{\alpha_3 x}{\beta_3 + x}
```

**Linearization:** Invert both sides:

```math
\frac{1}{y} = \frac{\beta_3 + x}{\alpha_3 x} = \frac{\beta_3}{\alpha_3} \cdot \frac{1}{x} + \frac{1}{\alpha_3}
```

This gives a linear form: $`\frac{1}{y} = \underbrace{\frac{\beta_3}{\alpha_3}}_{\text{slope}} \cdot \frac{1}{x} + \underbrace{\frac{1}{\alpha_3}}_{\text{intercept}}`$

> **Note:** **Graphical Effect.** A plot of $`y`$ vs $`x`$ shows a saturation curve. A plot of $`\frac{1}{y}`$ vs $`\frac{1}{x}`$ becomes a straight line with slope $`\frac{\beta_3}{\alpha_3}`$ and intercept $`\frac{1}{\alpha_3}`$.

| Model | Original Form | Transformed Variables | Slope | Intercept |
|---|---|---|---|---|
| Exponential | $`y = \alpha_1 e^{\beta_1 x}`$ | $`\ln y`$ vs $`x`$ | $`\beta_1`$ | $`\ln \alpha_1`$ |
| Power | $`y = \alpha_2 x^{\beta_2}`$ | $`\ln y`$ vs $`\ln x`$ | $`\beta_2`$ | $`\ln \alpha_2`$ |
| Saturation-Growth | $`y = \frac{\alpha_3 x}{\beta_3 + x}`$ | $`\frac{1}{y}`$ vs $`\frac{1}{x}`$ | $`\frac{\beta_3}{\alpha_3}`$ | $`\frac{1}{\alpha_3}`$ |

*Table 2. Transformed variables, slopes, and intercepts of the linearized models*

---

<br>

## 3. Polynomial Regression

### 3.1 Problem Formulation

When data shows a curvilinear pattern (e.g., parabolic), a linear line is not sufficient. We extend the linear regression idea to **polynomial regression**.

Confidence intervals for the parameters estimated by polynomial regression are obtained from the covariance matrix in Section 5.6.

Given a set of data $`\{x_i, y_i\}_{i=1}^{n}`$, find $`\beta_0, \beta_1, \beta_2`$ such that:

```math
\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2
```

The objective is to minimize the sum of squared errors:

```math
\min_{\beta_0, \beta_1, \beta_2} SSE
```

where:

```math
SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)^2
```

### 3.2 Deriving the Normal Equations

To determine $`\beta_0, \beta_1, \beta_2`$, take partial derivatives and set them to zero:

```math
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-1) = 0
```

```math
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i) = 0
```

```math
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i^2) = 0
```

After simplification, the **normal equations** become:

```math
\sum y_i = \beta_0 \cdot n + \beta_1 \sum x_i + \beta_2 \sum x_i^2
```

```math
\sum y_i x_i = \beta_0 \sum x_i + \beta_1 \sum x_i^2 + \beta_2 \sum x_i^3
```

```math
\sum y_i x_i^2 = \beta_0 \sum x_i^2 + \beta_1 \sum x_i^3 + \beta_2 \sum x_i^4
```

### 3.3 Normal Equations in Matrix Form

The normal equations can be written compactly as:

```math
\begin{bmatrix} \sum y_i \\ \sum y_i x_i \\ \sum y_i x_i^2 \end{bmatrix}
=
\begin{bmatrix} n & \sum x_i & \sum x_i^2 \\ \sum x_i & \sum x_i^2 & \sum x_i^3 \\ \sum x_i^2 & \sum x_i^3 & \sum x_i^4 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
```

> **Key Point:** This is a $`3 \times 3`$ linear system that can be solved with methods for linear algebraic equations such as Gauss elimination (Chapter 9) or LU factorization (Chapter 10).

### 3.4 Standard Error of the Estimate

For polynomial regression with 3 parameters ($`\beta_0, \beta_1, \beta_2`$):

```math
e_i = y_i - \hat{y}_i
```

```math
SSE = \sum e_i^2
```

```math
s_e = \sqrt{\frac{SSE}{n - 3}}
```

> **Note:** **Degrees of Freedom.** The denominator is $`n - 3`$ because three parameters ($`\beta_0, \beta_1, \beta_2`$) are used to construct $`\hat{y}`$. In general, for an $`m`$-th degree polynomial, the degrees of freedom are $`n - (m+1)`$.

---

<br>

## 4. Multiple Linear Regression

Multiple linear regression extends the concept from 1D to 2D (and beyond): fitting data with **two or more independent variables**.

### 4.1 Problem Formulation

The model fits a plane (or hyperplane) to the data:

```math
\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2
```

The objective is:

```math
SSE = \sum e_i^2 = \sum (y_i - \hat{y}_i)^2 = \sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})^2
```

### 4.2 Normal Equations

Seek $`\beta_0, \beta_1, \beta_2`$ such that $`\min_{\beta_0, \beta_1, \beta_2} SSE`$:

```math
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-1) = 0
```

```math
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{1i}) = 0
```

```math
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{2i}) = 0
```

In matrix form:

```math
\begin{bmatrix} \sum y_i \\ \sum y_i x_{1i} \\ \sum y_i x_{2i} \end{bmatrix}
=
\begin{bmatrix} n & \sum x_{1i} & \sum x_{2i} \\ \sum x_{1i} & \sum x_{1i}^2 & \sum x_{1i} x_{2i} \\ \sum x_{2i} & \sum x_{1i} x_{2i} & \sum x_{2i}^2 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
```

### 4.3 Standard Error of the Estimate

For multiple linear regression with two independent variables (3 parameters):

```math
s_e = \sqrt{\frac{SSE}{n - 3}}
```

### 4.4 Extension to Multivariable Power Models

Multiple linear regression can be extended to **multivariable power models**:

```math
y = \beta_0 x_1^{\beta_1} x_2^{\beta_2} \cdots x_m^{\beta_m} + \varepsilon
```

Ignoring the error term $`\varepsilon`$ and taking the natural logarithm:

```math
\ln y = \ln \beta_0 + \beta_1 \ln x_1 + \beta_2 \ln x_2 + \cdots + \beta_m \ln x_m
```

> **Key Point:** After the log transformation, this becomes a standard multiple linear regression problem in the transformed variables $`\ln y`$, $`\ln x_1`$, $`\ln x_2`$, ..., $`\ln x_m`$.

---

<br>

## 5. General Linear Least Square Regression

### 5.1 General Model Formulation

The general linear least square model unifies all previous cases. The model is linear in the parameters $`\beta_j`$, but the basis functions $`f_j`$ can be arbitrary (e.g., trigonometric):

**Example (Fourier-type model):**

```math
y = \beta_0 + \beta_1 \cos(\omega x) + \beta_2 \sin(\omega x)
```

**General form:**

```math
y = \beta_0 f_0(x_1, x_2, \ldots, x_p) + \beta_1 f_1(x_1, x_2, \ldots, x_p) + \cdots + \beta_m f_m(x_1, x_2, \ldots, x_p)
```

> **Key Point:** "Linear" here means linear in the **parameters** $`\beta_j`$, not necessarily in $`x`$. The functions $`f_j`$ can be nonlinear in $`x`$.

### 5.2 Matrix Formulation

Collect $`n`$ data points $`\{y_i, x_{1i}, x_{2i}, \ldots, x_{pi}\}_{i=1}^{n}`$.

For each data point $`i = 1, 2, \ldots, n`$:

```math
y_i = \beta_0 f_0(x_{1i}, x_{2i}, \ldots, x_{pi}) + \beta_1 f_1(x_{1i}, x_{2i}, \ldots, x_{pi}) + \cdots + \beta_m f_m(x_{1i}, x_{2i}, \ldots, x_{pi})
```

In matrix form:

```math
\underset{(n \times 1)}{\mathbf{y}}
=
\underset{(n \times (m+1))}{\mathbf{X}}
\underset{((m+1) \times 1)}{\boldsymbol{\beta}}
+ \underset{(n \times 1)}{\mathbf{e}}
```

where:

```math
\mathbf{y} = \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{bmatrix}, \quad
\mathbf{X} = \begin{bmatrix} f_{01} & f_{11} & f_{21} & \cdots & f_{m1} \\ f_{02} & f_{12} & f_{22} & \cdots & f_{m2} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ f_{0n} & f_{1n} & f_{2n} & \cdots & f_{mn} \end{bmatrix}, \quad
\boldsymbol{\beta} = \begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \\ \vdots \\ \beta_m \end{bmatrix}, \quad
\mathbf{e} = \begin{bmatrix} e_1 \\ e_2 \\ \vdots \\ e_n \end{bmatrix}
```

Here $`f_{ji}`$ denotes the value of basis function $`f_j`$ evaluated at data point $`i`$, and $`\mathbf{e}`$ is the residual (error) vector.

> **Note:** **Notation.** $`\mathbf{X}`$ has $`n`$ rows (data points) and $`m+1`$ columns (number of parameters including $`\beta_0`$).

### 5.3 Deriving the Normal Equations (Matrix Calculus)

The sum of squared errors in matrix form:

```math
SSE = \sum e_i^2 = \mathbf{e}^T \mathbf{e} = (\mathbf{y} - \hat{\mathbf{y}})^T (\mathbf{y} - \hat{\mathbf{y}})
```

To minimize, take the derivative with respect to $`\boldsymbol{\beta}`$ and set to zero:

```math
\frac{\partial SSE}{\partial \boldsymbol{\beta}} = \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{e}^T \mathbf{e})
= \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})
= (-2\mathbf{X})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta}) = \mathbf{0}
```

> **[Calculus]** The gradient of a multivariable function is the vector of its partial derivatives, and all partial derivatives of a differentiable function vanish at a local minimum. Since SSE is a quadratic function of $`\boldsymbol{\beta}`$, writing the componentwise partial derivatives of Section 3.2 as one vector gives the gradient above, and setting it to zero gives the normal equations.

This yields the **normal equation**:

```math
\mathbf{X}^T \mathbf{X} \boldsymbol{\beta} = \mathbf{X}^T \mathbf{y}
```

Solving for $`\boldsymbol{\beta}`$:

```math
\boxed{\boldsymbol{\beta} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}}
```

> **Key Point:** This is the most important result of the chapter. It provides a closed-form solution for the parameter vector $`\boldsymbol{\beta}`$ that minimizes the sum of squared errors, applicable to any general linear model.

### 5.4 Computing the Solution

Once $`\boldsymbol{\beta}`$ is computed:

| Quantity | Formula |
|---|---|
| Residual vector | $`\mathbf{e} = \mathbf{y} - \mathbf{X}\boldsymbol{\beta}`$ |
| Sum of squared errors | $`SSE = \mathbf{e}^T \mathbf{e}`$ |
| Predicted values | $`\hat{\mathbf{y}} = \mathbf{X}\boldsymbol{\beta}`$ |

*Table 3. Quantities computed after obtaining the parameter vector*

### 5.5 Goodness-of-Fit Statistics

```math
SST = \sum (y_i - \bar{y})^2
```

```math
SSR = \sum (\hat{y}_i - \bar{y})^2
```

```math
R^2 = \frac{SSR}{SST}
```

```math
SSE = SST - SSR
```

**Standard error of the estimate** (general case with $`m+1`$ parameters):

```math
s_e = \sqrt{\frac{SSE}{n - (m+1)}}
```

> **Note:** **Degrees of Freedom.** The denominator $`n - (m+1)`$ accounts for the $`m+1`$ parameters used to fit the model.

### 5.6 Covariance Matrix of Parameter Estimates

The covariance matrix of the estimated parameters:

```math
\text{cov}(\hat{\boldsymbol{\beta}}) = \sigma_e^2 (\mathbf{X}^T \mathbf{X})^{-1}
```

where $`\sigma_e^2`$ is estimated by $`s_e^2 = \frac{SSE}{n-(m+1)}`$.

> **Note:** **Practical Use.** The diagonal elements of $`\text{cov}(\hat{\boldsymbol{\beta}})`$ give the variance of each individual parameter estimate, and their square roots provide the standard errors of the parameters. These can be used to construct confidence intervals.

---

<br>

## 6. Nonlinear Regression (Gauss-Newton)

### 6.1 Motivation: The Antoine Equation

Some models are **inherently nonlinear** in the parameters and cannot be linearized by transformation. For example, the **Antoine Equation** for vapor pressure:

```math
\log_{10} P_v = A - \frac{B}{C + T}
```

where:
- $`P_v`$ = vapor pressure
- $`T`$ = temperature
- $`A, B, C`$ = adjustable parameters

Since $`C`$ appears inside a nonlinear expression ($`C + T`$ in the denominator), this equation **cannot** be expressed as a general linear model. We need nonlinear regression (optimization) to fit it.

> **Note:** **Gauss-Newton Method.** The Gauss-Newton method, a technique designed specifically for nonlinear regression, linearizes the model $`f(x_i; a_0, \ldots, a_m)`$ around the current parameter estimates with a first-order Taylor series. This gives a linear relation between the residual vector $`\mathbf{d}`$ ($`d_i = y_i - f(x_i)`$) and the parameter correction $`\Delta\mathbf{a}`$ through the matrix $`\mathbf{Z}`$ whose elements are the partial derivatives $`\partial f / \partial a_k`$, so the correction is found by solving the general linear least-squares normal equations $`\mathbf{Z}^T \mathbf{Z} \, \Delta\mathbf{a} = \mathbf{Z}^T \mathbf{d}`$. The parameters are corrected and the linearization is repeated until the correction becomes sufficiently small. The code in Section 6.2 instead minimizes the SSE directly with a general-purpose optimizer.

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

> **Note:** **Method.** The `scipy.optimize.minimize` function with `method='L-BFGS-B'` uses a quasi-Newton optimization algorithm. An initial guess is required, and the algorithm iteratively refines the parameters to minimize $`SSE`$.

> **Caution:** **Initial Guess Sensitivity.** Nonlinear regression methods depend on the initial guess. A poor initial guess may lead to convergence to a local minimum rather than the global minimum. In the code above, the parameters also differ greatly in scale, and L-BFGS-B declares convergence at a point where $`B`$ has stayed near its initial value of 1500. From the same initial guess, the Levenberg-Marquardt method, an improvement of the Gauss-Newton method and the default of `scipy.optimize.curve_fit`, reaches $`A \approx 8.01`$, $`B \approx 1697`$, $`C \approx 230.6`$ with a smaller SSE.

---

<br>

## Summary

| Concept | Key Summary |
|---|---|
| Linear regression | The model $`\hat{y} = \beta_0 + \beta_1 x`$ has 2 parameters and is solved in closed form (Chapter 14), and its standard error of the estimate is $`s_e = \sqrt{\frac{SSE}{n-2}}`$. |
| Polynomial regression | The model $`\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2`$ has 3 parameters and is solved with $`3 \times 3`$ normal equations, and $`s_e = \sqrt{\frac{SSE}{n-3}}`$. |
| Multiple linear regression | The model $`\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2`$ has 3 parameters and is solved with $`3 \times 3`$ normal equations, and $`s_e = \sqrt{\frac{SSE}{n-3}}`$. |
| General linear least squares | The model $`\hat{y} = \sum \beta_j f_j(\mathbf{x})`$ has $`m+1`$ parameters and is solved by $`\boldsymbol{\beta} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\mathbf{y}`$, and $`s_e = \sqrt{\frac{SSE}{n-(m+1)}}`$. |
| Linearized models | The exponential, power, and saturation-growth models have 2 parameters and are solved by transforming the variables and then applying linear regression, and their $`s_e`$ depends on the transform. |
| Nonlinear regression | Models such as the Antoine equation have a model-dependent number of parameters and are solved by iterative optimization (e.g., L-BFGS-B), and with $`p`$ parameters $`s_e = \sqrt{\frac{SSE}{n-p}}`$. |

---

<br>

## Self-Check Questions

1. **Standard error of regression and coefficient of determination:** For the example in Section 1.3 (heights = [170, 172, 168, 174, 176] cm, $`x = [1, 2, 3, 4, 5]`$), find the least-squares line, $`SSE`$, $`s_e`$, and $`R^2`$.

   <details>
   <summary>Answer</summary>

   With $`\bar{x} = 3`$, $`\bar{y} = 172`$, $`\sum (x_i - \bar{x})(y_i - \bar{y}) = 14`$, and $`\sum (x_i - \bar{x})^2 = 10`$, the slope is $`\beta_1 = 1.4`$ and the intercept is $`\beta_0 = 172 - 1.4 \cdot 3 = 167.8`$. The predicted values are 169.2, 170.6, 172.0, 173.4, 174.8 and the residuals are 0.8, 1.4, -4.0, 0.6, 1.2, so $`SSE = 20.4`$. Therefore $`s_e = \sqrt{20.4/3} \approx 2.61`$ cm. Since $`SST = 40`$, $`R^2 = 1 - 20.4/40 = 0.49`$, so the line explains only 49% of the variation in height.

   </details>

2. **Linearizing the saturation-growth-rate model:** A linear regression of $`1/y`$ on $`1/x`$ gives slope 0.5 and intercept 0.2. Find $`\alpha_3`$ and $`\beta_3`$ of the saturation-growth-rate model $`y = \frac{\alpha_3 x}{\beta_3 + x}`$.

   <details>
   <summary>Answer</summary>

   The intercept is $`1/\alpha_3 = 0.2`$, so $`\alpha_3 = 5`$. The slope is $`\beta_3/\alpha_3 = 0.5`$, so $`\beta_3 = 0.5 \cdot 5 = 2.5`$. The model is therefore $`y = \frac{5x}{2.5 + x}`$, and $`y`$ approaches the saturation value 5 as $`x`$ grows.

   </details>

3. **Degrees of freedom:** A second-degree polynomial fitted to 10 data points gives $`SSE = 3.5`$. Find the standard error of the estimate, and explain how the value would change if the same $`SSE`$ came from a straight-line fit.

   <details>
   <summary>Answer</summary>

   A second-degree polynomial has 3 parameters, so the degrees of freedom are $`n - 3 = 7`$ and $`s_e = \sqrt{3.5/7} = \sqrt{0.5} \approx 0.707`$. A straight line has 2 parameters, so the degrees of freedom are $`n - 2 = 8`$ and $`s_e = \sqrt{3.5/8} \approx 0.661`$. More parameters leave fewer degrees of freedom, so the same $`SSE`$ gives a larger $`s_e`$.

   </details>

4. **Linearity in the parameters:** Decide whether each model can be fitted directly by general linear least squares, fitted by linear regression after a transformation, or requires nonlinear regression. (a) $`y = \beta_0 + \beta_1 \cos(\omega x) + \beta_2 \sin(\omega x)`$ with known $`\omega`$ (b) $`y = \alpha_1 e^{\beta_1 x}`$ (c) $`\log_{10} P_v = A - \frac{B}{C + T}`$

   <details>
   <summary>Answer</summary>

   In (a) the parameters are the coefficients of the basis functions $`1`$, $`\cos(\omega x)`$, and $`\sin(\omega x)`$, and the model is linear in them, so it is fitted directly by general linear least squares. Model (b) is nonlinear in its parameters, but the transformation $`\ln y = \ln \alpha_1 + \beta_1 x`$ allows a linear regression fit. In (c), $`C`$ sits inside $`C + T`$ in the denominator, so $`\log_{10} P_v`$ cannot be written as a linear combination of the parameters, and nonlinear regression is required.

   </details>
