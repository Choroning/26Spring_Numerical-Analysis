# Chapter 15 — Nonlinear Regression

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-11-19
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 15

> **Learning Objectives:**
> 1. Linearize exponential, power, and saturation-growth-rate models and fit them by linear regression
> 2. Derive and apply the normal equations of polynomial and multiple linear regression
> 3. Derive the normal equations of the general linear least-squares model in matrix form and compute goodness-of-fit statistics
> 4. Fit models that are inherently nonlinear in their parameters by iterative optimization
> 5. Explain why the condition number of the normal equations is squared and solve least-squares problems by QR factorization
> 6. Derive the Gauss-Newton method, iterate it by hand, and perform nonlinear regression with SciPy's `curve_fit` and `least_squares`

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
  - [3.5 Polynomial Regression Example](#35-polynomial-regression-example)
- [4. Multiple Linear Regression](#4-multiple-linear-regression)
  - [4.1 Problem Formulation](#41-problem-formulation)
  - [4.2 Normal Equations](#42-normal-equations)
  - [4.3 Standard Error of the Estimate](#43-standard-error-of-the-estimate)
  - [4.4 Extension to Multivariable Power Models](#44-extension-to-multivariable-power-models)
  - [4.5 Multiple Linear Regression Example](#45-multiple-linear-regression-example)
- [5. General Linear Least Square Regression](#5-general-linear-least-square-regression)
  - [5.1 General Model Formulation](#51-general-model-formulation)
  - [5.2 Matrix Formulation](#52-matrix-formulation)
  - [5.3 Deriving the Normal Equations (Matrix Calculus)](#53-deriving-the-normal-equations-matrix-calculus)
  - [5.4 Computing the Solution](#54-computing-the-solution)
  - [5.5 Goodness-of-Fit Statistics](#55-goodness-of-fit-statistics)
  - [5.6 Covariance Matrix of Parameter Estimates](#56-covariance-matrix-of-parameter-estimates)
  - [5.7 Conditioning of the Normal Equations and QR Factorization](#57-conditioning-of-the-normal-equations-and-qr-factorization)
  - [5.8 Model Selection: Adjusted and Predicted Coefficients of Determination](#58-model-selection-adjusted-and-predicted-coefficients-of-determination)
  - [5.9 General Linear Least-Squares Example with Arbitrary Basis Functions](#59-general-linear-least-squares-example-with-arbitrary-basis-functions)
- [6. Nonlinear Regression (Gauss-Newton)](#6-nonlinear-regression-gauss-newton)
  - [6.1 Motivation: The Antoine Equation](#61-motivation-the-antoine-equation)
  - [6.2 Python Implementation with scipy.optimize](#62-python-implementation-with-scipyoptimize)
  - [6.3 Deriving the Gauss-Newton Method](#63-deriving-the-gauss-newton-method)
  - [6.4 Gauss-Newton Example](#64-gauss-newton-example)
  - [6.5 Linearized Versus Nonlinear Fitting](#65-linearized-versus-nonlinear-fitting)
  - [6.6 curve_fit and least_squares](#66-curve_fit-and-least_squares)
  - [6.7 Minimizing the Sum of Squares with minimize](#67-minimizing-the-sum-of-squares-with-minimize)
  - [6.8 Initial Estimates and Convergence Failure](#68-initial-estimates-and-convergence-failure)
- [7. Common Mistakes](#7-common-mistakes)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)
- [References](#references)

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

**Interpretation:** Individual height measurements differ from the mean by roughly 3.16 cm.

### 1.2 Standard Error of the Mean (SE)

**Standard Error** (SE) measures how precisely the sample mean estimates the true population mean.

```math
SE = \frac{s_y}{\sqrt{n}}
```

**Example (continued):**
- $`SE = \frac{3.16}{\sqrt{5}} = 1.41`$ cm

**Interpretation:** The average height is estimated to be 172 cm, with a standard error of $`\pm 1.41`$ cm.

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

**Interpretation:** The observed $`y`$-values deviate from the regression line by about 2.61 cm.

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

**Interpretation:** The relationship between the two variables explains 88% (approximately 90%) of the variation in the data.

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

**Graphical Effect:** A plot of $`y`$ vs $`x`$ shows an exponential curve. A plot of $`\ln y`$ vs $`x`$ becomes a straight line with slope $`\beta_1`$ and intercept $`\ln \alpha_1`$.

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

**Graphical Effect:** A plot of $`y`$ vs $`x`$ shows a power curve. A plot of $`\ln y`$ vs $`\ln x`$ becomes a straight line with slope $`\beta_2`$ and intercept $`\ln \alpha_2`$.

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

**Graphical Effect:** A plot of $`y`$ vs $`x`$ shows a saturation curve. A plot of $`\frac{1}{y}`$ vs $`\frac{1}{x}`$ becomes a straight line with slope $`\frac{\beta_3}{\alpha_3}`$ and intercept $`\frac{1}{\alpha_3}`$.

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

**Degrees of Freedom:** The denominator is $`n - 3`$ because three parameters ($`\beta_0, \beta_1, \beta_2`$) are used to construct $`\hat{y}`$. In general, for an $`m`$-th degree polynomial, the degrees of freedom are $`n - (m+1)`$.

### 3.5 Polynomial Regression Example

A quadratic $`y = \beta_0 + \beta_1x + \beta_2x^2`$ is fitted to the following six data points.

```math
x = (-2,\ -1,\ 0,\ 1,\ 2,\ 3), \qquad y = (5.2,\ 1.9,\ 1.1,\ 2.0,\ 5.1,\ 9.8)
```

The sums needed for the normal equations are $`n = 6`$, $`\sum x_i = 3`$, $`\sum x_i^2 = 19`$, $`\sum x_i^3 = 27`$, $`\sum x_i^4 = 115`$, $`\sum y_i = 25.1`$, $`\sum x_iy_i = 29.3`$, and $`\sum x_i^2y_i = 133.3`$. Substituting into the matrix form of Section 3.3 gives the following.

```math
\begin{bmatrix} 6 & 3 & 19 \\ 3 & 19 & 27 \\ 19 & 27 & 115 \end{bmatrix}\begin{Bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{Bmatrix} = \begin{Bmatrix} 25.1 \\ 29.3 \\ 133.3 \end{Bmatrix}
```

Solving this symmetric system by Gauss elimination (Chapter 9) or Cholesky factorization (Chapter 10) gives $`\beta_0 = 1.057143`$, $`\beta_1 = -0.035714`$, and $`\beta_2 = 0.992857`$, so the fitted curve is $`y = 1.0571 - 0.0357x + 0.9929x^2`$. The residuals $`(0.100,\ -0.186,\ 0.043,\ -0.014,\ 0.143,\ -0.086)`$ are small and have mixed signs. The goodness-of-fit statistics are the following.

```math
SSE = 0.074286, \quad SST = 52.908333, \quad R^2 = 1 - \frac{0.074286}{52.908333} = 0.998596, \quad s_e = \sqrt{\frac{0.074286}{6 - 3}} = 0.157359
```

The denominator of the standard error is $`n - (m + 1) = 6 - 3 = 3`$ because three coefficients have been estimated (Section 3.4). The entries of the coefficient matrix have the form $`\sum x_i^{j+k}`$, so as the degree rises their sizes differ more and more sharply, which is the source of the conditioning problem treated in Section 5.7.

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

### 4.5 Multiple Linear Regression Example

$`y = \beta_0 + \beta_1x_1 + \beta_2x_2`$ is fitted to the following data in two independent variables $`x_1`$ and $`x_2`$.

| $`x_1`$ | 1 | 2 | 3 | 4 | 5 | 6 |
|:-----|-----:|-----:|-----:|-----:|-----:|-----:|
| $`x_2`$ | 2 | 1 | 4 | 3 | 6 | 5 |
| $`y`$ | 4.10 | 5.30 | 9.15 | 10.40 | 14.05 | 15.50 |

*Table 3. Data for the multiple linear regression example*

The required sums are $`n = 6`$, $`\sum x_1 = 21`$, $`\sum x_2 = 21`$, $`\sum x_1^2 = 91`$, $`\sum x_2^2 = 91`$, $`\sum x_1x_2 = 88`$, $`\sum y = 58.5`$, $`\sum x_1y = 247`$, and $`\sum x_2y = 243.1`$. Substituting into the normal equations of Section 4.2 gives the following.

```math
\begin{bmatrix} 6 & 21 & 21 \\ 21 & 91 & 88 \\ 21 & 88 & 91 \end{bmatrix}\begin{Bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{Bmatrix} = \begin{Bmatrix} 58.5 \\ 247 \\ 243.1 \end{Bmatrix}
```

Solving gives $`\beta_0 = 0.934375`$, $`\beta_1 = 1.909375`$, and $`\beta_2 = 0.609375`$, so the fitted plane is $`y = 0.9344 + 1.9094x_1 + 0.6094x_2`$. $`\beta_1`$ means that $`y`$ increases by about 1.91 when $`x_1`$ increases by 1 with $`x_2`$ held fixed, and it is interpreted as a partial effect with the influence of the other variable removed. The residuals are $`(0.0375,\ -0.0625,\ 0.05,\ 0,\ -0.0875,\ 0.0625)`$, and the goodness of fit is the following.

```math
SSE = 0.019375, \quad SST = 104.06, \quad R^2 = 0.999814, \quad s_e = \sqrt{\frac{0.019375}{6 - 3}} = 0.080364
```

Taking square roots of the diagonal entries of the covariance matrix $`s_e^2(\mathbf{X}^T\mathbf{X})^{-1}`$ of Section 5.6 gives the standard errors of the coefficients, $`0.0776`$, $`0.0343`$, and $`0.0343`$. The standard errors are very small compared with the coefficients, so both variables have a meaningful effect on $`y`$. If the two variables are strongly correlated (for example, if $`x_2`$ is nearly a multiple of $`x_1`$), $`\mathbf{X}^T\mathbf{X}`$ becomes nearly singular, the standard errors grow, and the coefficients become hard to interpret separately.

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

**Notation:** $`\mathbf{X}`$ has $`n`$ rows (data points) and $`m+1`$ columns (number of parameters including $`\beta_0`$).

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

*Table 4. Quantities computed after obtaining the parameter vector*

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

**Degrees of Freedom:** The denominator $`n - (m+1)`$ accounts for the $`m+1`$ parameters used to fit the model.

### 5.6 Covariance Matrix of Parameter Estimates

The covariance matrix of the estimated parameters:

```math
\text{cov}(\hat{\boldsymbol{\beta}}) = \sigma_e^2 (\mathbf{X}^T \mathbf{X})^{-1}
```

where $`\sigma_e^2`$ is estimated by $`s_e^2 = \frac{SSE}{n-(m+1)}`$.

**Practical Use:** The diagonal elements of $`\text{cov}(\hat{\boldsymbol{\beta}})`$ give the variance of each individual parameter estimate, and their square roots provide the standard errors of the parameters. These can be used to construct confidence intervals.

### 5.7 Conditioning of the Normal Equations and QR Factorization

The normal equations $`\mathbf{X}^T\mathbf{X}\,\boldsymbol{\beta} = \mathbf{X}^T\mathbf{y}`$ are simple to derive but have a numerical weakness. If the singular values of $`\mathbf{X}`$ are $`\sigma_1 \geq \cdots \geq \sigma_n`$, the eigenvalues of $`\mathbf{X}^T\mathbf{X}`$ are $`\sigma_i^2`$, so the 2-norm condition number of Chapter 11 is squared.

```math
\kappa_2(\mathbf{X}^T\mathbf{X}) = \frac{\sigma_1^2}{\sigma_n^2} = \kappa_2(\mathbf{X})^2
```

By the precision rule of Chapter 11, a condition number of $`10^{k}`$ costs about $`k`$ digits, so forming the normal equations doubles the number of digits lost. The standard way to avoid this is **QR factorization**. Writing $`\mathbf{X}`$ ($`n \times (m+1)`$) as the product $`\mathbf{X} = \mathbf{Q}\mathbf{R}`$ of a matrix $`\mathbf{Q}`$ with orthonormal columns and an upper triangular matrix $`\mathbf{R}`$, we have $`\mathbf{Q}^T\mathbf{Q} = I`$, so the normal equations become the following.

```math
\mathbf{R}^T\mathbf{Q}^T\mathbf{Q}\mathbf{R}\,\boldsymbol{\beta} = \mathbf{R}^T\mathbf{Q}^T\mathbf{y} \quad\Longrightarrow\quad \mathbf{R}\,\boldsymbol{\beta} = \mathbf{Q}^T\mathbf{y}
```

The last system is solved by a single back substitution, and since $`\kappa_2(\mathbf{R}) = \kappa_2(\mathbf{X})`$, the condition number is not squared. The following code takes data generated at 20 points of $`x \in [0, 1]`$ from a degree-7 polynomial whose coefficients are all 1, refits it by three methods, and compares the results with the known true coefficients.

```python
import numpy as np

x = np.linspace(0, 1, 20)
X = np.vander(x, 8, increasing=True)     # columns 1, x, ..., x^7
a_true = np.ones(8)
y = X @ a_true                           # data from a known polynomial

print(f"cond(X) = {np.linalg.cond(X):.2e}, cond(X^T X) = {np.linalg.cond(X.T @ X):.2e}")

a_ne = np.linalg.solve(X.T @ X, X.T @ y)         # normal equations
Q, R = np.linalg.qr(X)                            # reduced QR: X = Q R
a_qr = np.linalg.solve(R, Q.T @ y)                # R a = Q^T y
a_ls = np.linalg.lstsq(X, y, rcond=None)[0]       # library least squares

for name, a in (("normal eq.", a_ne), ("QR", a_qr), ("lstsq", a_ls)):
    print(f"{name:10s} max coefficient error = {np.abs(a - a_true).max():.1e}")
```

```text
cond(X) = 1.05e+05, cond(X^T X) = 1.11e+10
normal eq. max coefficient error = 3.1e-07
QR         max coefficient error = 6.7e-12
lstsq      max coefficient error = 5.7e-13
```

$`\kappa(\mathbf{X}) \approx 10^{5}`$ while $`\kappa(\mathbf{X}^T\mathbf{X}) \approx 10^{10}`$, exactly the square. Even though the data are error-free, the coefficients from the normal equations carry an error of about $`3 \times 10^{-7}`$, losing about 10 of the 16 double-precision digits, whereas QR factorization and `np.linalg.lstsq` are accurate to the $`10^{-12}`$ level. `np.linalg.lstsq` takes $`\mathbf{X}`$ directly and returns the least-squares solution; it treats singular values that are small relative to `rcond` as zero in determining the coefficients and also returns the singular values [1]. If the polynomial degree is low and the data are well spread, the normal equations suffice, but if the degree is high or the columns are nearly linearly dependent, QR factorization or `lstsq` is used.

> **Note:** **Centering the variable.** Using $`x - \bar{x}`$ or a variable mapped to $`[-1, 1]`$ instead of $`x`$ reduces the linear dependence among the columns of $`\mathbf{X}`$ and makes the condition number much smaller. This is why data are centered and scaled first in polynomial regression.

### 5.8 Model Selection: Adjusted and Predicted Coefficients of Determination

Adding terms to a model, such as raising the polynomial degree, can only decrease $`SSE`$, so $`R^2`$ always increases. $`R^2`$ alone therefore cannot decide how many terms to use. The textbook advises judging a model by both performance and adequacy and introduces the following two performance criteria.

The **adjusted $`R^2`$** reflects the cost of adding terms by dividing $`SSE`$ and $`SST`$ by their degrees of freedom ($`m + 1`$ is the number of coefficients).

```math
R^2_{adj} = 1 - \frac{SSE/(n - (m + 1))}{SST/(n - 1)}
```

The **predicted $`R^2`$** measures performance by how well a model fitted without each point, in turn, predicts that point. The sum of squares of these left-out prediction residuals is called PRESS, and it can be computed without refitting $`n`$ times from the diagonal entries $`h_{ii}`$ of the hat matrix $`\mathbf{H} = \mathbf{X}(\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T`$.

```math
PRESS = \sum_{i=1}^{n}\left(\frac{e_i}{1 - h_{ii}}\right)^2, \qquad R^2_{pred} = 1 - \frac{PRESS}{SST}
```

An overfitted model predicts the left-out points poorly, so its $`R^2_{pred}`$ actually decreases. The following code fits polynomials of degree 1 through 5 to 10 data points and compares the three measures. To avoid the conditioning problem of Section 5.7, it uses the standardized variable $`z = (x - \bar{x})/s_x`$.

```python
import numpy as np

x = np.linspace(0, 3, 10)
y = np.array([1.05, 1.32, 1.98, 2.83, 4.05, 5.62, 7.31, 9.48, 11.80, 14.65])
n = len(x)
z = (x - x.mean()) / x.std(ddof=1)             # standardized variable
sst = np.sum((y - y.mean())**2)
for m in range(1, 6):                          # polynomial degree
    X = np.vander(z, m + 1, increasing=True)   # columns 1, z, ..., z^m
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ beta
    sse = np.sum(e**2)
    H = X @ np.linalg.solve(X.T @ X, X.T)      # hat matrix
    press = np.sum((e / (1 - np.diag(H)))**2)
    r2 = 1 - sse / sst
    r2adj = 1 - (sse / (n - (m + 1))) / (sst / (n - 1))
    r2pred = 1 - press / sst
    print(f"{m}  {r2:.5f}  {r2adj:.5f}  {r2pred:.5f}")
```

```text
1  0.93931  0.93173  0.88495
2  0.99990  0.99988  0.99975
3  0.99992  0.99987  0.99957
4  0.99994  0.99989  0.99959
5  0.99995  0.99988  0.99869
```

The columns are, in order, the degree, $`R^2`$, $`R^2_{adj}`$, and $`R^2_{pred}`$. $`R^2`$ grows slightly with each degree, but apart from the step from degree 1 to 2 the increase shows only in the fifth decimal place. $`R^2_{adj}`$ barely changes after degree 2, and $`R^2_{pred}`$ is largest at degree 2 and falls to 0.99869 at degree 5. The quadratic model is therefore the appropriate choice. As the textbook stresses, after choosing by performance measures, the adequacy of the model must also be checked with the residual plot of Chapter 14, Section 5.6.

### 5.9 General Linear Least-Squares Example with Arbitrary Basis Functions

As the textbook stresses, the "linear" in general linear least squares refers to linearity in the coefficients, so the columns of $`\mathbf{X}`$ can hold not only powers of $`x`$ but any functions $`z_j(x)`$. The model $`y = \beta_0 + \beta_1x + \beta_2e^{-x}`$, which has a constant term, a linear trend, and an exponentially vanishing initial component, is fitted to the following data.

```math
x = (0,\ 1,\ 2,\ 3,\ 4,\ 5), \qquad y = (5.050,\ 4.196,\ 4.301,\ 4.580,\ 5.077,\ 5.484)
```

The three columns of $`\mathbf{X}`$ are $`1`$, $`x_i`$, and $`e^{-x_i}`$, and the coefficient matrix and right-hand side of the normal equations are the following.

```math
\mathbf{X}^T\mathbf{X} = \begin{bmatrix} 6 & 15 & 1.5781 \\ 15 & 55 & 0.8949 \\ 1.5781 & 0.8949 & 1.1565 \end{bmatrix}, \qquad \mathbf{X}^T\mathbf{y} = \begin{Bmatrix} 28.688 \\ 74.266 \\ 7.5337 \end{Bmatrix}
```

```python
import numpy as np

x = np.array([0, 1, 2, 3, 4, 5], dtype=float)
y = np.array([5.050, 4.196, 4.301, 4.580, 5.077, 5.484])
X = np.column_stack([np.ones_like(x), x, np.exp(-x)])   # basis functions 1, x, e^(-x)
beta = np.linalg.solve(X.T @ X, X.T @ y)
e = y - X @ beta
n, p = X.shape
se2 = np.sum(e**2) / (n - p)
std = np.sqrt(np.diag(se2 * np.linalg.inv(X.T @ X)))
t = 3.182                                                 # t value for 95% with n - p = 3
for b, s in zip(beta, std):
    print(f"{b:.4f}  +/- {t * s:.4f}")
print(f"R2 = {1 - np.sum(e**2) / np.sum((y - y.mean())**2):.5f}")
```

```text
2.9780  +/- 0.2470
0.5046  +/- 0.0657
2.0602  +/- 0.3194
R2 = 0.99500
```

The code stacks the values of the basis functions as columns to build $`\mathbf{X}`$, solves the normal equations, obtains the standard errors from the covariance matrix of Section 5.6, and computes the half-widths of the 95% confidence intervals with $`t_{0.025,3} = 3.182`$ for $`n - p = 3`$ degrees of freedom. The result is $`y = 2.978 + 0.505x + 2.060e^{-x}`$, and the intervals of the three coefficients, $`[2.731,\ 3.225]`$, $`[0.439,\ 0.570]`$, and $`[1.741,\ 2.380]`$, all exclude 0, so all three terms are needed. These data were in fact generated from $`3 + 0.5x + 2e^{-x}`$ plus small errors, and all the true values lie inside the intervals. Even with nonlinear basis functions, a model that is linear in its coefficients is solved in one pass without iteration, which is the difference from the nonlinear regression of Section 6.

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

> **Note:** **The code of Section 6.2 and the Gauss-Newton method.** The code of Section 6.2 minimizes the SSE directly with a general-purpose optimizer instead of a technique dedicated to nonlinear regression. The Gauss-Newton method, which linearizes the model by a first-order Taylor series and repeats linear least squares, is derived in Section 6.3.

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

**Method:** The `scipy.optimize.minimize` function with `method='L-BFGS-B'` uses a quasi-Newton optimization algorithm. An initial guess is required, and the algorithm iteratively refines the parameters to minimize $`SSE`$.

> **Caution:** **Initial Guess Sensitivity.** Nonlinear regression methods depend on the initial guess. A poor initial guess may lead to convergence to a local minimum rather than the global minimum. In the code above, the parameters also differ greatly in scale, and L-BFGS-B declares convergence at a point where $`B`$ has stayed near its initial value of 1500. From the same initial guess, the Levenberg-Marquardt method, an improvement of the Gauss-Newton method and the default of `scipy.optimize.curve_fit`, reaches $`A \approx 8.01`$, $`B \approx 1697`$, $`C \approx 230.6`$ with a smaller SSE.

### 6.3 Deriving the Gauss-Newton Method

Fitting a model $`f(x; \mathbf{a})`$ that is nonlinear in its parameters $`a_0, \ldots, a_m`$ to data $`(x_i, y_i)`$ requires minimizing $`SSE(\mathbf{a}) = \sum_i\bigl(y_i - f(x_i; \mathbf{a})\bigr)^2`$. Since this cannot be solved in one pass like the normal equations, the model is linearized by a first-order Taylor series near the current estimate $`\mathbf{a}^{(k)}`$.

```math
f(x_i; \mathbf{a}^{(k)} + \Delta\mathbf{a}) \approx f(x_i; \mathbf{a}^{(k)}) + \sum_{j=0}^{m}\frac{\partial f(x_i; \mathbf{a}^{(k)})}{\partial a_j}\,\Delta a_j
```

With the matrix of partial derivatives $`\mathbf{Z}_k`$ (whose $`(i, j)`$ entry is $`\partial f(x_i)/\partial a_j`$) and the current residuals $`d_i = y_i - f(x_i; \mathbf{a}^{(k)})`$, the approximate residual is $`\mathbf{d} - \mathbf{Z}_k\Delta\mathbf{a}`$. This is a least-squares problem linear in $`\Delta\mathbf{a}`$, so the normal equations of Section 5.3 apply directly.

```math
\mathbf{Z}_k^T\mathbf{Z}_k\,\Delta\mathbf{a} = \mathbf{Z}_k^T\mathbf{d}, \qquad \mathbf{a}^{(k+1)} = \mathbf{a}^{(k)} + \Delta\mathbf{a}
```

The iteration is repeated until the correction $`\Delta\mathbf{a}`$ becomes small enough. If the model is linear in its parameters, $`\mathbf{Z}`$ does not depend on $`\mathbf{a}`$ and equals the $`\mathbf{X}`$ of Section 5, so a single iteration reaches the general linear least-squares solution of Section 5. That is, the Gauss-Newton method applies general linear least squares repeatedly, the same idea as the Newton-Raphson method for systems of nonlinear equations in Chapter 12. Near the solution, with small residuals, it converges quickly, but it can diverge if the initial estimate is poor or $`\mathbf{Z}^T\mathbf{Z}`$ is nearly singular. The Levenberg-Marquardt method remedies this weakness by adding $`\mu I`$ to the diagonal of the normal equations and solving $`(\mathbf{Z}^T\mathbf{Z} + \mu I)\Delta\mathbf{a} = \mathbf{Z}^T\mathbf{d}`$, adjusting $`\mu`$ to move between the Gauss-Newton method and steepest descent (Chapter 7).

### 6.4 Gauss-Newton Example

The exponential model $`f(x; a_0, a_1) = a_0e^{a_1x}`$ is fitted by nonlinear regression to the following data.

```math
x = (0,\ 1,\ 2,\ 3,\ 4), \qquad y = (1.9,\ 3.1,\ 4.8,\ 8.3,\ 12.6)
```

The partial derivatives are $`\partial f/\partial a_0 = e^{a_1x}`$ and $`\partial f/\partial a_1 = a_0xe^{a_1x}`$. At the initial estimate $`\mathbf{a}^{(0)} = (2,\ 0.4)`$, $`\mathbf{Z}_0`$ and the residuals are the following.

```math
\mathbf{Z}_0 = \begin{bmatrix} 1 & 0 \\ 1.4918 & 2.9836 \\ 2.2255 & 8.9022 \\ 3.3201 & 19.9207 \\ 4.9530 & 39.6243 \end{bmatrix}, \qquad \mathbf{d} = \begin{Bmatrix} -0.1000 \\ 0.1164 \\ 0.3489 \\ 1.6598 \\ 2.6939 \end{Bmatrix}
```

The normal equations are the following.

```math
\begin{bmatrix} 43.7343 & 286.6625 \\ 286.6625 & 2055.0670 \end{bmatrix}\begin{Bmatrix} \Delta a_0 \\ \Delta a_1 \end{Bmatrix} = \begin{Bmatrix} 19.7039 \\ 143.2622 \end{Bmatrix}
```

Solving gives $`\Delta\mathbf{a} = (-0.074676,\ 0.080128)`$, so $`\mathbf{a}^{(1)} = (1.925324,\ 0.480128)`$. The code that repeats the same process and its results are as follows.

```python
import numpy as np

x = np.array([0, 1, 2, 3, 4], dtype=float)
y = np.array([1.9, 3.1, 4.8, 8.3, 12.6])

def model(a, x):
    return a[0] * np.exp(a[1] * x)

def jacobian(a, x):
    # columns: partial derivatives with respect to a0 and a1
    return np.column_stack([np.exp(a[1] * x), a[0] * x * np.exp(a[1] * x)])

a = np.array([2.0, 0.4])                 # initial guess
for k in range(1, 6):
    Z = jacobian(a, x)
    d = y - model(a, x)                  # residuals at the current guess
    da = np.linalg.solve(Z.T @ Z, Z.T @ d)
    a = a + da
    sse = np.sum((y - model(a, x))**2)
    print(k, np.round(a, 6), f"{sse:.6f}", f"{np.abs(da).max():.1e}")
```

```text
1 [1.925324 0.480128] 0.373488 8.0e-02
2 [1.960323 0.467503] 0.165345 3.5e-02
3 [1.959738 0.467527] 0.165333 5.8e-04
4 [1.959741 0.467527] 0.165333 2.8e-06
5 [1.959741 0.467527] 0.165333 4.8e-08
```

| Iteration | $`a_0`$ | $`a_1`$ | $`SSE`$ | $`\max\vert\Delta a_j\vert`$ |
|-----:|-----:|-----:|-----:|-----:|
| 1 | 1.925324 | 0.480128 | 0.373488 | $`8.0 \times 10^{-2}`$ |
| 2 | 1.960323 | 0.467503 | 0.165345 | $`3.5 \times 10^{-2}`$ |
| 3 | 1.959738 | 0.467527 | 0.165333 | $`5.8 \times 10^{-4}`$ |
| 4 | 1.959741 | 0.467527 | 0.165333 | $`2.8 \times 10^{-6}`$ |

*Table 5. Gauss-Newton iterations for the exponential model*

The `jacobian` function builds the two columns of $`\mathbf{Z}_k`$, the `model` function builds the current predictions, and the loop performs the two formulas of Section 6.3 directly. The largest absolute correction shrinks as $`8.0 \times 10^{-2}`$, $`3.5 \times 10^{-2}`$, $`5.8 \times 10^{-4}`$, $`2.8 \times 10^{-6}`$, converging to six decimal places in four iterations. The final model is $`y = 1.9597e^{0.4675x}`$.

### 6.5 Linearized Versus Nonlinear Fitting

The same exponential model can also be linearized as in Section 2.1, $`\ln y = \ln a_0 + a_1x`$, and fitted by straight-line regression. Comparing the two methods by the $`SSE`$ computed in the original $`y`$ space gives the following.

| Method | $`a_0`$ | $`a_1`$ | $`SSE`$ in the original space |
|:-----|-----:|-----:|-----:|
| Linearized straight-line regression | 1.905333 | 0.476854 | 0.187935 |
| Gauss-Newton nonlinear regression | 1.959741 | 0.467527 | 0.165333 |

*Table 6. Coefficients and SSE of linearized and nonlinear fits*

The two results differ because they minimize different quantities. The linearized fit minimizes $`\sum(\ln y_i - \ln\hat{y}_i)^2`$, a quantity close to the relative error, so it weighs the errors of small values relatively heavily. The nonlinear fit directly minimizes $`\sum(y_i - \hat{y}_i)^2`$ in the original space, so by that criterion it always gives the smaller $`SSE`$. If the measurement error is constant regardless of the size of $`y`$, the nonlinear fit matches the assumptions of Chapter 14, Section 4.6 better; if the error is proportional to $`y`$, the linearized fit in log space does. In practice, the result of the linearized fit is often used as the initial estimate for nonlinear regression.

### 6.6 curve_fit and least_squares

SciPy provides functions dedicated to nonlinear least squares. `scipy.optimize.curve_fit` takes a model function `f(x, *params)` and the data and returns the optimal parameters `popt` and their estimated covariance matrix `pcov`. Without bounds the default method is Levenberg-Marquardt (`'lm'`), and the one-standard-deviation errors of the parameters are obtained as `np.sqrt(np.diag(pcov))` [2]. `scipy.optimize.least_squares` takes a residual function and minimizes $`\frac{1}{2}\sum f_i^2`$; its default method is the Trust Region Reflective method (`'trf'`), and the returned `cost` is the value of this objective function at the solution [3].

```python
import numpy as np
from scipy.optimize import curve_fit, least_squares

x = np.array([0, 1, 2, 3, 4], dtype=float)
y = np.array([1.9, 3.1, 4.8, 8.3, 12.6])

def f(x, a0, a1):
    return a0 * np.exp(a1 * x)

popt, pcov = curve_fit(f, x, y, p0=[2.0, 0.4])
print(np.round(popt, 6), np.round(np.sqrt(np.diag(pcov)), 6))

res = least_squares(lambda a: f(x, *a) - y, x0=[2.0, 0.4])
print(np.round(res.x, 6), round(2 * res.cost, 6))
```

```text
[1.959741 0.467527] [0.109695 0.015901]
[1.959741 0.467527] 0.165333
```

Both functions give the same parameters as the Gauss-Newton result of Section 6.4. The second array from `curve_fit` means that the standard errors of $`a_0`$ and $`a_1`$ are about 0.11 and 0.016, which corresponds to computing the covariance matrix $`s_e^2(\mathbf{X}^T\mathbf{X})^{-1}`$ of Section 5.6 with $`\mathbf{X}`$ replaced by the $`\mathbf{Z}`$ at the solution. The `cost` of `least_squares` is $`\frac{1}{2}SSE`$, so doubling it gives $`SSE = 0.165333`$. This covariance results from approximating the model as linear near the solution, so if the nonlinearity is strong, it may not represent the uncertainty accurately [2].

The standard errors from `curve_fit` also give confidence intervals for the parameters. With $`n - 2 = 3`$ degrees of freedom, multiplying by $`t_{0.025,3} = 3.182`$ gives $`a_0 = 1.9597 \pm 0.3491`$ and $`a_1 = 0.4675 \pm 0.0506`$. When evaluating the results of nonlinear regression, the textbook also examines the **coefficient of variation**, the standard error divided by the estimate. In this example it is $`0.1097/1.9597 = 5.6\%`$ and $`0.0159/0.4675 = 3.4\%`$, so both parameters are determined fairly accurately. These intervals, however, rely on a linear approximation near the solution, so with strong nonlinearity or few data they may differ from the actual uncertainty.

### 6.7 Minimizing the Sum of Squares with minimize

As in Section 6.2, nonlinear regression can also be done by passing $`SSE`$ as the objective function to a general-purpose optimizer. Without constraints or bounds the default method of `scipy.optimize.minimize` is the quasi-Newton method BFGS, and in the returned object `x` is the solution and `success` tells whether it succeeded [4]. Applying it to the exponential model of Section 6.4 gives the following.

```python
import numpy as np
from scipy.optimize import minimize

x = np.array([0, 1, 2, 3, 4], dtype=float)
y = np.array([1.9, 3.1, 4.8, 8.3, 12.6])

def sse(a):
    return np.sum((y - a[0] * np.exp(a[1] * x))**2)

res = minimize(sse, x0=[2.0, 0.4])               # no bounds: BFGS
print(res.success, np.round(res.x, 6), round(res.fun, 6))
res_nm = minimize(sse, x0=[2.0, 0.4], method="Nelder-Mead")
print(res_nm.success, np.round(res_nm.x, 6), round(res_nm.fun, 6))
```

```text
True [1.959742 0.467527] 0.165333
True [1.959789 0.46752 ] 0.165333
```

Both methods reach the same minimum $`SSE = 0.165333`$ as the Gauss-Newton result $`(1.959741,\ 0.467527)`$, but Nelder-Mead, which uses no derivatives, gives parameters that differ slightly in the fourth decimal place. A general-purpose optimizer needs only the objective function and works for any model, but it does not exploit the structure of $`SSE`$ as a sum of squared residuals, so it needs more iterations than the Gauss-Newton or Levenberg-Marquardt methods and does not directly give the standard errors of the parameters. When the parameters differ greatly in size, as in the Antoine equation of Section 6.2, the convergence test can trigger too early, so it is safer to scale the parameters to similar sizes or to use a dedicated function.

**Comparing goodness of fit.** Nonlinear models can also be compared by computing $`R^2 = 1 - SSE/SST`$ and $`s_e = \sqrt{SSE/(n - 2)}`$ from the $`SSE`$ in the original space. For the two exponential models of Section 6.5, $`SST = 75.412`$, so the linearized fit has $`R^2 = 0.99751`$ and $`s_e = 0.2503`$, and the nonlinear fit $`R^2 = 0.99781`$ and $`s_e = 0.2348`$. By the criterion of the original space the nonlinear fit is better, but the difference is small, and which one is more appropriate is decided by the nature of the measurement errors (Section 6.5).

### 6.8 Initial Estimates and Convergence Failure

The weakness of the Gauss-Newton method mentioned in Section 6.3 shows when the initial estimate is poor. Starting the exponential model of Section 6.4 from $`(a_0, a_1) = (20,\ -1)`$, far from the solution, gives the following.

```python
import numpy as np
from scipy.optimize import curve_fit

x = np.array([0, 1, 2, 3, 4], dtype=float)
y = np.array([1.9, 3.1, 4.8, 8.3, 12.6])
a = np.array([20.0, -1.0])                        # poor initial guess
with np.errstate(over="ignore", invalid="ignore"):
    for k in range(1, 10):
        Z = np.column_stack([np.exp(a[1] * x), a[0] * x * np.exp(a[1] * x)])
        d = y - a[0] * np.exp(a[1] * x)
        a = a + np.linalg.solve(Z.T @ Z, Z.T @ d)
        print(k, np.round(a, 4), f"{np.sum((y - a[0] * np.exp(a[1] * x))**2):.4g}")

popt, _ = curve_fit(lambda x, a0, a1: a0 * np.exp(a1 * x), x, y, p0=[20.0, -1.0])
print(np.round(popt, 6))
```

```text
1 [ 0.7987 -0.0566] 224.7
2 [0.6124 4.082 ] 5.71e+13
3 [2.000e-04 4.082e+00] 3.817e+06
4 [2.0000e-04 3.8335e+00] 5.118e+05
5 [3.0000e-04 3.3113e+00] 3.145e+04
6 [1.6000e-03 2.1329e+00] 108.8
7 [ 0.0531 -5.7555] 263.7
8 [1.90000000e+00 1.85750841e+04] inf
9 [nan nan] nan
[1.959741 0.467527]
```

In the second iteration $`a_1`$ jumps to 4.082 and $`SSE`$ becomes $`5.7 \times 10^{13}`$; afterward $`a_0`$ drops near zero, the parameters swing widely again, and in the eighth iteration the exponential overflows (`inf`) and the computation collapses. The first-order Taylor approximation holds only near the current point, so starting far away makes the correction $`\Delta\mathbf{a}`$ from the normal equations far too large. `errstate` is used only to hide the overflow warnings of this process. From the same initial values, `curve_fit`, which uses the Levenberg-Marquardt method, reaches the solution $`(1.959741,\ 0.467527)`$ because the $`\mu I`$ term shrinks the corrections. In practice, initial values are chosen from a linearized fit (Section 6.5), and methods with safeguards are used, such as halving the correction whenever $`SSE`$ does not decrease.

---

<br>

## 7. Common Mistakes

The following mistakes are common in general linear and nonlinear regression.

1. **Confusing linearity in the parameters with linearity in $`x`$.** $`y = \beta_0 + \beta_1x + \beta_2x^2`$ is nonlinear in $`x`$ but linear in its parameters, so it is solved in one pass by general linear least squares. Nonlinear regression is needed for models such as $`a_0e^{a_1x}`$, in which the parameters enter nonlinearly.
2. **Solving high-degree polynomials with the normal equations.** As in Section 5.7, the condition number is squared and much accuracy is lost. Use QR factorization or `np.linalg.lstsq` and center the variable.
3. **Reading the $`R^2`$ of a linearized fit as the goodness of fit of the original model.** The $`R^2`$ in log space is the fraction of the variation of the logarithms that is explained and differs from the $`SSE`$ or $`R^2`$ in the original space.
4. **Starting nonlinear regression without an initial estimate.** If `p0` is not given, `curve_fit` sets all initial values to 1 [2]. If the parameters differ greatly from 1 in size, the iteration may fail to converge or go to another minimum, so choose the initial values from a linearized fit or from physical meaning.
5. **Using the wrong degrees of freedom in the standard error.** With $`m + 1`$ parameters, $`s_e = \sqrt{SSE/(n - (m + 1))}`$. If the number of parameters approaches the number of data points, the model may be overfitted even when $`SSE`$ is small.

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
| Conditioning of the normal equations | $`\kappa(\mathbf{X}^T\mathbf{X}) = \kappa(\mathbf{X})^2`$, so for high-degree polynomial regression QR factorization solving $`\mathbf{R}\boldsymbol{\beta} = \mathbf{Q}^T\mathbf{y}`$ or `lstsq` is more accurate. |
| Gauss-Newton method | Linearize the model by a first-order Taylor series and repeatedly solve $`\mathbf{Z}_k^T\mathbf{Z}_k\Delta\mathbf{a} = \mathbf{Z}_k^T\mathbf{d}`$. |
| Linearized versus nonlinear fitting | The linearized fit minimizes the error in the transformed space and the nonlinear fit the error in the original space, so their results differ. |
| `curve_fit`, `least_squares` | `curve_fit` returns the parameters and their covariance, and `least_squares` minimizes $`\frac{1}{2}\sum f_i^2`$ for a residual function. |
| Interpreting multiple regression coefficients | The normal equations fit a plane, and each coefficient is a partial effect with the other variables held fixed. |
| Model selection | $`R^2`$ always grows with more terms, so choose the degree with $`R^2_{adj}`$ and the $`R^2_{pred}`$ computed from PRESS, then check adequacy with the residuals. |
| Regression with `minimize` | Minimizing $`SSE`$ as an objective works for any model but is slower than dedicated methods and gives no standard errors. |

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

5. **One Gauss-Newton step:** For the model $`f(x) = a_0e^{a_1x}`$ and the data $`(0,\ 2)`$, $`(1,\ 4)`$, find $`\mathbf{Z}`$ and the residuals $`\mathbf{d}`$ at the initial estimate $`(a_0, a_1) = (2,\ 0.5)`$.

   <details>
   <summary>Answer</summary>

   Since $`\partial f/\partial a_0 = e^{a_1x}`$ and $`\partial f/\partial a_1 = a_0xe^{a_1x}`$, the rows are $`(1,\ 0)`$ at $`x = 0`$ and $`(e^{0.5},\ 2e^{0.5}) = (1.6487,\ 3.2974)`$ at $`x = 1`$. The predictions are $`2`$ and $`2e^{0.5} = 3.2974`$, so $`\mathbf{d} = (0,\ 0.7026)`$. With two data points and two parameters, $`\mathbf{Z}`$ is a square nonsingular matrix, so the normal equations reduce to $`\mathbf{Z}\Delta\mathbf{a} = \mathbf{d}`$.

   </details>

6. **Squared condition number:** If the 2-norm condition number of the design matrix $`\mathbf{X}`$ is $`10^{6}`$, roughly how many digits can be lost in double precision when solving with the normal equations and with QR factorization?

   <details>
   <summary>Answer</summary>

   The coefficient matrix of the normal equations has condition number $`10^{12}`$, so about 12 digits can be lost, leaving only about 4. With QR factorization, $`\mathbf{R}`$ has condition number $`10^{6}`$, so about 6 digits are lost and about 10 remain.

   </details>

7. **Adjusted coefficient of determination:** A model with 4 coefficients fitted to 12 data points gave $`R^2 = 0.95`$. Find $`R^2_{adj}`$.

   <details>
   <summary>Answer</summary>

   $`R^2_{adj} = 1 - (1 - R^2)\frac{n - 1}{n - (m + 1)} = 1 - 0.05 \times \frac{11}{8} = 1 - 0.06875 = 0.93125`$.

   </details>

8. **Interpreting multiple regression coefficients:** For $`y = 0.9344 + 1.9094x_1 + 0.6094x_2`$, find the prediction at $`x_1 = 3`$, $`x_2 = 4`$, and state how much the prediction changes when only $`x_2`$ increases by 1.

   <details>
   <summary>Answer</summary>

   $`y = 0.9344 + 1.9094 \times 3 + 0.6094 \times 4 = 0.9344 + 5.7282 + 2.4376 = 9.1002`$. With $`x_1`$ held fixed, increasing only $`x_2`$ by 1 increases the prediction by the coefficient $`\beta_2 = 0.6094`$.

   </details>

9. **Meaning of PRESS:** A point has residual $`e_i = 0.2`$ and hat-matrix diagonal entry $`h_{ii} = 0.6`$. Find the residual of predicting this point with the model fitted without it, and explain why points with $`h_{ii}`$ close to 1 are dangerous.

   <details>
   <summary>Answer</summary>

   The left-out prediction residual is $`e_i/(1 - h_{ii}) = 0.2/0.4 = 0.5`$, 2.5 times the ordinary residual. If $`h_{ii}`$ is close to 1, the point pulls the fit strongly toward itself, so its ordinary residual looks small, but the prediction goes badly wrong when the point is left out. Such a point is a high-leverage point with a large influence on the model and should be checked separately.

   </details>

---

<br>

## References

1. NumPy Developers, "numpy.linalg.lstsq," *NumPy Documentation*. [Online]. Available: https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
2. SciPy Developers, "scipy.optimize.curve_fit," *SciPy Documentation*. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
3. SciPy Developers, "scipy.optimize.least_squares," *SciPy Documentation*. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
4. SciPy Developers, "scipy.optimize.minimize," *SciPy Documentation*. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
