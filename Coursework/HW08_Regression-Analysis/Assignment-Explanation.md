# HW8 Assignment Explanation — Regression Analysis

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 14, 15

> **Learning Objectives:**
> 1. Complete a function that computes the slope and intercept of simple linear regression under the least-squares criterion
> 2. Evaluate the goodness of fit of a regression model with SST, SSE, SSR, and the coefficient of determination
> 3. Linearize the power model by a log transformation to obtain its parameters and compare the fit in the original space
> 4. Perform polynomial regression with the design matrix and the normal equations and interpret how the fit statistics change with the degree
> 5. Estimate the parameters of a model that cannot be linearized with `curve_fit` and interpret the result

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Coefficients of Simple Linear Regression](#3-problem-1-coefficients-of-simple-linear-regression)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Goodness-of-Fit Statistics](#4-problem-2-goodness-of-fit-statistics)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: Linearizing the Power Model](#5-problem-3-linearizing-the-power-model)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Problem 4: Normal Equations for Polynomial Regression](#6-problem-4-normal-equations-for-polynomial-regression)
  - [6.1 Approach](#61-approach)
  - [6.2 Code Written](#62-code-written)
  - [6.3 Results and Verification](#63-results-and-verification)
- [7. Problem 5: Nonlinear Regression of the Antoine Equation](#7-problem-5-nonlinear-regression-of-the-antoine-equation)
  - [7.1 Approach](#71-approach)
  - [7.2 Code Written](#72-code-written)
  - [7.3 Results and Verification](#73-results-and-verification)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the five blanks of the lab notebook `Lab09-Regression.ipynb` distributed by the instructor and run it. The notebook has six parts: why regression is needed, simple linear regression, linearizable nonlinear regression, general linear least squares, nonlinear regression with `curve_fit`, and model comparison and selection, and the blanks lie in Sections 2-5. The functions completed in earlier blanks are used again in later cells. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab09-Regression.ipynb](Original/Lab09-Regression.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab09-Regression_Completed.ipynb](Solution/Lab09-Regression_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 2.1 Least-squares formulas | `a1` and `a0` in the function `linreg` |
| Problem 2 | 2.2 Goodness-of-fit statistics | `SST`, `SSE`, `SSR`, and `R2` in the function `goodness_of_fit` |
| Problem 3 | 3.1 Power model | `X`, `Y`, `alpha`, and `beta` in the function `power_fit` |
| Problem 4 | 4.1 Normal equations | `X[:, j]` and `c` in the function `poly_lstsq` |
| Problem 5 | 5.1 Models that cannot be linearized | `popt_ant` and `pcov_ant` in the Antoine-equation cell |

*Table 1. The five blanks that make up the assignment*

The least-squares line and the quantification of the fit are covered in Sections 4 and 5 of the [Chapter 14 lecture note](../../Lecture-Notes/Ch14_Linear-Regression/Concepts_Lecture.md), and linearization, polynomial regression, general linear least-squares regression, and nonlinear regression are covered in Sections 2, 3, 5, and 6 of the [Chapter 15 lecture note](../../Lecture-Notes/Ch15_Nonlinear-Regression/Concepts_Lecture.md). This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell stays as distributed, and only the blanks are filled in.

1. **Problem 1:** `linreg` must compute the slope `a1` and the intercept `a0` with the least-squares formulas of Section 2.1 of the notebook and return them in the order `(a0, a1)`. The function is also used for the power-model and exponential-model fits in Section 3.
2. **Problem 2:** `goodness_of_fit` must compute `SST`, `SSE`, `SSR`, and `R2` as defined in Table 1, and the code after the blanks uses these values to compute `se` and `R2adj`.
3. **Problem 3:** `power_fit` must transform the data to common-logarithm space and recover `alpha` and `beta` from the log-space intercept and slope obtained with `linreg`.
4. **Problem 4:** `poly_lstsq` must build the design matrix whose $`j`$th column is $`x^j`$, solve the normal equations $`\mathbf{X}^T\mathbf{X}\mathbf{c} = \mathbf{X}^T\mathbf{y}`$, and return the coefficients `c`. The function is also used in the model comparison of Section 6.
5. **Problem 5:** The optimal parameters of the Antoine equation and their covariance matrix must be obtained with `curve_fit`, starting from the given initial guess `p0_ant`.

---

<br>

## 3. Problem 1: Coefficients of Simple Linear Regression

### 3.1 Approach

The sums $`\sum`$ in the formulas of Section 2.1 of the notebook become `np.sum`, and the means $`\bar{x}`$ and $`\bar{y}`$ become `np.mean`. The intercept formula uses the slope, so `a1` is computed first.

$$
a_1 = \frac{n\sum x_i y_i - \sum x_i \sum y_i}{n\sum x_i^2 - \left(\sum x_i\right)^2}, \qquad a_0 = \bar{y} - a_1 \bar{x}
$$

The second form given in the notebook is the same expression with the numerator and denominator divided by $`n`$. The intercept formula also means that the least-squares line always passes through the point $`(\bar{x}, \bar{y})`$. The derivation of the normal equations and the solution for the coefficients are covered in Sections 4.3 and 4.4 of the [Chapter 14 lecture note](../../Lecture-Notes/Ch14_Linear-Regression/Concepts_Lecture.md).

### 3.2 Code Written

```python
    a1 = (n*np.sum(xi*yi) - np.sum(xi)*np.sum(yi)) / (n*np.sum(xi**2) - np.sum(xi)**2)
    a0 = np.mean(yi) - a1 * np.mean(xi)
```

### 3.3 Results and Verification

```text
Slope     a1 = 0.210357 mm/kN
Intercept a0 = 0.000000 mm
Stiffness  k = 1/a1 = 4.7538 kN/mm
```

The output matches the output in the distributed notebook. For the beam-deflection data, $`n = 8`$, $`\sum x_i = 280`$, $`\sum y_i = 58.9`$, $`\sum x_i y_i = 2945`$, and $`\sum x_i^2 = 14000`$, so the hand calculation is as follows.

$$
a_1 = \frac{8 \times 2945 - 280 \times 58.9}{8 \times 14000 - 280^2} = \frac{7068}{33600} = 0.210357
$$

With $`\bar{x} = 35`$ and $`\bar{y} = 7.3625`$, the equality $`7068 \times 35 / 33600 = 7.3625`$ holds exactly, so $`a_0 = \bar{y} - a_1 \bar{x} = 0`$. For this data the least-squares line through $`(\bar{x}, \bar{y})`$ also passes through the origin, which agrees with the interpretation that the deflection is proportional to the load. The stiffness is $`k = 1/a_1 = 4.7538`$ kN/mm. The comparison cell that follows also fits the same data with `scipy.stats.linregress` and `np.polyfit`.

```text
=== Method Comparison ===
Method                        a1 (slope)    a0 (intercept)          R2
------------------------------------------------------------------------
Our linreg                    0.21035714        0.00000000  0.99974064
scipy.stats.linregress        0.21035714        0.00000000  0.99974064
np.polyfit(deg=1)             0.21035714       -0.00000000
```

The slopes of the three methods agree to eight decimal places. The `-0.00000000` from `np.polyfit` is a rounding error in which the intercept comes out as a negative number close to 0.

---

<br>

## 4. Problem 2: Goodness-of-Fit Statistics

### 4.1 Approach

The definitions of Table 1 are translated directly. `SST` is the sum of squared differences between the observations and their mean, `SSE` is the sum of squared differences between the observations and the predictions (the residuals), `SSR` is their difference, and `R2` is $`1 - S_r/S_T`$. The code after the blanks uses these values to compute the standard error of the estimate $`s_e = \sqrt{S_r/(n - p)}`$ and the adjusted coefficient of determination (adjusted $`R^2`$), where $`p`$ is the number of parameters `n_params`. The decomposition into SST, SSE, and SSR and the coefficient of determination are covered in Section 5 of the [Chapter 14 lecture note](../../Lecture-Notes/Ch14_Linear-Regression/Concepts_Lecture.md).

### 4.2 Code Written

```python
    SST = np.sum((yi - y_bar)**2)   # total sum of squares
    SSE = np.sum((yi - yi_hat)**2)  # residual sum of squares
    SSR = SST - SSE                 # regression sum of squares
    R2  = 1 - SSE / SST             # coefficient of determination
```

### 4.3 Results and Verification

```text
=== Goodness of Fit: Beam Deflection (Linear) ===
SST    = 185.898750   (total variation)
SSE    = 0.048214   (residual/unexplained)
SSR    = 185.850536   (explained by model)
R2     = 0.99974064
Adj R2 = 0.99969742
se     = 0.089642 mm
SST = SSR + SSE? True
```

The output matches the output in the distributed notebook. `SST` equals `SST (corrected) = 185.898750` printed by the sample-statistics cell of Section 1.3 of the notebook, which agrees with $`(n-1)s_y^2 = 7 \times 5.1533^2`$ computed from the sample standard deviation 5.1533. The remaining values are confirmed as follows.

- $`R^2 = 1 - 0.048214 / 185.898750 = 0.99974064`$
- $`s_e = \sqrt{0.048214 / (8 - 2)} = 0.089642`$
- Adjusted coefficient of determination $`= 1 - (0.048214 / 6) / (185.898750 / 7) = 0.99969742`$

> **Note:** The last line `SST = SSR + SSE? True` is always true because the code defines `SSR` as `SST - SSE`. The substantive identity is that the explained sum of squares $`\sum(\hat{y}_i - \bar{y})^2`$ equals $`S_T - S_r`$, which holds for a least-squares fit with an intercept. Computing $`\sum(\hat{y}_i - \bar{y})^2`$ separately for the beam-deflection line also gives 185.850536, the same as `SSR`.

---

<br>

## 5. Problem 3: Linearizing the Power Model

### 5.1 Approach

Transforming the data with `np.log10` into $`X = \log_{10} x`$ and $`Y = \log_{10} y`$ lets the call `linreg(X, Y)` between the blanks fit the log-space line $`Y = A_0 + A_1 X`$. By the relation in Section 3.1 of the notebook, the slope is the exponent and the intercept is the common logarithm of the coefficient, so the parameters are recovered as $`\beta = A_1`$ and $`\alpha = 10^{A_0}`$. Linearization of the power model is covered in Section 2.2 of the [Chapter 15 lecture note](../../Lecture-Notes/Ch15_Nonlinear-Regression/Concepts_Lecture.md).

### 5.2 Code Written

```python
    X = np.log10(xi)
    Y = np.log10(yi)
```

```python
    alpha = 10**A0
    beta  = A1
```

### 5.3 Results and Verification

```text
Power law fit:  y = 0.2741 * x^1.9842
Exponent beta ~ 1.98  (response ~ stress squared)

Linear:    R2=0.88049,  se=189.79
Power law: R2=0.80882,  se=240.04
```

The output matches the output in the distributed notebook. Repeating the computation gives the log-space intercept and slope $`A_0 = -0.56203`$ and $`A_1 = 1.98418`$, and $`10^{-0.56203} = 0.27414`$ matches $`\alpha`$ in the output.

A point worth noting is that $`R^2`$ computed in the original $`y`$ space is larger for the straight line (0.88049) than for the power model (0.80882). The coefficient of determination of the same line computed in log space is a high 0.948. Log-linearization minimizes the residual sum of squares of $`\log y`$, so it reduces relative errors evenly, and as a result the absolute errors at large responses such as 1220, 830, and 1450 remain large in the SSE of the original space. The comparison cell of Section 5 of the notebook shows this difference by fitting the same model in the original space with `curve_fit`.

```text
=== Power Law: Log-linearization vs. curve_fit ===
                               alpha          beta          R2
--------------------------------------------------------------
Log-linearization           0.274137      1.984176    0.808818
curve_fit                   2.538367      1.435856    0.876898
curve_fit std err           3.620273      0.340430

Note: curve_fit minimizes SSE in the original y-space;
      log-linearization minimizes SSE in log(y) space.
```

The $`R^2`$ of `curve_fit` (0.876898), which minimizes the SSE in the original space directly, is larger than that of log-linearization. It is still slightly smaller than that of the straight line, because the line $`a_0 + a_1 x`$ has an intercept whereas $`\alpha x^\beta`$ must pass through the origin, so the two models are different families of functions.

---

<br>

## 6. Problem 4: Normal Equations for Polynomial Regression

### 6.1 Approach

The basis functions of a polynomial of degree $`m`$ are $`z_j(x) = x^j`$, so the $`j`$th column of the design matrix is `xi**j`. For $`j = 0`$ every entry is 1, and that column carries the constant term. The coefficients are obtained by solving the normal equations with `np.linalg.solve`, without computing the inverse $`(\mathbf{X}^T\mathbf{X})^{-1}`$ separately. The normal equations of polynomial regression and general linear least-squares regression are covered in Sections 3 and 5 of the [Chapter 15 lecture note](../../Lecture-Notes/Ch15_Nonlinear-Regression/Concepts_Lecture.md).

### 6.2 Code Written

```python
        X[:, j] = xi**j  # (j-th column = x to the power j)
```

```python
    c = np.linalg.solve(X.T @ X, X.T @ yi)
```

### 6.3 Results and Verification

```text
Quadratic: y = 2.4786 + 2.3593t + 1.8607t²
```

The output matches the output in the distributed notebook. Applying `np.polyfit(x_vel, y_vel, 2)` to the same data returns $`[1.86071, 2.35929, 2.47857]`$ in descending order, so the coefficients agree. The legend of the figure drawn by the cell shows $`R^2`$ and $`s_e`$ for each degree, with the following values.

| Degree | $`R^2`$ | $`s_e`$ |
|:-----|:-----|:-----|
| 1 | 0.94708 | 5.766 |
| 2 | 0.99851 | 1.118 |
| 3 | 0.99866 | 1.299 |

*Table 2. Goodness of fit by degree for polynomials fitted to the velocity data*

$`R^2`$ grows with every increase in degree, but at degree 3 the value $`s_e = \sqrt{S_r/(n - m - 1)}`$ actually grows. This means the extra parameter did not reduce the SSE enough, so the quadratic model is adequate, as the comment of the cell anticipates.

`poly_lstsq` is also used in the model-comparison cell of Section 6 of the notebook. This cell fits polynomials of degrees 1-8 to 13 viscosity points using standardized temperature.

```text
 Degree          R2        Adj R2          se
---------------------------------------------
      1   0.8514521     0.8379478    0.190883
      2   0.9837378     0.9804854    0.066240
      3   0.9985234     0.9980312    0.021040
      4   0.9998335     0.9997503    0.007493
      5   0.9999585     0.9999289    0.003999
      6   0.9999905     0.9999810    0.002065
      7   0.9999993     0.9999984    0.000601
      8   0.9999997     0.9999990    0.000474

Observation: Adj. R² peaks around degree 4-5, then decreases -> overfitting!
```

> **Note:** The sentence printed at the end of the cell says that the adjusted coefficient of determination peaks around degree 4-5 and then decreases, but the adjusted coefficient of determination and $`s_e`$ in the table keep improving through degree 8. For this data the two metrics alone cannot detect overfitting, and the sentence prints an expectation written in advance rather than a computed result.

---

<br>

## 7. Problem 5: Nonlinear Regression of the Antoine Equation

### 7.1 Approach

Following the form given in Section 5.1 of the notebook, `curve_fit(model function, x data, y data, p0=initial guess)` is called. The model function `antoine_model(T, A, B, C)` takes the independent variable as its first argument and the parameters to estimate as the rest, which matches the form that `curve_fit` requires. Nonlinear least squares is iterative and needs an initial guess, so the given `p0_ant = [8.0, 1500.0, 200.0]` is passed. The returned `popt_ant` holds the optimal parameters, and `pcov_ant` is their covariance matrix. Nonlinear regression and the Antoine equation are covered in Section 6 of the [Chapter 15 lecture note](../../Lecture-Notes/Ch15_Nonlinear-Regression/Concepts_Lecture.md).

### 7.2 Code Written

```python
popt_ant, pcov_ant = curve_fit(antoine_model, T_ant, P_ant, p0=p0_ant)
```

### 7.3 Results and Verification

```text
=== Antoine Equation Fit (water vapor pressure) ===
A = 7.997574 +/- 0.007401
B = 1687.0509 +/- 4.6041
C = 229.7113 +/- 0.4237
R2 = 0.99999999
Literature: A=8.07131, B=1730.63, C=233.426
```

The output matches the output in the distributed notebook. The vapor pressures computed with the fitted parameters are 55.28, 92.51, 149.40, 233.70, 355.14, 525.75, and 760.02 mmHg, each within 0.05 mmHg of the data.

The fitted parameters differ from the literature values, yet both explain the data well. The vapor pressures computed with the literature values are 55.19, 92.30, 149.04, 233.17, 354.53, 525.27, and 760.09 mmHg, with $`R^2 = 0.9999973`$. The SSE is 1.07 for the literature values and 0.0050 for the fit, so `curve_fit` found parameters closer to these seven points. The correlation coefficients between the parameters computed from `pcov_ant` are 0.99964 for $`A`$ and $`B`$, 0.99847 for $`A`$ and $`C`$, and 0.99959 for $`B`$ and $`C`$, all close to 1. A change in one parameter is almost canceled by the others, so different parameter sets give nearly the same curve over this temperature range.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Simple linear regression | The slope comes from the least-squares formula built from sums and the intercept from $`\bar{y} - a_1 \bar{x}`$, and the line for the beam-deflection data has slope 0.210357 and passes through the origin. |
| Goodness-of-fit statistics | $`R^2 = 1 - S_r/S_T`$ and $`s_e`$ evaluate the fit, and when the code defines SSR as $`S_T - S_r`$ the decomposition check is always true. |
| Log-linearization | The power model is fitted as a line in common-logarithm space and recovered as $`\beta = A_1`$ and $`\alpha = 10^{A_0}`$, and because it minimizes the log-space SSE its $`R^2`$ in the original space can be lower than that of `curve_fit`. |
| Normal equations | Building the design matrix with $`x^j`$ as its $`j`$th column and solving $`\mathbf{X}^T\mathbf{X}\mathbf{c} = \mathbf{X}^T\mathbf{y}`$ gives the polynomial-regression coefficients, which agree with `np.polyfit`. |
| Degree selection | $`R^2`$ grows with the degree, so the degree is chosen with $`s_e`$ and the adjusted coefficient of determination, which account for the number of parameters. |
| Nonlinear regression | `curve_fit` starts from an initial guess and minimizes the SSE in the original space directly, and when the parameters are strongly correlated, different parameter sets can give similar curves. |

---

<br>

## Self-Check Questions

1. **A line through the origin:** Explain with a property of the least-squares line why the intercept in Problem 1 came out exactly 0.

   <details>
   <summary>Answer</summary>

   The least-squares line has $`a_0 = \bar{y} - a_1 \bar{x}`$, so it always passes through $`(\bar{x}, \bar{y}) = (35, 7.3625)`$. For this data $`a_1 \bar{x} = 7068 \times 35 / 33600 = 7.3625 = \bar{y}`$ holds exactly, so $`a_0 = 0`$ and the line also passes through the origin.

   </details>

2. **A check that is always true:** Why is `SST = SSR + SSE? True` in the output of Problem 2 true regardless of how good the fit is?

   <details>
   <summary>Answer</summary>

   The code defines `SSR = SST - SSE`, so `SSR + SSE` equals `SST` by definition. A meaningful check computes $`\sum(\hat{y}_i - \bar{y})^2`$ separately and compares it with $`S_T - S_r`$, and this equality holds for a least-squares fit with an intercept.

   </details>

3. **Low R² of the power model:** Explain why the $`R^2`$ of the power model in Problem 3 came out lower than that of the straight line.

   <details>
   <summary>Answer</summary>

   Log-linearization minimizes the residual sum of squares of $`\log_{10} y`$, so it does not minimize the SSE in the original space. The absolute errors at large responses remain large in the original-space SSE, which lowers $`R^2`$ to 0.80882. Fitting the same model in the original space with `curve_fit` raises $`R^2`$ to 0.876898.

   </details>

4. **s_e of the cubic:** For the velocity data of Problem 4, why did $`s_e`$ grow for the cubic polynomial even though its $`R^2`$ is larger than that of the quadratic?

   <details>
   <summary>Answer</summary>

   In $`s_e = \sqrt{S_r/(n - m - 1)}`$ the denominator shrinks as the number of parameters grows. Adding the cubic term reduced $`S_r`$ by only about 10%, from 3.75 to 3.37, while the denominator shrank from 3 to 2, so $`s_e`$ grew from 1.118 to 1.299. The added term did not improve the fit meaningfully.

   </details>

5. **Parameters different from the literature:** Why is $`R^2`$ close to 1 in Problem 5 even though the fitted parameters differ from the literature values?

   <details>
   <summary>Answer</summary>

   The correlation coefficients among $`A`$, $`B`$, and $`C`$ are 0.998 or higher, close to 1, so a change in one parameter is almost canceled by the others. Different parameter sets, such as the literature values and the fitted values, therefore give nearly the same curve over 40-100 °C. The literature values also give $`R^2 = 0.9999973`$, and `curve_fit` found parameters that reduce the SSE of these seven points further, from 1.07 to 0.0050.

   </details>
