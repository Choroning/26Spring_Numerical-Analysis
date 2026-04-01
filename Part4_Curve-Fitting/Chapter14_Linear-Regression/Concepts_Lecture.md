# Chapter 14 Lecture -- Linear Regression

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 14

> **Prerequisites**: [Statistics] Basic statistics. [Calculus] Least squares concepts (Ch 1-13).
>
> **Learning Objectives**:
> 1. Derive and apply simple linear regression
> 2. Assess regression quality using R-squared and standard error
> 3. Apply linearization for nonlinear models

---

<br>

## Table of Contents

- [1. Overview](#1-overview)
- [2. Statistics Review](#2-statistics-review)
  - [2.1 Descriptive Statistics](#21-descriptive-statistics)
  - [2.2 Standard Deviation and Variance](#22-standard-deviation-and-variance)
  - [2.3 Coefficient of Variation](#23-coefficient-of-variation)
  - [2.4 Median Absolute Deviation (MAD)](#24-median-absolute-deviation-mad)
- [3. Normal Distribution](#3-normal-distribution)
  - [3.1 Central Limit Theorem](#31-central-limit-theorem)
  - [3.2 Standard Normal Distribution](#32-standard-normal-distribution)
  - [3.3 Symmetry Properties](#33-symmetry-properties)
  - [3.4 General Normal Distribution](#34-general-normal-distribution)
  - [3.5 Standardization](#35-standardization)
  - [3.6 Normal CDF and PDF](#36-normal-cdf-and-pdf)
- [4. Straight Line Least Square Regression](#4-straight-line-least-square-regression)
  - [4.1 Problem Setup](#41-problem-setup)
  - [4.2 Alternative Error Criteria and Why Least Squares](#42-alternative-error-criteria-and-why-least-squares)
  - [4.3 Deriving the Normal Equations](#43-deriving-the-normal-equations)
  - [4.4 Solving for Coefficients](#44-solving-for-coefficients)
- [5. Quantifying the Fit Quality](#5-quantifying-the-fit-quality)
  - [5.1 SST, SSE, and SSR Decomposition](#51-sst-sse-and-ssr-decomposition)
  - [5.2 Proof that SST = SSE + SSR](#52-proof-that-sst--sse--ssr)
  - [5.3 Coefficient of Determination (r^2)](#53-coefficient-of-determination-r2)
  - [5.4 Standard Error of the Estimate](#54-standard-error-of-the-estimate)
- [6. Summary Table](#6-summary-table)

<br>

---

<br>

## 1. Overview

Linear regression is a fundamental technique in **curve fitting** (Part IV of numerical methods). The goal is to derive a **best-fit straight line** through a set of data points using the **least-squares** criterion.

**Learning objectives:**

- Review basic statistics (mean, standard deviation, distributions)
- Compute the slope and intercept of a best-fit straight line
- Compute the standard error of the estimate and analyze residual error

<br>

---

<br>

## 2. Statistics Review

### 2.1 Descriptive Statistics

**Statistics** is the discipline that deals with data from which meaningful information is inferred.

**Sampling** is the process of choosing a smaller group of items from a larger population to gain insights about the entire group.

#### Sample Average (Mean)

If you have a sample of $n$ observations $y_1, y_2, \ldots, y_n$, the **sample average** is:

$$\bar{y} = \frac{1}{n} \sum_{j=1}^{n} y_j$$

#### Median

The **median** $\tilde{y}$ is the middle value in a dataset when sorted in order.

- Odd count: Data = [3, 5, 7] --> Median = 5
- Even count: Data = [2, 4, 6, 8] --> Median = (4 + 6) / 2 = 5

#### Mode

The **mode** is the value that occurs most frequently in a dataset.

- Data = [2, 4, 4, 4, 5] --> Mode = 4

#### Outlier

An **outlier** is a data point that is significantly different from the rest.

- Data = [2, 4, 4, 4, 100] --> Outlier = 100

<br>

### 2.2 Standard Deviation and Variance

#### Population Standard Deviation

The **standard deviation** measures how much data points spread out from the mean:

$$\sigma = \sqrt{\frac{1}{n} \sum_{j=1}^{n} (y_j - \mu)^2}$$

where $\mu$ is the **true mean** of the entire population, and $n$ is the total number of the population.

> **[Statistics]** The population standard deviation uses $n$ in the denominator because all members of the population are accounted for, so no degree of freedom is lost.

#### Sample Standard Deviation

$$s_y = \sqrt{\frac{1}{n-1} \sum_{j=1}^{n} (y_j - \bar{y})^2}$$

where $\bar{y}$ is the **sample mean** and $n$ is the total number of data points in a sample.

The quantity inside the square root's numerator is called the **Total Corrected Sum of Squares (SST)**:

$$SST = \sum_{j=1}^{n} (y_j - \bar{y})^2$$

> **[Statistics]** Why $n - 1$ instead of $n$? Using the sample mean $\bar{y}$ instead of the true mean $\mu$ **underestimates the variability**. Since $\bar{y}$ is a random variable that already used 1 degree of freedom, we divide by $n - 1$ (Bessel's correction) to get an **unbiased** estimate of the population variance.

<br>

### 2.3 Coefficient of Variation

The **Coefficient of Variation (c.v.)** is the ratio of the standard deviation to the mean, expressed as a percentage:

$$c.v. = \frac{s_y}{\bar{y}} \cdot 100\%$$

> **[Statistics]** The c.v. normalizes the spread, so you can **compare datasets even if their scales are different**. A low c.v. implies data is tightly clustered around the mean; a high c.v. implies significant variability relative to the mean.

<br>

### 2.4 Median Absolute Deviation (MAD)

The **Median Absolute Deviation (MAD)** measures how far values deviate from the median:

$$MAD = \text{median}(|y_i - \tilde{y}|)$$

**Example:**

- Data: [2, 4, 6, 8, 100] --> median = 6
- Absolute deviations: |2-6|=4, |4-6|=2, |6-6|=0, |8-6|=2, |100-6|=94
- Sorted deviations: [0, 2, 2, 4, 94] --> MAD = 2

> **[Statistics]** MAD is **robust to outliers**. In the example above, the huge outlier 100 produces a deviation of 94, but since MAD takes the median of deviations, it completely ignores this extreme value. This makes MAD a more reliable measure of spread when outliers are present, compared to standard deviation.

<br>

---

<br>

## 3. Normal Distribution

### 3.1 Central Limit Theorem

Take many independent random variables (r.v.s), each with finite mean and variance, and add them up. Then the sum of a large number of i.i.d. random variables **tends to follow a normal distribution**, no matter what the distribution of the individual random variables is.

> **[Statistics]** The Central Limit Theorem (CLT) is the foundational reason why the normal distribution is so prevalent. Measurement errors, which are sums of many small independent effects, are typically normally distributed.

<br>

### 3.2 Standard Normal Distribution

The **standard normal distribution** is very popular due to the Central Limit Theorem.

**PDF (Probability Density Function):**

$$\varphi(z) = \frac{1}{\sqrt{2\pi}} e^{-\frac{z^2}{2}}, \quad \text{for } -\infty < z < \infty$$

**CDF (Cumulative Distribution Function):**

$$\Phi(z) = \int_{-\infty}^{z} \varphi(t) \, dt, \quad \text{for } -\infty < z < \infty$$

We denote this by $Z \sim \mathcal{N}(0, 1)$, meaning **mean 0 and variance 1**.

<br>

### 3.3 Symmetry Properties

**1. Symmetry of the PDF:**

$$\varphi(z) = \varphi(-z)$$

The bell curve is symmetric about $z = 0$.

**2. Symmetry of the tail area:**

$$P(Z \le -2) = P(Z \ge 2)$$

which is equivalent to:

$$\Phi(-z) = 1 - \Phi(z)$$

**3. Symmetry of $Z$ and $-Z$:**

If $Z \sim \mathcal{N}(0, 1)$, then $-Z \sim \mathcal{N}(0, 1)$ as well.

$$P(-Z \le z) = P(Z \ge -z) = 1 - \Phi(-z) = \Phi(z)$$

<br>

### 3.4 General Normal Distribution

If $Z \sim \mathcal{N}(0, 1)$, then:

$$X = \mu + \sigma Z$$

has the normal distribution with **mean** $\mu$ **and variance** $\sigma^2$, denoted $X \sim \mathcal{N}(\mu, \sigma^2)$.

**Proof of mean and variance:**

$$E(\mu + \sigma Z) = E(\mu) + E(\sigma Z) = \mu$$

$$Var(\mu + \sigma Z) = Var(\sigma Z) = \sigma^2 Var(Z) = \sigma^2$$

> **[Statistics]** The linearity of expectation ($E[aX + b] = aE[X] + b$) and the scaling property of variance ($Var(aX + b) = a^2 Var(X)$) are the key properties used here. Adding a constant shifts the mean but does not change the variance.

<br>

### 3.5 Standardization

For $X \sim \mathcal{N}(\mu, \sigma^2)$, the **standardization** of $X$ is:

$$Z = \frac{X - \mu}{\sigma} \sim \mathcal{N}(0, 1)$$

This transforms any normal random variable into a standard normal, enabling the use of standard normal tables.

<br>

### 3.6 Normal CDF and PDF

Let $X \sim \mathcal{N}(\mu, \sigma^2)$.

**CDF of $X$:**

$$F(x) = P(X \le x) = P\left(\frac{X - \mu}{\sigma} \le \frac{x - \mu}{\sigma}\right) = \Phi\left(\frac{x - \mu}{\sigma}\right)$$

**PDF of $X$:**

$$f(x) = F'(x) = \varphi\left(\frac{x - \mu}{\sigma}\right) \cdot \frac{1}{\sigma} = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mu)^2}{2\sigma^2}}$$

> **[Statistics]** The $\frac{1}{\sigma}$ factor in the PDF comes from the chain rule when differentiating the CDF. This ensures the total area under the PDF remains 1 regardless of $\sigma$.

<br>

---

<br>

## 4. Straight Line Least Square Regression

### 4.1 Problem Setup

Given data points $(x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)$, we want to find a **linear model**:

$$\hat{y} = a_0 + a_1 x$$

where:
- $a_0$ = **y-intercept**
- $a_1$ = **slope**

Each observed value can be written as:

$$y_i = \hat{y}_i + e_i = a_0 + a_1 x_i + e_i$$

where $e_i = y_i - \hat{y}_i$ is the **residual (error)** -- the difference between the observed value and the predicted value.

<br>

### 4.2 Alternative Error Criteria and Why Least Squares

We want the "best" line through the data by minimizing some measure of error. Several criteria can be considered:

**Criterion 1 -- Minimize the sum of errors:**

$$\min_{a_0, a_1} \sum e_i = \min_{a_0, a_1} \sum (y_i - a_0 - a_1 x_i)$$

Problem: Positive and negative errors **cancel each other out**, giving a misleading result.

**Criterion 2 -- Minimize the sum of absolute errors:**

$$\min_{a_0, a_1} \sum |e_i| = \min_{a_0, a_1} \sum |y_i - a_0 - a_1 x_i|$$

Problems: The absolute value function is **non-differentiable** at zero, and the solution is **not unique**.

**Criterion 3 -- Minimax (minimize the maximum error):**

$$\min_{a_0, a_1} \max_i |y_i - a_0 - a_1 x_i|$$

Problems: **Non-differentiable**, and heavily **sensitive to outliers** since it emphasizes the single largest error.

**Criterion 4 -- Minimize the sum of squared errors (Least Squares):**

$$\min_{a_0, a_1} \sum e_i^2 = \min_{a_0, a_1} \sum (y_i - a_0 - a_1 x_i)^2$$

This is the **least squares** approach, which yields the **Sum of Squared Errors (SSE)**. It is differentiable and yields a **unique** line for a given set of data.

> **[Statistics]** The least squares criterion is preferred because: (1) squaring removes sign issues, (2) the objective function is smooth and differentiable, (3) it has a unique closed-form solution, and (4) under the assumption of normally distributed errors, it corresponds to the maximum likelihood estimate.

<br>

### 4.3 Deriving the Normal Equations

The objective function to minimize is:

$$SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i)^2$$

To find the minimum, differentiate SSE with respect to each unknown coefficient and set the derivatives to zero:

$$\frac{\partial SSE}{\partial a_0} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i) = 0$$

$$\frac{\partial SSE}{\partial a_1} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i) x_i = 0$$

Since $a_0$ and $a_1$ are constants (not dependent on $i$), we can factor them out of the sums:

- $\sum a_0 = a_0 \cdot n$
- $\sum a_0 x_i = a_0 \sum x_i$
- $\sum a_1 x_i = a_1 \sum x_i$
- $\sum a_1 x_i^2 = a_1 \sum x_i^2$

Setting the partial derivatives to zero yields the **normal equations**:

$$\sum y_i = n \cdot a_0 + a_1 \sum x_i$$

$$\sum y_i x_i = a_0 \sum x_i + a_1 \sum x_i^2$$

**Matrix form:**

$$\begin{bmatrix} \sum y_i \\ \sum y_i x_i \end{bmatrix} = \begin{bmatrix} n & \sum x_i \\ \sum x_i & \sum x_i^2 \end{bmatrix} \begin{bmatrix} a_0 \\ a_1 \end{bmatrix} \quad (\ast)$$

> **[Linear Algebra]** This is a $2 \times 2$ linear system. The coefficient matrix is symmetric and positive definite (assuming the $x_i$ values are not all identical), guaranteeing a unique solution.

<br>

### 4.4 Solving for Coefficients

#### Using Cramer's Rule

The determinant of the coefficient matrix:

$$D = n \sum x_i^2 - \left(\sum x_i\right)^2$$

The intercept:

$$a_0 = \frac{1}{D} \left( \sum y_i \cdot \sum x_i^2 - \sum x_i \cdot \sum y_i x_i \right)$$

The slope:

$$a_1 = \frac{1}{D} \left( n \sum y_i x_i - \sum x_i \cdot \sum y_i \right)$$

#### Alternative Form Using Sample Means

From the first row of $(\ast)$:

$$\sum y_i = n \cdot a_0 + \left(\sum x_i\right) a_1$$

Dividing both sides by $n$:

$$\frac{\sum y_i}{n} = a_0 + \frac{\sum x_i}{n} \cdot a_1$$

Since $\bar{y} = \frac{\sum y_i}{n}$ and $\bar{x} = \frac{\sum x_i}{n}$ are the sample means:

$$\boxed{a_0 = \bar{y} - \bar{x} \cdot a_1} \quad (\ast\ast)$$

> **[Statistics]** This result has an elegant geometric interpretation: the regression line **always passes through the point $(\bar{x}, \bar{y})$**, the centroid of the data.

<br>

---

<br>

## 5. Quantifying the Fit Quality

### 5.1 SST, SSE, and SSR Decomposition

Three key quantities measure different aspects of variability:

| Quantity | Full Name | Formula | Interpretation |
|----------|-----------|---------|----------------|
| **SST** | Sum of Squared Total (total corrected sum of squares) | $\displaystyle \sum_{i=1}^{n}(y_i - \bar{y})^2$ | Total deviation of data points from the mean |
| **SSE** | Sum of Squared Errors (residual sum of squares) | $\displaystyle \sum_{i=1}^{n}(y_i - \hat{y}_i)^2$ | Unexplained variability (error between observed and predicted) |
| **SSR** | Sum of Squared Regression (explained sum of squares) | $\displaystyle \sum_{i=1}^{n}(\hat{y}_i - \bar{y})^2$ | Explained variability (how well $\hat{y}$ fits the data) |

The fundamental decomposition:

$$\boxed{SST = SSE + SSR}$$

$$\text{Total Variability} = \text{Unexplained Variability} + \text{Explained Variability}$$

> **[Statistics]** If $SSR > SSE$, then the regression model explains more variability than it leaves unexplained, meaning the model is worthwhile.

<br>

### 5.2 Proof that SST = SSE + SSR

**Step 1:** Express SST by adding and subtracting $\hat{y}_i$:

$$SST = \sum (y_i - \bar{y})^2 = \sum (y_i - \hat{y}_i + \hat{y}_i - \bar{y})^2$$

**Step 2:** Expand the square:

$$= \sum \left[(y_i - \hat{y}_i)^2 + 2(y_i - \hat{y}_i)(\hat{y}_i - \bar{y}) + (\hat{y}_i - \bar{y})^2\right]$$

**Step 3:** The cross term vanishes. We need to show:

$$\sum (y_i - \hat{y}_i)(\hat{y}_i - \bar{y}) = 0$$

**Proof that the cross term is zero:**

Let $e_i = y_i - \hat{y}_i$. Then:

$$\sum e_i (\hat{y}_i - \bar{y}) = \sum e_i (a_0 + a_1 x_i - \bar{y})$$

$$= (a_0 - \bar{y}) \sum e_i + a_1 \sum e_i x_i$$

We need two properties from the normal equations:

**Property 1:** $\sum e_i = 0$

$$\sum e_i = \sum (y_i - \hat{y}_i) = \sum (y_i - a_0 - a_1 x_i)$$

$$= \sum y_i - a_0 n - a_1 \sum x_i = n\left(\bar{y} - a_0 - a_1 \bar{x}\right) = 0 \quad \text{(by } (\ast\ast)\text{)}$$

**Property 2:** $\sum e_i x_i = 0$

$$\sum e_i x_i = \sum (y_i - \hat{y}_i) x_i = \sum y_i x_i - \sum (a_0 + a_1 x_i) x_i$$

$$= \sum y_i x_i - a_0 \sum x_i - a_1 \sum x_i^2 = 0 \quad \text{(by the second normal equation } (\ast)\text{)}$$

Therefore:

$$(a_0 - \bar{y}) \underbrace{\sum e_i}_{= 0} + a_1 \underbrace{\sum e_i x_i}_{= 0} = 0$$

**Conclusion:**

$$SST = \sum (y_i - \hat{y}_i)^2 + \sum (\hat{y}_i - \bar{y})^2 = SSE + SSR$$

<br>

### 5.3 Coefficient of Determination ($r^2$)

The **coefficient of determination** quantifies how well the regression line fits the data:

$$r^2 = \frac{SSR}{SST} = 1 - \frac{SSE}{SST}$$

| $r^2$ Value | Interpretation |
|-------------|----------------|
| $r^2 = 1$ | Perfect fit; all points lie on the line ($SSE = 0$) |
| $r^2 = 0$ | The model explains none of the variability ($SSR = 0$) |
| $r^2 \approx 1$ | Strong linear relationship |
| $r^2 \approx 0$ | Weak or no linear relationship |

> **[Statistics]** The value $r = \sqrt{r^2}$ (with appropriate sign) is the **Pearson correlation coefficient**. It measures the strength and direction of the linear relationship between $x$ and $y$. The sign of $r$ matches the sign of the slope $a_1$.

<br>

### 5.4 Standard Error of the Estimate

The **standard error of the estimate** quantifies the spread of data points around the regression line:

$$s_{y/x} = \sqrt{\frac{SSE}{n - 2}}$$

> **[Statistics]** We divide by $n - 2$ (not $n$) because we have estimated **two parameters** ($a_0$ and $a_1$) from the data, consuming 2 degrees of freedom. This is analogous to dividing by $n - 1$ in the sample standard deviation where 1 parameter (the mean) is estimated.

The standard error $s_{y/x}$ quantifies the typical size of a residual. A smaller $s_{y/x}$ indicates a tighter fit of the regression line to the data.

<br>

---

<br>

## 6. Summary Table

| Concept | Formula | Purpose |
|---------|---------|---------|
| Linear model | $\hat{y} = a_0 + a_1 x$ | Best-fit straight line |
| SSE (objective) | $\displaystyle \sum_{i=1}^{n}(y_i - a_0 - a_1 x_i)^2$ | Quantity to minimize |
| Normal equations (matrix) | $\begin{bmatrix} n & \sum x_i \\\\ \sum x_i & \sum x_i^2 \end{bmatrix} \begin{bmatrix} a_0 \\\\ a_1 \end{bmatrix} = \begin{bmatrix} \sum y_i \\\\ \sum y_i x_i \end{bmatrix}$ | System of equations for coefficients |
| Slope ($a_1$) | $\displaystyle \frac{n\sum x_i y_i - \sum x_i \sum y_i}{n\sum x_i^2 - (\sum x_i)^2}$ | Rate of change |
| Intercept ($a_0$) | $\bar{y} - a_1 \bar{x}$ | y-value when $x = 0$ |
| SST decomposition | $SST = SSE + SSR$ | Total = Unexplained + Explained |
| Coefficient of determination | $r^2 = SSR / SST = 1 - SSE / SST$ | Goodness of fit (0 to 1) |
| Standard error | $s_{y/x} = \sqrt{SSE / (n - 2)}$ | Spread of residuals |
| Sample std. dev. | $s_y = \sqrt{SST / (n - 1)}$ | Spread of data around mean |
| Coefficient of variation | $c.v. = (s_y / \bar{y}) \cdot 100\%$ | Normalized spread |
| MAD | $\text{median}(\|y_i - \tilde{y}\|)$ | Robust spread measure |
| Standard normal PDF | $\varphi(z) = \frac{1}{\sqrt{2\pi}} e^{-z^2/2}$ | Bell curve density |
| General normal PDF | $f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-(x-\mu)^2/(2\sigma^2)}$ | Normal density with mean $\mu$, variance $\sigma^2$ |
| Standardization | $Z = (X - \mu) / \sigma$ | Transform to $\mathcal{N}(0,1)$ |
