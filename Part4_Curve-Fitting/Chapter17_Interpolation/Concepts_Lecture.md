# Chapter 17 Lecture -- Polynomial Interpolation

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 17

> **Prerequisites**: [Calculus] Polynomials (Ch 14-16).
>
> **Learning Objectives**:
> 1. Apply Newton and Lagrange interpolation polynomials
> 2. Analyze interpolation error and Runge's phenomenon
> 3. Implement inverse interpolation for root finding

---

<br>

## Table of Contents

- [1. Polynomial Interpolation](#1-polynomial-interpolation)
  - [1.1 What Is Interpolation?](#11-what-is-interpolation)
  - [1.2 Polynomial Form and Unknowns](#12-polynomial-form-and-unknowns)
  - [1.3 Determining Coefficients](#13-determining-coefficients)
  - [1.4 Vandermonde Matrix](#14-vandermonde-matrix)
  - [1.5 Standardization](#15-standardization)
- [2. Newton Interpolating Polynomial](#2-newton-interpolating-polynomial)
  - [2.1 Linear Interpolation (1st Order)](#21-linear-interpolation-1st-order)
  - [2.2 Quadratic Interpolation (2nd Order)](#22-quadratic-interpolation-2nd-order)
  - [2.3 Finite Divided Differences](#23-finite-divided-differences)
  - [2.4 Recursive Nature of Finite Differences](#24-recursive-nature-of-finite-differences)
  - [2.5 General Form of Newton's Interpolating Polynomial](#25-general-form-of-newtons-interpolating-polynomial)
  - [2.6 Example 17.4 -- Estimating ln 2](#26-example-174----estimating-ln-2)
- [3. Lagrange Interpolating Polynomial](#3-lagrange-interpolating-polynomial)
  - [3.1 Linear Lagrange Interpolation (1st Order)](#31-linear-lagrange-interpolation-1st-order)
  - [3.2 Basis Functions N1(x) and N2(x)](#32-basis-functions-n1x-and-n2x)
  - [3.3 2nd-Order Lagrange Interpolating Polynomial](#33-2nd-order-lagrange-interpolating-polynomial)
  - [3.4 General (n-1)th-Order Lagrange Polynomial](#34-general-n-1th-order-lagrange-polynomial)
- [4. Inverse Interpolation](#4-inverse-interpolation)
  - [4.1 Concept](#41-concept)
  - [4.2 Example](#42-example)
- [5. Extrapolation and Oscillations](#5-extrapolation-and-oscillations)
  - [5.1 Extrapolation](#51-extrapolation)
  - [5.2 Oscillations (Runge Phenomenon)](#52-oscillations-runge-phenomenon)
- [Summary Table](#summary-table)

---

<br>

## 1. Polynomial Interpolation

### 1.1 What Is Interpolation?

Interpolation is the technique of estimating the value of a function at an intermediate point, given a set of discrete data points. The core question is: **How do we fill in the values between two known data points?**

The answer: use **polynomial interpolation**. The idea is to pass a polynomial through the known data points and then evaluate the polynomial at the desired intermediate location.

There are three common types of polynomial interpolation, depending on how many data points are used:

| Type | Data Points Required | Polynomial Order |
|------|---------------------|-----------------|
| Linear interpolation | 2 | 1st order |
| Quadratic interpolation | 3 | 2nd order |
| Cubic interpolation | 4 | 3rd order |

> **Motivating Example:** Suppose you are given a table of temperature $T$ [C] and density $\rho$ [kg/m$^3$] of air:
>
> | $T$ [C] | $\rho$ [kg/m$^3$] |
> |---------|-------------------|
> | $-40$   | $1.52$            |
> | $0$     | $1.29$            |
> | $20$    | $1.20$            |
> | $50$    | $1.09$            |
> | $\vdots$ | $\vdots$         |
> | $500$   | $0.46$            |
>
> What is the density at $T = 30$ C? Interpolation answers this question.

### 1.2 Polynomial Form and Unknowns

An $n$-th order polynomial has the general form:

$$f(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_n x^n$$

This polynomial has **$n + 1$ unknowns** (the coefficients $a_0, a_1, \ldots, a_n$), so we need **$n + 1$ data points** to uniquely determine them.

> **Key Principle:** To fit an $n$-th order polynomial, we need exactly $n + 1$ data points.

### 1.3 Determining Coefficients

**Example 17.1:** Suppose we want a 2nd-order polynomial through three data points:

$$f(x) = p_0 x^2 + p_1 x + p_2 \qquad (\ast)$$

Given the data:
- $f(300) = 0.616$
- $f(400) = 0.525$
- $f(500) = 0.457$

Substitute each data point into $(\ast)$:

$$0.616 = p_0 (300)^2 + p_1 (300) + p_2$$

$$0.525 = p_0 (400)^2 + p_1 (400) + p_2$$

$$0.457 = p_0 (500)^2 + p_1 (500) + p_2$$

This is a system of 3 equations with 3 unknowns ($p_0, p_1, p_2$).

### 1.4 Vandermonde Matrix

The system from Section 1.3 can be rewritten in **vector/matrix form**:

$$\begin{pmatrix} 0.616 \\ 0.525 \\ 0.457 \end{pmatrix} = \begin{pmatrix} 300^2 & 300 & 1 \\ 400^2 & 400 & 1 \\ 500^2 & 500 & 1 \end{pmatrix} \begin{pmatrix} p_0 \\ p_1 \\ p_2 \end{pmatrix}$$

In general, for a 2nd-order polynomial with data at $x_1, x_2, x_3$:

$$\begin{pmatrix} f(x_1) \\ f(x_2) \\ f(x_3) \end{pmatrix} = \underbrace{\begin{pmatrix} x_1^2 & x_1 & 1 \\ x_2^2 & x_2 & 1 \\ x_3^2 & x_3 & 1 \end{pmatrix}}_{V \text{ (Vandermonde matrix)}} \begin{pmatrix} p_0 \\ p_1 \\ p_2 \end{pmatrix}$$

The coefficient matrix $V$ is known as the **Vandermonde matrix**. Solve $V \mathbf{p} = \mathbf{f}$ to find the polynomial coefficients.

> **Problem:** The Vandermonde matrix is often **ill-conditioned** when the $x$-values are large. For the example above, $\text{cond}(V) \sim 10^6$, which causes numerical instability.

### 1.5 Standardization

To improve conditioning, define a **new standardized variable**:

$$z := \frac{x - 400}{100}$$

This maps the original $x$-values $\{300, 400, 500\}$ to $\{-1, 0, 1\}$, producing a new Vandermonde matrix:

$$V' = \begin{pmatrix} 1 & -1 & 1 \\ 0 & 0 & 1 \\ 1 & 1 & 1 \end{pmatrix}$$

Now $\text{cond}(V') = 3.26$, which is dramatically better than $10^6$.

> **Standardization** is necessary to improve the conditioning of the Vandermonde matrix and stabilize the computation.

---

<br>

## 2. Newton Interpolating Polynomial

Newton's approach builds the interpolating polynomial **incrementally**, adding one term at a time. Each new term introduces a higher-order correction (curvature).

### 2.1 Linear Interpolation (1st Order)

Given two data points $(x_1, f(x_1))$ and $(x_2, f(x_2))$, the **Newton linear interpolation formula** is:

$$\hat{f}(x) = f(x_1) + \frac{f(x_2) - f(x_1)}{x_2 - x_1}(x - x_1)$$

This is simply a straight line through the two points (slope-intercept form using the slope of the secant line).

> **Accuracy** increases as the interval $|x_2 - x_1|$ decreases and as $x$ is closer to the data points.

### 2.2 Quadratic Interpolation (2nd Order)

Given three data points $f(x_1), f(x_2), f(x_3)$, construct a quadratic interpolation formula:

$$\hat{f}(x) = b_1 + b_2(x - x_1) + b_3(x - x_1)(x - x_2)$$

**Determining the coefficients:**

**Step 1:** Set $\hat{f}(x_1) = f(x_1)$:

$$f(x_1) = b_1 \quad \Longrightarrow \quad b_1 = f(x_1)$$

**Step 2:** Set $\hat{f}(x_2) = f(x_2)$:

$$f(x_2) = f(x_1) + b_2(x_2 - x_1)$$

$$\therefore \quad b_2 = \frac{f(x_2) - f(x_1)}{x_2 - x_1} =: f[x_2, x_1]$$

This is the **first-order finite divided difference** (the slope).

**Step 3:** Set $\hat{f}(x_3) = f(x_3)$ and solve for $b_3$:

$$f(x_3) = f(x_1) + \frac{f(x_2) - f(x_1)}{x_2 - x_1}(x_3 - x_1) + b_3(x_3 - x_1)(x_3 - x_2)$$

After algebraic manipulation:

$$b_3 = \frac{\dfrac{f(x_3) - f(x_2)}{x_3 - x_2} - \dfrac{f(x_2) - f(x_1)}{x_2 - x_1}}{x_3 - x_1} = \frac{f[x_3, x_2] - f[x_2, x_1]}{x_3 - x_1} =: f[x_3, x_2, x_1]$$

This is the **second-order finite divided difference**.

### 2.3 Finite Divided Differences

Finite divided differences are defined recursively:

**Zeroth-order** (function values):

$$f[x_i] = f(x_i)$$

**First-order** (slopes):

$$f[x_i, x_j] = \frac{f(x_i) - f(x_j)}{x_i - x_j}$$

**Second-order** (curvature):

$$f[x_i, x_j, x_k] = \frac{f[x_i, x_j] - f[x_j, x_k]}{x_i - x_k}$$

**$n$-th order** (general):

$$f[x_n, x_{n-1}, \ldots, x_1] = \frac{f[x_n, x_{n-1}, \ldots, x_2] - f[x_{n-1}, x_{n-2}, \ldots, x_1]}{x_n - x_1}$$

### 2.4 Recursive Nature of Finite Differences

The divided differences can be organized in a **triangular table**:

| $x_i$ | $f(x_i)$ | 1st Divided Diff. | 2nd Divided Diff. | 3rd Divided Diff. |
|--------|-----------|-------------------|--------------------|---------------------|
| $x_1$ | $f(x_1) \; {}^{\prime\prime}b_1{}^{\prime\prime}$ | | | |
| | | $f[x_2, x_1] \; {}^{\prime\prime}b_2{}^{\prime\prime}$ | | |
| $x_2$ | $f(x_2)$ | | $f[x_3, x_2, x_1] \; {}^{\prime\prime}b_3{}^{\prime\prime}$ | |
| | | $f[x_3, x_2]$ | | $f[x_4, x_3, x_2, x_1] \; {}^{\prime\prime}b_4{}^{\prime\prime}$ |
| $x_3$ | $f(x_3)$ | | $f[x_4, x_3, x_2]$ | |
| | | $f[x_4, x_3]$ | | |
| $x_4$ | $f(x_4)$ | | | |

The **top diagonal** of the table gives the coefficients $b_1, b_2, b_3, b_4, \ldots$ of the Newton polynomial.

### 2.5 General Form of Newton's Interpolating Polynomial

For $n$ data points, the **(n-1)th-order** Newton interpolating polynomial is:

$$\hat{f}_{n-1}(x) = f(x_1) + f[x_2, x_1](x - x_1) + f[x_3, x_2, x_1](x - x_1)(x - x_2) + \cdots$$

$$+ f[x_n, x_{n-1}, \ldots, x_1](x - x_1)(x - x_2) \cdots (x - x_{n-1})$$

In compact notation:

$$\hat{f}_{n-1}(x) = \sum_{k=1}^{n} b_k \prod_{j=1}^{k-1}(x - x_j)$$

where

- $b_1 = f(x_1)$
- $b_2 = f[x_2, x_1]$
- $b_3 = f[x_3, x_2, x_1]$
- $b_k = f[x_k, x_{k-1}, \ldots, x_1]$

> **Advantage:** Newton's form is **incremental** -- adding a new data point only requires computing one additional divided difference and appending one more term. The previously computed coefficients remain unchanged.

### 2.6 Example 17.4 -- Estimating ln 2

**Problem:** Estimate $\ln 2$ using a 3rd-order Newton interpolating polynomial.

**Data points** (using $f(x) = \ln(x)$):

| $i$ | $x_i$ | $f(x_i) = \ln(x_i)$ |
|-----|--------|----------------------|
| 1   | 1      | 0                    |
| 2   | 4      | 1.386294             |
| 3   | 6      | 1.791759             |
| 4   | 5      | 1.609438             |

**Step 1:** Build the divided difference table:

| | $f(x_i)$ | 1st | 2nd | 3rd |
|---|----------|-----|-----|-----|
| $x_1 = 1$ | $0 = b_1$ | | | |
| | | $f[x_2, x_1]$ | | |
| $x_2 = 4$ | $1.386294$ | | $f[x_3, x_2, x_1]$ | |
| | | $f[x_3, x_2]$ | | $f[x_4, x_3, x_2, x_1] = b_4$ |
| $x_3 = 6$ | $1.791759$ | | $f[x_4, x_3, x_2]$ | |
| | | $f[x_4, x_3]$ | | |
| $x_4 = 5$ | $1.609438$ | | | |

**Step 2:** Compute each divided difference:

$$f[x_2, x_1] = \frac{1.386294 - 0}{4 - 1} = 0.462098 = b_2$$

$$f[x_3, x_2] = \frac{1.791759 - 1.386294}{6 - 4} = 0.202733$$

$$f[x_4, x_3] = \frac{1.609438 - 1.791759}{5 - 6} = 0.182321$$

$$f[x_3, x_2, x_1] = \frac{0.202733 - 0.462098}{6 - 1} = -0.051873 = b_3$$

$$f[x_4, x_3, x_2] = \frac{0.182321 - 0.202733}{5 - 4} = -0.020412$$

$$f[x_4, x_3, x_2, x_1] = \frac{-0.020412 - (-0.051873)}{5 - 1} = 0.007865 = b_4$$

**Step 3:** Construct the 3rd-order polynomial:

$$\hat{f}_3(x) = 0 + 0.462098(x - 1) + (-0.051873)(x - 1)(x - 4) + 0.007865(x - 1)(x - 4)(x - 6)$$

**Step 4:** Evaluate at $x = 2$:

$$\hat{f}_3(2) = 0.462098(1) + (-0.051873)(1)(-2) + 0.007865(1)(-2)(-4)$$

$$= 0.462098 + 0.103746 + 0.062920 = 0.628764$$

The exact value is $\ln 2 = 0.693147$, so the estimate has some error due to the spacing of data points.

```python
import numpy as np

# Data points
x = np.array([1, 4, 6, 5])
f = np.array([0, 1.386294, 1.791759, 1.609438])

# Divided differences
def divided_diff(x, f):
    n = len(x)
    table = np.zeros((n, n))
    table[:, 0] = f
    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (x[i+j] - x[i])
    return table[0]  # top row = coefficients b1, b2, b3, ...

# Newton evaluation
def newton_eval(coeffs, x_data, x_eval):
    n = len(coeffs)
    result = coeffs[-1]
    for k in range(n - 2, -1, -1):
        result = result * (x_eval - x_data[k]) + coeffs[k]
    return result

coeffs = divided_diff(x, f)
print(f"Coefficients: {coeffs}")
print(f"f(2) = {newton_eval(coeffs, x, 2):.6f}")
print(f"Exact ln(2) = {np.log(2):.6f}")
```

---

<br>

## 3. Lagrange Interpolating Polynomial

Lagrange's method provides an **alternative** formulation of the interpolating polynomial that does not require solving a system of equations or building a divided difference table. Instead, it directly constructs the polynomial using **basis functions** (also called **shape functions**).

### 3.1 Linear Lagrange Interpolation (1st Order)

Given two data points $(x_1, f(x_1))$ and $(x_2, f(x_2))$, the straight line through them is:

$$y - f(x_1) = \frac{f(x_2) - f(x_1)}{x_2 - x_1}(x - x_1)$$

Rearranging:

$$y = \frac{x_2 - x}{x_2 - x_1} f(x_1) + \frac{x - x_1}{x_2 - x_1} f(x_2)$$

Define the **basis functions** (shape functions):

$$N_1(x) = \frac{x_2 - x}{x_2 - x_1}, \qquad N_2(x) = \frac{x - x_1}{x_2 - x_1}$$

Then the **linear Lagrange interpolating polynomial** is:

$$\hat{f}_1(x) = N_1(x) \, f(x_1) + N_2(x) \, f(x_2)$$

### 3.2 Basis Functions N1(x) and N2(x)

The basis functions have the important **partition of unity** property:

- $N_1(x_1) = 1, \quad N_1(x_2) = 0$
- $N_2(x_1) = 0, \quad N_2(x_2) = 1$

At $x = x_1$: the interpolation gives $f(x_1)$ exactly (weighted 100% by $N_1$).
At $x = x_2$: the interpolation gives $f(x_2)$ exactly (weighted 100% by $N_2$).
Between $x_1$ and $x_2$: the result is a **weighted average** of $f(x_1)$ and $f(x_2)$.

> **Geometric Interpretation:** $N_1(x)$ and $N_2(x)$ are linear functions that cross at an intermediate point. The interpolated value $y = N_1(x) f(x_1) + N_2(x) f(x_2)$ is the weighted combination of function values at the data points.

### 3.3 2nd-Order Lagrange Interpolating Polynomial

For three data points $(x_1, f(x_1)), (x_2, f(x_2)), (x_3, f(x_3))$, define three Lagrange basis functions:

$$L_1(x) = \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)}$$

$$L_2(x) = \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)}$$

$$L_3(x) = \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)}$$

The 2nd-order Lagrange interpolating polynomial is:

$$\hat{f}_2(x) = L_1(x) \, f(x_1) + L_2(x) \, f(x_2) + L_3(x) \, f(x_3) = \sum_{i=1}^{3} L_i(x) \, f(x_i)$$

Each $L_i(x)$ satisfies:

$$L_i(x_j) = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{if } i \neq j \end{cases}$$

> **Note:** 3 data points are needed for a 2nd-order polynomial.

### 3.4 General (n-1)th-Order Lagrange Polynomial

For $n$ data points $x_1, x_2, \ldots, x_n$, the **(n-1)th-order Lagrange interpolating polynomial** is:

$$\hat{f}_{n-1}(x) = \sum_{i=1}^{n} L_i(x) \, f(x_i)$$

where the Lagrange basis functions are:

$$L_i(x) = \prod_{\substack{j=1 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}$$

Here:
- $n$ is the **number of data points**
- $\prod$ denotes the **product** over all $j$ except $j = i$

```python
import numpy as np

def lagrange_interp(x_data, f_data, x_eval):
    """Lagrange interpolation at a single point x_eval."""
    n = len(x_data)
    result = 0.0
    for i in range(n):
        # Compute L_i(x_eval)
        Li = 1.0
        for j in range(n):
            if j != i:
                Li *= (x_eval - x_data[j]) / (x_data[i] - x_data[j])
        result += Li * f_data[i]
    return result

# Example: Estimate ln(2) using the same 4 data points
x_data = np.array([1, 4, 6, 5], dtype=float)
f_data = np.log(x_data)

estimate = lagrange_interp(x_data, f_data, 2.0)
print(f"Lagrange f(2) = {estimate:.6f}")
print(f"Exact ln(2) = {np.log(2):.6f}")
```

> **Newton vs. Lagrange:** Both methods produce the **same** interpolating polynomial -- they are algebraically equivalent. Newton's form is more efficient for **incrementally adding** new data points. Lagrange's form is conceptually simpler and more straightforward to implement.

---

<br>

## 4. Inverse Interpolation

### 4.1 Concept

In standard interpolation, we are given $x$ and estimate $f(x)$. In **inverse interpolation**, we are given a target value $y$ and want to find $x$ such that $f(x) = y$.

Given data points $(x_0, f(x_0)), (x_1, f(x_1)), (x_2, f(x_2)), \ldots$, we want to find $x$ such that $f(x) = y$.

**Key idea:** No explicit inverse function $f^{-1}(x)$ is needed. Instead, **reverse the roles** of $x$ and $f(x)$ -- treat $f(x_i)$ as the independent variable and $x_i$ as the dependent variable, then apply standard interpolation.

$$x(y) = \sum_{i=0}^{n} x_i \, L_i(y)$$

where $L_i(y)$ are Lagrange basis functions constructed using $f(x_i)$ values as the "x-data."

### 4.2 Example

Given three data points:

$$
(x_0, f(x_0)) = (2, 1), \quad
(x_1, f(x_1)) = (4, 2), \quad
(x_2, f(x_2)) = (8, 3)
$$

**Goal:** Find $x$ such that $f(x) = 5$.

**Step 1:** Reverse the pairs:

$$(f(x_0), x_0) = (1, 2), \quad (f(x_1), x_1) = (2, 4), \quad (f(x_2), x_2) = (3, 8)$$

**Step 2:** Apply Lagrange interpolation with $y$-values as the independent variable:

$$x(y) = \sum_{i=0}^{2} x_i \, L_i(y)$$

**Step 3:** Evaluate at $y = 5$:

$$x(5) = \sum_{i=0}^{2} x_i \, L_i(5)$$

Compute each basis function at $y = 5$:

$$L_0(5) = \frac{(5-2)(5-3)}{(1-2)(1-3)} = \frac{(3)(2)}{(-1)(-2)} = 3$$

$$L_1(5) = \frac{(5-1)(5-3)}{(2-1)(2-3)} = \frac{(4)(2)}{(1)(-1)} = -8$$

$$L_2(5) = \frac{(5-1)(5-2)}{(3-1)(3-2)} = \frac{(4)(3)}{(2)(1)} = 6$$

$$x(5) = 2 \cdot 3 + 4 \cdot (-8) + 8 \cdot 6 = 6 - 32 + 48 = 22$$

```python
import numpy as np

# Inverse interpolation example
y_data = np.array([1, 2, 3], dtype=float)   # f(x) values as independent variable
x_data = np.array([2, 4, 8], dtype=float)   # x values as dependent variable

def lagrange_interp(ind, dep, target):
    n = len(ind)
    result = 0.0
    for i in range(n):
        Li = 1.0
        for j in range(n):
            if j != i:
                Li *= (target - ind[j]) / (ind[i] - ind[j])
        result += Li * dep[i]
    return result

x_at_y5 = lagrange_interp(y_data, x_data, 5.0)
print(f"x such that f(x)=5: {x_at_y5:.2f}")
```

---

<br>

## 5. Extrapolation and Oscillations

### 5.1 Extrapolation

**Extrapolation** means estimating a value of $f(x)$ that lies **outside** the range of known data points $[x_1, x_n]$.

- **Interpolation:** $x \in [x_1, x_n]$ (between known points)
- **Extrapolation:** $x < x_1$ or $x > x_n$ (outside known points)

> **Warning:** Extrapolation using polynomial interpolation is generally **unreliable**. The polynomial may diverge rapidly outside the data range because there is no data to constrain its behavior in those regions. The error can grow dramatically.

### 5.2 Oscillations (Runge Phenomenon)

Using a **high-order polynomial** for interpolation can introduce **spurious oscillations**, especially near the edges of the data interval. This is known as **Runge's phenomenon**.

**Example 17.7 -- Runge Function:**

$$f(x) = \frac{1}{1 + 25x^2}, \qquad x \in [-1, 1]$$

When this function is interpolated with a high-order polynomial (e.g., 4th order or higher) using equally spaced points, the polynomial develops large oscillations near the boundaries $x = \pm 1$, even though it passes exactly through all the data points.

> **Runge's Phenomenon:** As the polynomial order increases with equally spaced nodes, the interpolation error can actually **increase** near the endpoints, rather than decrease. This is a fundamental limitation of high-order polynomial interpolation.

**Remedies for oscillation:**
1. Use **lower-order** piecewise polynomials (splines) instead of a single high-order polynomial
2. Use **non-equally spaced** nodes (e.g., Chebyshev nodes) that cluster near the boundaries
3. Keep the polynomial order reasonably low and use more intervals

```python
import numpy as np
import matplotlib.pyplot as plt

# Runge function demonstration
f_runge = lambda x: 1 / (1 + 25 * x**2)

x_fine = np.linspace(-1, 1, 500)
y_true = f_runge(x_fine)

for n in [5, 9, 13]:
    x_nodes = np.linspace(-1, 1, n)
    y_nodes = f_runge(x_nodes)
    # Lagrange interpolation
    y_interp = np.zeros_like(x_fine)
    for i in range(n):
        Li = np.ones_like(x_fine)
        for j in range(n):
            if j != i:
                Li *= (x_fine - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        y_interp += Li * y_nodes[i]
    plt.plot(x_fine, y_interp, label=f"Order {n-1}")

plt.plot(x_fine, y_true, 'k--', linewidth=2, label="True f(x)")
plt.ylim(-0.5, 1.5)
plt.legend()
plt.title("Runge Phenomenon")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.show()
```

---

<br>

## Summary Table

| Topic | Key Formula / Concept | Points Needed | Order |
|-------|----------------------|---------------|-------|
| **Polynomial Interpolation** | $f(x) = a_0 + a_1 x + \cdots + a_n x^n$ | $n+1$ | $n$ |
| **Vandermonde Matrix** | $V \mathbf{p} = \mathbf{f}$ | Solves for coefficients | Often ill-conditioned |
| **Standardization** | $z = \frac{x - \bar{x}}{\Delta x}$ | Improves $\text{cond}(V)$ | Same order |
| **Newton Linear** | $\hat{f}(x) = f(x_1) + f[x_2,x_1](x - x_1)$ | 2 | 1 |
| **Newton Quadratic** | $+ \, f[x_3,x_2,x_1](x-x_1)(x-x_2)$ | 3 | 2 |
| **Newton General** | $\hat{f}_{n-1}(x) = \sum b_k \prod_{j=1}^{k-1}(x - x_j)$ | $n$ | $n-1$ |
| **Divided Differences** | $f[x_i,x_j] = \frac{f(x_i)-f(x_j)}{x_i - x_j}$ | Recursive | Builds table |
| **Lagrange Linear** | $\hat{f}_1 = N_1(x)f(x_1) + N_2(x)f(x_2)$ | 2 | 1 |
| **Lagrange General** | $\hat{f}_{n-1}(x) = \sum L_i(x) f(x_i)$ | $n$ | $n-1$ |
| **Lagrange Basis** | $L_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}$ | Kronecker delta property | -- |
| **Inverse Interpolation** | Swap roles of $x$ and $f(x)$, then interpolate | Same as forward | Same |
| **Extrapolation** | Evaluation outside $[x_1, x_n]$ | -- | Unreliable |
| **Oscillations (Runge)** | $f(x) = \frac{1}{1+25x^2}$ | High order + equal spacing | Diverges at edges |
