# Chapter 4 Lecture — Truncation Errors and the Taylor Series

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 4

> **Prerequisites**: [Calculus] Taylor series expansion (Ch 3).
>
> **Learning Objectives**:
> 1. Estimate truncation errors using Taylor series
> 2. Apply the remainder term for error bounds
> 3. Analyze the order of approximation methods

---

<br>

## Table of Contents

- [1. Taylor Series Expansion](#1-taylor-series-expansion)
  - [1.1 The Taylor Series Formula](#11-the-taylor-series-formula)
  - [1.2 The Remainder Term (Truncation Error)](#12-the-remainder-term-truncation-error)
- [2. Order of Approximation](#2-order-of-approximation)
  - [2.1 Zeroth Through nth Order](#21-zeroth-through-nth-order)
  - [2.2 Big-O Notation and Its Meaning](#22-big-o-notation-and-its-meaning)
  - [2.3 Example: Taylor Series of $e^x$](#23-example-taylor-series-of-ex)
- [3. Numerical Differentiation Using Taylor Series](#3-numerical-differentiation-using-taylor-series)
  - [3.1 Forward Difference (First Derivative)](#31-forward-difference-first-derivative)
  - [3.2 Backward Difference (First Derivative)](#32-backward-difference-first-derivative)
  - [3.3 Central Difference (First Derivative)](#33-central-difference-first-derivative)
  - [3.4 Second Derivative (Central Difference)](#34-second-derivative-central-difference)
  - [3.5 Higher-Order Forward Difference](#35-higher-order-forward-difference)
- [4. Detailed Example: Numerical Differentiation](#4-detailed-example-numerical-differentiation)
  - [4.1 Problem Setup](#41-problem-setup)
  - [4.2 Results with $h = 0.5$](#42-results-with-h--05)
  - [4.3 Effect of Reducing $h$](#43-effect-of-reducing-h)
- [5. Total Numerical Error](#5-total-numerical-error)
  - [5.1 Truncation Error vs. Round-Off Error](#51-truncation-error-vs-round-off-error)
  - [5.2 Optimal Step Size](#52-optimal-step-size)
- [6. Stability and Condition Number](#6-stability-and-condition-number)
- [7. Lab: Finite Differences in Python](#7-lab-finite-differences-in-python)
  - [7.1 Comparing Forward and Central Differences](#71-comparing-forward-and-central-differences)
  - [7.2 Step Size Analysis and Optimal $h$](#72-step-size-analysis-and-optimal-h)
- [Summary](#summary)

---

<br>

## 1. Taylor Series Expansion

### 1.1 The Taylor Series Formula

The **Taylor series** is one of the most fundamental tools in numerical analysis. It allows us to represent any sufficiently smooth function as an infinite sum of polynomial terms, each constructed from the function's derivatives evaluated at a single point.

Given a function $f(x)$ that is $(n+1)$ times continuously differentiable, the value of $f$ at a nearby point $x_{i+1}$ can be expressed in terms of $f$ and its derivatives at $x_i$:

$$f(x_{i+1}) = f(x_i) + f'(x_i)(x_{i+1} - x_i) + \frac{f''(x_i)}{2!}(x_{i+1} - x_i)^2 + \frac{f'''(x_i)}{3!}(x_{i+1} - x_i)^3 + \cdots + R_n$$

Defining the **step size** $h = x_{i+1} - x_i$, this becomes the more compact form:

$$f(x_i + h) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots + \frac{f^{(n)}(x_i)}{n!}h^n + R_n$$

Each successive term in the series includes a higher-order derivative of $f$ at $x_i$, divided by the corresponding factorial, and multiplied by a higher power of $h$. The more terms we include, the better the approximation — provided $h$ is small enough.

> **[Calculus]** Taylor's theorem is one of the most fundamental results in calculus. It states that any sufficiently smooth function can be approximated by a polynomial. The remainder term $R_n$ quantifies exactly how good this approximation is. The unknown $\xi$ makes $R_n$ impractical to compute exactly, but we can bound it using the maximum value of $f^{(n+1)}$ on the interval.

### 1.2 The Remainder Term (Truncation Error)

When we truncate the Taylor series after $n$ terms, the omitted terms are collectively represented by the **remainder** $R_n$:

$$R_n = \frac{f^{(n+1)}(\xi)}{(n+1)!} h^{n+1}$$

where $\xi$ is some (unknown) value satisfying $x_i \leq \xi \leq x_{i+1}$.

Key properties of the remainder:

- $R_n$ is the **truncation error** — the error introduced by truncating an infinite series to a finite number of terms
- The exact value of $\xi$ is generally unknown, which means we cannot compute $R_n$ exactly
- However, we can **bound** $R_n$ by finding the maximum of $|f^{(n+1)}(x)|$ on the interval $[x_i, x_{i+1}]$:
  $$|R_n| \leq \frac{\max_{x_i \leq x \leq x_{i+1}} |f^{(n+1)}(x)|}{(n+1)!} |h|^{n+1}$$
- For small $h$, the dominant behavior of $R_n$ is proportional to $h^{n+1}$, which we write as $R_n = O(h^{n+1})$

---

<br>

## 2. Order of Approximation

### 2.1 Zeroth Through nth Order

Truncating the Taylor series at different points yields approximations of different orders. Each additional term improves the approximation:

| Order | Approximation | Truncation Error |
|:------|:-------------|:-----------------|
| 0th | $f(x_{i+1}) \approx f(x_i)$ (constant) | $O(h)$ |
| 1st | $f(x_{i+1}) \approx f(x_i) + f'(x_i)h$ (linear) | $O(h^2)$ |
| 2nd | adds $\frac{f''(x_i)}{2!}h^2$ term (quadratic) | $O(h^3)$ |
| 3rd | adds $\frac{f'''(x_i)}{3!}h^3$ term (cubic) | $O(h^4)$ |
| nth | includes up to $\frac{f^{(n)}(x_i)}{n!}h^n$ term | $O(h^{n+1})$ |

### 2.2 Big-O Notation and Its Meaning

The **Big-O notation** $O(h^{n+1})$ describes how the truncation error scales with the step size $h$. Specifically, $O(h^{n+1})$ means that for sufficiently small $h$, the error is bounded by a constant times $h^{n+1}$:

$$|\text{error}| \leq C \cdot h^{n+1} \quad \text{for some constant } C$$

The practical consequence is:

- **0th-order** ($O(h)$): halving $h$ halves the error
- **1st-order** ($O(h^2)$): halving $h$ reduces the error by a factor of 4
- **2nd-order** ($O(h^3)$): halving $h$ reduces the error by a factor of 8
- **nth-order** ($O(h^{n+1})$): halving $h$ reduces the error by a factor of $2^{n+1}$

> **Note:** The Big-O notation $O(h^{n+1})$ means the error is proportional to $h^{n+1}$ for small $h$. Halving $h$ in a first-order method halves the error; halving $h$ in a second-order method quarters the error.

### 2.3 Example: Taylor Series of $e^x$

Consider approximating $e^1$ using the Taylor series expansion of $f(x) = e^x$ about $x_i = 0$ (i.e., $h = 1$):

Since $f(x) = e^x$, all derivatives are also $e^x$, so $f(0) = f'(0) = f''(0) = \cdots = 1$.

The Taylor series becomes:

$$e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \frac{x^4}{4!} + \cdots$$

Evaluating at $x = 1$ ($h = 1$) and adding terms one at a time:

| Order | Approximation | Value | True Error | $\varepsilon_t$ |
|:------|:-------------|:------|:-----------|:----------------|
| 0th | $f \approx 1$ | 1.000000 | 1.71828 | 63.2% |
| 1st | $f \approx 1 + 1$ | 2.000000 | 0.71828 | 26.4% |
| 2nd | $f \approx 2 + 0.5$ | 2.500000 | 0.21828 | 8.0% |
| 3rd | $f \approx 2.5 + 1/6$ | 2.666667 | 0.05161 | 1.9% |
| 4th | $f \approx 2.667 + 1/24$ | 2.708333 | 0.00995 | 0.366% |
| 5th | $f \approx 2.708 + 1/120$ | 2.716667 | 0.00161 | 0.0594% |
| 6th | $f \approx 2.717 + 1/720$ | 2.718056 | 0.00023 | 0.00842% |

The true value is $e = 2.718282\ldots$. Each additional term significantly reduces the error, demonstrating the convergence of the Taylor series. By the 6th-order approximation, the relative error is already below $0.01\%$.

---

<br>

## 3. Numerical Differentiation Using Taylor Series

This section derives the finite difference formulas that are used throughout numerical analysis. These formulas approximate derivatives using discrete function values — a necessity when working with tabulated data or when analytical differentiation is impractical.

### 3.1 Forward Difference (First Derivative)

Starting from the Taylor series expansion:

$$f(x_{i+1}) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots$$

Solving for $f'(x_i)$:

$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \frac{f''(x_i)}{2!}h - \frac{f'''(x_i)}{3!}h^2 - \cdots$$

Truncating the higher-order terms gives the **first-order forward difference** formula:

$$\boxed{f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} + O(h)}$$

The leading error term is $-\frac{f''(x_i)}{2}h$, so the truncation error is $O(h)$ — this is a **first-order accurate** approximation.

### 3.2 Backward Difference (First Derivative)

Expanding $f(x_{i-1})$ in a Taylor series about $x_i$ (substituting $-h$ for $h$):

$$f(x_{i-1}) = f(x_i) - f'(x_i)h + \frac{f''(x_i)}{2!}h^2 - \frac{f'''(x_i)}{3!}h^3 + \cdots$$

Solving for $f'(x_i)$:

$$\boxed{f'(x_i) = \frac{f(x_i) - f(x_{i-1})}{h} + O(h)}$$

Like the forward difference, this is also **first-order accurate** with $O(h)$ truncation error. The backward difference uses the point behind $x_i$ rather than ahead.

### 3.3 Central Difference (First Derivative)

The central difference is derived by subtracting the backward Taylor expansion from the forward Taylor expansion.

**Forward expansion:**

$$f(x_{i+1}) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots$$

**Backward expansion:**

$$f(x_{i-1}) = f(x_i) - f'(x_i)h + \frac{f''(x_i)}{2!}h^2 - \frac{f'''(x_i)}{3!}h^3 + \cdots$$

**Subtracting** backward from forward:

$$f(x_{i+1}) - f(x_{i-1}) = 2f'(x_i)h + \frac{2f'''(x_i)}{3!}h^3 + \cdots$$

Note that all even-powered terms ($h^2, h^4, \ldots$) cancel. Solving for $f'(x_i)$:

$$\boxed{f'(x_i) = \frac{f(x_{i+1}) - f(x_{i-1})}{2h} + O(h^2)}$$

The truncation error is $O(h^2)$ — this is **second-order accurate**, which is significantly more accurate than either the forward or backward difference.

> **[Calculus]** The central difference is more accurate because the even-powered error terms cancel when we subtract the backward expansion from the forward expansion. This is a general principle: symmetric formulas tend to be more accurate.

### 3.4 Second Derivative (Central Difference)

The central difference formula for the second derivative is derived by **adding** the forward and backward Taylor expansions instead of subtracting them.

**Adding** forward and backward:

$$f(x_{i+1}) + f(x_{i-1}) = 2f(x_i) + 2\frac{f''(x_i)}{2!}h^2 + 2\frac{f^{(4)}(x_i)}{4!}h^4 + \cdots$$

Note that all odd-powered terms ($h, h^3, h^5, \ldots$) cancel. Solving for $f''(x_i)$:

$$f''(x_i) = \frac{f(x_{i+1}) + f(x_{i-1}) - 2f(x_i)}{h^2} - \frac{f^{(4)}(x_i)}{12}h^2 - \cdots$$

Truncating:

$$\boxed{f''(x_i) = \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2} + O(h^2)}$$

This is a **second-order accurate** approximation of the second derivative, requiring three function evaluations at equally spaced points.

### 3.5 Higher-Order Forward Difference

It is possible to achieve higher-order accuracy using only forward points (useful at the boundary of a domain where backward points are unavailable).

Using three forward points $f(x_i)$, $f(x_{i+1})$, and $f(x_{i+2})$, we can derive the **second-order forward difference** for $f'(x_i)$:

$$\boxed{f'(x_i) = \frac{-f(x_{i+2}) + 4f(x_{i+1}) - 3f(x_i)}{2h} + O(h^2)}$$

This achieves $O(h^2)$ accuracy using only forward-direction points — the same order as the central difference, but without requiring points behind $x_i$.

---

<br>

## 4. Detailed Example: Numerical Differentiation

### 4.1 Problem Setup

Consider the polynomial:

$$f(x) = -0.1x^4 - 0.15x^3 - 0.5x^2 - 0.25x + 1.2$$

The analytically computed (true) derivative is:

$$f'(x) = -0.4x^3 - 0.45x^2 - x - 0.25$$

At $x = 0.5$:

$$f'(0.5) = -0.4(0.125) - 0.45(0.25) - 0.5 - 0.25 = -0.9125$$

We also need the function values:
- $f(0.0) = 1.2$
- $f(0.25) = 1.10351563$
- $f(0.5) = 0.925$
- $f(0.75) = 0.63632813$
- $f(1.0) = 0.2$

### 4.2 Results with $h = 0.5$

Using step size $h = 0.5$ to approximate $f'(0.5)$:

**Forward difference:**

$$f'(0.5) \approx \frac{f(1.0) - f(0.5)}{0.5} = \frac{0.2 - 0.925}{0.5} = -1.155$$

$$\varepsilon_t = \left|\frac{-0.9125 - (-1.155)}{-0.9125}\right| \times 100\% = 26.6\%$$

**Backward difference:**

$$f'(0.5) \approx \frac{f(0.5) - f(0.0)}{0.5} = \frac{0.925 - 1.2}{0.5} = -0.714$$

$$\varepsilon_t = \left|\frac{-0.9125 - (-0.714)}{-0.9125}\right| \times 100\% = 21.7\%$$

**Central difference:**

$$f'(0.5) \approx \frac{f(1.0) - f(0.0)}{2 \times 0.5} = \frac{0.2 - 1.2}{1.0} = -0.934$$

$$\varepsilon_t = \left|\frac{-0.9125 - (-0.934)}{-0.9125}\right| \times 100\% = 2.4\%$$

The central difference is dramatically more accurate than both forward and backward differences, consistent with its $O(h^2)$ error versus $O(h)$.

### 4.3 Effect of Reducing $h$

Using a smaller step size $h = 0.25$ with the central difference:

$$f'(0.5) \approx \frac{f(0.75) - f(0.25)}{2 \times 0.25} = \frac{0.63632813 - 1.10351563}{0.5} = -0.934375$$

$$\varepsilon_t \approx 0.6\%$$

Reducing $h$ from $0.5$ to $0.25$ (a factor of 2) reduced the central difference error by roughly a factor of 4, confirming the $O(h^2)$ behavior.

| Method | $h = 0.5$ | $h = 0.25$ | Error Reduction |
|:-------|:----------|:-----------|:----------------|
| Forward | $\varepsilon_t = 26.6\%$ | $\varepsilon_t \approx 13\%$ | $\sim 2\times$ (confirms $O(h)$) |
| Backward | $\varepsilon_t = 21.7\%$ | $\varepsilon_t \approx 11\%$ | $\sim 2\times$ (confirms $O(h)$) |
| Central | $\varepsilon_t = 2.4\%$ | $\varepsilon_t \approx 0.6\%$ | $\sim 4\times$ (confirms $O(h^2)$) |

---

<br>

## 5. Total Numerical Error

### 5.1 Truncation Error vs. Round-Off Error

When performing numerical differentiation on a computer, the total error has two components:

$$\text{Total error} = \text{Truncation error} + \text{Round-off error}$$

These two error sources have **opposite** dependencies on the step size $h$:

| Component | As $h$ decreases | Cause |
|:----------|:-----------------|:------|
| Truncation error | **Decreases** | Finite difference better approximates the derivative |
| Round-off error | **Increases** | Subtracting nearly equal floating-point numbers (catastrophic cancellation) |

When $h$ is large, the truncation error dominates. When $h$ is very small, the function values $f(x+h)$ and $f(x)$ become nearly identical, and their difference is contaminated by floating-point round-off. The total error curve has a characteristic **V-shape** on a log-log plot.

> **Note:** This is one of the most important practical insights in numerical analysis. You cannot make $h$ arbitrarily small to improve accuracy — at some point, round-off error takes over and the result actually gets WORSE.

### 5.2 Optimal Step Size

There exists an **optimal step size** $h^*$ that minimizes the total error. The optimal value depends on the finite difference formula being used and the machine epsilon $\varepsilon_{mach}$:

**For forward difference** ($O(h)$ truncation error):

$$h^* \approx \sqrt{\varepsilon_{mach} \cdot \left|\frac{f(x)}{f''(x)}\right|}$$

Since $\varepsilon_{mach} \approx 2.2 \times 10^{-16}$ for 64-bit double precision:

$$h^* \sim \sqrt{\varepsilon_{mach}} \approx 1.5 \times 10^{-8}$$

**For central difference** ($O(h^2)$ truncation error):

$$h^* \sim \varepsilon_{mach}^{1/3} \approx 6 \times 10^{-6}$$

The central difference has a larger optimal $h$ because its truncation error decreases faster ($O(h^2)$ vs. $O(h)$), so the balance point between truncation and round-off shifts to a larger step size.

| Formula | Truncation Error | Optimal $h^*$ | Minimum Total Error |
|:--------|:----------------|:---------------|:--------------------|
| Forward difference | $O(h)$ | $\sim \sqrt{\varepsilon_{mach}} \approx 10^{-8}$ | $\sim \sqrt{\varepsilon_{mach}} \approx 10^{-8}$ |
| Central difference | $O(h^2)$ | $\sim \varepsilon_{mach}^{1/3} \approx 10^{-6}$ | $\sim \varepsilon_{mach}^{2/3} \approx 10^{-11}$ |

---

<br>

## 6. Stability and Condition Number

The **condition number** of a function quantifies how sensitive the output is to small perturbations in the input. For a function $f(x)$, the condition number is:

$$\text{Cond} = \left|\frac{x \cdot f'(x)}{f(x)}\right|$$

The condition number measures the ratio of the relative change in the output to the relative change in the input:

- If $\text{Cond} \gg 1$: the problem is **ill-conditioned** — small input errors produce large output errors. Numerical results may be unreliable regardless of the algorithm used.
- If $\text{Cond} \approx 1$: the problem is **well-conditioned** — input errors are not amplified significantly.

For example, consider $f(x) = \frac{1}{x - 1}$ near $x = 1.00001$:

$$\text{Cond} = \left|\frac{x \cdot (-1/(x-1)^2)}{1/(x-1)}\right| = \left|\frac{x}{x-1}\right|$$

At $x = 1.00001$, $\text{Cond} \approx 10^5$ — extremely ill-conditioned. A tiny perturbation in $x$ causes a massive change in $f(x)$.

> **[Linear Algebra]** The condition number concept extends to matrices ($\kappa(A) = \|A\| \cdot \|A^{-1}\|$), which we'll study in Chapter 11. Both measure sensitivity to perturbations.

---

<br>

## 7. Lab: Finite Differences in Python

### 7.1 Comparing Forward and Central Differences

The following code compares the error of forward and central difference formulas across a range of step sizes, demonstrating both the accuracy advantage of central differences and the effect of round-off error at very small step sizes:

```python
import numpy as np
import matplotlib.pyplot as plt

# Test function
f = lambda x: -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x + 1.2
df_true = lambda x: -0.4*x**3 - 0.45*x**2 - x - 0.25

x = 0.5
true_val = df_true(x)

# Compare finite differences at various step sizes
h_values = np.logspace(-1, -15, 15)
E_fwd = np.zeros(len(h_values))
E_cen = np.zeros(len(h_values))

for i, h in enumerate(h_values):
    fwd = (f(x + h) - f(x)) / h
    cen = (f(x + h) - f(x - h)) / (2 * h)
    E_fwd[i] = abs(true_val - fwd)
    E_cen[i] = abs(true_val - cen)

plt.loglog(h_values, E_fwd, 'bo-', label='Forward O(h)')
plt.loglog(h_values, E_cen, 'rs-', label='Central O(h²)')
plt.xlabel('Step size h')
plt.ylabel('|Error|')
plt.title('Finite Difference Error vs Step Size')
plt.legend()
plt.grid(True)
plt.show()
```

The resulting log-log plot shows two key features:

1. **For large $h$**: both errors decrease as $h$ decreases, with central difference error decreasing faster (slope $\approx -2$ vs. $-1$)
2. **For very small $h$**: both errors increase due to round-off, forming the characteristic V-shape. The minimum error for central difference is much lower and occurs at a larger $h$.

### 7.2 Step Size Analysis and Optimal $h$

The following code systematically finds the optimal step size for the forward difference:

```python
def step_size_analysis(func, dfunc, x, n_steps=15):
    """Analyze error vs step size for forward difference."""
    h_vals = np.logspace(-1, -15, n_steps)
    errors = np.zeros(n_steps)
    true_deriv = dfunc(x)

    for i, h in enumerate(h_vals):
        approx = (func(x + h) - func(x)) / h
        errors[i] = abs(true_deriv - approx)

    # Find optimal h
    i_opt = np.argmin(errors)
    print(f'Optimal h ≈ {h_vals[i_opt]:.2e}, min error ≈ {errors[i_opt]:.2e}')
    # Expected: h* ≈ sqrt(eps_mach) ≈ 1.5e-8 for forward difference

    return h_vals, errors

h_vals, errors = step_size_analysis(f, df_true, 0.5)
```

> **Note:** The optimal step size for forward difference is approximately $h^* \approx \sqrt{\varepsilon_{mach}} \approx 1.5 \times 10^{-8}$, and for central difference it is approximately $h^* \approx \varepsilon_{mach}^{1/3} \approx 6 \times 10^{-6}$. Central difference has a larger optimal $h$ because its truncation error decreases faster ($O(h^2)$ vs. $O(h)$), so the balance point with round-off error occurs at a larger step size.

---

<br>

## Summary

| Topic | Key Formula / Concept |
|:------|:---------------------|
| Taylor Series | $f(x+h) = f(x) + f'(x)h + f''(x)h^2/2! + \cdots$ |
| Remainder | $R_n = f^{(n+1)}(\xi) h^{n+1} / (n+1)!$ |
| Forward Difference | $f'(x) = [f(x+h) - f(x)] / h + O(h)$ |
| Backward Difference | $f'(x) = [f(x) - f(x-h)] / h + O(h)$ |
| Central Difference | $f'(x) = [f(x+h) - f(x-h)] / (2h) + O(h^2)$ |
| 2nd Derivative | $f''(x) = [f(x+h) - 2f(x) + f(x-h)] / h^2 + O(h^2)$ |
| Higher-Order Forward | $f'(x) = [-f(x+2h) + 4f(x+h) - 3f(x)] / (2h) + O(h^2)$ |
| Optimal $h$ (Forward) | $h^* \sim \sqrt{\varepsilon_{mach}} \approx 1.5 \times 10^{-8}$ |
| Optimal $h$ (Central) | $h^* \sim \varepsilon_{mach}^{1/3} \approx 6 \times 10^{-6}$ |
| Condition Number | $\text{Cond} = \|x f'(x) / f(x)\|$; $\gg 1$ means ill-conditioned |

---
