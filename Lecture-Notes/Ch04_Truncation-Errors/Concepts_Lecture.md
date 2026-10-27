# Chapter 4 — Truncation Errors

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-27
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 4

> **Learning Objectives:**
> 1. Estimate truncation errors using Taylor series
> 2. Bound the truncation error with the remainder term
> 3. Analyze the order of accuracy of approximation methods and finite difference formulas
> 4. Explain the optimal step size as the balance between truncation and roundoff error
> 5. Decide from the condition number of a function whether a problem is ill-conditioned

---

## Table of Contents

- [1. Taylor Series Expansion](#1-taylor-series-expansion)
  - [1.1 The Taylor Series Formula](#11-the-taylor-series-formula)
  - [1.2 The Remainder Term (Truncation Error)](#12-the-remainder-term-truncation-error)
- [2. Order of Approximation](#2-order-of-approximation)
  - [2.1 Zeroth Through nth Order](#21-zeroth-through-nth-order)
  - [2.2 Big-O Notation and Its Meaning](#22-big-o-notation-and-its-meaning)
  - [2.3 Example: Taylor Series of $`e^x`$](#23-example-taylor-series-of-ex)
- [3. Numerical Differentiation Using Taylor Series](#3-numerical-differentiation-using-taylor-series)
  - [3.1 Forward Difference (First Derivative)](#31-forward-difference-first-derivative)
  - [3.2 Backward Difference (First Derivative)](#32-backward-difference-first-derivative)
  - [3.3 Central Difference (First Derivative)](#33-central-difference-first-derivative)
  - [3.4 Second Derivative (Central Difference)](#34-second-derivative-central-difference)
  - [3.5 Higher-Order Forward Difference](#35-higher-order-forward-difference)
- [4. Detailed Example: Numerical Differentiation](#4-detailed-example-numerical-differentiation)
  - [4.1 Problem Setup](#41-problem-setup)
  - [4.2 Results with $`h = 0.5`$](#42-results-with-h--05)
  - [4.3 Effect of Reducing $`h`$](#43-effect-of-reducing-h)
- [5. Total Numerical Error](#5-total-numerical-error)
  - [5.1 Truncation Error vs. Roundoff Error](#51-truncation-error-vs-roundoff-error)
  - [5.2 Optimal Step Size](#52-optimal-step-size)
- [6. Stability and Condition Number](#6-stability-and-condition-number)
- [7. Lab: Finite Differences in Python](#7-lab-finite-differences-in-python)
  - [7.1 Comparing Forward and Central Differences](#71-comparing-forward-and-central-differences)
  - [7.2 Step Size Analysis and Optimal $`h`$](#72-step-size-analysis-and-optimal-h)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Taylor Series Expansion

### 1.1 The Taylor Series Formula

The **Taylor series** is one of the most fundamental tools in numerical analysis. It allows us to represent a sufficiently smooth function as a sum of polynomial terms, each constructed from the function's derivatives evaluated at a single point, plus a remainder term; for an analytic function, the remainder tends to zero within the radius of convergence and the function equals the infinite sum.

Given a function $`f(x)`$ that is $`(n+1)`$ times continuously differentiable, the value of $`f`$ at a nearby point $`x_{i+1}`$ can be expressed in terms of $`f`$ and its derivatives at $`x_i`$ as follows.

```math
f(x_{i+1}) = f(x_i) + f'(x_i)(x_{i+1} - x_i) + \frac{f''(x_i)}{2!}(x_{i+1} - x_i)^2 + \frac{f'''(x_i)}{3!}(x_{i+1} - x_i)^3 + \cdots + R_n
```

Defining the **step size** $`h = x_{i+1} - x_i`$, this becomes the more compact form:

```math
f(x_i + h) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots + \frac{f^{(n)}(x_i)}{n!}h^n + R_n
```

Each successive term in the series includes a higher-order derivative of $`f`$ at $`x_i`$, divided by the corresponding factorial, and multiplied by a higher power of $`h`$. The more terms we include, the better the approximation, provided $`h`$ is small enough.

> **[Calculus]** Taylor's theorem is the calculus result that a function that is $`(n+1)`$ times continuously differentiable can be written exactly, near a point, as an $`n`$th-degree polynomial plus a remainder term. The truncation error analysis in this chapter starts by reading the remainder term of this theorem (Section 1.2) as the error of the polynomial approximation.

### 1.2 The Remainder Term (Truncation Error)

When we truncate the Taylor series after the $`n`$th-order term, the omitted terms are collectively represented by the **remainder** $`R_n`$:

```math
R_n = \frac{f^{(n+1)}(\xi)}{(n+1)!} h^{n+1}
```

where $`\xi`$ is some (unknown) value satisfying $`x_i \leq \xi \leq x_{i+1}`$.

Key properties of the remainder are as follows.

- $`R_n`$ is the **truncation error**, the error introduced by truncating an infinite series to a finite number of terms.
- The exact value of $`\xi`$ is generally unknown, which means we cannot compute $`R_n`$ exactly.
- However, we can **bound** $`R_n`$ by finding the maximum of $`|f^{(n+1)}(x)|`$ on the interval $`[x_i, x_{i+1}]`$:
  ```math
  |R_n| \leq \frac{\max_{x_i \leq x \leq x_{i+1}} |f^{(n+1)}(x)|}{(n+1)!} |h|^{n+1}
  ```
- For small $`h`$, the dominant behavior of $`R_n`$ is proportional to $`h^{n+1}`$, which we write as $`R_n = O(h^{n+1})`$.

---

<br>

## 2. Order of Approximation

### 2.1 Zeroth Through nth Order

Truncating the Taylor series at different points yields approximations of different orders. Each additional term raises the order of the truncation error by one, so the approximation improves when $`h`$ is small enough.

| Order | Approximation | Truncation Error |
|:------|:-------------|:-----------------|
| 0th | $`f(x_{i+1}) \approx f(x_i)`$ (constant) | $`O(h)`$ |
| 1st | $`f(x_{i+1}) \approx f(x_i) + f'(x_i)h`$ (linear) | $`O(h^2)`$ |
| 2nd | adds $`\frac{f''(x_i)}{2!}h^2`$ term (quadratic) | $`O(h^3)`$ |
| 3rd | adds $`\frac{f'''(x_i)}{3!}h^3`$ term (cubic) | $`O(h^4)`$ |
| nth | includes up to $`\frac{f^{(n)}(x_i)}{n!}h^n`$ term | $`O(h^{n+1})`$ |

*Table 1. Taylor approximations and truncation errors by truncation order*

### 2.2 Big-O Notation and Its Meaning

The **Big-O notation** $`O(h^{n+1})`$ describes how the truncation error scales with the step size $`h`$. Specifically, $`O(h^{n+1})`$ means that for sufficiently small $`h`$, the error is bounded by a constant times $`h^{n+1}`$:

```math
|\text{error}| \leq C \cdot h^{n+1} \quad \text{for some constant } C
```

The practical consequences are as follows.

- **0th-order** ($`O(h)`$): Halving $`h`$ halves the error.
- **1st-order** ($`O(h^2)`$): Halving $`h`$ reduces the error by a factor of 4.
- **2nd-order** ($`O(h^3)`$): Halving $`h`$ reduces the error by a factor of 8.
- **nth-order** ($`O(h^{n+1})`$): Halving $`h`$ reduces the error by a factor of $`2^{n+1}`$.

> **Note:** The Big-O notation $`O(h^{n+1})`$ means the error is proportional to $`h^{n+1}`$ for small $`h`$. The orders in the list above are the order $`n`$ of the Taylor approximation, whereas the order of accuracy of a numerical method is named by the exponent of its error. Halving $`h`$ in a first-order accurate method (error $`O(h)`$) halves the error; halving $`h`$ in a second-order accurate method (error $`O(h^2)`$) quarters the error.

### 2.3 Example: Taylor Series of $`e^x`$

Consider approximating $`e^1`$ using the Taylor series expansion of $`f(x) = e^x`$ about $`x_i = 0`$ (evaluating at $`x = 1`$, so $`h = 1`$).

Since $`f(x) = e^x`$, all derivatives are also $`e^x`$, so $`f(0) = f'(0) = f''(0) = \cdots = 1`$.

The Taylor series becomes:

```math
e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \frac{x^4}{4!} + \cdots
```

Evaluating at $`x = 1`$ ($`h = 1`$) and adding terms one at a time gives the following.

| Order | Approximation | Value | True Error | $`\varepsilon_t`$ |
|:------|:-------------|:------|:-----------|:----------------|
| 0th | $`f \approx 1`$ | 1.000000 | 1.71828 | 63.2% |
| 1st | $`f \approx 1 + 1`$ | 2.000000 | 0.71828 | 26.4% |
| 2nd | $`f \approx 2 + 0.5`$ | 2.500000 | 0.21828 | 8.0% |
| 3rd | $`f \approx 2.5 + 1/6`$ | 2.666667 | 0.05162 | 1.9% |
| 4th | $`f \approx 2.667 + 1/24`$ | 2.708333 | 0.00995 | 0.366% |
| 5th | $`f \approx 2.708 + 1/120`$ | 2.716667 | 0.00162 | 0.0594% |
| 6th | $`f \approx 2.717 + 1/720`$ | 2.718056 | 0.00023 | 0.00832% |

*Table 2. Taylor approximations of e and their errors as terms are added*

The true value is $`e = 2.718282\ldots`$. Each additional term significantly reduces the error, demonstrating the convergence of the Taylor series. By the 6th-order approximation, the relative error is already below $`0.01\%`$.

---

<br>

## 3. Numerical Differentiation Using Taylor Series

This section derives the finite difference formulas that are used throughout numerical analysis. These formulas approximate derivatives using discrete function values, which is a necessity when working with tabulated data or when analytical differentiation is impractical.

```text
x_{i-2}   x_{i-1}    x_i     x_{i+1}   x_{i+2}
   |         |        |         |         |
                      ●─────────●              forward difference      O(h)
             ●────────●                        backward difference     O(h)
             ●────────○─────────●              central difference      O(h^2)
             ●────────●─────────●              second derivative       O(h^2)
                      ●─────────●─────────●    higher-order forward    O(h^2)

● point used by the formula    ○ point not used by the formula
```

*Figure 1. Grid points used by each finite difference formula*

### 3.1 Forward Difference (First Derivative)

We start from the Taylor series expansion.

```math
f(x_{i+1}) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots
```

Solving for $`f'(x_i)`$ gives the following.

```math
f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \frac{f''(x_i)}{2!}h - \frac{f'''(x_i)}{3!}h^2 - \cdots
```

Truncating the higher-order terms gives the **first-order forward difference** formula.

```math
\boxed{f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} + O(h)}
```

The leading error term is $`-\frac{f''(x_i)}{2}h`$, so the truncation error is $`O(h)`$. This is a **first-order accurate** approximation.

### 3.2 Backward Difference (First Derivative)

Expand $`f(x_{i-1})`$ in a Taylor series about $`x_i`$ (substituting $`-h`$ for $`h`$).

```math
f(x_{i-1}) = f(x_i) - f'(x_i)h + \frac{f''(x_i)}{2!}h^2 - \frac{f'''(x_i)}{3!}h^3 + \cdots
```

Solving for $`f'(x_i)`$ gives the following.

```math
\boxed{f'(x_i) = \frac{f(x_i) - f(x_{i-1})}{h} + O(h)}
```

Like the forward difference, this is also **first-order accurate** with $`O(h)`$ truncation error. The backward difference uses the point behind $`x_i`$ rather than ahead.

### 3.3 Central Difference (First Derivative)

The central difference is derived by subtracting the backward Taylor expansion from the forward Taylor expansion.

**Forward expansion:**

```math
f(x_{i+1}) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots
```

**Backward expansion:**

```math
f(x_{i-1}) = f(x_i) - f'(x_i)h + \frac{f''(x_i)}{2!}h^2 - \frac{f'''(x_i)}{3!}h^3 + \cdots
```

**Subtracting** backward from forward gives the following.

```math
f(x_{i+1}) - f(x_{i-1}) = 2f'(x_i)h + \frac{2f'''(x_i)}{3!}h^3 + \cdots
```

Note that all even-powered terms, including $`f(x_i)`$ and the $`h^2, h^4, \ldots`$ terms, cancel. Solving for $`f'(x_i)`$ gives the following.

```math
\boxed{f'(x_i) = \frac{f(x_{i+1}) - f(x_{i-1})}{2h} + O(h^2)}
```

The truncation error is $`O(h^2)`$. This is **second-order accurate**, which is significantly more accurate than either the forward or backward difference.

> **Note:** **Accuracy of symmetric formulas.** The central difference is more accurate because the even-powered error terms cancel when we subtract the backward expansion from the forward expansion. This is a general principle: symmetric formulas tend to be more accurate.

### 3.4 Second Derivative (Central Difference)

The central difference formula for the second derivative is derived by **adding** the forward and backward Taylor expansions instead of subtracting them.

**Adding** forward and backward gives the following.

```math
f(x_{i+1}) + f(x_{i-1}) = 2f(x_i) + 2\frac{f''(x_i)}{2!}h^2 + 2\frac{f^{(4)}(x_i)}{4!}h^4 + \cdots
```

Note that all odd-powered terms ($`h, h^3, h^5, \ldots`$) cancel. Solving for $`f''(x_i)`$ gives the following.

```math
f''(x_i) = \frac{f(x_{i+1}) + f(x_{i-1}) - 2f(x_i)}{h^2} - \frac{f^{(4)}(x_i)}{12}h^2 - \cdots
```

Truncating gives the following.

```math
\boxed{f''(x_i) = \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2} + O(h^2)}
```

This is a **second-order accurate** approximation of the second derivative, requiring three function evaluations at equally spaced points.

### 3.5 Higher-Order Forward Difference

It is possible to achieve higher-order accuracy using only forward points. This is useful at the boundary of a domain where backward points are unavailable.

Using three forward points $`f(x_i)`$, $`f(x_{i+1})`$, and $`f(x_{i+2})`$, we can derive the **second-order forward difference** for $`f'(x_i)`$.

```math
\boxed{f'(x_i) = \frac{-f(x_{i+2}) + 4f(x_{i+1}) - 3f(x_i)}{2h} + O(h^2)}
```

This achieves $`O(h^2)`$ accuracy using only forward-direction points. It has the same order as the central difference but does not require points behind $`x_i`$.

---

<br>

## 4. Detailed Example: Numerical Differentiation

### 4.1 Problem Setup

Consider the following polynomial.

```math
f(x) = -0.1x^4 - 0.15x^3 - 0.5x^2 - 0.25x + 1.2
```

The analytically computed (true) derivative is:

```math
f'(x) = -0.4x^3 - 0.45x^2 - x - 0.25
```

At $`x = 0.5`$:

```math
f'(0.5) = -0.4(0.125) - 0.45(0.25) - 0.5 - 0.25 = -0.9125
```

We also need the function values:
- $`f(0.0) = 1.2`$
- $`f(0.25) = 1.10351563`$
- $`f(0.5) = 0.925`$
- $`f(0.75) = 0.63632813`$
- $`f(1.0) = 0.2`$

### 4.2 Results with $`h = 0.5`$

We use step size $`h = 0.5`$ to approximate $`f'(0.5)`$.

**Forward difference:**

```math
f'(0.5) \approx \frac{f(1.0) - f(0.5)}{0.5} = \frac{0.2 - 0.925}{0.5} = -1.45
```

```math
\varepsilon_t = \left|\frac{-0.9125 - (-1.45)}{-0.9125}\right| \times 100\% = 58.9\%
```

**Backward difference:**

```math
f'(0.5) \approx \frac{f(0.5) - f(0.0)}{0.5} = \frac{0.925 - 1.2}{0.5} = -0.55
```

```math
\varepsilon_t = \left|\frac{-0.9125 - (-0.55)}{-0.9125}\right| \times 100\% = 39.7\%
```

**Central difference:**

```math
f'(0.5) \approx \frac{f(1.0) - f(0.0)}{2 \times 0.5} = \frac{0.2 - 1.2}{1.0} = -1.0
```

```math
\varepsilon_t = \left|\frac{-0.9125 - (-1.0)}{-0.9125}\right| \times 100\% = 9.6\%
```

The central difference is dramatically more accurate than both forward and backward differences, consistent with its $`O(h^2)`$ error versus $`O(h)`$.

### 4.3 Effect of Reducing $`h`$

We use a smaller step size $`h = 0.25`$ with the central difference.

```math
f'(0.5) \approx \frac{f(0.75) - f(0.25)}{2 \times 0.25} = \frac{0.63632813 - 1.10351563}{0.5} = -0.934375
```

```math
\varepsilon_t = \left|\frac{-0.9125 - (-0.934375)}{-0.9125}\right| \times 100\% \approx 2.4\%
```

With the same $`h = 0.25`$, the forward difference gives $`(0.63632813 - 0.925)/0.25 = -1.1546875`$ ($`\varepsilon_t \approx 26.5\%`$), and the backward difference gives $`(0.925 - 1.10351563)/0.25 = -0.7140625`$ ($`\varepsilon_t \approx 21.7\%`$).

Halving $`h`$ from $`0.5`$ to $`0.25`$ reduced the central difference error from $`9.6\%`$ to $`2.4\%`$, roughly a factor of 4, confirming the $`O(h^2)`$ behavior. The forward and backward difference errors fell by roughly a factor of 2, consistent with $`O(h)`$ behavior.

| Method | $`h = 0.5`$ | $`h = 0.25`$ | Error Reduction |
|:-------|:----------|:-----------|:----------------|
| Forward | $`\varepsilon_t = 58.9\%`$ | $`\varepsilon_t = 26.5\%`$ | $`\sim 2.2\times`$ (confirms $`O(h)`$) |
| Backward | $`\varepsilon_t = 39.7\%`$ | $`\varepsilon_t = 21.7\%`$ | $`\sim 1.8\times`$ (confirms $`O(h)`$) |
| Central | $`\varepsilon_t = 9.6\%`$ | $`\varepsilon_t = 2.4\%`$ | $`\sim 4\times`$ (confirms $`O(h^2)`$) |

*Table 3. Error of the three difference formulas when the step size is halved*

---

<br>

## 5. Total Numerical Error

### 5.1 Truncation Error vs. Roundoff Error

When performing numerical differentiation on a computer, the total error has two components.

```math
\text{Total error} = \text{Truncation error} + \text{Roundoff error}
```

These two error sources have **opposite** dependencies on the step size $`h`$.

| Component | As $`h`$ decreases | Cause |
|:----------|:-----------------|:------|
| Truncation error | **Decreases** | Finite difference better approximates the derivative |
| Roundoff error | **Increases** | Subtracting nearly equal floating-point numbers (subtractive cancellation, Chapter 3) |

*Table 4. Change in truncation and roundoff error as the step size decreases*

When $`h`$ is large, the truncation error dominates. When $`h`$ is very small, the function values $`f(x+h)`$ and $`f(x)`$ become nearly identical, and their difference is contaminated by floating-point roundoff. The total error curve has a characteristic **V-shape** on a log-log plot.

> **Note:** This is one of the most important practical insights in numerical analysis. We cannot make $`h`$ arbitrarily small to improve accuracy. At some point, roundoff error takes over and the result actually gets worse.

### 5.2 Optimal Step Size

There exists an **optimal step size** $`h^*`$ that minimizes the total error. The optimal value depends on the finite difference formula being used and the machine epsilon $`\varepsilon_{mach}`$ (Section 5.4 of Chapter 3).

**For forward difference** ($`O(h)`$ truncation error):

```math
h^* \approx \sqrt{\varepsilon_{mach} \cdot \left|\frac{f(x)}{f''(x)}\right|}
```

Since $`\varepsilon_{mach} \approx 2.2 \times 10^{-16}`$ for 64-bit double precision:

```math
h^* \sim \sqrt{\varepsilon_{mach}} \approx 1.5 \times 10^{-8}
```

**For central difference** ($`O(h^2)`$ truncation error):

```math
h^* \sim \varepsilon_{mach}^{1/3} \approx 6 \times 10^{-6}
```

The central difference has a larger optimal $`h`$ because its truncation error decreases faster ($`O(h^2)`$ vs. $`O(h)`$), so the balance point between truncation and roundoff shifts to a larger step size.

| Formula | Truncation Error | Optimal $`h^*`$ | Minimum Total Error |
|:--------|:----------------|:---------------|:--------------------|
| Forward difference | $`O(h)`$ | $`\sim \sqrt{\varepsilon_{mach}} \approx 10^{-8}`$ | $`\sim \sqrt{\varepsilon_{mach}} \approx 10^{-8}`$ |
| Central difference | $`O(h^2)`$ | $`\sim \varepsilon_{mach}^{1/3} \approx 10^{-6}`$ | $`\sim \varepsilon_{mach}^{2/3} \approx 10^{-11}`$ |

*Table 5. Optimal step size and minimum total error of the forward and central differences*

---

<br>

## 6. Stability and Condition Number

The **condition number** of a function quantifies how sensitive the output is to small perturbations in the input. For a function $`f(x)`$, the condition number is:

```math
\text{Cond} = \left|\frac{x \cdot f'(x)}{f(x)}\right|
```

The condition number measures the ratio of the relative change in the output to the relative change in the input.

- If $`\text{Cond} \gg 1`$, the problem is **ill-conditioned**. Small input errors produce large output errors, and numerical results may be unreliable regardless of the algorithm used.
- If $`\text{Cond}`$ is about 1 or smaller, the problem is **well-conditioned**. Input errors are not amplified significantly.

For example, consider $`f(x) = \frac{1}{x - 1}`$ near $`x = 1.00001`$:

```math
\text{Cond} = \left|\frac{x \cdot (-1/(x-1)^2)}{1/(x-1)}\right| = \left|\frac{x}{x-1}\right|
```

At $`x = 1.00001`$, $`\text{Cond} \approx 10^5`$, so the problem is extremely ill-conditioned. A tiny perturbation in $`x`$ causes a massive relative change in $`f(x)`$.

> **Note:** The matrix condition number $`\kappa(A) = \|A\| \cdot \|A^{-1}\|`$, which extends the same sensitivity concept to matrices, is covered in Chapter 11. Just as the condition number of a function measures the sensitivity of the output to input perturbations, the matrix condition number measures the sensitivity of the solution of a linear system to perturbations in the coefficients and the right-hand side.

---

<br>

## 7. Lab: Finite Differences in Python

### 7.1 Comparing Forward and Central Differences

The following code compares the error of forward and central difference formulas across a range of step sizes, demonstrating both the accuracy advantage of central differences and the effect of roundoff error at very small step sizes.

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

The resulting log-log plot shows two key features.

1. **For large $`h`$:** Both errors decrease as $`h`$ decreases, with the central difference error decreasing faster (the slopes on the plot are about 2 versus 1). For example, reducing $`h`$ from $`10^{-1}`$ to $`10^{-4}`$ reduces the forward difference error from $`9.1 \times 10^{-2}`$ to $`8.8 \times 10^{-5}`$ and the central difference error from $`3.5 \times 10^{-3}`$ to $`3.5 \times 10^{-9}`$.
2. **For very small $`h`$:** Both errors increase again due to roundoff, forming the characteristic V-shape. The minimum error for the central difference (about $`5.4 \times 10^{-12}`$ at $`h = 10^{-6}`$) is much lower than that of the forward difference (about $`8.9 \times 10^{-9}`$ at $`h = 10^{-8}`$) and occurs at a larger $`h`$.

### 7.2 Step Size Analysis and Optimal $`h`$

The following code systematically finds the optimal step size for the forward difference.

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

Running this code prints `Optimal h ≈ 1.00e-08, min error ≈ 8.89e-09`. Among the grid values $`10^{-1}, 10^{-2}, \ldots, 10^{-15}`$, the value $`10^{-8}`$, which is closest to the expected $`\sqrt{\varepsilon_{mach}} \approx 1.5 \times 10^{-8}`$, is selected as optimal.

> **Note:** As analyzed in Section 5.2, the optimal step size is about $`h^* \approx \sqrt{\varepsilon_{mach}} \approx 1.5 \times 10^{-8}`$ for the forward difference and about $`h^* \approx \varepsilon_{mach}^{1/3} \approx 6 \times 10^{-6}`$ for the central difference. In the experiment of Section 7.1, the central difference error is also smallest at $`h = 10^{-6}`$, in agreement with this prediction.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Taylor Series | $`f(x+h) = f(x) + f'(x)h + f''(x)h^2/2! + \cdots`$ expresses the function value at a nearby point through the derivatives at one point. |
| Remainder | $`R_n = f^{(n+1)}(\xi) h^{n+1} / (n+1)!`$ is the truncation error of truncating at order $`n`$ and is $`O(h^{n+1})`$. |
| Forward Difference | $`f'(x) = [f(x+h) - f(x)] / h + O(h)`$ is first-order accurate. |
| Backward Difference | $`f'(x) = [f(x) - f(x-h)] / h + O(h)`$ is first-order accurate. |
| Central Difference | $`f'(x) = [f(x+h) - f(x-h)] / (2h) + O(h^2)`$ is second-order accurate. |
| Central Difference for the 2nd Derivative | $`f''(x) = [f(x+h) - 2f(x) + f(x-h)] / h^2 + O(h^2)`$ is second-order accurate. |
| Higher-Order Forward Difference | $`f'(x) = [-f(x+2h) + 4f(x+h) - 3f(x)] / (2h) + O(h^2)`$ reaches second-order accuracy without points behind $`x`$. |
| Optimal Step Size (Forward) | The balance point between truncation and roundoff error is $`h^* \sim \sqrt{\varepsilon_{mach}} \approx 1.5 \times 10^{-8}`$. |
| Optimal Step Size (Central) | The truncation error shrinks faster, so the balance point grows to $`h^* \sim \varepsilon_{mach}^{1/3} \approx 6 \times 10^{-6}`$. |
| Condition Number | It is $`\text{Cond} = \vert x f'(x) / f(x) \vert`$, and the problem is ill-conditioned when $`\text{Cond} \gg 1`$. |

---

<br>

## Self-Check Questions

1. **Bounding the remainder:** Approximating $`e^{0.1}`$ by expanding $`f(x) = e^x`$ about $`x_i = 0`$ up to the second-order term, find the approximation and an upper bound for $`R_2`$.

   <details>
   <summary>Answer</summary>

   The approximation is $`1 + 0.1 + 0.1^2/2 = 1.105`$. Since $`R_2 = e^{\xi} \cdot 0.1^3 / 3!`$ and $`e^{\xi} \leq e^{0.1} \approx 1.10517`$ for $`0 \leq \xi \leq 0.1`$, we have $`|R_2| \leq 1.10517 \times 0.001 / 6 \approx 1.84 \times 10^{-4}`$. The actual error $`e^{0.1} - 1.105 \approx 1.71 \times 10^{-4}`$ lies within this bound.

   </details>

2. **Applying difference formulas:** For $`f(x) = x^3`$, $`x = 2`$, and $`h = 0.1`$, approximate $`f'(2)`$ with the forward, backward, and central differences and find each error.

   <details>
   <summary>Answer</summary>

   The true value is $`f'(2) = 12`$. The forward difference gives $`(2.1^3 - 2^3)/0.1 = (9.261 - 8)/0.1 = 12.61`$, an error of $`0.61`$. The backward difference gives $`(8 - 1.9^3)/0.1 = (8 - 6.859)/0.1 = 11.41`$, an error of $`0.59`$. The central difference gives $`(9.261 - 6.859)/0.2 = 12.01`$, an error of $`0.01`$, which matches the leading error term $`f'''(x)h^2/6 = 6 \times 0.01 / 6 = 0.01`$.

   </details>

3. **Using the order of accuracy:** The central difference error was $`4.0 \times 10^{-3}`$ at $`h = 0.2`$. Ignoring roundoff, what error do you expect at $`h = 0.05`$?

   <details>
   <summary>Answer</summary>

   The central difference error is $`O(h^2)`$, so reducing $`h`$ to a quarter reduces the error to about $`(1/4)^2 = 1/16`$. The expected error is therefore $`4.0 \times 10^{-3} / 16 = 2.5 \times 10^{-4}`$.

   </details>

4. **Optimal step size:** Using an error model, explain why the optimal $`h`$ is about $`\sqrt{\varepsilon_{mach}}`$ for the forward difference and about $`\varepsilon_{mach}^{1/3}`$ for the central difference.

   <details>
   <summary>Answer</summary>

   The roundoff error divides a difference of two function values by $`h`$, so it is proportional to about $`\varepsilon_{mach}/h`$. The total error of the forward difference has the form $`C_1 h + C_2 \varepsilon_{mach}/h`$, which is smallest where the two terms balance, at $`h \sim \sqrt{\varepsilon_{mach}}`$. The central difference has the form $`C_1 h^2 + C_2 \varepsilon_{mach}/h`$, which is smallest at $`h^3 \sim \varepsilon_{mach}`$, that is, $`h \sim \varepsilon_{mach}^{1/3}`$, where the minimum error is of order $`\varepsilon_{mach}^{2/3}`$.

   </details>

5. **Condition number:** Find the condition number of $`f(x) = \sqrt{x}`$ for $`x > 0`$ and decide whether the problem is ill-conditioned.

   <details>
   <summary>Answer</summary>

   Since $`f'(x) = 1/(2\sqrt{x})`$, $`\text{Cond} = |x \cdot (1/(2\sqrt{x})) / \sqrt{x}| = 1/2`$. For every $`x > 0`$, the relative input error reaches the output halved, so the problem is well-conditioned.

   </details>
