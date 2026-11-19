# Chapter 6 — Roots: Open Methods

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-11-19
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 6

> **Learning Objectives:**
> 1. Apply the Newton-Raphson method to root finding
> 2. Implement the secant and modified secant methods
> 3. Compare the convergence of open and bracketing methods

---

## Table of Contents

- [1. Objectives](#1-objectives)
- [2. Open Methods vs Bracketing Methods](#2-open-methods-vs-bracketing-methods)
- [3. Fixed-Point Iteration](#3-fixed-point-iteration)
  - [3.1 Concept and Definition](#31-concept-and-definition)
  - [3.2 Algorithm](#32-algorithm)
  - [3.3 Fixed Point and Root Relationship](#33-fixed-point-and-root-relationship)
  - [3.4 Convergence Analysis](#34-convergence-analysis)
  - [3.5 Convergence Proof via Mean Value Theorem](#35-convergence-proof-via-mean-value-theorem)
  - [3.6 Example 6.1: Convergent Case](#36-example-61-convergent-case)
  - [3.7 Example 6.1: Divergent Case](#37-example-61-divergent-case)
  - [3.8 Example 6.1: Cycling Case](#38-example-61-cycling-case)
  - [3.9 Graphical Interpretation of Fixed-Point Iteration](#39-graphical-interpretation-of-fixed-point-iteration)
  - [3.10 Rate of Convergence and Error Estimation](#310-rate-of-convergence-and-error-estimation)
- [4. The Wegstein Method](#4-the-wegstein-method)
  - [4.1 Motivation from False Position](#41-motivation-from-false-position)
  - [4.2 Geometric Idea](#42-geometric-idea)
  - [4.3 Derivation](#43-derivation)
  - [4.4 General Formula](#44-general-formula)
- [5. Newton-Raphson Method](#5-newton-raphson-method)
  - [5.1 Geometric Derivation](#51-geometric-derivation)
  - [5.2 Two-Step Approach](#52-two-step-approach)
  - [5.3 Derivation from Taylor Expansion](#53-derivation-from-taylor-expansion)
  - [5.4 Quadratic Convergence](#54-quadratic-convergence)
  - [5.5 Proof of Quadratic Convergence](#55-proof-of-quadratic-convergence)
  - [5.6 Example 6.5: Slowly Converging Function](#56-example-65-slowly-converging-function)
  - [5.7 Numerical Confirmation of Quadratic Convergence](#57-numerical-confirmation-of-quadratic-convergence)
  - [5.8 When the Newton-Raphson Method Fails](#58-when-the-newton-raphson-method-fails)
  - [5.9 Slow Convergence at Multiple Roots](#59-slow-convergence-at-multiple-roots)
  - [5.10 Newton-Raphson Method for Systems of Nonlinear Equations](#510-newton-raphson-method-for-systems-of-nonlinear-equations)
- [6. Secant Method](#6-secant-method)
  - [6.1 Standard Secant Method](#61-standard-secant-method)
  - [6.2 Modified Secant Method](#62-modified-secant-method)
  - [6.3 Order of Convergence of the Secant Method](#63-order-of-convergence-of-the-secant-method)
  - [6.4 Comparing the Convergence of Open and Bracketing Methods](#64-comparing-the-convergence-of-open-and-bracketing-methods)
- [7. Inverse Quadratic Interpolation](#7-inverse-quadratic-interpolation)
  - [7.1 Concept](#71-concept)
  - [7.2 Lagrange Polynomial Formulation](#72-lagrange-polynomial-formulation)
  - [7.3 Comparison with Secant Method](#73-comparison-with-secant-method)
- [8. Brent's Method](#8-brents-method)
  - [8.1 Overview](#81-overview)
  - [8.2 Python Implementation with SciPy](#82-python-implementation-with-scipy)
  - [8.3 Steps of Brent's Method](#83-steps-of-brents-method)
  - [8.4 Comparing SciPy Root-Finding Functions](#84-comparing-scipy-root-finding-functions)
- [9. Summary Table](#9-summary-table)
- [10. Roots of Polynomials](#10-roots-of-polynomials)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)
- [References](#references)

---

<br>

## 1. Objectives

- Understand the **fixed-point iteration**.
- Implement the **Wegstein method**.
- Understand the quadratically convergent **Newton-Raphson method**.
- Understand the **secant method**.
- Understand **Brent's method**.

> **Note:** **Connection to Chapter 5.** This chapter builds on Chapter 5 (Roots: Bracketing Methods), which covered bisection and false position. The concepts of bracketing, error tolerance, and iterative convergence from Chapter 5 are directly extended here.

---

<br>

## 2. Open Methods vs Bracketing Methods

Unlike **bracketing methods** (bisection, false position), open methods do **not** require two initial guesses that bracket a root. Their key characteristics are as follows.

| Property | Bracketing Methods | Open Methods |
|---|---|---|
| Initial guesses | Two points $`x_l, x_u`$ with $`f(x_l) \cdot f(x_u) < 0`$ | One or two initial guesses (no sign-change requirement) |
| Convergence guarantee | Always converge (guaranteed) | May **diverge** or move away from the root |
| Speed | Slower (linear convergence) | Faster convergence **when they converge** |
| Interval behavior | Interval shrinks each iteration | Single point(s) iterate toward root |

*Table 1. Comparison of bracketing and open methods*

> **Key Point:** Open methods trade convergence guarantees for speed. If they converge, they approach the root much more quickly than bracketing methods.

---

<br>

## 3. Fixed-Point Iteration

### 3.1 Concept and Definition

It is also known by the following names.
- **Simple one-point iteration**
- **Method of successive substitutions**

Given an equation $`f(x) = 0`$ where $`f: [a, b] \to \mathbb{R}`$, we rearrange it into the following form.

```math
f(x) = x - g(x) = 0 \quad \Longrightarrow \quad x = g(x)
```

The iteration $`x_{i+1} = g(x_i)`$ is called **fixed-point iteration**.

> **Definition:** A **fixed point** of a function $`g`$ is a real number $`x^*`$ such that:
>
> ```math
> g(x^*) = x^*
> ```

> **Note:** **Difference from bracketing methods.** While bracketing methods (Chapter 5) work directly with $`f(x) = 0`$, fixed-point iteration reformulates the problem as $`x = g(x)`$, which can be viewed as finding the intersection of $`y = g(x)`$ and $`y = x`$.

### 3.2 Algorithm

To find the root $`x^* \in (a, b)`$, proceed as follows.

1. Choose a starting point $`x_0 \in (a, b)`$.
2. Sequentially calculate $`x_1, x_2, \ldots, x_{i+1}`$ using the following formula.

```math
x_{i+1} = g(x_i)
```

3. Repeat until $`x_{i+1} \approx x_i`$, that is, until convergence is achieved.

This is a **circular calculation**: the output $`x_{i+1}`$ from one step becomes the input $`x_i`$ for the next.

```mermaid
flowchart LR
    I["input x_i"] --> G["g"]
    G --> O["output x_{i+1}"]
    O --> I
```

*Figure 1. The circular calculation of fixed-point iteration*

### 3.3 Fixed Point and Root Relationship

Let $`g`$ be a continuous function. Let the sequence $`\{x_i : i \geq 0\}`$ be generated by the fixed-point iteration $`x_{i+1} = g(x_i)`$.

If $`\displaystyle\lim_{i \to \infty} x_i = x^*`$, then $`x^*`$ is a fixed point of $`g`$.

Thus $`x^*`$ is the **root** of $`f`$:

```math
f(x^*) = x^* - g(x^*) = 0
```

### 3.4 Convergence Analysis

Convergence depends on the **slope of $`g`$** at the fixed point. Using the mean value theorem gives the following.

```math
\text{err}_{i+1} = \text{err}_i \cdot |g'(\xi)|
```

where $`\xi`$ is some point between $`x_i`$ and $`x^*`$.

**Convergence criterion:**

| Condition | Result |
|---|---|
| $`\vert g'(\xi)\vert < 1`$ | The iteration **converges** |
| $`\vert g'(\xi)\vert > 1`$ | The iteration **diverges** |
| $`\vert g'(\xi)\vert = 1`$ | Inconclusive (may cycle) |

*Table 2. Convergence of fixed-point iteration by the magnitude of the slope of g*

### 3.5 Convergence Proof via Mean Value Theorem

Let $`x^* \in (a, b)`$ be a fixed point of $`g`$ and let $`x_i \in (a, b)`$.

By the **Mean Value Theorem**, there exists $`\xi_i \in (x_i, x^*)`$ or $`(x^*, x_i)`$ such that:

```math
g(x_i) - g(x^*) = g'(\xi_i)(x_i - x^*)
```

> **[Calculus]** The Mean Value Theorem states that for a function $`g`$ continuous on $`[a, b]`$ and differentiable on $`(a, b)`$ there exists $`\xi \in (a, b)`$ with $`g(b) - g(a) = g'(\xi)(b - a)`$. The convergence analysis of fixed-point iteration uses it to turn the difference of function values at two points into the slope $`g'(\xi_i)`$ times the distance, showing that the error is multiplied by $`|g'(\xi_i)|`$ at each iteration. The analysis therefore holds under the assumption that $`g`$ is differentiable.

Since $`x_i = g(x_{i-1})`$ and $`x^* = g(x^*)`$, we obtain the following.

```math
|x_i - x^*| = |g(x_{i-1}) - g(x^*)| = |g'(\xi_{i-1})| \cdot |x_{i-1} - x^*|
```

If $`|g'(\xi)| \leq q < 1`$ for every $`\xi \in (a, b)`$ (assuming that all iterates stay in $`(a, b)`$), then:

```math
|x_i - x^*| \leq q \cdot |x_{i-1} - x^*| \leq q^2 \cdot |x_{i-2} - x^*| \leq \cdots \leq q^i \cdot |x_0 - x^*|
```

Since $`q < 1`$, we have $`q^i \to 0`$ as $`i \to \infty`$, so:

```math
|x_i - x^*| \to 0
```

> **Key Point:** The fixed-point iteration converges if the **contraction constant** $`q = \max|g'(\xi)| < 1`$ on the interval containing the root and iterates.

### 3.6 Example 6.1: Convergent Case

**Problem:** Estimate the root of $`f(x) = x - e^{-x}`$.

**Setup:** From $`f(x) = 0`$, we set the following.

```math
x = e^{-x} = g(x)
```

We start with $`x_0 = 0`$.

| Iteration | $`x_i`$ | $`g(x_i) = e^{-x_i}`$ |
|---|---|---|
| 0 | 0 | $`e^{0} = 1.0000`$ |
| 1 | 1.0000 | $`e^{-1} = 0.3679`$ |
| 2 | 0.3679 | $`e^{-0.3679} = 0.6922`$ |
| 3 | 0.6922 | $`\cdots`$ |
| $`\vdots`$ | $`\vdots`$ | converges to $`x^* \approx 0.5671`$ |

*Table 3. The first steps of a convergent fixed-point iteration*

**Why it converges:** $`g(x) = e^{-x}`$, so $`g'(x) = -e^{-x}`$. At the root $`x^* \approx 0.5671`$, $`|g'(x^*)| = e^{-0.5671} \approx 0.5671 < 1`$.

Graphically, this is a **spiral convergence**. Since $`g'(x^*) < 0`$, the iterates alternate between the two sides of the fixed point and spiral inward toward the intersection of $`y = g(x) = e^{-x}`$ and $`y = x`$.

```python
import numpy as np

def fixed_point_ex1(x0, tol=1e-6, max_iter=100):
    """Fixed-point iteration for f(x) = x - e^(-x)"""
    g = lambda x: np.exp(-x)
    x = x0
    for i in range(max_iter):
        x_new = g(x)
        if abs(x_new - x) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.6f}")
            return x_new
        x = x_new
    print("Did not converge")
    return x
```

Running `fixed_point_ex1(0)` converges to `x* = 0.567143` after 26 iterations.

### 3.7 Example 6.1: Divergent Case

**Problem:** Estimate the root of $`f(x) = x + \ln x`$.

**Setup:** From $`f(x) = 0`$, we set the following.

```math
x = -\ln(x) = g(x)
```

**Result:** The iteration **diverges** (spiral outward). For example, starting from $`x_0 = 0.5`$ gives $`0.693, 0.367, 1.004, -0.0037`$, whose oscillation grows until the value becomes negative and the logarithm can no longer be computed.

**Why it diverges:** $`g(x) = -\ln(x)`$, so $`g'(x) = -1/x`$. At the root $`x^* \approx 0.5671`$, $`|g'(x^*)| = 1/0.5671 \approx 1.763 > 1`$.

> **Note:** **Choice of rearrangement.** The same root can be found by one rearrangement of $`f(x) = 0`$ but not another. The choice of $`g(x)`$ critically determines convergence.

### 3.8 Example 6.1: Cycling Case

Consider $`g(x) = \dfrac{1}{x}`$ with $`x_0 = 4`$.

| Iteration | $`x_i`$ | $`g(x_i) = 1/x_i`$ |
|---|---|---|
| 0 | 4 | $`1/4 = 0.25`$ |
| 1 | 0.25 | $`1/0.25 = 4`$ |
| 2 | 4 | $`1/4 = 0.25`$ |

*Table 4. A fixed-point iteration cycling between two points*

The iteration **cycles between two points** $`x_0`$ and $`x_1`$ indefinitely, never converging to the fixed point $`x^* = 1`$.

**Why:** $`g'(x) = -1/x^2`$, so $`|g'(1)| = 1`$ (boundary case).

### 3.9 Graphical Interpretation of Fixed-Point Iteration

Fixed-point iteration can be pictured as a staircase or cobweb path on the graphs of $`y = g(x)`$ and $`y = x`$. Going vertically from $`x_i`$ to the curve $`y = g(x)`$ gives the height $`x_{i+1} = g(x_i)`$, and going horizontally to the line $`y = x`$ makes the horizontal coordinate $`x_{i+1}`$. If the path traced by repeating this approaches the intersection of the two graphs (the fixed point), the iteration converges; if it moves away, it diverges.

![Figure 2. Paths of a converging and a diverging fixed-point iteration](../../images/Ch06_Fixed-Point-Cobweb.png)

*Figure 2. Paths of a converging and a diverging fixed-point iteration*

The left panel of Figure 2 is $`g(x) = e^{-x}`$ from Section 3.6, and the right panel is $`g(x) = -\ln x`$ from Section 3.7. Both curves meet $`y = x`$ at the same fixed point $`x^* \approx 0.5671`$, but on the left the magnitude of the slope near the intersection is $`0.567 < 1`$, so the path winds inward, while on the right it is $`1.763 > 1`$, so the path unwinds outward. The shape of the path is set by the sign of $`g'(x^*)`$.

| Range of $`g'(x^*)`$ | Shape of the Path | Result |
|:-----|:-----|:-----|
| $`0 < g'(x^*) < 1`$ | Staircase approaching from one side | Monotone convergence |
| $`-1 < g'(x^*) < 0`$ | Cobweb approaching from alternating sides | Oscillatory convergence |
| $`g'(x^*) > 1`$ | Staircase moving away to one side | Monotone divergence |
| $`g'(x^*) < -1`$ | Cobweb moving away on alternating sides | Oscillatory divergence |

*Table 5. Shape of the fixed-point iteration path by the slope of g*

### 3.10 Rate of Convergence and Error Estimation

From the result of Section 3.5, $`|x_{i+1} - x^*| = |g'(\xi_i)| \cdot |x_i - x^*|`$, if $`x_i \to x^*`$ then $`\xi_i \to x^*`$, so the ratio of successive errors converges to $`|g'(x^*)|`$. A fixed-point iteration with $`0 < |g'(x^*)| < 1`$ therefore converges linearly with asymptotic error constant $`C = |g'(x^*)|`$ (Section 7.4 of Chapter 5). Applying $`g(x) = e^{-x}`$ of Section 3.6 with $`x_0 = 0`$ gives the following table.

| $`i`$ | $`x_i`$ | $`\vert x_i - x^* \vert`$ | Error Ratio | $`\varepsilon_a`$ (%) |
|-----:|-----:|-----:|-----:|-----:|
| 1 | 1.000000 | $`4.329 \times 10^{-1}`$ | - | 100 |
| 2 | 0.367879 | $`1.993 \times 10^{-1}`$ | 0.4603 | 172 |
| 3 | 0.692201 | $`1.251 \times 10^{-1}`$ | 0.6276 | 46.9 |
| 4 | 0.500474 | $`6.667 \times 10^{-2}`$ | 0.5331 | 38.3 |
| 5 | 0.606244 | $`3.910 \times 10^{-2}`$ | 0.5865 | 17.4 |
| 6 | 0.545396 | $`2.175 \times 10^{-2}`$ | 0.5562 | 11.2 |
| 8 | 0.560115 | $`7.028 \times 10^{-3}`$ | 0.5636 | 3.48 |
| 10 | 0.564879 | $`2.264 \times 10^{-3}`$ | 0.5660 | 1.11 |
| 12 | 0.566415 | $`7.286 \times 10^{-4}`$ | 0.5668 | 0.356 |

*Table 6. Fixed-point iteration of g(x) = exp(-x) and the ratio of errors*

The error ratio approaches $`|g'(x^*)| = e^{-x^*} = x^* \approx 0.5671`$, and $`x_i`$ alternates between the two sides of the fixed point because $`g'(x^*) < 0`$. The error is multiplied by about 0.567 in each iteration, so each significant digit costs about $`1/\log_{10}(1/0.567) \approx 4`$ iterations. This is why the code of Section 3.6 took 26 iterations for a tolerance of $`10^{-6}`$.

Even when the true value is unknown, the error can be estimated from the difference of two successive approximations. With contraction constant $`q`$, $`|x_{i+1} - x^*| \leq q|x_i - x^*|`$, and by the triangle inequality $`|x_i - x^*| \leq |x_i - x_{i+1}| + |x_{i+1} - x^*|`$; rearranging gives:

```math
|x_{i+1} - x^*| \leq \frac{q}{1 - q}\,|x_{i+1} - x_i|
```

With $`q = 0.567`$, $`q/(1 - q) \approx 1.31`$, so if the difference of two successive approximations is below the tolerance, the true error is at most 1.31 times that. If $`q`$ is close to 1, however, for example $`q = 0.99`$, the factor is 99, and the true error can be up to 99 times larger even when the approximation barely changes. Stopping a slowly converging iteration on $`|x_{i+1} - x_i|`$ alone overestimates its accuracy.

---

<br>

## 4. The Wegstein Method

### 4.1 Motivation from False Position

Recall that in the **false position method** (Chapter 5), a straight line between two guesses on $`f(x)`$ crossed the x-axis to produce the next estimate.

The Wegstein method applies a similar linear interpolation idea, but to the function $`g(x)`$ instead of $`f(x)`$.

### 4.2 Geometric Idea

In the Wegstein method, a straight line between two guesses on the curve $`y = g(x)`$ **intersects with the line $`y = x`$** to determine the next estimate of the solution.

Given two points $`(x_0, g(x_0))`$ and $`(x_1, g(x_1))`$ on the curve $`g(x)`$, proceed as follows.
- Draw a straight line through these two points.
- Find where this line intersects $`y = x`$.
- That intersection gives the next estimate $`x_2`$.

### 4.3 Derivation

The slope of the line through $`(x_0, g(x_0))`$ and $`(x_1, g(x_1))`$ equals the slope through the two points $`(x_1, g(x_1))`$ and $`(x_2, g(x_2))`$ on the same line. Here $`g(x_2)`$ denotes not the value of the curve $`g`$ but the height of this line at $`x_2`$.

```math
\text{slope} = \frac{g(x_1) - g(x_0)}{x_1 - x_0} = \frac{g(x_2) - g(x_1)}{x_2 - x_1}
```

Since the intersection is on the line $`y = x`$, this height is $`g(x_2) = x_2`$. Therefore:

```math
\frac{g(x_1) - g(x_0)}{x_1 - x_0} = \frac{x_2 - g(x_1)}{x_2 - x_1}
```

Cross-multiplying and solving for $`x_2`$ gives the following.

```math
x_2(g(x_1) - g(x_0)) - x_1(g(x_1) - g(x_0)) = (x_1 - x_0)x_2 - g(x_1)(x_1 - x_0)
```

```math
x_2 = \frac{x_1 g(x_0) - x_0 g(x_1)}{x_1 - x_0 - g(x_1) + g(x_0)}
```

### 4.4 General Formula

The general Wegstein iteration formula is:

```math
\boxed{x_{i+1} = \frac{x_i \, g(x_{i-1}) - x_{i-1} \, g(x_i)}{x_i - x_{i-1} - g(x_i) + g(x_{i-1})}}
```

> **Key Point:** The last two estimates are used to determine the next estimate. This requires **two initial guesses** (like the secant method), but they do not need to bracket the root.

> **Note:** **Relation to false position and the secant method.** The Wegstein method is essentially the false position idea (Chapter 5) applied to $`g(x)`$ intersecting $`y = x`$ rather than $`f(x)`$ intersecting $`y = 0`$. Since it keeps no bracket and always uses the last two estimates, however, its formula reduces to the secant method (Section 6.1) applied to $`F(x) = x - g(x)`$.

```python
def wegstein(g, x0, x1, tol=1e-6, max_iter=100):
    """Wegstein method for solving x = g(x)"""
    for i in range(max_iter):
        g0, g1 = g(x0), g(x1)
        denom = x1 - x0 - g1 + g0
        if abs(denom) < 1e-15:
            print("Denominator near zero")
            return x1
        x_new = (x1 * g0 - x0 * g1) / denom
        if abs(x_new - x1) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.6f}")
            return x_new
        x0, x1 = x1, x_new
    print("Did not converge")
    return x1
```

For example, applying it to $`g(x) = e^{-x}`$ with $`x_0 = 0`$ and $`x_1 = 1`$ converges to `x* = 0.567143` in 5 iterations, and applying it with $`x_0 = 0.5`$ and $`x_1 = 0.6`$ to $`g(x) = -\ln x`$, for which fixed-point iteration diverged, converges to the same root in 4 iterations.

---

<br>

## 5. Newton-Raphson Method

### 5.1 Geometric Derivation

Given $`f(x) = 0`$, we use the **first-order derivative** (tangent line) at the current estimate.

At the point $`(x_i, f(x_i))`$, the slope of the tangent line is:

```math
\text{slope} = f'(x_i) = \frac{0 - f(x_i)}{x_{i+1} - x_i}
```

Solving for $`x_{i+1}`$ gives the following.

```math
\boxed{x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}}
```

### 5.2 Two-Step Approach

The Newton-Raphson method can be viewed as a two-step process.

**Step 1:** Compute the correction (Newton step).

```math
\delta x_i = -\frac{f(x_i)}{f'(x_i)}
```

**Step 2:** Update the estimate.

```math
x_{i+1} = x_i + \delta x_i
```

> **Note:** **Comparison with bracketing methods.** Newton-Raphson is an open method requiring only **one** initial guess and the derivative $`f'(x)`$. Unlike bracketing methods, there is no guarantee of convergence, but when it converges to a root with $`f'(x^*) \neq 0`$, it does so quadratically.

### 5.3 Derivation from Taylor Expansion

We can derive the Newton-Raphson formula from the **Taylor expansion**.

```math
f(x_{i+1}) = f(x_i + \delta x_i) = f(x_i) + f'(x_i)\,\delta x_i + \frac{f''(x_i)}{2}(\delta x_i)^2 + \cdots
```

Truncating after the first-order term gives the following.

```math
f(x_{i+1}) \approx f(x_i) + f'(x_i)\,\delta x_i
```

We **hope** $`f(x_{i+1}) = 0`$ (i.e., $`x_{i+1}`$ is the root $`x^*`$), so we set the above to zero.

```math
0 = f(x_i) + f'(x_i)\,\delta x_i
```

```math
\delta x_i = -\frac{f(x_i)}{f'(x_i)}
```

This recovers the Newton-Raphson formula.

### 5.4 Quadratic Convergence

The Newton-Raphson method has **quadratic convergence**, meaning:

```math
|x_{i+1} - x^*| \lesssim |x_i - x^*|^2
```

**Intuitive meaning:** When the constant $`M`$ below is of order 1, if the current error is $`\sim 10^{-2}`$, the next error will be $`\sim 10^{-4}`$. The number of correct digits roughly **doubles** each iteration.

**Formal statement:** Assume that $`f`$ is twice continuously differentiable on $`(a, b)`$. Suppose there exists $`x^* \in (a, b)`$ with $`f(x^*) = 0`$ and $`f'(x^*) \neq 0`$. Define the Newton-Raphson iteration as follows.

```math
x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}, \quad i = 1, 2, \ldots
```

Assume that $`x_i \to x^*`$ as $`i \to \infty`$. Then for sufficiently large $`i`$:

```math
|x_{i+1} - x^*| \leq M |x_i - x^*|^2 \quad \text{where} \quad M > \frac{|f''(x^*)|}{2|f'(x^*)|}
```

### 5.5 Proof of Quadratic Convergence

**Proof:**

Let $`e_i = x_i - x^*`$ (the error at step $`i`$), so $`x^* = x_i - e_i`$.

Expanding $`f`$ in a Taylor series about $`x_i`$ and evaluating at $`x_i - e_i`$ gives the following.

```math
f(x_i - e_i) = f(x_i) - f'(x_i)\,e_i + f''(\xi_i)\frac{e_i^2}{2}
```

where $`\xi_i`$ is between $`x_i`$ and $`x^*`$ (by Taylor's theorem with the Lagrange form of the remainder).

Since $`x^* = x_i - e_i`$ and $`f(x^*) = 0`$:

```math
0 = f(x_i) - f'(x_i)\,e_i + f''(\xi_i)\frac{e_i^2}{2}
```

Since $`f'`$ is continuous and $`f'(x^*) \neq 0`$, we have $`f'(x_i) \neq 0`$ as long as $`x_i`$ is close enough to $`x^*`$. Dividing by $`f'(x_i)`$ gives the following.

```math
0 = \frac{f(x_i)}{f'(x_i)} - e_i + \frac{f''(\xi_i)}{f'(x_i)}\frac{e_i^2}{2}
```

By definition of Newton's method, $`x_{i+1} = x_i - \dfrac{f(x_i)}{f'(x_i)}`$, so:

```math
x_{i+1} - x^* = x_i - x^* - \frac{f(x_i)}{f'(x_i)} = e_i - \frac{f(x_i)}{f'(x_i)}
```

From the equation above:

```math
0 = (x^* - x_{i+1}) + \frac{f''(\xi_i)}{f'(x_i)} \cdot \frac{e_i^2}{2}
```

Therefore:

```math
x_{i+1} - x^* = \frac{f''(\xi_i)}{f'(x_i)} \cdot \frac{1}{2}(x_i - x^*)^2
```

Taking absolute values:

```math
|x_{i+1} - x^*| = \frac{|f''(\xi_i)|}{2|f'(x_i)|} \cdot |x_i - x^*|^2
```

When $`x_i`$ is close enough to $`x^*`$, there is a constant $`M`$ with $`\dfrac{|f''(\xi_i)|}{2|f'(x_i)|} \leq M`$, so:

```math
\boxed{|x_{i+1} - x^*| \leq M |x_i - x^*|^2}
```

This proves **quadratic convergence**. $`\blacksquare`$

```python
def newton_raphson(f, df, x0, tol=1e-6, max_iter=100):
    """Newton-Raphson method"""
    x = x0
    for i in range(max_iter):
        fx = f(x)
        dfx = df(x)
        if abs(dfx) < 1e-15:
            print("Derivative near zero, method fails")
            return x
        delta_x = -fx / dfx            # Step 1: Newton step
        x = x + delta_x                # Step 2: Update
        if abs(delta_x) < tol:
            print(f"Converged in {i+1} iterations: x* = {x:.10f}")
            return x
    print("Did not converge")
    return x
```

### 5.6 Example 6.5: Slowly Converging Function

**Problem:** $`f(x) = x^{60} - 1`$

This is a **slowly converging function**: for $`0 < x < 1`$, away from the root, the slope is nearly zero, so the Newton step becomes huge, and after jumping far out to $`x > 1`$ each iteration reduces $`x`$ only by a factor of about $`59/60`$. The reason is that $`x - (x^{60} - 1)/(60x^{59}) \approx (59/60)x`$ for large $`x`$.

- With $`x_0 = 0.5`$, the slope is tiny ($`f'(0.5) = 60 \cdot 0.5^{59} \approx 1.04 \times 10^{-16}`$), so the first step jumps to $`x_1 \approx 9.6 \times 10^{15}`$. In exact arithmetic, each later iteration shrinks $`x`$ by about $`59/60`$, so the root is reached only after about 2,200 iterations, and in double precision $`x^{60}`$ overflows and the computation fails. The `newton_raphson` function above stops immediately at its check that the derivative is smaller than $`10^{-15}`$.
- With $`x_0 = 0.9`$, the guess is also far from the root $`x^* = 1`$. The first step jumps to $`x_1 \approx 9.23`$, so the function above prints "Did not converge" with its default setting (`max_iter=100`) and converges after 138 iterations when the iteration limit is raised. With a guess close enough to the root, such as $`x_0 = 0.99`$, it converges in 5 iterations.

> **Note:** **When Newton-Raphson struggles.** Newton-Raphson can fail or converge slowly in the following cases.
> - $`f'(x_i) \approx 0`$: a nearly horizontal tangent causes overshooting.
> - The initial guess is far from the root.
> - The function is very flat near the root.

### 5.7 Numerical Confirmation of Quadratic Convergence

Apply the Newton-Raphson method to the same root as Section 3.6 with $`f(x) = e^{-x} - x`$. Since $`f'(x) = -e^{-x} - 1`$, the update is:

```math
x_{i+1} = x_i - \frac{e^{-x_i} - x_i}{-e^{-x_i} - 1}
```

Starting from $`x_0 = 0`$ gives $`x_1 = 0 - (1 - 0)/(-1 - 1) = 0.5`$ and $`x_2 = 0.5 - (0.606531 - 0.5)/(-1.606531) = 0.566311`$. Continuing gives the following table.

| $`i`$ | $`x_i`$ | $`\vert e_i \vert = \vert x_i - x^* \vert`$ | $`\vert e_i \vert / \vert e_{i-1} \vert^2`$ |
|-----:|-----:|-----:|-----:|
| 0 | 0 | $`5.671 \times 10^{-1}`$ | - |
| 1 | 0.5 | $`6.714 \times 10^{-2}`$ | 0.2087 |
| 2 | 0.566311003197218 | $`8.323 \times 10^{-4}`$ | 0.1846 |
| 3 | 0.567143165034862 | $`1.254 \times 10^{-7}`$ | 0.1810 |
| 4 | 0.567143290409781 | $`2.776 \times 10^{-15}`$ | 0.1766 |

*Table 7. Errors of the Newton-Raphson method for f(x) = exp(-x) - x*

The exponents of the errors, $`-1, -2, -4, -7, -15`$, roughly double in each iteration, and four iterations reach the limit of double precision. The ratio in the last column approaches the constant of Section 5.4, $`|f''(x^*)|/(2|f'(x^*)|) = e^{-x^*}/(2(e^{-x^*} + 1)) = 0.5671/3.1343 \approx 0.1809`$, confirming quadratic convergence. For the same root, fixed-point iteration (Table 6) reached only an error of $`7.3 \times 10^{-4}`$ after 12 iterations.

### 5.8 When the Newton-Raphson Method Fails

The proof of Section 5.5 assumes that the initial guess is close enough to the root. When this assumption fails, the Newton-Raphson method can fail in the following ways.

- **Cycling.** For $`f(x) = x^3 - 2x + 2`$ with $`x_0 = 0`$, $`x_1 = 0 - 2/(-2) = 1`$ and $`x_2 = 1 - 1/1 = 0`$, so the iteration alternates between 0 and 1 forever. The tangents at the two points point to each other.
- **Divergence.** The root of $`f(x) = \tan^{-1} x`$ is 0, but starting from $`x_0 = 1.5`$ the iterates $`-1.694, 2.321, -5.114, 32.30, -1575, \ldots`$ grow in amplitude and diverge. Far from the root the tangent is nearly horizontal, so the step is too large. With $`x_0 = 1.3`$ the same function converges as $`-1.162, 0.859, -0.374, 0.034, \ldots`$, so convergence is sensitive to the initial guess.
- **Near an extremum.** If $`f'(x_i) = 0`$ the step cannot be computed, and if $`f'(x_i) \approx 0`$ the iterate jumps far away as in Section 5.6. This often happens for functions with an extremum near the root.
- **Near an inflection point.** If there is an inflection point ($`f''(x) = 0`$) near the root, the iteration can oscillate around the root and move away. The root 0 of $`\tan^{-1} x`$ is exactly an inflection point.

In practice these failures are prevented by setting a maximum number of iterations, falling back to bisection when the step does not shrink or $`|f|`$ grows, and, where possible, first obtaining an initial guess close enough to the root with a plot or a bracketing method. Brent's method in Section 8 automates this strategy.

### 5.9 Slow Convergence at Multiple Roots

The proof of Section 5.5 assumes $`f'(x^*) \neq 0`$. At a **multiple root**, that is, a root with $`f(x^*) = f'(x^*) = 0`$, this assumption fails. A root $`x^*`$ has **multiplicity** $`m`$ if $`f(x) = (x - x^*)^m h(x)`$ with $`h(x^*) \neq 0`$. Then $`f(x)/f'(x) \approx (x - x^*)/m`$, so the Newton step gives:

```math
e_{i+1} \approx e_i - \frac{e_i}{m} = \frac{m - 1}{m}\,e_i
```

At a multiple root the Newton-Raphson method therefore drops to linear convergence with asymptotic error constant $`(m - 1)/m`$. For a double root ($`m = 2`$) the error is halved in each iteration, the same speed as bisection. $`f(x) = (x - 3)(x - 1)^2`$ has a double root at $`x = 1`$, and starting from $`x_0 = 0`$ proceeds as follows.

| $`i`$ | Standard Newton $`x_i`$ | Error | Error Ratio | Modified Newton ($`m = 2`$) $`x_i`$ | Error |
|-----:|-----:|-----:|-----:|-----:|-----:|
| 1 | 0.428571 | $`5.714 \times 10^{-1}`$ | 0.5714 | 0.8571428571 | $`1.429 \times 10^{-1}`$ |
| 2 | 0.685714 | $`3.143 \times 10^{-1}`$ | 0.5500 | 0.9953917051 | $`4.608 \times 10^{-3}`$ |
| 3 | 0.832865 | $`1.671 \times 10^{-1}`$ | 0.5318 | 0.9999947092 | $`5.291 \times 10^{-6}`$ |
| 4 | 0.913330 | $`8.667 \times 10^{-2}`$ | 0.5186 | 1.0000000000 | $`6.998 \times 10^{-12}`$ |
| 5 | 0.955783 | $`4.422 \times 10^{-2}`$ | 0.5102 | - | - |
| 6 | 0.977655 | $`2.234 \times 10^{-2}`$ | 0.5053 | - | - |

*Table 8. Standard and modified Newton-Raphson methods at a double root*

The error ratio of the standard Newton-Raphson method approaches the theoretical value $`(2 - 1)/2 = 0.5`$. If the multiplicity $`m`$ is known, the **modified Newton-Raphson method**, which multiplies the step by $`m`$, $`x_{i+1} = x_i - m f(x_i)/f'(x_i)`$, recovers quadratic convergence; in the table the exponents of the errors double as $`-1, -3, -6, -12`$. If the multiplicity is unknown, apply the Newton-Raphson method to $`u(x) = f(x)/f'(x)`$. Since $`u`$ has a simple root (multiplicity 1) at $`x^*`$, the following iteration converges quadratically.

```math
x_{i+1} = x_i - \frac{u(x_i)}{u'(x_i)} = x_i - \frac{f(x_i) f'(x_i)}{[f'(x_i)]^2 - f(x_i) f''(x_i)}
```

In the same example this method converges as $`1.105263, 1.003082, 1.0000024, 1.0000000000014`$, at the cost of a second derivative in each iteration. At a multiple root $`f`$ and $`f'`$ approach zero together, so both the numerator and the denominator are strongly affected by roundoff, and the attainable accuracy itself is lower than for a simple root. At a root of even multiplicity, $`f`$ also does not change sign, so bracketing methods cannot be used (Section 7.6 of Chapter 5).

### 5.10 Newton-Raphson Method for Systems of Nonlinear Equations

The Newton-Raphson method extends to systems of nonlinear equations $`\mathbf{F}(\mathbf{x}) = \mathbf{0}`$ in several unknowns. Two things change from the single-variable case. The derivative $`f'`$ is replaced by the **Jacobian matrix** $`\mathbf{J}`$, which collects the partial derivatives $`\partial F_i/\partial x_j`$, and instead of the division $`f/f'`$, each iteration solves the linear system $`\mathbf{J}(\mathbf{x}^{(k)})\,\Delta\mathbf{x} = -\mathbf{F}(\mathbf{x}^{(k)})`$ and updates $`\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \Delta\mathbf{x}`$. Quadratic convergence near a root and the need for an initial guess close to the root carry over from one variable, and when the number of unknowns $`n`$ is large, computing the Jacobian matrix and solving an $`n \times n`$ linear system in every iteration becomes costly. The derivation and an example are given in Section 9 of Chapter 12.

---

<br>

## 6. Secant Method

### 6.1 Standard Secant Method

The secant method replaces the derivative in the Newton-Raphson method with a **finite difference approximation**. It is useful when **analytic derivatives are not available**.

We start from Newton-Raphson.

```math
x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}
```

Replace $`f'(x_i)`$ with a finite difference approximation using two previous points (the same form as the backward difference of Chapter 4).

```math
f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}}
```

Substituting gives the following.

```math
\boxed{x_{i+1} = x_i - \frac{f(x_i)(x_i - x_{i-1})}{f(x_i) - f(x_{i-1})}}
```

**Properties:**
- It requires **two initial guesses** $`x_0`$ and $`x_1`$ (like Wegstein), but they need not bracket the root.
- It is similar to the Wegstein method in using two guesses.
- Its convergence order is **superlinear** ($`\approx 1.618`$, the golden ratio), faster than linear but slower than quadratic.

> **Note:** **Difference from false position.** The secant method formula is identical in form to the false position formula (Chapter 5), but the secant method does **not** maintain a bracket. It always uses the two most recent points, whereas false position retains the bracket around the root.

```python
def secant_method(f, x0, x1, tol=1e-6, max_iter=100):
    """Standard secant method"""
    for i in range(max_iter):
        f0, f1 = f(x0), f(x1)
        if abs(f1 - f0) < 1e-15:
            print("Division by zero in secant method")
            return x1
        x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x_new - x1) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.10f}")
            return x_new
        x0, x1 = x1, x_new
    print("Did not converge")
    return x1
```

For example, `secant_method(lambda x: x - np.exp(-x), 0, 1)` converges to `x* = 0.5671432904` in 5 iterations.

### 6.2 Modified Secant Method

If the two points $`x_{i-1}`$ and $`x_i`$ are **not close enough**, we can use the **modified secant method**.

Instead of two separate points, we use a point obtained by perturbing $`x_i`$ by a small fraction. Writing the perturbation as $`\delta x_i`$ gives the following.

```math
f'(x_i) \approx \frac{f(x_i + \delta x_i) - f(x_i)}{\delta x_i}
```

Substituting into the Newton-Raphson formula gives the following.

```math
\boxed{x_{i+1} = x_i - \frac{f(x_i) \cdot \delta x_i}{f(x_i + \delta x_i) - f(x_i)}}
```

where $`\delta x_i`$ is a small perturbation (e.g., $`\delta x_i = \epsilon \cdot x_i`$ for some small $`\epsilon`$).

> **Note:** **Advantage.** It requires only **one** initial guess (like Newton-Raphson) and avoids computing the analytic derivative. The perturbation $`\delta x_i`$ controls the finite difference step size.

```python
def modified_secant(f, x0, delta=1e-4, tol=1e-6, max_iter=100):
    """Modified secant method using perturbation"""
    x = x0
    for i in range(max_iter):
        fx = f(x)
        dx = delta * x if abs(x) > 1e-10 else delta
        f_perturbed = f(x + dx)
        if abs(f_perturbed - fx) < 1e-15:
            print("Denominator near zero")
            return x
        x_new = x - fx * dx / (f_perturbed - fx)
        if abs(x_new - x) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.10f}")
            return x_new
        x = x_new
    print("Did not converge")
    return x
```

For example, `modified_secant(lambda x: x - np.exp(-x), 1.0)` converges to the same root in 4 iterations.

### 6.3 Order of Convergence of the Secant Method

The order of convergence 1.618 of the secant method is derived as follows. Approximating $`f`$ near the root by its second-order Taylor polynomial and substituting into the secant formula gives the following relation between the errors.

```math
e_{i+1} \approx \frac{f''(x^*)}{2f'(x^*)}\,e_i\,e_{i-1}
```

This is the $`e_{i+1} \approx M e_i^2`$ of the Newton-Raphson method with $`e_i e_{i-1}`$ in place of $`e_i^2`$. The previous error $`e_{i-1}`$ is larger than $`e_i`$, so the method is slower than Newton-Raphson. Assuming $`|e_{i+1}| \approx C|e_i|^p`$ gives $`|e_i| \approx C|e_{i-1}|^p`$, that is, $`|e_{i-1}| \approx (|e_i|/C)^{1/p}`$; substituting this into the relation above and comparing the exponents of $`|e_i|`$ on both sides gives $`p = 1 + 1/p`$, that is, $`p^2 - p - 1 = 0`$. The positive root is the golden ratio $`p = (1 + \sqrt{5})/2 \approx 1.618`$.

Applying the secant method to the same function as Section 5.7 with $`x_0 = 0`$ and $`x_1 = 1`$ gives the following table. The last column is the order estimated from three successive errors, $`p \approx \ln(|e_{i}|/|e_{i-1}|)/\ln(|e_{i-1}|/|e_{i-2}|)`$.

| $`i`$ | $`x_i`$ | $`\vert e_i \vert`$ | Estimated Order $`p`$ |
|-----:|-----:|-----:|-----:|
| 0 | 0 | $`5.671 \times 10^{-1}`$ | - |
| 1 | 1 | $`4.329 \times 10^{-1}`$ | - |
| 2 | 0.612699836780282 | $`4.556 \times 10^{-2}`$ | - |
| 3 | 0.563838389161074 | $`3.305 \times 10^{-3}`$ | 1.17 |
| 4 | 0.567170358419745 | $`2.707 \times 10^{-5}`$ | 1.83 |
| 5 | 0.567143306604963 | $`1.620 \times 10^{-8}`$ | 1.54 |
| 6 | 0.567143290409705 | $`7.927 \times 10^{-14}`$ | 1.65 |

*Table 9. Errors and estimated order of the secant method for f(x) = exp(-x) - x*

The estimated order fluctuates at first and then gathers near 1.618, and the exponents of the errors grow as $`-2, -3, -5, -8, -13`$, roughly the sum of the previous two exponents. Taking logarithms of $`e_{i+1} \propto e_i e_{i-1}`$ makes the exponents add like a Fibonacci sequence, and this is the same principle by which the ratio of successive Fibonacci numbers converges to the golden ratio.

### 6.4 Comparing the Convergence of Open and Bracketing Methods

Figure 3 plots the true error against the number of iterations when the root $`x^* \approx 0.567143`$ is found by five methods. The bracketing methods start from $`[0, 1]`$, fixed-point iteration from $`g(x) = e^{-x}`$ and $`x_0 = 0`$, the Newton-Raphson method from $`x_0 = 0`$, and the secant method from $`x_0 = 0`$ and $`x_1 = 1`$.

![Figure 3. Errors of five methods for the same root (semilog axes)](../../images/Ch06_Convergence-Comparison.png)

*Figure 3. Errors of five methods for the same root (semilog axes)*

On semilog axes a linearly converging method appears as a straight line whose slope reflects the asymptotic error constant $`C`$. Fixed-point iteration ($`C \approx 0.567`$) is the shallowest, bisection ($`C = 1/2`$, fluctuating) is similar, and false position, which for this nearly linear function has $`C \approx 0.11`$, is much steeper. The superlinear secant method and the quadratic Newton-Raphson method trace curves that grow steeper and bring the error below $`10^{-13}`$ in 5 and 4 iterations, respectively.

Compared by the number of iterations alone, the Newton-Raphson method is the fastest, but cost must be counted in function evaluations. The Newton-Raphson method evaluates $`f`$ and $`f'`$ once each per iteration, while the secant method evaluates only $`f`$ once. If computing $`f'`$ costs the same as $`f`$, the exponent of the error grows by a factor of $`2^{1/2} \approx 1.414`$ per function evaluation for Newton-Raphson and $`1.618`$ for the secant method. When the derivative is expensive or no analytical derivative exists, the secant method can therefore be faster in practice. Bracketing methods, on the other hand, are slow but guaranteed to converge, so combining the two families, as Brent's method in Section 8 does, is the practical choice.

---

<br>

## 7. Inverse Quadratic Interpolation

### 7.1 Concept

Suppose we have **three points** $`(x_{i-2}, y_{i-2})`$, $`(x_{i-1}, y_{i-1})`$, $`(x_i, y_i)`$ where $`y_k = f(x_k)`$.

The key idea is as follows.
- A forward quadratic $`y = p(x)`$ (a quadratic in $`x`$) through three points may **not cross the x-axis**.
- But the **inverse function** $`x = g(y)`$ through these three points is a quadratic in $`y`$, so it always **crosses** the x-axis. That is, we can evaluate $`g(0)`$ to find the root estimate.

We fit a **quadratic function of $`y`$**.

```math
x = ay^2 + by + c
```

When the three $`y`$-values are distinct, the 3 unknowns $`(a, b, c)`$ and 3 equations (one for each point) determine the system uniquely.

### 7.2 Lagrange Polynomial Formulation

Using **Lagrange polynomials**, the inverse quadratic interpolation is:

```math
g(y) = \frac{(y - y_{i-1})(y - y_i)}{(y_{i-2} - y_{i-1})(y_{i-2} - y_i)} x_{i-2} + \frac{(y - y_{i-2})(y - y_i)}{(y_{i-1} - y_{i-2})(y_{i-1} - y_i)} x_{i-1} + \frac{(y - y_{i-2})(y - y_{i-1})}{(y_i - y_{i-2})(y_i - y_{i-1})} x_i
```

This can be written compactly as:

```math
g(y) = \ell_{i-2}(y) \cdot x_{i-2} + \ell_{i-1}(y) \cdot x_{i-1} + \ell_i(y) \cdot x_i
```

where $`\ell_k(y)`$ are the Lagrange basis polynomials.

The next root estimate is $`x_{i+1} = g(0)`$, i.e., we set $`y = 0`$ in the formula above.

> **Note:** **Falling back to the secant method.** If two $`y`$-values coincide (e.g., $`y_{i-2} = y_{i-1}`$), a denominator of the Lagrange polynomial becomes zero and the formula degenerates. In this case, fall back to the **secant method** using two points instead.

### 7.3 Comparison with Secant Method

| Feature | Secant Method | Inverse Quadratic Interpolation |
|---|---|---|
| Points used | 2 | 3 |
| Interpolation | Linear (straight line through 2 points on $`f`$) | Quadratic (parabola through 3 points on inverse) |
| Convergence order | $`\approx 1.618`$ (golden ratio) | $`\approx 1.839`$ |
| Robustness | More robust | Can fail if $`y`$-values coincide |

*Table 10. Comparison of the secant method and inverse quadratic interpolation*

---

<br>

## 8. Brent's Method

### 8.1 Overview

**Brent's method** (Richard Brent, 1973) is a hybrid root-finding algorithm that combines the following methods.
- **Bisection method:** a bracketing method with guaranteed convergence.
- **Secant method:** a fast open method.
- **Inverse quadratic interpolation:** even faster when applicable.

```math
\text{Brent's method} = \text{Bracketing methods} + \text{Open methods}
```

**Strategy:** Use the faster open methods (secant or IQI) when they produce a result within the current bracket, and fall back to bisection when they fail or produce results outside the bracket.

```mermaid
flowchart TD
    A["Bracket [a, b] around the root"] --> B["Compute a candidate by IQI or secant"]
    B --> C{"Is the candidate inside the bracket?"}
    C -- yes --> D["Accept the candidate"]
    C -- no --> E["Take the bisection midpoint"]
    D --> F["Update the bracket by sign change"]
    E --> F
    F --> G{"Tolerance met?"}
    G -- no --> B
    G -- yes --> H["Output the root"]
```

*Figure 4. The strategy of Brent's method*

> **Key Point:** Brent's method inherits the **guaranteed convergence** of bisection while achieving the **speed** of open methods in favorable cases. It is the default root-finding algorithm in many numerical libraries.

### 8.2 Python Implementation with SciPy

SciPy provides Brent's method through `scipy.optimize.brentq`.

```python
from scipy import optimize

def f(x):
    return (x**2 - 1)

# Find root in [-2, 0]
root = optimize.brentq(f, -2, 0)
print(root)  # -1.0

# Find root in [0, 2]
root = optimize.brentq(f, 0, 2)
print(root)  # 1.0
```

**Function signature:** `brentq(f, a, b)`
- `f`: the function of the equation $`f(x) = 0`$ to be solved.
- `a`: the lower (left) end of the bracket, $`x_l`$.
- `b`: the upper (right) end of the bracket, $`x_u`$.
- As in bisection, $`f(x_l) \cdot f(x_u) < 0`$ is required.

> **Note:** **Practical view.** In real-world applications where a sign-changing bracket is available, Brent's method is usually preferred over pure bisection, Newton-Raphson, or secant methods because of its combination of reliability and speed.

### 8.3 Steps of Brent's Method

The detailed rules of Brent's method differ slightly between implementations, but the skeleton is as follows.

1. **Keeping three points.** Keep two points $`b`$ and $`c`$ that bracket the root, where $`b`$ is the best estimate so far ($`|f(b)| \leq |f(c)|`$). Remember the previous $`b`$ as $`a`$.
2. **Computing a candidate.** If the function values at the three points $`a`$, $`b`$, and $`c`$ are all different, compute a candidate by inverse quadratic interpolation (Section 7); otherwise use the secant method (Section 6.1).
3. **Testing the candidate.** Accept the candidate only if it lies within a prescribed range between $`b`$ and the opposite endpoint $`c`$ and the step has shrunk sufficiently compared with the previous steps. Otherwise use the bisection midpoint.
4. **Updating the bracket.** Update the bracket according to the sign of the function value at the new point, and swap points by comparing the magnitudes of the function values so that $`b`$ is again the best estimate.
5. **Stopping.** Stop when half the bracket, $`|c - b|/2`$, is below the tolerance (the sum of a relative and an absolute tolerance) or $`f(b) = 0`$.

Thanks to the test in step 3, the method converges superlinearly near a root where interpolation works well, and when interpolation fails or the step does not shrink enough, the bisection midpoint is used instead, so the bracket around the root is kept and the convergence guarantee of bisection is not lost. For ordinary smooth functions it therefore converges in about as many steps as the secant method, yet it does not diverge on functions for which interpolation fits poorly.

### 8.4 Comparing SciPy Root-Finding Functions

`scipy.optimize.newton` uses the Newton-Raphson method when the derivative `fprime` is given and the secant method otherwise. `root_scalar` uses bracketing methods when a bracket (`bracket`) is given and open methods when initial guesses (`x0`, `x1`) are given, and the method can be named with `method` [1]. The following code finds the root of Section 5.7 with four methods and compares the numbers of iterations and function calls.

```python
import numpy as np
from scipy import optimize

f = lambda x: np.exp(-x) - x
fprime = lambda x: -np.exp(-x) - 1

print(optimize.newton(f, 0.0, fprime=fprime))     # Newton-Raphson
print(optimize.newton(f, 0.0))                    # secant (no derivative)
for method, kw in [('bisect', dict(bracket=[0, 1])),
                   ('brentq', dict(bracket=[0, 1])),
                   ('secant', dict(x0=0.0, x1=1.0)),
                   ('newton', dict(x0=0.0, fprime=fprime))]:
    sol = optimize.root_scalar(f, method=method, xtol=1e-12, **kw)
    print(f"{method:7s} {sol.root:.15f} {sol.iterations:3d} {sol.function_calls:3d}")
```

```text
0.567143290409784
0.5671432904097838
bisect  0.567143290410058  40  42
brentq  0.567143290409784   6   7
secant  0.567143290409784   6   7
newton  0.567143290409784   5  10
```

For the same tolerance of $`10^{-12}`$, bisection needs 40 iterations and 42 function calls. Brent's method and the secant method each finish in 6 iterations and 7 function calls. The Newton-Raphson method needs the fewest iterations, 5, but evaluates $`f`$ and $`f'`$ in every iteration, so it makes 10 function calls. The efficiency comparison of Section 6.4 is confirmed by the actual call counts. When a bracket is known, `brentq`, which combines guaranteed convergence with speed, is the default choice; when no bracket is known but a good initial guess and a cheap derivative are available, the Newton-Raphson method is suitable.

---

<br>

## 9. Summary Table

| Method | Initial Guesses | Derivative Required | Convergence Order | Convergence Guarantee | Key Formula |
|---|---|---|---|---|---|
| **Fixed-Point Iteration** | 1 | No | Linear ($`\sim q^i`$) | Only if $`\vert g'(x^*)\vert < 1`$ | $`x_{i+1} = g(x_i)`$ |
| **Wegstein** | 2 | No | Superlinear | No | $`x_{i+1} = \frac{x_i g(x_{i-1}) - x_{i-1} g(x_i)}{x_i - x_{i-1} - g(x_i) + g(x_{i-1})}`$ |
| **Newton-Raphson** | 1 | Yes ($`f'`$) | Quadratic ($`\sim e_i^2`$) | No | $`x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}`$ |
| **Secant** | 2 | No | Superlinear ($`\approx 1.618`$) | No | $`x_{i+1} = x_i - \frac{f(x_i)(x_i - x_{i-1})}{f(x_i) - f(x_{i-1})}`$ |
| **Modified Secant** | 1 | No | Superlinear | No | $`x_{i+1} = x_i - \frac{f(x_i)\,\delta x_i}{f(x_i + \delta x_i) - f(x_i)}`$ |
| **Inverse Quadratic Interp.** | 3 | No | $`\approx 1.839`$ | No | Lagrange inverse polynomial at $`y = 0`$ |
| **Brent's Method** | 2 (bracket) | No | Superlinear (hybrid) | Yes (via bisection fallback) | Adaptive: bisection + secant + IQI |

*Table 11. Comparison of the open methods and Brent's method*

---

<br>

## 10. Roots of Polynomials

Polynomials are a special case of nonlinear equations. A polynomial of degree $`n`$ is written as:

```math
f_n(x) = a_1 x^n + a_2 x^{n-1} + \cdots + a_n x + a_{n+1}
```

By the fundamental theorem of algebra, a polynomial of degree $`n`$ has exactly $`n`$ roots counted with multiplicity, and if all coefficients are real, the roots are real or occur in complex conjugate pairs. Polynomials appear often in curve fitting (Chapters 14-18) and in the characteristic equations of dynamic systems described by linear differential equations, and these applications often need all the roots (including complex ones) rather than a single root. The methods of this chapter are designed to find one real root, so all the roots of a polynomial are found with one of the following two families of methods.

**One root at a time with deflation.** Find one root, divide the polynomial by the corresponding linear factor to reduce its degree by one (**deflation**), and find the next root in the reduced polynomial. For example, after finding the root 3 of $`x^3 - 6x^2 + 11x - 6`$, dividing by $`x - 3`$ gives $`x^2 - 3x + 2`$, and the remaining two roots follow from the quadratic formula. Repeated deflation accumulates the errors of the roots already found into the coefficients of the next polynomial, so roots obtained from deflated polynomials are polished with a few Newton-Raphson steps on the original polynomial. Finding complex roots as well requires methods that work in complex arithmetic (such as Muller's method, which uses quadratic interpolation).

**Conversion to an eigenvalue problem.** Dividing by $`a_1 \neq 0`$ gives $`x^n = -\frac{a_2}{a_1}x^{n-1} - \cdots - \frac{a_{n+1}}{a_1}`$. The $`n \times n`$ matrix with $`-a_2/a_1, -a_3/a_1, \ldots, -a_{n+1}/a_1`$ in its first row, ones on the subdiagonal below it, and zeros elsewhere is called the **companion matrix**, and its eigenvalues equal the roots of the polynomial. All the roots can therefore be found at once with an eigenvalue algorithm (Chapter 13). NumPy's `np.roots` uses this approach [2].

```python
import numpy as np

p = [1, -6, 11, -6]                 # x^3 - 6x^2 + 11x - 6
print(np.roots(p))
print(np.roots([1, 2, 5]))          # x^2 + 2x + 5
print(np.poly([1, 2, 3]))           # coefficients from roots
print(np.polyval(p, 2.5))           # p(2.5)

C = np.array([[6.0, -11.0, 6.0],    # companion matrix of p
              [1.0,   0.0, 0.0],
              [0.0,   1.0, 0.0]])
print(np.linalg.eigvals(C))
```

```text
[3. 2. 1.]
[-1.+2.j -1.-2.j]
[ 1. -6. 11. -6.]
-0.375
[3. 2. 1.]
```

- The `p` of `np.roots(p)` holds the coefficients from the highest-degree term to the constant term. It returns the roots 3, 2, 1 of $`x^3 - 6x^2 + 11x - 6 = (x - 1)(x - 2)(x - 3)`$.
- $`x^2 + 2x + 5`$ has discriminant $`4 - 20 < 0`$, so it has the complex conjugate roots $`-1 \pm 2i`$, and the result is a complex array.
- `np.poly` conversely builds coefficients from roots, and `np.polyval(p, x)` evaluates the polynomial: $`p(2.5) = 15.625 - 37.5 + 27.5 - 6 = -0.375`$.
- `C` is the companion matrix of this polynomial, whose first row is $`-(-6)/1, -11/1, -(-6)/1`$. That `np.linalg.eigvals(C)` gives the same values as `np.roots` confirms that the two approaches agree.

> **Caution:** The roots of a polynomial can be very sensitive to small changes in the coefficients. Especially when roots are close together or multiple, roundoff error in the coefficients alone can move the roots considerably, and a real multiple root may even be computed as a nearby pair of complex roots. When finding the roots of high-degree polynomials, keep this ill-conditioning in mind and check the results by substituting them into the original polynomial.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Open Methods | They start from one or two initial guesses without a bracket, are not guaranteed to converge, and are fast when they converge. |
| Fixed-Point Iteration | It rewrites the equation as $`x = g(x)`$ and repeats $`x_{i+1} = g(x_i)`$, converging linearly if $`\vert g'\vert < 1`$ near the fixed point. |
| Rate of Fixed-Point Iteration | It converges linearly with asymptotic error constant $`\vert g'(x^*) \vert`$, and the true error is at most $`q/(1-q)`$ times the difference of successive approximations. |
| Wegstein Method | It uses the intersection of $`y = x`$ with the line through two points on $`y = g(x)`$ as the next estimate, which equals the secant method applied to $`x - g(x)`$. |
| Newton-Raphson Method | It repeats the tangent's x-intercept $`x_{i+1} = x_i - f(x_i)/f'(x_i)`$ and converges quadratically near a root with $`f'(x^*) \neq 0`$. |
| Multiple Roots | At a root of multiplicity $`m`$ the Newton-Raphson method drops to linear convergence with $`C = (m-1)/m`$, and multiplying the step by $`m`$ or applying it to $`f/f'`$ restores quadratic convergence. |
| Newton-Raphson for Systems | Each iteration solves $`\mathbf{J}\Delta\mathbf{x} = -\mathbf{F}`$ and updates, converging quadratically near a root. |
| Secant Method | It approximates the derivative by the difference of the two most recent points and converges superlinearly with order about 1.618. |
| Computational Efficiency | Per function evaluation, the exponent of the error grows by about 1.414 for Newton-Raphson and about 1.618 for the secant method. |
| Modified Secant Method | It approximates the derivative by a difference with a small perturbation of one point, so it works with a single initial guess. |
| Inverse Quadratic Interpolation | It uses $`g(0)`$ of the quadratic $`x = g(y)`$ through three points as the next estimate and converges with order about 1.839. |
| Brent's Method | It uses inverse quadratic interpolation or the secant step inside the bracket and falls back to bisection when needed, combining guaranteed convergence with speed. |
| Roots of Polynomials | A polynomial of degree $`n`$ has $`n`$ roots, found one at a time with deflation and Newton-Raphson or all at once as eigenvalues of the companion matrix (`np.roots`). |

---

<br>

## Self-Check Questions

1. **Convergence of fixed-point iteration:** To find the root $`x^* = 3`$ of $`x^2 - 2x - 3 = 0`$, consider the three rearrangements $`g_1(x) = \sqrt{2x + 3}`$, $`g_2(x) = 3/(x - 2)`$, and $`g_3(x) = (x^2 - 3)/2`$. Which rearrangements converge near $`x^* = 3`$?

   <details>
   <summary>Answer</summary>

   Since $`g_1'(x) = 1/\sqrt{2x + 3}`$, $`|g_1'(3)| = 1/3 < 1`$, so it converges (from $`x_0 = 4`$ it goes 3.317, 3.104, 3.034, ...). Since $`g_2'(x) = -3/(x - 2)^2`$, $`|g_2'(3)| = 3 > 1`$, and since $`g_3'(x) = x`$, $`|g_3'(3)| = 3 > 1`$, so these two rearrangements move away from $`x^* = 3`$.

   </details>

2. **One Newton-Raphson step:** For $`f(x) = x^2 - 2`$ and $`x_0 = 1`$, find $`x_1`$ and $`x_2`$ by the Newton-Raphson method.

   <details>
   <summary>Answer</summary>

   Since $`f'(x) = 2x`$, $`x_1 = 1 - (-1)/2 = 1.5`$ and $`x_2 = 1.5 - 0.25/3 \approx 1.416667`$. This iteration gives the same values as the Babylonian method of the Chapter 4 lab.

   </details>

3. **One secant step:** For $`f(x) = x^2 - 2`$ with $`x_0 = 1`$ and $`x_1 = 2`$, find $`x_2`$ by the secant method.

   <details>
   <summary>Answer</summary>

   Since $`f(1) = -1`$ and $`f(2) = 2`$, $`x_2 = 2 - \frac{2 \cdot (2 - 1)}{2 - (-1)} = 2 - 2/3 \approx 1.333333`$.

   </details>

4. **Error under quadratic convergence:** In the Newton-Raphson method with $`M = 0.5`$ and a current error $`|e_i| = 10^{-3}`$, what are the error bounds for the next two iterations?

   <details>
   <summary>Answer</summary>

   $`|e_{i+1}| \leq 0.5 \times (10^{-3})^2 = 5 \times 10^{-7}`$ and $`|e_{i+2}| \leq 0.5 \times (5 \times 10^{-7})^2 = 1.25 \times 10^{-13}`$. The number of correct digits roughly doubles each iteration.

   </details>

5. **Cause of slow convergence:** Applying the Newton-Raphson method to $`f(x) = x^{60} - 1`$ with $`x_0 = 0.9`$, find $`x_1`$ and explain why convergence is slow afterward.

   <details>
   <summary>Answer</summary>

   Since $`f(0.9) \approx -0.9982`$ and $`f'(0.9) = 60 \cdot 0.9^{59} \approx 0.120`$, $`x_1 = 0.9 + 0.9982/0.120 \approx 9.23`$. For large $`x`$, $`x_{i+1} \approx (59/60)x_i`$, so each iteration reduces $`x`$ by only about 1.7%, and quadratic convergence begins only near the root. That is why 138 iterations are needed to reach the tolerance $`10^{-6}`$.

   </details>

6. **Multiple roots and Newton-Raphson:** Apply the Newton-Raphson method to $`f(x) = (x - 2)^3`$ with $`x_0 = 3`$ and find $`x_1`$, $`x_2`$, and the error ratio. What happens with the modified Newton-Raphson method?

   <details>
   <summary>Answer</summary>

   $`f/f' = (x - 2)^3/(3(x - 2)^2) = (x - 2)/3`$, so $`x_{i+1} = x_i - (x_i - 2)/3`$ and the error is multiplied by exactly $`2/3`$ in each iteration. $`x_1 = 3 - 1/3 \approx 2.6667`$ and $`x_2 = 2.6667 - 0.6667/3 \approx 2.4444`$. This equals the theoretical value $`(m - 1)/m = 2/3`$ for multiplicity $`m = 3`$. The modified Newton-Raphson method $`x_{i+1} = x_i - 3f/f'`$ gives $`x_1 = 3 - 3 \times 1/3 = 2`$, reaching the root in one step.

   </details>

7. **One Newton step for a system:** Update $`x^2 + y^2 = 2`$, $`x - y = 0`$ once by the Newton-Raphson method from $`\mathbf{x}^{(0)} = (2, 1)`$.

   <details>
   <summary>Answer</summary>

   $`\mathbf{F} = (x^2 + y^2 - 2,\ x - y) = (3, 1)`$ and $`\mathbf{J} = \begin{bmatrix} 2x & 2y \\ 1 & -1 \end{bmatrix} = \begin{bmatrix} 4 & 2 \\ 1 & -1 \end{bmatrix}`$. Solving $`4\Delta x + 2\Delta y = -3`$ and $`\Delta x - \Delta y = -1`$ gives $`\Delta y = 1/6`$ and $`\Delta x = -5/6`$, so $`\mathbf{x}^{(1)} = (2 - 5/6,\ 1 + 1/6) = (7/6, 7/6) \approx (1.1667, 1.1667)`$. It has moved closer to the true root $`(1, 1)`$.

   </details>

8. **Error estimate for fixed-point iteration:** In a fixed-point iteration with contraction constant $`q = 0.9`$, the difference of two successive approximations was $`10^{-4}`$. Find an upper bound on the true error, and find how small the successive difference must be to guarantee a true error of at most $`10^{-4}`$.

   <details>
   <summary>Answer</summary>

   By the formula of Section 3.10, the true error is at most $`q/(1 - q) \times 10^{-4} = 9 \times 10^{-4}`$. Guaranteeing a true error of at most $`10^{-4}`$ requires $`9|x_{i+1} - x_i| \leq 10^{-4}`$, that is, a successive difference of at most about $`1.11 \times 10^{-5}`$. The closer $`q`$ is to 1, the smaller than the tolerance the successive difference must be.

   </details>

---

<br>

## References

1. SciPy Developers, "scipy.optimize.root_scalar," *SciPy Reference Guide*. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html
2. NumPy Developers, "numpy.roots," *NumPy Reference*. [Online]. Available: https://numpy.org/doc/stable/reference/generated/numpy.roots.html
