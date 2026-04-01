# Chapter 5 Lecture -- Roots: Bracketing Methods

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 5

> **Prerequisites**: [Calculus] Continuity, intermediate value theorem (Ch 1-4).
>
> **Learning Objectives**:
> 1. Apply bisection method for root finding
> 2. Implement false position method
> 3. Analyze convergence rates of bracketing methods

---

<br>

## Table of Contents

1. [Motivation: The Quadratic Formula](#1-motivation-the-quadratic-formula)
2. [What Are Roots Problems?](#2-what-are-roots-problems)
3. [Graphical Methods](#3-graphical-methods)
   - 3.1 [Savings Account Example](#31-savings-account-example)
   - 3.2 [Graphical Root Estimation](#32-graphical-root-estimation)
4. [Bracketing vs. Open Methods](#4-bracketing-vs-open-methods)
5. [Incremental Search](#5-incremental-search)
   - 5.1 [Sign-Change Principle](#51-sign-change-principle)
   - 5.2 [Algorithm and Drawbacks](#52-algorithm-and-drawbacks)
6. [Bisection Method (Half-Interval Method)](#6-bisection-method-half-interval-method)
   - 6.1 [Algorithm](#61-algorithm)
   - 6.2 [Step-by-Step Procedure](#62-step-by-step-procedure)
7. [Convergence of the Bisection Method](#7-convergence-of-the-bisection-method)
   - 7.1 [Error Bound Derivation](#71-error-bound-derivation)
   - 7.2 [A Priori Error Estimation](#72-a-priori-error-estimation)
   - 7.3 [Required Number of Iterations](#73-required-number-of-iterations)
8. [False Position Method (Regula Falsi)](#8-false-position-method-regula-falsi)
   - 8.1 [Derivation of the Formula](#81-derivation-of-the-formula)
   - 8.2 [Convergence Properties](#82-convergence-properties)
   - 8.3 [One-Sidedness Problem](#83-one-sidedness-problem)
9. [Modified False Position Method](#9-modified-false-position-method)
   - 9.1 [Introducing Scaling Parameters](#91-introducing-scaling-parameters)
   - 9.2 [Two Consecutive Points on Same Side](#92-two-consecutive-points-on-same-side)
10. [Summary Table](#10-summary-table)

---

<br>

## 1. Motivation: The Quadratic Formula

For a quadratic equation:

$$f(x) = ax^2 + bx + c = 0$$

we can derive the roots analytically through completing the square:

> **[Calculus]** Completing the square is a fundamental algebraic technique that transforms a quadratic into vertex form, providing a direct path to the quadratic formula.

**Step 1.** Factor out $a$:

$$a\left(x^2 + \frac{b}{a}x\right) + c = 0$$

**Step 2.** Complete the square:

$$a\left(x + \frac{b}{2a}\right)^2 - \frac{b^2}{4a} + c = 0$$

**Step 3.** Solve for $x$ (assuming $a \neq 0$):

$$\left(x + \frac{b}{2a}\right)^2 = \frac{b^2 - 4ac}{4a^2}$$

$$x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$$

The roots are the values of $x$ where $f(x) = 0$. Geometrically, the roots are the points where the curve **crosses the x-axis**.

However, for higher-degree polynomials and transcendental equations, **no closed-form formula** generally exists. We need **numerical methods** to find roots:

- **Chapter 5**: Bracketing methods (guaranteed convergence, slower)
- **Chapter 6**: Open methods (may diverge, faster)

---

<br>

## 2. What Are Roots Problems?

"Roots" problems occur when some function $f$ can be written in terms of one or more dependent variables $x$, where the solutions to $f(x) = 0$ yield the solution to the problem.

These problems often occur when a **design problem** presents an **implicit equation** for a required parameter.

> **[Calculus]** The Intermediate Value Theorem (IVT) is the theoretical foundation for all bracketing methods: if $f$ is continuous on $[a, b]$ and $f(a) \cdot f(b) < 0$, then there exists at least one $x^* \in (a, b)$ such that $f(x^*) = 0$.

---

<br>

## 3. Graphical Methods

### 3.1 Savings Account Example

**Problem Setup:** Shin saves money by making regular monthly deposits of $p$ to his account. The annual interest rate is $r$. The total amount $A$ after $n$ payments is:

$$A = p\left(1 + \frac{r}{12}\right) + p\left(1 + \frac{r}{12}\right)^2 + \cdots + p\left(1 + \frac{r}{12}\right)^n$$

> **[Calculus]** This is a geometric series with first term $p(1 + r/12)$ and common ratio $(1 + r/12)$.

Applying the geometric series formula:

$$A = \sum_{i=1}^{n} p\left(1 + \frac{r}{12}\right)^i = p \cdot \frac{\left(1 + \frac{r}{12}\right)^n - 1}{\left(1 + \frac{r}{12}\right) - 1} = \frac{12p}{r}\left[\left(1 + \frac{r}{12}\right)^n - 1\right]$$

**Root-Finding Formulation:** Suppose $A$ is the target amount. Shin wants to have $A$ within $n$ months. What is the interest rate $r$ to achieve this goal?

$$\text{Find } r: \quad f(r) = \frac{12p}{r}\left[\left(1 + \frac{r}{12}\right)^n - 1\right] - A = 0$$

This equation **cannot be solved analytically** for $r$ -- it is a root-finding problem.

### 3.2 Graphical Root Estimation

We can use the graphical approach to estimate the root by plotting $g(x) = f(r)$ and visually identifying where the curve crosses the x-axis ($g(x) = 0$).

```python
# ex
def func(interest, duration, deposit, target):
    ratio = 1. + interest / 12.
    return 12.0 * deposit / interest * (ratio**duration - 1.0) - target

from functools import partial

duration = 24
deposit = 60.
target = 1500.
g = partial(func, duration=duration, deposit=deposit, target=target)
```

By plotting over progressively narrower intervals:

| Zoom Level | Interest Range | Observation |
|:---:|:---:|:---|
| 1st | 1% -- 10% | Root is near ~4% |
| 2nd | 4.0% -- 5.0% | Root is near ~4.2% |
| 3rd | 4.20% -- 4.30% | Root is approximately **4.23%** |

> **[Numerical Analysis]** Graphical methods provide a rough initial estimate but are imprecise. They are typically used as a starting point for more systematic numerical methods.

---

<br>

## 4. Bracketing vs. Open Methods

Root-finding methods are classified into two families:

### Bracketing Methods

- Based on **two initial guesses** that **bracket** (enclose) the root
- The root lies within the interval $[x_l, x_u]$
- **Guaranteed to converge** to a solution
- Slower convergence rate

### Open Methods

- Can involve **one or more initial guesses**, but there is **no bracket**
- **Can diverge** (not guaranteed to converge)
- **Faster** convergence rate when they do converge

> **[Numerical Analysis]** The trade-off between guaranteed convergence and speed is a recurring theme. In practice, a hybrid strategy often starts with a bracketing method to get close, then switches to an open method for rapid refinement.

---

<br>

## 5. Incremental Search

### 5.1 Sign-Change Principle

If $f$ is **real and continuous** in the interval from $x_l$ to $x_u$, and $f$ changes sign on opposite sides of the root, then:

$$f(x_l) \cdot f(x_u) < 0$$

This means a root exists somewhere between $x_l$ (lower bound) and $x_u$ (upper bound).

> **[Calculus]** This is a direct application of the Intermediate Value Theorem. Continuity of $f$ on $[x_l, x_u]$ is a necessary condition.

### 5.2 Algorithm and Drawbacks

**Algorithm:**

1. Divide the interval $[x_l, x_u]$ into smaller subintervals
2. Iterate over each subinterval, checking if a root exists within it (via the sign-change test: $f(x_i) \cdot f(x_{i+1}) < 0$)
3. If a sign change is detected, the root is in that subinterval

**Drawback:** It is **sensitive to the increment length** (step size).

- If the step size is **too large**, root brackets may be **missed** because an even number of roots could fall between two consecutive points
- If the step size is **too small**, the method becomes **computationally expensive**

> **[Numerical Analysis]** A function like $\sin(x^2)$ can have many closely-spaced roots. If the increment is larger than the spacing between adjacent roots, pairs of roots will be missed entirely because the function returns to the same sign.

---

<br>

## 6. Bisection Method (Half-Interval Method)

The bisection method is the simplest bracketing method. It systematically halves the search interval at each iteration.

### 6.1 Algorithm

**The idea:**

1. Split the search space in half at each iteration by choosing the **midpoint** of the interval
2. Determine which of the two halves contains the root
3. Repeat steps 1 and 2 until the search interval is smaller than a tolerance $\varepsilon$

### 6.2 Step-by-Step Procedure

**Step 0 (Initialization):** Start with an interval $[x_l, x_u]$ such that $f(x_l) \cdot f(x_u) < 0$.

Set $a_0 = x_l$, $b_0 = x_u$.

Compute the midpoint:

$$c_0 = \frac{x_l + x_u}{2} = \frac{a_0 + b_0}{2}$$

**Step 1 (Determine which half):** Check sign changes:

- If $f(a_0) \cdot f(c_0) < 0$ : the root is in the **left** interval $[a_0, c_0]$
  - Set $a_1 = a_0$, $b_1 = c_0$
- If $f(c_0) \cdot f(b_0) < 0$ : the root is in the **right** interval $[c_0, b_0]$
  - Set $a_1 = c_0$, $b_1 = b_0$

**Step 2 (Convergence check):** If $(b_1 - a_1) < \varepsilon$, then the output is:

$$x^* \approx \frac{a_1 + b_1}{2}$$

Otherwise, continue the search with the interval $[a_1, b_1]$.

**General iteration:** At iteration $n$:

$$c_n = \frac{a_n + b_n}{2}$$

Then update the interval based on the sign-change test and repeat.

```python
def bisection(f, a, b, tol=1e-6, max_iter=100):
    """
    Bisection method for finding a root of f in [a, b].
    Assumes f(a) * f(b) < 0.
    """
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs")

    for n in range(max_iter):
        c = (a + b) / 2.0
        if (b - a) < tol:
            return c, n
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2.0, max_iter
```

---

<br>

## 7. Convergence of the Bisection Method

### 7.1 Error Bound Derivation

At each iteration, the error (distance between the true root $x^*$ and the midpoint estimate $c_n$) is bounded by half the current interval length.

**After the 1st iteration:**

$$\text{err}_1 = |x^* - c_0| \leq \frac{1}{2}(b_0 - a_0)$$

**After the 2nd iteration:**

$$\text{err}_2 = |x^* - c_1| \leq \frac{1}{2}(b_1 - a_1)$$

**After the 3rd iteration:**

$$\text{err}_3 = |x^* - c_2| \leq \frac{1}{2}(b_2 - a_2)$$

**After the 4th iteration:**

$$\text{err}_4 = |x^* - c_3| \leq \frac{1}{2}(b_3 - a_3)$$

Since each interval is half the previous one:

$$\frac{1}{2}(b_3 - a_3) = \frac{1}{2^2}(b_2 - a_2) = \frac{1}{2^3}(b_1 - a_1) = \frac{1}{2^4}(b_0 - a_0)$$

### 7.2 A Priori Error Estimation

After $n$ iterations, the error is bounded above by:

$$\boxed{\text{err}_n = |x^* - c_{n-1}| \leq \frac{\Delta x}{2^n}}$$

where $\Delta x = b_0 - a_0$ is the **initial interval width**.

> **[Numerical Analysis]** This is called an **a priori error estimation** because we can predict the maximum error *before* running the algorithm, based solely on the initial interval and the number of iterations.

**Practical examples:**
- After **10 iterations**: the original interval uncertainty is reduced by $2^{10} = 1024$ (about 3 decimal digits of accuracy)
- After **20 iterations**: reduced by $2^{20} \approx 10^6$ (about 6 decimal digits -- "a million-fold" improvement)

> **[Calculus]** The convergence is **linear** -- each iteration adds roughly one bit of accuracy. This is in contrast to higher-order methods (e.g., Newton's method) that converge quadratically.

### 7.3 Required Number of Iterations

**Question:** How many iterations are needed for a given error tolerance $\varepsilon_s$?

We need:

$$\text{err}_n \leq \frac{\Delta x}{2^n} < \varepsilon_s$$

Rearranging:

$$2^n > \frac{\Delta x}{\varepsilon_s}$$

Taking $\log_2$ of both sides:

$$\boxed{n > \log_2\left(\frac{\Delta x}{\varepsilon_s}\right)}$$

**Example:** If $\Delta x = 1$ and $\varepsilon_s = 10^{-5}$:

$$n > \log_2\left(\frac{1}{10^{-5}}\right) = \log_2(10^5) \approx 16.6$$

Therefore, **17 iterations** will guarantee the tolerance.

```python
import math

def bisection_iterations_needed(delta_x, tolerance):
    """Calculate minimum number of bisection iterations needed."""
    return math.ceil(math.log2(delta_x / tolerance))

# Example
n = bisection_iterations_needed(1.0, 1e-5)
print(f"Iterations needed: {n}")  # Output: 17
```

---

<br>

## 8. False Position Method (Regula Falsi)

### 8.1 Derivation of the Formula

The False Position (FP) method is similar to the bisection method, but instead of the **midpoint**, it proposes a new estimation of the root using **linear interpolation**.

It takes the point $c$ where a straight line between $(x_l, f(x_l))$ and $(x_u, f(x_u))$ crosses the x-axis.

> **[Linear Algebra]** The method constructs the equation of a line (secant line) passing through two points and finds its x-intercept. This is the simplest form of polynomial interpolation.

**Derivation:** The line through $(a, f(a))$ and $(b, f(b))$ is:

$$y - f(b) = \frac{f(b) - f(a)}{b - a}(x - b)$$

To find the x-intercept, set $y = 0$:

$$-f(b) \cdot \frac{b - a}{f(b) - f(a)} = x - b$$

$$x = b + \frac{-b \cdot f(b) + a \cdot f(b)}{f(b) - f(a)}$$

$$x = \frac{a \cdot f(b) - b \cdot f(a)}{f(b) - f(a)}$$

Therefore, the **False Position formula** at iteration $n$ is:

$$\boxed{c_n = \frac{a_n \cdot f(b_n) - b_n \cdot f(a_n)}{f(b_n) - f(a_n)}}$$

The algorithm then proceeds like bisection: check which subinterval $[a_n, c_n]$ or $[c_n, b_n]$ contains the root, update the interval, and repeat until $|f(c_n)| < \varepsilon$.

```python
def false_position(f, a, b, tol=1e-6, max_iter=100):
    """
    False Position (Regula Falsi) method.
    Assumes f(a) * f(b) < 0.
    """
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs")

    for n in range(max_iter):
        fa, fb = f(a), f(b)
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)

        if abs(fc) < tol:
            return c, n

        if fa * fc < 0:
            b = c
        else:
            a = c
    return c, max_iter
```

### 8.2 Convergence Properties

- False position **may or may not** converge more quickly than bisection
- It **does not** provide a guarantee of meeting a criterion within a given number of iterations
- We **cannot predict the error** using the number of iterations (unlike bisection's a priori bound)

> **[Numerical Analysis]** Unlike bisection, there is no simple formula $\text{err}_n \leq \Delta x / 2^n$ for false position. The convergence rate depends heavily on the shape (curvature) of $f$ near the root.

### 8.3 One-Sidedness Problem

A key weakness of the false position method is the **one-sidedness problem**: one of the bracket points remains **stationary** (does not move). This causes the false position method to exhibit **poor convergence**.

**Example:** $f(x) = 2x^3 - 4x^2 + 3x$ on $[-1, 1]$.

When the function has significant curvature, one endpoint gets "stuck" while the other slowly approaches the root from one side. The interval does not shrink symmetrically.

---

<br>

## 9. Modified False Position Method

### 9.1 Introducing Scaling Parameters

To guarantee the length of the interval containing a root tends to zero, the Regular-Falsi method can be **modified** by introducing scaling parameters $\alpha$ and $\beta$ (with $\alpha \cdot \beta = 1$):

**Original formula:**

$$c_n = \frac{a_n \cdot f(b_n) - b_n \cdot f(a_n)}{f(b_n) - f(a_n)}$$

**Modified formula:**

$$c_n = \frac{a_n \cdot \beta \cdot f(b_n) - b_n \cdot \alpha \cdot f(a_n)}{\beta \cdot f(b_n) - \alpha \cdot f(a_n)}$$

### 9.2 Two Consecutive Points on Same Side

The modification is triggered when **two consecutive midpoint estimates** $c_0$ and $c_1$ appear on the **same side** of the root (i.e., $f(c_0) \cdot f(c_1) > 0$).

**Case 1: Two consecutive points to the RIGHT of the root** (i.e., $f(c_0) \cdot f(c_1) > 0$ and the stationary endpoint is $a_n$):

Set $\alpha = \frac{1}{2}$, which halves the function value at the **stationary** endpoint:

$$c = \frac{a \cdot f(b) - b \cdot \frac{1}{2} f(a)}{f(b) - \frac{1}{2} f(a)}$$

By reducing $f(a)$ by half, the secant line tilts, moving the next estimate closer to the stationary side and thus closer to the actual root.

**Case 2: Two consecutive points to the LEFT of the root** (stationary endpoint is $b_n$):

Set $\beta = \frac{1}{2}$, which halves the function value at the stationary endpoint:

$$c = \frac{a \cdot \frac{1}{2} f(b) - b \cdot f(a)}{\frac{1}{2} f(b) - f(a)}$$

> **[Numerical Analysis]** The modified false position method guarantees that the interval width shrinks to zero, overcoming the one-sidedness problem of standard Regula Falsi. The Illinois algorithm is a well-known variant that uses this halving strategy.

```python
def modified_false_position(f, a, b, tol=1e-6, max_iter=100):
    """
    Modified False Position method (Illinois algorithm variant).
    Halves the function value at the stationary endpoint when
    two consecutive estimates appear on the same side.
    """
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs")

    fa, fb = f(a), f(b)
    side = 0  # Track which side was updated last

    for n in range(max_iter):
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)

        if abs(fc) < tol:
            return c, n

        if fa * fc < 0:
            # Root is in [a, c]
            if side == -1:
                # Two consecutive on the right side: halve f(a)
                fa *= 0.5
            b = c
            fb = fc
            side = -1
        else:
            # Root is in [c, b]
            if side == 1:
                # Two consecutive on the left side: halve f(b)
                fb *= 0.5
            a = c
            fa = fc
            side = 1

    return c, max_iter
```

---

<br>

## 10. Summary Table

| Property | Incremental Search | Bisection | False Position | Modified False Position |
|:---|:---:|:---:|:---:|:---:|
| **Type** | Bracketing | Bracketing | Bracketing | Bracketing |
| **Initial Requirement** | Interval $[x_l, x_u]$ | $f(x_l) \cdot f(x_u) < 0$ | $f(x_l) \cdot f(x_u) < 0$ | $f(x_l) \cdot f(x_u) < 0$ |
| **Midpoint Formula** | Uniform step | $c = \frac{a + b}{2}$ | $c = \frac{a f(b) - b f(a)}{f(b) - f(a)}$ | Scaled version of FP |
| **Convergence Guarantee** | Depends on step size | Always (if continuous) | Always (if continuous) | Always (if continuous) |
| **A Priori Error Bound** | None | $\frac{\Delta x}{2^n}$ | None | None |
| **Iterations for Tolerance** | Unknown | $n > \log_2(\Delta x / \varepsilon_s)$ | Unpredictable | Unpredictable |
| **One-Sidedness Issue** | N/A | No | Yes | No (corrected) |
| **Convergence Rate** | N/A | Linear | Varies (can stall) | Generally better than FP |
| **Key Advantage** | Simple to implement | Predictable, robust | Often faster than bisection | Fixes FP stalling |
| **Key Disadvantage** | Sensitive to step size | Slow (linear convergence) | Can stall (one-sided) | Slightly more complex |
