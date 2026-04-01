# Chapter 6 Lecture — Roots: Open Methods

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 6

> **Prerequisites**: [Calculus] Derivatives (Ch 1-5).
>
> **Learning Objectives**:
> 1. Apply Newton-Raphson method for root finding
> 2. Implement secant and modified Newton methods
> 3. Compare convergence of open vs bracketing methods

---

<br>

## Table of Contents

1. [Objectives](#1-objectives)
2. [Open Methods vs Bracketing Methods](#2-open-methods-vs-bracketing-methods)
3. [Fixed-Point Iteration](#3-fixed-point-iteration)
   - 3.1 [Concept and Definition](#31-concept-and-definition)
   - 3.2 [Algorithm](#32-algorithm)
   - 3.3 [Fixed Point and Root Relationship](#33-fixed-point-and-root-relationship)
   - 3.4 [Convergence Analysis](#34-convergence-analysis)
   - 3.5 [Convergence Proof via Mean Value Theorem](#35-convergence-proof-via-mean-value-theorem)
   - 3.6 [Example 6.1 — Convergent Case](#36-example-61--convergent-case)
   - 3.7 [Example 6.1 — Divergent Case](#37-example-61--divergent-case)
   - 3.8 [Example 6.1 — Cycling Case](#38-example-61--cycling-case)
4. [The Wegstein Method](#4-the-wegstein-method)
   - 4.1 [Motivation from False Position](#41-motivation-from-false-position)
   - 4.2 [Geometric Idea](#42-geometric-idea)
   - 4.3 [Derivation](#43-derivation)
   - 4.4 [General Formula](#44-general-formula)
5. [Newton-Raphson Method](#5-newton-raphson-method)
   - 5.1 [Geometric Derivation](#51-geometric-derivation)
   - 5.2 [Two-Step Approach](#52-two-step-approach)
   - 5.3 [Derivation from Taylor Expansion](#53-derivation-from-taylor-expansion)
   - 5.4 [Quadratic Convergence](#54-quadratic-convergence)
   - 5.5 [Proof of Quadratic Convergence](#55-proof-of-quadratic-convergence)
   - 5.6 [Example 6.5 — Slowly Converging Function](#56-example-65--slowly-converging-function)
6. [Secant Method](#6-secant-method)
   - 6.1 [Standard Secant Method](#61-standard-secant-method)
   - 6.2 [Modified Secant Method](#62-modified-secant-method)
7. [Inverse Quadratic Interpolation](#7-inverse-quadratic-interpolation)
   - 7.1 [Concept](#71-concept)
   - 7.2 [Lagrange Polynomial Formulation](#72-lagrange-polynomial-formulation)
   - 7.3 [Comparison with Secant Method](#73-comparison-with-secant-method)
8. [Brent's Method](#8-brents-method)
   - 8.1 [Overview](#81-overview)
   - 8.2 [Python Implementation with SciPy](#82-python-implementation-with-scipy)
9. [Summary Table](#9-summary-table)

<br>

---

<br>

## 1. Objectives

- Understand the **fixed-point iteration**
- Implement the **Wegstein method**
- Understand the quadratically convergent **Newton-Raphson method**
- Understand the **Secant method**
- Understand **Brent's method**

> **Prerequisite:** Chapter 5 — Roots: Bracketing Methods (bisection, false position). The concepts of bracketing, error tolerance, and iterative convergence from Chapter 5 are directly extended here.

<br>

---

<br>

## 2. Open Methods vs Bracketing Methods

Unlike **bracketing methods** (bisection, false position), open methods do **not** require two initial guesses that bracket a root. Key characteristics:

| Property | Bracketing Methods | Open Methods |
|---|---|---|
| Initial guesses | Two points $x_l, x_u$ with $f(x_l) \cdot f(x_u) < 0$ | One or two initial guesses (no sign-change requirement) |
| Convergence guarantee | Always converge (guaranteed) | May **diverge** or move away from the root |
| Speed | Slower (linear convergence) | Faster convergence **when they converge** |
| Interval behavior | Interval shrinks each iteration | Single point(s) iterate toward root |

> **Key Insight:** Open methods trade convergence guarantees for speed. If they converge, they approach the root much more quickly than bracketing methods.

<br>

---

<br>

## 3. Fixed-Point Iteration

### 3.1 Concept and Definition

Also known as:
- **Simple one-point iteration**
- **Method of successive substitutions**

Given an equation $f(x) = 0$ where $f: [a, b] \to \mathbb{R}$, we rearrange it into the form:

$$f(x) = x - g(x) = 0 \quad \Longrightarrow \quad x = g(x)$$

The iteration $x_{i+1} = g(x_i)$ is called **fixed-point iteration**.

**Definition:** A **fixed point** of a function $g$ is a real number $x^*$ such that:

$$g(x^*) = x^*$$

> **Prerequisite Connection:** While bracketing methods (Ch.5) work directly with $f(x) = 0$, fixed-point iteration reformulates the problem as $x = g(x)$, which can be viewed as finding the intersection of $y = g(x)$ and $y = x$.

### 3.2 Algorithm

To find the root $x^* \in (a, b)$:

1. Choose a starting point $x_0 \in (a, b)$
2. Sequentially calculate $x_1, x_2, \ldots, x_{i+1}$ using:

$$x_{i+1} = g(x_i)$$

3. Repeat until $x_{i+1} \approx x_i$ (i.e., convergence is achieved)

This is a **circular calculation**: the output $x_{i+1}$ from one step becomes the input $x_i$ for the next.

```
input x_i  -->  [ g ]  -->  output x_{i+1}
     ^                            |
     |____________________________|
```

### 3.3 Fixed Point and Root Relationship

Let $g$ be a continuous function. Let the sequence $\{x_i : i \geq 0\}$ be generated by the fixed-point iteration $x_{i+1} = g(x_i)$.

If $\displaystyle\lim_{i \to \infty} x_i = x^*$, then $x^*$ is a fixed point of $g$.

Thus $x^*$ is the **root** of $f$:

$$f(x^*) = x^* - g(x^*) = 0$$

### 3.4 Convergence Analysis

Convergence depends on the **slope of $g$** at the fixed point. Using the mean value theorem:

$$\text{err}_{i+1} = \text{err}_i \cdot |g'(\xi)|$$

where $\xi$ is some point between $x_i$ and $x^*$.

**Convergence criterion:**

| Condition | Result |
|---|---|
| $\|g'(\xi)\| < 1$ | $x_{i+1}$ **converges** |
| $\|g'(\xi)\| > 1$ | $x_{i+1}$ **diverges** |
| $\|g'(\xi)\| = 1$ | Inconclusive (may cycle) |

### 3.5 Convergence Proof via Mean Value Theorem

Let $x^* \in (a, b)$ be a fixed point of $g$ and let $x_i \in (a, b)$.

By the **Mean Value Theorem**, there exists $\xi_i \in (x_i, x^*)$ or $(x^*, x_i)$ such that:

$$g(x_i) - g(x^*) = g'(\xi_i)(x_i - x^*)$$

Since $x_i = g(x_{i-1})$ and $x^* = g(x^*)$:

$$|x_i - x^*| = |g(x_{i-1}) - g(x^*)| = |g'(\xi_{i-1})| \cdot |x_{i-1} - x^*|$$

If for any $\xi \in (a, b)$ we have $|g'(\xi)| \leq q < 1$, then:

$$|x_i - x^*| \leq q \cdot |x_{i-1} - x^*| \leq q^2 \cdot |x_{i-2} - x^*| \leq \cdots \leq q^i \cdot |x_0 - x^*|$$

Since $q < 1$, we have $q^i \to 0$ as $i \to \infty$, so:

$$|x_i - x^*| \to 0$$

> **Key Result:** The fixed-point iteration converges if the **contraction constant** $q = \max|g'(\xi)| < 1$ on the interval containing the root and iterates.

### 3.6 Example 6.1 — Convergent Case

**Problem:** Estimate the root of $f(x) = x - e^{-x}$.

**Setup:** From $f(x) = 0$, we view:

$$x = e^{-x} = g(x)$$

Starting with $x_0 = 0$:

| Iteration | $x_i$ | $g(x_i) = e^{-x_i}$ |
|---|---|---|
| 0 | 0 | $e^{0} = 1.0000$ |
| 1 | 1.0000 | $e^{-1} = 0.3679$ |
| 2 | 0.3679 | $e^{-0.3679} = 0.6922$ |
| 3 | 0.6922 | $\cdots$ |
| $\vdots$ | $\vdots$ | converges to $x^* \approx 0.5671$ |

**Why it converges:** $g(x) = e^{-x}$, so $g'(x) = -e^{-x}$. At the root $x^* \approx 0.5671$, $|g'(x^*)| = e^{-0.5671} \approx 0.5671 < 1$.

Graphically, this is a **spiral convergence** — the iterates spiral inward toward the intersection of $y = g(x) = e^{-x}$ and $y = x$.

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

### 3.7 Example 6.1 — Divergent Case

**Problem:** Estimate the root of $f(x) = x + \ln x$.

**Setup:** From $f(x) = 0$:

$$x = -\ln(x) = g(x)$$

**Result:** The iteration **diverges** (spiral outward).

**Why it diverges:** $g(x) = -\ln(x)$, so $g'(x) = -1/x$. At the root $x^* \approx 0.5671$, $|g'(x^*)| = 1/0.5671 \approx 1.763 > 1$.

> **Lesson:** The same root can be found by one rearrangement of $f(x) = 0$ but not another. The choice of $g(x)$ critically determines convergence.

### 3.8 Example 6.1 — Cycling Case

Consider $g(x) = \dfrac{1}{x}$ with $x_0 = 4$:

| Iteration | $x_i$ | $g(x_i) = 1/x_i$ |
|---|---|---|
| 0 | 4 | $1/4 = 0.25$ |
| 1 | 0.25 | $1/0.25 = 4$ |
| 2 | 4 | $1/4 = 0.25$ |

The iteration **cycles between two points** $x_0$ and $x_1$ indefinitely, never converging to the fixed point $x^* = 1$.

**Why:** $g'(x) = -1/x^2$, so $|g'(1)| = 1$ (boundary case).

<br>

---

<br>

## 4. The Wegstein Method

### 4.1 Motivation from False Position

Recall that in the **false position method** (Chapter 5), a straight line between two guesses on $f(x)$ crossed the x-axis to produce the next estimate.

The Wegstein method applies a similar linear interpolation idea, but to the function $g(x)$ instead of $f(x)$.

### 4.2 Geometric Idea

In the Wegstein method, a straight line between two guesses on the curve $y = g(x)$ **intersects with the line $y = x$** to determine the next estimate of the solution.

Given two points $(x_0, g(x_0))$ and $(x_1, g(x_1))$ on the curve $g(x)$:
- Draw a straight line through these two points
- Find where this line intersects $y = x$
- That intersection gives the next estimate $x_2$

### 4.3 Derivation

The slope of the line through $(x_0, g(x_0))$ and $(x_1, g(x_1))$ equals the slope through $(x_1, g(x_1))$ and $(x_2, g(x_2))$:

$$\text{slope} = \frac{g(x_1) - g(x_0)}{x_1 - x_0} = \frac{g(x_2) - g(x_1)}{x_2 - x_1}$$

Since the intersection is on the line $y = x$, we have $g(x_2) = x_2$. Therefore:

$$\frac{g(x_1) - g(x_0)}{x_1 - x_0} = \frac{x_2 - g(x_1)}{x_2 - x_1}$$

Cross-multiplying and solving for $x_2$:

$$x_2(g(x_1) - g(x_0)) - x_1(g(x_1) - g(x_0)) = (x_1 - x_0)x_2 - g(x_1)(x_1 - x_0)$$

$$x_2 = \frac{x_1 g(x_0) - x_0 g(x_1)}{x_1 - x_0 - g(x_1) + g(x_0)}$$

### 4.4 General Formula

The general Wegstein iteration formula is:

$$\boxed{x_{i+1} = \frac{x_i \, g(x_{i-1}) - x_{i-1} \, g(x_i)}{x_i - x_{i-1} - g(x_i) + g(x_{i-1})}}$$

> **Key Property:** The last two estimates are used to determine the next estimate. This requires **two initial guesses** (like the secant method), but they do not need to bracket the root.

> **Prerequisite Connection:** The Wegstein method is essentially the false position idea (Ch.5) applied to $g(x)$ intersecting $y = x$ rather than $f(x)$ intersecting $y = 0$.

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

<br>

---

<br>

## 5. Newton-Raphson Method

### 5.1 Geometric Derivation

Given $f(x) = 0$, we use the **first-order derivative** (tangent line) at the current estimate.

At the point $(x_i, f(x_i))$, the slope of the tangent line is:

$$\text{slope} = f'(x_i) = \frac{0 - f(x_i)}{x_{i+1} - x_i}$$

Solving for $x_{i+1}$:

$$\boxed{x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}}$$

### 5.2 Two-Step Approach

The Newton-Raphson method can be viewed as a two-step process:

**Step 1:** Compute the correction (Newton step):

$$\delta x_i = -\frac{f(x_i)}{f'(x_i)}$$

**Step 2:** Update the estimate:

$$x_{i+1} = x_i + \delta x_i$$

> **Prerequisite Connection:** Newton-Raphson is an open method requiring only **one** initial guess and the derivative $f'(x)$. Unlike bracketing methods, there is no guarantee of convergence, but when it converges, it does so quadratically.

### 5.3 Derivation from Taylor Expansion

We can derive the Newton-Raphson formula from the **Taylor expansion**:

$$f(x_{i+1}) = f(x_i + \delta x_i) = f(x_i) + f'(x_i)\,\delta x_i + \frac{f''(x_i)}{2}(\delta x_i)^2 + \cdots$$

Truncating after the first-order term:

$$f(x_{i+1}) \approx f(x_i) + f'(x_i)\,\delta x_i$$

We **hope** $f(x_{i+1}) = 0$ (i.e., $x_{i+1}$ is the root $x^*$). Setting the above to zero:

$$0 = f(x_i) + f'(x_i)\,\delta x_i$$

$$\delta x_i = -\frac{f(x_i)}{f'(x_i)}$$

This recovers the Newton-Raphson formula.

### 5.4 Quadratic Convergence

Newton-Raphson method has **quadratic convergence**, meaning:

$$|x_{i+1} - x^*| \lesssim |x_i - x^*|^2$$

**Intuitive meaning:** If the current error is $\sim 10^{-2}$, the next error will be $\sim 10^{-4}$. The number of correct digits roughly **doubles** each iteration.

**Formal statement:** Assume that $f$ is twice continuously differentiable on $(a, b)$. Suppose there exists $x^* \in (a, b)$ with $f'(x^*) \neq 0$. Define Newton-Raphson iteration:

$$x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}, \quad i = 1, 2, \ldots$$

Assume that $x_i \to x^*$ as $i \to \infty$. Then for sufficiently large $i$:

$$|x_{i+1} - x^*| \leq M |x_i - x^*|^2 \quad \text{where} \quad M \geq \frac{|f''(x^*)|}{2|f'(x^*)|}$$

### 5.5 Proof of Quadratic Convergence

**Proof:**

Let $e_i = x_i - x^*$ (the error at step $i$), so $x^* = x_i - e_i$.

Consider a Taylor expansion of $f$ at $x_i - e_i$:

$$f(x_i - e_i) = f(x_i) - f'(x_i)\,e_i + f''(\xi_i)\frac{e_i^2}{2}$$

where $\xi_i$ is between $x_i$ and $x^*$ (by the mean value theorem applied to $f''$).

Since $x^* = x_i - e_i$ and $f(x^*) = 0$:

$$0 = f(x_i) - f'(x_i)\,e_i + f''(\xi_i)\frac{e_i^2}{2}$$

Since $f'(x_i) \neq 0$ (as long as $x_i$ is close to $x^*$), divide by $f'(x_i)$:

$$0 = \frac{f(x_i)}{f'(x_i)} - e_i + \frac{f''(\xi_i)}{f'(x_i)}\frac{e_i^2}{2}$$

By definition of Newton's method, $x_{i+1} = x_i - \dfrac{f(x_i)}{f'(x_i)}$, so:

$$x_{i+1} - x^* = x_i - x^* - \frac{f(x_i)}{f'(x_i)} = e_i - \frac{f(x_i)}{f'(x_i)}$$

From the equation above:

$$0 = (x^* - x_{i+1}) + \frac{f''(\xi_i)}{f'(x_i)} \cdot \frac{e_i^2}{2}$$

Therefore:

$$x_{i+1} - x^* = \frac{f''(\xi_i)}{f'(x_i)} \cdot \frac{1}{2}(x_i - x^*)^2$$

Taking absolute values:

$$|x_{i+1} - x^*| = \frac{|f''(\xi_i)|}{2|f'(x_i)|} \cdot |x_i - x^*|^2$$

Since $\dfrac{|f''(\xi_i)|}{2|f'(x_i)|} \leq M$ for some constant $M$:

$$\boxed{|x_{i+1} - x^*| \leq M |x_i - x^*|^2}$$

This proves **quadratic convergence**. $\blacksquare$

```python
def newton_raphson(f, df, x0, tol=1e-6, max_iter=100):
    """Newton-Raphson method"""
    x = x0
    for i in range(max_iter):
        fx = f(x)
        dfx = df(x)
        if abs(dfx) < 1e-15:
            print("Derivative near zero — method fails")
            return x
        delta_x = -fx / dfx            # Step 1: Newton step
        x = x + delta_x                # Step 2: Update
        if abs(delta_x) < tol:
            print(f"Converged in {i+1} iterations: x* = {x:.10f}")
            return x
    print("Did not converge")
    return x
```

### 5.6 Example 6.5 — Slowly Converging Function

**Problem:** $f(x) = x^{60} - 1$

This is a **slowly converging function** because the slope is nearly zero for $x$ values away from the root.

- With $x_0 = 0.5$: the slope at $x = 0.5$ is very small ($f'(0.5) = 60 \cdot 0.5^{59} \approx 0$), causing huge Newton steps and potential divergence.
- With $x_0 = 0.9$: the method converges within **3 iterations** because $x_0$ is close enough to the root $x^* = 1$.

> **Lesson:** Newton-Raphson can fail or converge slowly if:
> - $f'(x_i) \approx 0$ (nearly horizontal tangent causes overshooting)
> - The initial guess is far from the root
> - The function is very flat near the root

<br>

---

<br>

## 6. Secant Method

### 6.1 Standard Secant Method

Secant methods are a **finite difference approximation** of the Newton-Raphson method. They are useful when **analytic derivatives are not available**.

Starting from Newton-Raphson:

$$x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}$$

Replace $f'(x_i)$ with a finite difference approximation using two previous points:

$$f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}}$$

Substituting:

$$\boxed{x_{i+1} = x_i - \frac{f(x_i)(x_i - x_{i-1})}{f(x_i) - f(x_{i-1})}}$$

**Properties:**
- Requires **two initial guesses** $x_0$ and $x_1$ (like Wegstein), but they need not bracket the root
- Similar to Wegstein method in using two guesses
- Convergence order is **superlinear** ($\approx 1.618$, the golden ratio), faster than linear but slower than quadratic

> **Prerequisite Connection:** The secant method formula is identical in form to the false position formula (Ch.5), but the secant method does **not** maintain a bracket. It always uses the two most recent points, whereas false position retains the bracket.

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

### 6.2 Modified Secant Method

If the two points $x_{i-1}$ and $x_i$ are **not close enough**, we can use the **modified secant method**.

Instead of using two separate points, we perturb $x_i$ by a small fraction $\delta x_i$:

$$f'(x_i) \approx \frac{f(x_i + \delta x_i) - f(x_i)}{\delta x_i}$$

Substituting into Newton-Raphson:

$$\boxed{x_{i+1} = x_i - \frac{f(x_i) \cdot \delta x_i}{f(x_i + \delta x_i) - f(x_i)}}$$

where $\delta x_i$ is a small perturbation (e.g., $\delta x_i = \epsilon \cdot x_i$ for some small $\epsilon$).

> **Advantage:** Only requires **one** initial guess (like Newton-Raphson) and avoids computing the analytic derivative. The perturbation $\delta x_i$ controls the finite difference step size.

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

<br>

---

<br>

## 7. Inverse Quadratic Interpolation

### 7.1 Concept

Suppose we have **three points** $(x_{i-2}, y_{i-2})$, $(x_{i-1}, y_{i-1})$, $(x_i, y_i)$ where $y_k = f(x_k)$.

The key idea:
- A forward quadratic $y = f(x)$ through three points may **not cross the x-axis**
- But the **inverse function** $x = g(y)$ through these three points **does cross** the x-axis (i.e., we can evaluate $g(0)$ to find the root estimate)

We fit a **quadratic function of $y$**:

$$x = ay^2 + by + c$$

With 3 unknowns $(a, b, c)$ and 3 equations (one for each point), the system is uniquely determined.

### 7.2 Lagrange Polynomial Formulation

Using **Lagrange polynomials**, the inverse quadratic interpolation is:

$$g(y) = \frac{(y - y_{i-1})(y - y_i)}{(y_{i-2} - y_{i-1})(y_{i-2} - y_i)} x_{i-2} + \frac{(y - y_{i-2})(y - y_i)}{(y_{i-1} - y_{i-2})(y_{i-1} - y_i)} x_{i-1} + \frac{(y - y_{i-2})(y - y_{i-1})}{(y_i - y_{i-2})(y_i - y_{i-1})} x_i$$

This can be written compactly as:

$$g(y) = \ell_{i-2}(y) \cdot x_{i-2} + \ell_{i-1}(y) \cdot x_{i-1} + \ell_i(y) \cdot x_i$$

where $\ell_k(y)$ are the Lagrange basis polynomials.

The next root estimate is $x_{i+1} = g(0)$, i.e., we set $y = 0$ in the formula above.

> **Fallback:** If $y_{i-2} = y_{i-1}$ (two $y$-values coincide), the Lagrange polynomial degenerates. In this case, fall back to the **secant method** using two points instead.

### 7.3 Comparison with Secant Method

| Feature | Secant Method | Inverse Quadratic Interpolation |
|---|---|---|
| Points used | 2 | 3 |
| Interpolation | Linear (straight line through 2 points on $f$) | Quadratic (parabola through 3 points on inverse) |
| Convergence order | $\approx 1.618$ (golden ratio) | $\approx 1.839$ |
| Robustness | More robust | Can fail if $y$-values coincide |

<br>

---

<br>

## 8. Brent's Method

### 8.1 Overview

**Brent's method** (Richard Brent, 1973) is a hybrid root-finding algorithm that combines:
- **Bisection method** (bracketing, guaranteed convergence)
- **Secant method** (fast open method)
- **Inverse quadratic interpolation** (even faster when applicable)

$$\text{Brent's method} = \text{Bracketing methods} + \text{Open methods}$$

**Strategy:** Use the faster open methods (secant or IQI) when they produce a result within the current bracket, and fall back to bisection when they fail or produce results outside the bracket.

> **Key Advantage:** Brent's method inherits the **guaranteed convergence** of bisection while achieving the **speed** of open methods in favorable cases. It is the default root-finding algorithm in many numerical libraries.

### 8.2 Python Implementation with SciPy

SciPy provides Brent's method through `scipy.optimize.brentq`:

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

**Function signature:** `brentq(f, xl, xu)`
- `f` : function name for equation to be solved
- `xl` : lower (left) bracket value
- `xu` : upper (right) bracket value
- Requires $f(x_l) \cdot f(x_u) < 0$ (sign change, like bisection)

> **Practical Note:** In real-world applications, Brent's method is almost always preferred over pure bisection, Newton-Raphson, or secant methods because of its combination of reliability and speed.

<br>

---

<br>

## 9. Summary Table

| Method | Initial Guesses | Derivative Required | Convergence Order | Convergence Guarantee | Key Formula |
|---|---|---|---|---|---|
| **Fixed-Point Iteration** | 1 | No | Linear ($\sim q^i$) | Only if $\|g'(x^*)\| < 1$ | $x_{i+1} = g(x_i)$ |
| **Wegstein** | 2 | No | Superlinear | No | $x_{i+1} = \frac{x_i g(x_{i-1}) - x_{i-1} g(x_i)}{x_i - x_{i-1} - g(x_i) + g(x_{i-1})}$ |
| **Newton-Raphson** | 1 | Yes ($f'$) | Quadratic ($\sim e_i^2$) | No | $x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}$ |
| **Secant** | 2 | No | Superlinear ($\approx 1.618$) | No | $x_{i+1} = x_i - \frac{f(x_i)(x_i - x_{i-1})}{f(x_i) - f(x_{i-1})}$ |
| **Modified Secant** | 1 | No | Superlinear | No | $x_{i+1} = x_i - \frac{f(x_i)\,\delta x_i}{f(x_i + \delta x_i) - f(x_i)}$ |
| **Inverse Quadratic Interp.** | 3 | No | $\approx 1.839$ | No | Lagrange inverse polynomial at $y = 0$ |
| **Brent's Method** | 2 (bracket) | No | Superlinear (hybrid) | Yes (via bisection fallback) | Adaptive: bisection + secant + IQI |
