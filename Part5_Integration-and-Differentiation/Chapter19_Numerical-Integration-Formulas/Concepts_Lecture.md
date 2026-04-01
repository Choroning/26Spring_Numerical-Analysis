# Chapter 19 Lecture — Numerical Integration Formulas

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 19

> **Prerequisites**: [Calculus] Integration (Ch 1-18).
>
> **Learning Objectives**:
> 1. Apply trapezoidal rule and Simpson's rules
> 2. Derive error estimates for numerical integration
> 3. Implement Romberg integration and Gauss quadrature

---

<br>

## Table of Contents

- [1. Introduction to Numerical Integration](#1-introduction-to-numerical-integration)
  - [1.1 Motivation — Why Numerical Integration?](#11-motivation--why-numerical-integration)
  - [1.2 The General Idea: Polynomial Approximation](#12-the-general-idea-polynomial-approximation)
  - [1.3 Lagrange Polynomial Review](#13-lagrange-polynomial-review)
  - [1.4 Quadrature Formulas and Newton-Cotes](#14-quadrature-formulas-and-newton-cotes)
- [2. Trapezoidal Rule (n = 1)](#2-trapezoidal-rule-n--1)
  - [2.1 Derivation from Linear Interpolation](#21-derivation-from-linear-interpolation)
  - [2.2 Final Formula](#22-final-formula)
  - [2.3 Truncation Error](#23-truncation-error)
- [3. Simpson's 1/3 Rule (n = 2)](#3-simpsons-13-rule-n--2)
  - [3.1 Derivation from Quadratic Interpolation](#31-derivation-from-quadratic-interpolation)
  - [3.2 Computing the Basis Integrals via Change of Variable](#32-computing-the-basis-integrals-via-change-of-variable)
  - [3.3 Final Formula](#33-final-formula)
  - [3.4 Truncation Error](#34-truncation-error)
- [4. Simpson's 3/8 Rule (n = 3)](#4-simpsons-38-rule-n--3)
  - [4.1 Derivation from Cubic Interpolation](#41-derivation-from-cubic-interpolation)
  - [4.2 Final Formula](#42-final-formula)
  - [4.3 Truncation Error](#43-truncation-error)
- [5. Composite Trapezoidal Rule](#5-composite-trapezoidal-rule)
  - [5.1 Idea: Subdivide the Interval](#51-idea-subdivide-the-interval)
  - [5.2 Derivation](#52-derivation)
  - [5.3 Final Formula](#53-final-formula)
  - [5.4 Error Analysis](#54-error-analysis)
- [6. Composite Simpson's 1/3 Rule](#6-composite-simpsons-13-rule)
  - [6.1 Derivation](#61-derivation)
  - [6.2 Final Formula](#62-final-formula)
  - [6.3 Error Analysis](#63-error-analysis)
- [7. Composite Simpson's 3/8 Rule](#7-composite-simpsons-38-rule)
  - [7.1 Error Analysis](#71-error-analysis)
- [8. Combining Rules — Practical Application](#8-combining-rules--practical-application)
- [9. Python Implementations](#9-python-implementations)
- [Summary](#summary)

---

<br>

## 1. Introduction to Numerical Integration

### 1.1 Motivation — Why Numerical Integration?

Numerical integration deals with the problem of approximating the definite integral of a function over an interval:

$$I = \int_a^b f(x)\,dx \approx I_n$$

Many real-world integrals **cannot be evaluated analytically**. A classic motivating example from the lecture involves a semiconductor wafer manufacturer.

**Example — Wafer Quality Control:**

A manufacturer produces 2-inch diameter wafers. The actual diameter is normally distributed with mean $\mu = 2$ inches and standard deviation $\sigma = 0.01$ inches. Specifications require the diameter to be between 1.985 and 2.02 inches. What is the probability that a wafer is acceptable?

The normal probability density function is:

$$p(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\!\left(\frac{-(x - \mu)^2}{2\sigma^2}\right)$$

The probability of an acceptable wafer is:

$$P(1.985 \leq x \leq 2.02) = \int_{1.985}^{2.02} p(x)\,dx = \frac{1}{\sigma\sqrt{2\pi}} \int_{1.985}^{2.02} \exp\!\left(\frac{-(x - \mu)^2}{2\sigma^2}\right) dx$$

This integral has **no closed-form antiderivative**, so we need numerical integration methods.

### 1.2 The General Idea: Polynomial Approximation

The strategy is to **replace the integrand with a polynomial** that is easy to integrate:

$$p(x) \approx p_n(x) = a_0 + a_1 x + \cdots + a_n x^n$$

Using Lagrange interpolation, this polynomial can be written as:

$$p_n(x) = \sum_{j=0}^{n} p_n(x_j)\,\ell_j(x)$$

where $\ell_j(x)$ are the Lagrange basis polynomials and $x_0, x_1, \ldots, x_n$ are the interpolation nodes.

### 1.3 Lagrange Polynomial Review

For the first-order case with two points $x_1, x_2$, the Lagrange basis functions are:

$$N_1(x) = \frac{x - x_2}{x_1 - x_2}, \qquad N_2(x) = \frac{x_1 - x}{x_1 - x_2}$$

In general, the $n$-th order Lagrange basis polynomial is:

$$\ell_j(x) = \prod_{\substack{i=0 \\ j \neq i}}^{n} \frac{x - x_i}{x_j - x_i}, \qquad j = 0, 1, \ldots, n$$

### 1.4 Quadrature Formulas and Newton-Cotes

If $p(x) = p_n(x) + e_n(x)$ where $e_n(x)$ is the interpolation error, then:

$$\int_a^b p(x)\,dx = \int_a^b p_n(x)\,dx + \int_a^b e_n(x)\,dx$$

The polynomial part can be expanded:

$$\int_a^b p_n(x)\,dx = \int_a^b \sum_{j=0}^{n} p_n(x_j)\,\ell_j(x)\,dx = \sum_{j=0}^{n} p_n(x_j) \underbrace{\int_a^b \ell_j(x)\,dx}_{\text{a number (weight)}}$$

> **[Calculus]** The key insight is that $\int_a^b \ell_j(x)\,dx$ does **not** depend on $p_n(x)$ — it is a fixed numerical weight determined only by the node positions and the interval $[a,b]$. This is why these are called **quadrature formulas**.

The approach of approximating $I = \int_a^b f(x)\,dx$ by $I_n = \int_a^b f_n(x)\,dx$ where $f_n(x)$ is an $n$-th order polynomial is the foundation of **Newton-Cotes formulas**. The choice of $n$ determines the specific rule:

| $n$ | Rule Name | Polynomial Order |
|-----|-----------|-----------------|
| 1 | Trapezoidal Rule | Linear |
| 2 | Simpson's 1/3 Rule | Quadratic |
| 3 | Simpson's 3/8 Rule | Cubic |

---

<br>

## 2. Trapezoidal Rule (n = 1)

### 2.1 Derivation from Linear Interpolation

Approximate the integrand $f(x)$ by the first-order polynomial $p_1(x)$ that passes through the endpoints $(a, f(a))$ and $(b, f(b))$:

$$p_1(x) = f(a)\,N_1(x) + f(b)\,N_2(x) = f(a)\!\left(\frac{x - b}{a - b}\right) + f(b)\!\left(\frac{x - a}{b - a}\right)$$

Simplifying:

$$p_1(x) = \frac{f(b) - f(a)}{b - a}\,x + \frac{b\,f(a) - a\,f(b)}{b - a}$$

Now integrate:

$$\int_a^b f(x)\,dx \approx \int_a^b p_1(x)\,dx = \frac{f(b) - f(a)}{b - a} \int_a^b x\,dx + \frac{b\,f(a) - a\,f(b)}{b - a} \int_a^b 1\,dx$$

Evaluating the elementary integrals:

$$\int_a^b x\,dx = \frac{1}{2}(b^2 - a^2), \qquad \int_a^b 1\,dx = b - a$$

After algebraic simplification:

$$= \frac{f(b) - f(a)}{2}(a + b) + b\,f(a) - a\,f(b)$$

$$= \frac{1}{2}\bigl[a\,f(b) - a\,f(a) + b\,f(b) - b\,f(a)\bigr]$$

$$= \frac{1}{2}(b - a)\,f(b) + \frac{1}{2}(b - a)\,f(a)$$

### 2.2 Final Formula

$$\boxed{I \approx I_1 = \frac{b - a}{2}\bigl[f(a) + f(b)\bigr] = (b - a) \cdot \frac{f(a) + f(b)}{2}}$$

> **[Calculus]** The trapezoidal rule computes the area of a **trapezoid** formed by the function values at the endpoints. The factor $(b - a)$ is the **width** and $\frac{f(a) + f(b)}{2}$ is the **average height**.

### 2.3 Truncation Error

The error of the trapezoidal rule is:

$$\boxed{E_1 = I - I_1 = -\frac{1}{12}\,f''(\xi)\,(b - a)^3, \qquad a \leq \xi \leq b}$$

**Important consequence:** If $f(x) = mx + n$ (a linear function), then $f''(x) = 0$, so $E_1 = 0$. The **trapezoidal rule is exact for linear functions** (polynomials of degree $\leq 1$).

> **[Calculus]** The error is proportional to $(b-a)^3$ and depends on the second derivative $f''$. This tells us: (1) narrower intervals give smaller errors, and (2) functions with small curvature are well approximated by this rule.

---

<br>

## 3. Simpson's 1/3 Rule (n = 2)

### 3.1 Derivation from Quadratic Interpolation

Approximate $f(x)$ by the second-order polynomial $p_2(x)$ passing through three equally spaced points: $x_1 = a$, $x_2 = \frac{a+b}{2}$, $x_3 = b$, with step size $h = \frac{b - a}{2}$.

Define:

$$f_1 = f(a), \qquad f_2 = f\!\left(\frac{a+b}{2}\right), \qquad f_3 = f(b)$$

The interpolating polynomial is:

$$p_2(x) = f_1\,\ell_1(x) + f_2\,\ell_2(x) + f_3\,\ell_3(x)$$

where the Lagrange basis polynomials are:

$$\ell_1(x) = \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)}, \quad \ell_2(x) = \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)}, \quad \ell_3(x) = \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)}$$

So the integral becomes:

$$\int_a^b p_2(x)\,dx = f_1 \int_a^b \ell_1(x)\,dx + f_2 \int_a^b \ell_2(x)\,dx + f_3 \int_a^b \ell_3(x)\,dx$$

### 3.2 Computing the Basis Integrals via Change of Variable

To evaluate $\int_a^b \ell_j(x)\,dx$, apply the change of variable:

$$x = x(\zeta) = \frac{b - a}{2}\,\zeta + \frac{b + a}{2}, \qquad dx = \frac{b - a}{2}\,d\zeta$$

This maps $[a, b] \to [-1, 1]$, and the three nodes $x_1 = a,\; x_2 = \frac{a+b}{2},\; x_3 = b$ map to $\zeta = -1, 0, 1$ respectively.

**Basis functions in $\zeta$-space:**

$$\ell_1(\zeta) = \frac{(\zeta - 0)(\zeta - 1)}{(-1 - 0)(-1 - 1)} = \frac{\zeta(\zeta - 1)}{2}$$

$$\ell_2(\zeta) = \frac{(\zeta - (-1))(\zeta - 1)}{(0 - (-1))(0 - 1)} = (1 - \zeta)(1 + \zeta) = 1 - \zeta^2$$

$$\ell_3(\zeta) = \frac{(\zeta - (-1))(\zeta - 0)}{(1 - (-1))(1 - 0)} = \frac{\zeta(1 + \zeta)}{2}$$

**Computing the integrals:**

For $\ell_1$:

$$\int_a^b \ell_1(x)\,dx = \int_{-1}^{1} \frac{\zeta(\zeta - 1)}{2} \cdot \frac{b - a}{2}\,d\zeta = \frac{b - a}{4} \int_{-1}^{1} (\zeta^2 - \zeta)\,d\zeta$$

Since $\zeta$ is an **odd function** on $[-1, 1]$, its integral vanishes:

$$= \frac{b - a}{4} \cdot \frac{2}{3} = \frac{b - a}{6}$$

By symmetry, $\int_a^b \ell_3(x)\,dx = \frac{b - a}{6}$ as well.

For $\ell_2$:

$$\int_a^b \ell_2(x)\,dx = \int_{-1}^{1} (1 - \zeta^2) \cdot \frac{b - a}{2}\,d\zeta = \frac{b - a}{2}\left(\zeta - \frac{\zeta^3}{3}\right)\Bigg|_{-1}^{1} = \frac{b - a}{2} \cdot \frac{4}{3} = \frac{4(b - a)}{6}$$

### 3.3 Final Formula

Combining the weights:

$$\int_a^b f(x)\,dx \approx \frac{b - a}{6}\bigl(f_1 + 4f_2 + f_3\bigr)$$

With $h = \frac{b - a}{2}$:

$$\boxed{I \approx I_2 = \frac{h}{3}(f_1 + 4f_2 + f_3) = \frac{b - a}{6}\bigl[f(a) + 4f\!\left(\tfrac{a+b}{2}\right) + f(b)\bigr]}$$

> **[Calculus]** The name "1/3 rule" comes from the factor $\frac{h}{3}$ in the formula. The middle point receives **four times** the weight of the endpoints, reflecting the parabola's better fit near the center of the interval.

### 3.4 Truncation Error

The error for Simpson's 1/3 rule is:

$$E_2 = I - I_2 = -\frac{1}{90}\,h^5\,f^{(4)}(\xi), \qquad a \leq \xi \leq b$$

Since $h = \frac{b - a}{2}$:

$$\boxed{E_2 = -\frac{(b - a)^5}{2880}\,f^{(4)}(\xi)}$$

**Key insight:** Although Simpson's 1/3 rule uses a **second-order** polynomial, the error depends on the **fourth derivative** $f^{(4)}$. This means Simpson's 1/3 rule is **exact for polynomials up to degree 3** (since $f^{(4)} = 0$ for cubics).

> **[Calculus]** This "bonus order" of accuracy occurs because the odd-power error terms cancel out due to the symmetry of the integration nodes around the midpoint.

---

<br>

## 4. Simpson's 3/8 Rule (n = 3)

### 4.1 Derivation from Cubic Interpolation

Approximate $f(x)$ by a third-order polynomial $p_3(x)$ through four equally spaced points with step size $h = \frac{b - a}{3}$:

$$x_0 = a, \quad x_1 = a + h, \quad x_2 = a + 2h, \quad x_3 = b$$

The interpolating polynomial is:

$$p_3(x) = f_1\,\ell_1(x) + f_2\,\ell_2(x) + f_3\,\ell_3(x) + f_4\,\ell_4(x)$$

Integrating:

$$\int_a^b p_3(x)\,dx = f_1 \int_a^b \ell_1(x)\,dx + f_2 \int_a^b \ell_2(x)\,dx + f_3 \int_a^b \ell_3(x)\,dx + f_4 \int_a^b \ell_4(x)\,dx$$

### 4.2 Final Formula

$$\boxed{I \approx I_3 = \frac{3h}{8}(f_0 + 3f_1 + 3f_2 + f_3), \qquad h = \frac{b - a}{3}}$$

> **[Calculus]** The name "3/8 rule" comes from the coefficient $\frac{3}{8}$ multiplied by $h$. The weights follow the pattern $1 : 3 : 3 : 1$, reflecting the symmetric cubic interpolation.

### 4.3 Truncation Error

$$E_3 = I - I_3 = -\frac{3}{80}\,h^5\,f^{(4)}(\xi), \qquad a \leq \xi \leq b$$

Since $h = \frac{b - a}{3}$:

$$\boxed{E_3 = -\frac{(b - a)^5}{6480}\,f^{(4)}(\xi)}$$

Like Simpson's 1/3 rule, the 3/8 rule is also **exact for polynomials up to degree 3**, since the error depends on $f^{(4)}$.

> **[Calculus]** Comparing single-application errors: $|E_3| = \frac{(b-a)^5}{6480}|f^{(4)}(\xi)|$ vs. $|E_2| = \frac{(b-a)^5}{2880}|f^{(4)}(\xi)|$. The 3/8 rule is **more accurate** per single application ($6480 > 2880$), but requires one additional function evaluation.

---

<br>

## 5. Composite Trapezoidal Rule

### 5.1 Idea: Subdivide the Interval

**Q: How can we reduce the error of the trapezoidal rule?**

**A:** Divide the interval $[a, b]$ into $n$ subintervals (segments) and apply the trapezoidal rule to each segment. Smaller subintervals mean the linear approximation is more accurate on each piece.

With $n$ segments:

$$h = \frac{b - a}{n}, \qquad x_i = a + ih \;\;(i = 0, 1, \ldots, n)$$

where $x_0 = a$ and $x_n = b$, and $f_i = f(x_i)$.

### 5.2 Derivation

Split the integral into $n$ subintegrals and apply the trapezoidal rule to each:

$$I = \int_a^b f(x)\,dx = \int_{x_0}^{x_1} f(x)\,dx + \int_{x_1}^{x_2} f(x)\,dx + \cdots + \int_{x_{n-1}}^{x_n} f(x)\,dx$$

$$\approx (x_1 - x_0)\frac{f_0 + f_1}{2} + (x_2 - x_1)\frac{f_1 + f_2}{2} + \cdots + (x_n - x_{n-1})\frac{f_{n-1} + f_n}{2}$$

Since all subintervals have equal width $h$:

$$= \frac{h}{2}\bigl(f_0 + f_1 + f_1 + f_2 + f_2 + \cdots + f_{n-1} + f_{n-1} + f_n\bigr)$$

Collecting terms (interior points appear twice):

$$= \frac{h}{2}\left(f_0 + 2f_1 + 2f_2 + \cdots + 2f_{n-1} + f_n\right)$$

### 5.3 Final Formula

$$\boxed{I \approx \frac{h}{2}\left(f_0 + 2\sum_{i=1}^{n-1} f_i + f_n\right) = (b - a)\,\frac{f_0 + 2\displaystyle\sum_{i=1}^{n-1} f_i + f_n}{2n}}$$

This is the **Composite Trapezoidal Rule**. It has the same interpretation: **width** $\times$ **average height**, where the average height uses a weighted average with endpoints counted once and interior points counted twice.

### 5.4 Error Analysis

The local error on each subinterval $[x_{i-1}, x_i]$ is:

$$-\frac{1}{12}\,f''(\xi_i)\,h^3, \qquad x_{i-1} \leq \xi_i \leq x_i$$

Summing all $n$ local errors:

$$E_{1t} = -\frac{1}{12}\,h^3 \sum_{i=1}^{n} f''(\xi_i) = -\frac{1}{12}\,h^3 \cdot n \cdot \frac{\displaystyle\sum_{i=1}^{n} f''(\xi_i)}{n}$$

The average of $f''(\xi_i)$ values approximates $\overline{f''}$, and $n = \frac{b - a}{h}$:

$$\boxed{E_{1a} \approx -\frac{(b - a)^3}{12\,n^2}\,\overline{f''} = -\frac{(b - a)}{12}\,\overline{f''}\,h^2}$$

The error is **proportional to $h^2$**. Therefore:

> **If $n$ is doubled, the truncation error is reduced to one-fourth!** (Since $h^2 \to (h/2)^2 = h^2/4$.)

---

<br>

## 6. Composite Simpson's 1/3 Rule

### 6.1 Derivation

Divide $[a, b]$ into $n$ subintervals (**$n$ must be even**) with $h = \frac{b - a}{n}$. Apply Simpson's 1/3 rule to each consecutive pair of subintervals:

$$I = \int_{x_0}^{x_2} f(x)\,dx + \int_{x_2}^{x_4} f(x)\,dx + \cdots + \int_{x_{n-2}}^{x_n} f(x)\,dx$$

$$\approx \frac{h}{3}(f_0 + 4f_1 + f_2) + \frac{h}{3}(f_2 + 4f_3 + f_4) + \frac{h}{3}(f_4 + 4f_5 + f_6) + \cdots + \frac{h}{3}(f_{n-2} + 4f_{n-1} + f_n)$$

Collecting terms:

$$= \frac{h}{3}\,f_0 + \frac{4h}{3}(f_1 + f_3 + f_5 + \cdots + f_{n-1}) + \frac{2h}{3}(f_2 + f_4 + f_6 + \cdots + f_{n-2}) + \frac{h}{3}\,f_n$$

### 6.2 Final Formula

$$\boxed{I \approx \frac{h}{3}\left(f_0 + 4\!\!\sum_{\substack{i=1,3,5,\ldots}}^{n-1}\!\! f_i + 2\!\!\sum_{\substack{i=2,4,6,\ldots}}^{n-2}\!\! f_i + f_n\right)}$$

Equivalently:

$$= (b - a)\,\frac{f_0 + 4\displaystyle\sum_{\substack{i\,\text{odd}}} f_i + 2\displaystyle\sum_{\substack{i\,\text{even}}} f_i + f_n}{3n}$$

> **[Calculus]** The pattern of weights is $1, 4, 2, 4, 2, 4, \ldots, 2, 4, 1$. Odd-indexed interior points get weight 4, even-indexed interior points get weight 2, and endpoints get weight 1. This requires $n$ to be **even**.

### 6.3 Error Analysis

The local error for each pair of subintervals is:

$$-\frac{h^5}{90}\,f^{(4)}(\xi_k), \qquad x_{2k-2} \leq \xi_k \leq x_{2k}$$

Summing over all $n/2$ applications:

$$E_{2t} = -\frac{h^5}{90} \sum_{k=1}^{n/2} f^{(4)}(\xi_k) = -\frac{h^5}{90} \cdot \frac{n}{2} \cdot \overline{f^{(4)}}$$

Since $h = \frac{b-a}{n}$:

$$\boxed{E_{2a} \approx -\frac{(b - a)^5}{180\,n^4}\,\overline{f^{(4)}} = -\frac{(b - a)}{180}\,\overline{f^{(4)}}\,h^4}$$

The error is **proportional to $h^4$**. Therefore:

> **If $n$ is doubled, the truncation error is reduced by a factor of 16!** (Since $h^4 \to (h/2)^4 = h^4/16$.)

---

<br>

## 7. Composite Simpson's 3/8 Rule

Apply Simpson's 3/8 rule repeatedly over groups of 3 subintervals (**$n$ must be a multiple of 3**). The derivation follows the same pattern as the composite 1/3 rule.

### 7.1 Error Analysis

The local error for each group of 3 subintervals is:

$$-\frac{3}{80}\,h^5\,f^{(4)}(\xi_k)$$

Summing over all $n/3$ applications:

$$E_{3t} = -\frac{3}{80}\,h^5 \sum_{k=1}^{n/3} f^{(4)}(\xi_k) \approx -\frac{1}{80}\,h^5 \cdot n \cdot \overline{f^{(4)}}$$

$$\boxed{E_{3a} \approx -\frac{(b - a)}{80}\,\overline{f^{(4)}}\,h^4}$$

Like the composite 1/3 rule, the error is proportional to $h^4$, but the constant is slightly different.

> **[Calculus]** Although the 3/8 rule is more accurate than the 1/3 rule for a **single application**, the composite 3/8 rule has a larger error constant ($\frac{1}{80}$ vs. $\frac{1}{180}$) per unit of $h^4$ when applied in composite form. This is because the 3/8 rule uses more subintervals per application.

---

<br>

## 8. Combining Rules — Practical Application

In practice, when the number of segments $n$ is **not even** (required for composite Simpson's 1/3) and **not a multiple of 3** (required for composite Simpson's 3/8), we can **combine** the two rules:

- Apply Simpson's **1/3 rule** to a portion of the subintervals (uses pairs of segments)
- Apply Simpson's **3/8 rule** to the remaining subintervals (uses triples of segments)

For example, with 5 segments: use the 1/3 rule on the first 2 segments and the 3/8 rule on the last 3 segments (or vice versa).

---

<br>

## 9. Python Implementations

### Trapezoidal Rule (Single Application)

```python
def trapezoidal(f, a, b):
    """Single-application trapezoidal rule."""
    return (b - a) / 2 * (f(a) + f(b))
```

### Composite Trapezoidal Rule

```python
def composite_trapezoidal(f, a, b, n):
    """Composite trapezoidal rule with n subintervals."""
    h = (b - a) / n
    x = [a + i * h for i in range(n + 1)]
    result = f(x[0]) + f(x[n])
    for i in range(1, n):
        result += 2 * f(x[i])
    return h / 2 * result
```

### Simpson's 1/3 Rule (Single Application)

```python
def simpson_13(f, a, b):
    """Single-application Simpson's 1/3 rule."""
    h = (b - a) / 2
    return h / 3 * (f(a) + 4 * f(a + h) + f(b))
```

### Composite Simpson's 1/3 Rule

```python
def composite_simpson_13(f, a, b, n):
    """Composite Simpson's 1/3 rule with n subintervals (n must be even)."""
    if n % 2 != 0:
        raise ValueError("n must be even for Simpson's 1/3 rule")
    h = (b - a) / n
    x = [a + i * h for i in range(n + 1)]
    result = f(x[0]) + f(x[n])
    for i in range(1, n, 2):      # odd indices: weight 4
        result += 4 * f(x[i])
    for i in range(2, n - 1, 2):  # even indices: weight 2
        result += 2 * f(x[i])
    return h / 3 * result
```

### Simpson's 3/8 Rule (Single Application)

```python
def simpson_38(f, a, b):
    """Single-application Simpson's 3/8 rule."""
    h = (b - a) / 3
    return 3 * h / 8 * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(b))
```

### Complete Example — Comparing All Methods

```python
import math

def f(x):
    """Example: integrate e^x from 0 to 1 (exact answer = e - 1)."""
    return math.exp(x)

a, b = 0, 1
exact = math.e - 1

print(f"Exact value: {exact:.10f}")
print(f"Trapezoidal (n=1):       {trapezoidal(f, a, b):.10f}")
print(f"Simpson's 1/3 (n=2):     {simpson_13(f, a, b):.10f}")
print(f"Simpson's 3/8 (n=3):     {simpson_38(f, a, b):.10f}")
print(f"Composite Trap (n=10):   {composite_trapezoidal(f, a, b, 10):.10f}")
print(f"Composite S-1/3 (n=10):  {composite_simpson_13(f, a, b, 10):.10f}")
```

---

<br>

## Summary

| Rule | Formula | Points | Error (Single) | Error (Composite) | Exact For |
|------|---------|--------|----------------|-------------------|-----------|
| **Trapezoidal** | $(b-a)\,\dfrac{f(a)+f(b)}{2}$ | 2 | $O\bigl((b-a)^3\bigr)$ | $O(h^2)$ | Degree $\leq 1$ |
| **Simpson's 1/3** | $\dfrac{h}{3}(f_0 + 4f_1 + f_2)$ | 3 | $O\bigl((b-a)^5\bigr)$ | $O(h^4)$ | Degree $\leq 3$ |
| **Simpson's 3/8** | $\dfrac{3h}{8}(f_0 + 3f_1 + 3f_2 + f_3)$ | 4 | $O\bigl((b-a)^5\bigr)$ | $O(h^4)$ | Degree $\leq 3$ |

**Key takeaways from the lecture:**

1. All Newton-Cotes formulas work by **replacing the integrand with a polynomial** and integrating the polynomial exactly
2. The **trapezoidal rule** uses a 1st-order polynomial; its composite error decreases as $O(h^2)$ --- doubling $n$ reduces error by $4\times$
3. **Simpson's 1/3 rule** uses a 2nd-order polynomial but is exact up to degree 3 due to symmetry; its composite error decreases as $O(h^4)$ --- doubling $n$ reduces error by $16\times$
4. **Simpson's 3/8 rule** uses a 3rd-order polynomial; also exact up to degree 3 with $O(h^4)$ composite error; slightly more accurate per single application than the 1/3 rule
5. The **composite versions** subdivide the interval to reduce error, and different rules can be **combined** when $n$ does not match the divisibility requirements of a single rule
