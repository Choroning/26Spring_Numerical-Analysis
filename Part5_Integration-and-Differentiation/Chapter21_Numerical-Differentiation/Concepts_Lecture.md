# Chapter 21 Lecture — Numerical Differentiation

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 21

> **Prerequisites**: [Calculus] Derivatives, Taylor series (Ch 1-20).
>
> **Learning Objectives**:
> 1. Apply finite difference formulas for differentiation
> 2. Analyze truncation error in numerical derivatives
> 3. Implement Richardson extrapolation for improved accuracy

---

<br>

## Table of Contents

- [1. Introduction to Numerical Differentiation](#1-introduction-to-numerical-differentiation)
  - [1.1 The Derivative as a Limit](#11-the-derivative-as-a-limit)
  - [1.2 Second Derivative](#12-second-derivative)
  - [1.3 Partial Derivatives](#13-partial-derivatives)
  - [1.4 Physical Application: Heat Flux](#14-physical-application-heat-flux)
- [2. Finite Difference Approximation for the First Derivative](#2-finite-difference-approximation-for-the-first-derivative)
  - [2.1 Forward Difference (First-Order)](#21-forward-difference-first-order)
  - [2.2 Backward Difference (First-Order)](#22-backward-difference-first-order)
  - [2.3 Central Difference (Second-Order)](#23-central-difference-second-order)
- [3. Finite Difference Approximation for the Second Derivative](#3-finite-difference-approximation-for-the-second-derivative)
  - [3.1 Central Difference via Taylor Series](#31-central-difference-via-taylor-series)
  - [3.2 Forward Difference for the Second Derivative](#32-forward-difference-for-the-second-derivative)
- [4. Improving Accuracy: Higher-Order Forward Difference](#4-improving-accuracy-higher-order-forward-difference)
- [5. Richardson Extrapolation for Differentiation](#5-richardson-extrapolation-for-differentiation)
  - [5.1 Basic Idea](#51-basic-idea)
  - [5.2 Revisiting the Central Difference Error Structure](#52-revisiting-the-central-difference-error-structure)
  - [5.3 First Extrapolation Level](#53-first-extrapolation-level)
  - [5.4 Second Extrapolation Level](#54-second-extrapolation-level)
- [6. Derivatives of Unequally Spaced Data](#6-derivatives-of-unequally-spaced-data)
  - [6.1 Lagrange Polynomial Approach](#61-lagrange-polynomial-approach)
  - [6.2 Derivatives of the Lagrange Basis Functions](#62-derivatives-of-the-lagrange-basis-functions)
- [7. Numerical Partial Derivatives](#7-numerical-partial-derivatives)
  - [7.1 First-Order Partial Derivatives (Central Difference)](#71-first-order-partial-derivatives-central-difference)
  - [7.2 Mixed Partial Derivative](#72-mixed-partial-derivative)
- [8. Method of Undetermined Coefficients](#8-method-of-undetermined-coefficients)
  - [8.1 Problem Setup](#81-problem-setup)
  - [8.2 Taylor Expansion and Matching](#82-taylor-expansion-and-matching)
  - [8.3 Solving for Coefficients](#83-solving-for-coefficients)
  - [8.4 Truncation Error](#84-truncation-error)
- [9. Differentiable Programming (Brief Overview)](#9-differentiable-programming-brief-overview)
- [10. Python Implementation](#10-python-implementation)
- [Summary](#summary)

---

<br>

## 1. Introduction to Numerical Differentiation

Numerical differentiation is the process of estimating derivatives from discrete data or function evaluations, rather than computing them analytically. This is essential when:

- The function is known only at discrete points (experimental/field data)
- The analytical derivative is too complex or impossible to obtain
- Derivatives are needed within a larger numerical algorithm (e.g., solving ODEs/PDEs)

### 1.1 The Derivative as a Limit

The finite difference ratio provides a natural starting point:

$$\frac{\Delta y}{\Delta x} = \frac{f(x + \Delta x) - f(x)}{\Delta x}$$

Taking the limit as $\Delta x \to 0$:

$$\frac{dy}{dx} = \lim_{\Delta x \to 0} \frac{f(x + \Delta x) - f(x)}{\Delta x} =: f'(x)$$

This is **the first derivative of $y$ with respect to $x$**. Geometrically, the derivative is the slope of the tangent line to the curve at point $x_i$.

> **[Calculus]** The three diagrams from the lecture illustrate how as $\Delta x$ shrinks, the secant line (connecting two points on the curve) approaches the tangent line at $x_i$, and the slope of the secant converges to $f'(x_i)$.

### 1.2 Second Derivative

The second derivative is the derivative of the first derivative:

$$\frac{d^2 y}{dx^2} = \frac{d}{dx}\left(\frac{dy}{dx}\right)$$

It measures **how fast the slope is changing** -- i.e., the curvature or concavity of the function.

### 1.3 Partial Derivatives

For a function of two variables $f(x, y)$:

$$\frac{\partial f}{\partial x} = \lim_{\Delta x \to 0} \frac{f(x + \Delta x, y) - f(x, y)}{\Delta x}$$

$$\frac{\partial f}{\partial y} = \lim_{\Delta y \to 0} \frac{f(x, y + \Delta y) - f(x, y)}{\Delta y}$$

These are the **partial derivatives of $f$** with respect to $x$ and $y$, respectively. Each partial derivative measures the rate of change of $f$ while holding the other variable fixed.

### 1.4 Physical Application: Heat Flux

A classic engineering application is Fourier's law of heat conduction. Given a temperature distribution $T(x)$, the heat flux is:

$$\vec{q} = -k \frac{\partial T}{\partial x}$$

where $k$ is the thermal conductivity. The negative sign indicates that heat flows from high temperature to low temperature:

- A **positive slope** ($\partial T / \partial x > 0$) drives a **negative heat flow** (heat flows in the $-x$ direction)
- A **negative slope** ($\partial T / \partial x < 0$) drives a **positive heat flow** (heat flows in the $+x$ direction)

> **[Calculus]** The direction of heat flow is always opposite to the temperature gradient. This is a direct physical application of the derivative concept.

---

<br>

## 2. Finite Difference Approximation for the First Derivative

All finite difference formulas are derived from Taylor series expansions. Let $h$ denote the step size between equally spaced points: $x_{i+1} = x_i + h$.

### 2.1 Forward Difference (First-Order)

**Taylor expansion** about $x_i$ in the forward direction:

$$f(x_{i+1}) = f(x_i) + f'(x_i)\,h + f''(x_i)\frac{h^2}{2} + f'''(x_i)\frac{h^3}{6} + \cdots \quad (\ast)$$

Solving equation $(\ast)$ for $f'(x_i)$:

$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - f''(x_i)\frac{h}{2} + O(h^2)$$

Dropping the truncation error terms yields the **forward difference** approximation:

$$\boxed{f'(x_i) \approx \frac{f(x_{i+1}) - f(x_i)}{h}}$$

- **Truncation error**: $O(h)$ -- first-order accurate
- **Leading error term**: $-f''(x_i)\dfrac{h}{2}$
- Uses points $x_i$ and $x_{i+1}$ (the point itself and one point ahead)

### 2.2 Backward Difference (First-Order)

Similarly, expanding $f(x_{i-1})$ about $x_i$:

$$f(x_{i-1}) = f(x_i) - f'(x_i)\,h + f''(x_i)\frac{h^2}{2} - f'''(x_i)\frac{h^3}{6} + \cdots \quad (\ast\ast)$$

Solving for $f'(x_i)$:

$$\boxed{f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{h}}$$

- **Truncation error**: $O(h)$ -- first-order accurate
- Uses points $x_{i-1}$ and $x_i$

### 2.3 Central Difference (Second-Order)

Subtracting $(\ast\ast)$ from $(\ast)$:

$$f(x_{i+1}) - f(x_{i-1}) = 2f'(x_i)\,h + 2f'''(x_i)\frac{h^3}{6} + \cdots$$

Solving for $f'(x_i)$:

$$\boxed{f'(x_i) \approx \frac{f(x_{i+1}) - f(x_{i-1})}{2h}}$$

- **Truncation error**: $O(h^2)$ -- second-order accurate
- **Leading error term**: $-f'''(x_i)\dfrac{h^2}{6}$
- Uses points $x_{i-1}$ and $x_{i+1}$ (symmetric about $x_i$)

> **[Calculus]** The central difference achieves higher accuracy ($O(h^2)$ vs. $O(h)$) compared to forward/backward differences because the odd-order error terms cancel due to symmetry. This is a key advantage whenever both neighboring points are available.

---

<br>

## 3. Finite Difference Approximation for the Second Derivative

### 3.1 Central Difference via Taylor Series

Adding $(\ast)$ and $(\ast\ast)$ together:

$$f(x_{i+1}) + f(x_{i-1}) = 2f(x_i) + f''(x_i)\,h^2 + O(h^4)$$

Solving for $f''(x_i)$:

$$\boxed{f''(x_i) = \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2} + O(h^2)}$$

This is the **central difference approximation for the second derivative**:

- **Truncation error**: $O(h^2)$ -- second-order accurate
- Uses three equally spaced points: $x_{i-1}$, $x_i$, $x_{i+1}$

> **[Calculus]** This formula has a beautiful physical interpretation: it measures how much the function value at $x_i$ deviates from the average of its two neighbors. A large positive $f''$ means the function curves upward (concave up), while a large negative $f''$ means it curves downward (concave down).

### 3.2 Forward Difference for the Second Derivative

To derive a **forward** formula for $f''(x_i)$ using only forward points, expand $f(x_{i+2})$:

$$f(x_{i+2}) = f(x_i) + f'(x_i)\,(2h) + f''(x_i)\frac{(2h)^2}{2} + f'''(x_i)\frac{(2h)^3}{6} + O(h^4) \quad (\ast\ast\ast)$$

Computing $(\ast\ast\ast) - 2(\ast)$:

$$f(x_{i+2}) - 2f(x_{i+1}) = -f(x_i) + f''(x_i)\,h^2 + O(h^3)$$

Solving for $f''(x_i)$:

$$\boxed{f''(x_i) = \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{h^2} + O(h)}$$

- **Truncation error**: $O(h)$ -- first-order accurate (less accurate than the central formula)
- Uses three forward points: $x_i$, $x_{i+1}$, $x_{i+2}$

---

<br>

## 4. Improving Accuracy: Higher-Order Forward Difference

We can improve the forward difference for $f'(x_i)$ by incorporating the second derivative correction term. Recall:

$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \underbrace{f''(x_i)}_{\text{plug in FD approx.}}\frac{h}{2} + O(h^2)$$

Substituting the forward difference approximation for $f''(x_i)$:

$$f''(x_i) \approx \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{h^2}$$

into the error term:

$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{2h} + O(h^2)$$

Combining:

$$\boxed{f'(x_i) = \frac{-f(x_{i+2}) + 4f(x_{i+1}) - 3f(x_i)}{2h} + O(h^2)}$$

This is the **second-order forward difference** formula for $f'(x_i)$.

- **Accuracy improved** from $O(h)$ to $O(h^2)$ by including one additional point
- Adding the diffusion (second derivative) term enhances the accuracy

> **[Calculus]** The idea of substituting a numerical approximation for the error's leading term back into the formula is a general technique for increasing accuracy. It is closely related to Richardson extrapolation (discussed next).

---

<br>

## 5. Richardson Extrapolation for Differentiation

### 5.1 Basic Idea

Recall from numerical integration that Richardson extrapolation uses two estimates computed with different step sizes to eliminate leading error terms. For integration:

$$I = I(h_2) + \frac{1}{(h_1/h_2)^2 - 1}\bigl(I(h_2) - I(h_1)\bigr)$$

When $h_1/h_2 = 2$:

$$I = \frac{4}{3}I(h_2) - \frac{1}{3}I(h_1)$$

Similarly, for differentiation:

$$\boxed{D = \frac{4}{3}D(h_2) - \frac{1}{3}D(h_1)}$$

where $D(h_1)$ and $D(h_2)$ are derivative estimates computed with step sizes $h_1$ (larger) and $h_2 = h_1/2$ (smaller).

### 5.2 Revisiting the Central Difference Error Structure

Start with the Taylor expansions:

$$f(x+h) = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}h^k$$

$$f(x-h) = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}(-h)^k = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}(-1)^k h^k$$

Computing $f(x+h) - f(x-h)$:

$$f(x+h) - f(x-h) = 2h\,f'(x) + 2f'''(x)\frac{h^3}{3!} + 2f^{(5)}(x)\frac{h^5}{5!} + \cdots$$

Dividing by $2h$:

$$f'(x) = \frac{f(x+h) - f(x-h)}{2h} - f'''(x)\frac{h^2}{3!} - f^{(5)}(x)\frac{h^4}{5!} - \cdots$$

Define the central difference operator and error coefficients:

$$D(h) = \frac{f(x+h) - f(x-h)}{2h}$$

$$e_2 = -\frac{f^{(3)}(x)}{3!}, \qquad e_4 = -\frac{f^{(5)}(x)}{5!}$$

Then:

$$f'(x) = D(h) + e_2\,h^2 + e_4\,h^4 + \cdots$$

where the error $E(h) = e_2\,h^2 + e_4\,h^4 + \cdots$.

> **[Calculus]** A crucial observation is that the error coefficients $e_2, e_4, \ldots$ are **independent of $h$** -- they depend only on the function and the evaluation point $x$. This is exactly the property that makes Richardson extrapolation possible.

### 5.3 First Extrapolation Level

**Goal**: Eliminate the $h^2$ term.

Write the approximation at step sizes $h$ and $2h$:

$$f'(x) = D(h) + e_2\,h^2 + e_4\,h^4 + \cdots \quad \text{(multiply by 4)}$$

$$f'(x) = D(2h) + e_2\,(2h)^2 + e_4\,(2h)^4 + \cdots$$

Compute $4 \times (\text{first}) - (\text{second})$:

$$3f'(x) = 4D(h) - D(2h) - 12h^4 e_4 + O(h^6)$$

$$\boxed{f'(x) = \underbrace{\frac{4}{3}D(h) - \frac{1}{3}D(2h)}_{=:\,\bar{D}(h)} - 4h^4 e_4 + O(h^6)}$$

The combined estimate $\bar{D}(h)$ is now **$O(h^4)$** accurate -- the $h^2$ error has been eliminated.

### 5.4 Second Extrapolation Level

Continue the same process to eliminate the $h^4$ term. At step sizes $h$ and $2h$:

$$f'(x) = \bar{D}(h) - 4h^4 e_4 + O(h^6)$$

$$f'(x) = \bar{D}(2h) - 4(2h)^4 e_4 + O(h^6)$$

Compute $16 \times (\text{first}) - (\text{second})$:

$$15f'(x) = 16\bar{D}(h) - \bar{D}(2h) + O(h^6)$$

$$\boxed{f'(x) = \frac{16}{15}\bar{D}(h) - \frac{1}{15}\bar{D}(2h) + O(h^6)}$$

This is now **$O(h^6)$** accurate. The process can be continued indefinitely, each level eliminating the next even power of $h$.

> **[Calculus]** Richardson extrapolation is a systematic way to accelerate convergence. At each level, the leading error term is cancelled by combining estimates at two different step sizes. The general pattern for the central difference (where errors are in even powers of $h$) uses the ratio $4^n$ at the $n$th level: the weights are $\frac{4^n}{4^n - 1}$ and $\frac{-1}{4^n - 1}$.

---

<br>

## 6. Derivatives of Unequally Spaced Data

In practice, data from experiments or field studies is often collected at **unequal intervals**. The finite difference formulas derived above assume equally spaced points and cannot be applied directly.

### 6.1 Lagrange Polynomial Approach

One way to handle non-equispaced data is to fit a **Lagrange polynomial** and then differentiate it.

Given three data points $(x_0, y_0)$, $(x_1, y_1)$, $(x_2, y_2)$, construct a 2nd-order Lagrange polynomial:

$$f(x) = y_0\,N_0(x) + y_1\,N_1(x) + y_2\,N_2(x)$$

where $N_k(x)$ are the Lagrange basis functions. Differentiating:

$$f'(x) = y_0\,N_0'(x) + y_1\,N_1'(x) + y_2\,N_2'(x)$$

### 6.2 Derivatives of the Lagrange Basis Functions

For three points, the basis functions and their derivatives are:

**Basis function $N_0(x)$**:

$$N_0(x) = \frac{(x - x_1)(x - x_2)}{(x_0 - x_1)(x_0 - x_2)} = \frac{x^2 - (x_1 + x_2)x + x_1 x_2}{(x_0 - x_1)(x_0 - x_2)}$$

$$N_0'(x) = \frac{2x - x_1 - x_2}{(x_0 - x_1)(x_0 - x_2)}$$

**Basis function $N_1(x)$**:

$$N_1(x) = \frac{(x - x_0)(x - x_2)}{(x_1 - x_0)(x_1 - x_2)} = \frac{x^2 - (x_0 + x_2)x + x_0 x_2}{(x_1 - x_0)(x_1 - x_2)}$$

$$N_1'(x) = \frac{2x - x_0 - x_2}{(x_1 - x_0)(x_1 - x_2)}$$

**Basis function $N_2(x)$**:

$$N_2(x) = \frac{(x - x_0)(x - x_1)}{(x_2 - x_0)(x_2 - x_1)} = \frac{x^2 - (x_0 + x_1)x + x_0 x_1}{(x_2 - x_0)(x_2 - x_1)}$$

$$N_2'(x) = \frac{2x - x_0 - x_1}{(x_2 - x_0)(x_2 - x_1)}$$

The derivative at any point $x$ is then:

$$f'(x) = y_0 \cdot \frac{2x - x_1 - x_2}{(x_0 - x_1)(x_0 - x_2)} + y_1 \cdot \frac{2x - x_0 - x_2}{(x_1 - x_0)(x_1 - x_2)} + y_2 \cdot \frac{2x - x_0 - x_1}{(x_2 - x_0)(x_2 - x_1)}$$

> **[Calculus]** This approach works for any spacing of the data points. When the points happen to be equally spaced (i.e., $x_1 - x_0 = x_2 - x_1 = h$), substituting $x = x_1$ into the formula above recovers the standard central difference formula $f'(x_1) \approx (y_2 - y_0)/(2h)$.

---

<br>

## 7. Numerical Partial Derivatives

### 7.1 First-Order Partial Derivatives (Central Difference)

For a function $f(x, y)$, the centered finite difference approximations for partial derivatives are:

$$\boxed{\frac{\partial f}{\partial x} \approx \frac{f(x + \Delta x,\, y) - f(x - \Delta x,\, y)}{2\,\Delta x}}$$

$$\boxed{\frac{\partial f}{\partial y} \approx \frac{f(x,\, y + \Delta y) - f(x,\, y - \Delta y)}{2\,\Delta y}}$$

Both are $O(\Delta x^2)$ and $O(\Delta y^2)$ accurate, respectively.

### 7.2 Mixed Partial Derivative

The mixed second partial derivative $\dfrac{\partial^2 f}{\partial x \,\partial y}$ is computed by applying the central difference twice:

$$\frac{\partial^2 f}{\partial x \,\partial y} = \frac{\partial}{\partial x}\left(\frac{\partial f}{\partial y}\right)$$

$$\approx \frac{1}{2\,\Delta x}\left(\left.\frac{\partial f}{\partial y}\right|_{x+\Delta x,\,y} - \left.\frac{\partial f}{\partial y}\right|_{x-\Delta x,\,y}\right)$$

Applying the central difference for $\partial f / \partial y$ at each point:

$$\approx \frac{1}{2\,\Delta x}\left(\frac{f(x+\Delta x,\, y+\Delta y) - f(x+\Delta x,\, y-\Delta y)}{2\,\Delta y} - \frac{f(x-\Delta x,\, y+\Delta y) - f(x-\Delta x,\, y-\Delta y)}{2\,\Delta y}\right)$$

Simplifying:

$$\boxed{\frac{\partial^2 f}{\partial x \,\partial y} \approx \frac{f(x+\Delta x,\, y+\Delta y) - f(x+\Delta x,\, y-\Delta y) - f(x-\Delta x,\, y+\Delta y) + f(x-\Delta x,\, y-\Delta y)}{4\,\Delta x\,\Delta y}}$$

This requires four function evaluations at the "diagonal" points around $(x, y)$.

---

<br>

## 8. Method of Undetermined Coefficients

This is a systematic technique for deriving finite difference formulas when you want to approximate a derivative using a specified set of points.

### 8.1 Problem Setup

Seek an approximation of the second derivative $f''(x)$ based on three equally spaced points $f(x-h)$, $f(x)$, $f(x+h)$:

$$f''(x) \approx a\,f(x-h) + b\,f(x) + c\,f(x+h)$$

**Question**: How can we determine $a$, $b$, $c$?

### 8.2 Taylor Expansion and Matching

Expand $f(x-h)$ and $f(x+h)$ in Taylor series (with remainders involving $\xi_-$ and $\xi_+$ in $[x-h, x+h]$):

$$f(x-h) = f(x) - f'(x)\,h + f''(x)\frac{h^2}{2} - f'''(x)\frac{h^3}{6} + f^{(4)}(\xi_-)\frac{h^4}{24}$$

$$f(x+h) = f(x) + f'(x)\,h + f''(x)\frac{h^2}{2} + f'''(x)\frac{h^3}{6} + f^{(4)}(\xi_+)\frac{h^4}{24}$$

Substituting into the approximation and grouping by derivatives:

$$a\,f(x-h) + b\,f(x) + c\,f(x+h) = (a + b + c)\,f(x) + (-a + c)\,h\,f'(x) + (a + c)\frac{h^2}{2}\,f''(x)$$
$$\quad + (c - a)\frac{h^3}{6}\,f'''(x) + \frac{h^4}{24}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)$$

### 8.3 Solving for Coefficients

For this expression to equal $f''(x)$, we need:

| Condition | Equation |
|:----------|:---------|
| Coefficient of $f(x) = 0$ | $a + b + c = 0$ |
| Coefficient of $f'(x) = 0$ | $-a + c = 0$ (equivalently, $c = a$) |
| Coefficient of $f''(x) = 1$ | $(a + c)\dfrac{h^2}{2} = 1$ |

From $c = a$ and $(a + c)\dfrac{h^2}{2} = 1$:

$$2a \cdot \frac{h^2}{2} = 1 \implies a = \frac{1}{h^2}$$

Therefore:

$$\boxed{a = \frac{1}{h^2}, \quad c = \frac{1}{h^2}, \quad b = -\frac{2}{h^2}}$$

### 8.4 Truncation Error

Plugging back in, the third-order term vanishes since $c - a = 0$. The truncation error comes from the fourth-order remainder:

$$f''(x) = \frac{1}{h^2}\bigl(f(x-h) - 2f(x) + f(x+h)\bigr) - \frac{h^4}{24}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)$$

By the **Intermediate Value Theorem**, since $f^{(4)}$ is continuous and $f^{(4)}(\xi_-)$ and $f^{(4)}(\xi_+)$ are two values of $f^{(4)}$, there exists some $\xi \in [x-h, x+h]$ such that:

$$f^{(4)}(\xi) = \frac{1}{2}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)$$

Therefore:

$$\boxed{f''(x) = \frac{f(x-h) - 2f(x) + f(x+h)}{h^2} - \frac{h^2}{12}\,f^{(4)}(\xi)}$$

This confirms the central difference for $f''$ is $O(h^2)$ accurate, with the explicit truncation error term $-\dfrac{h^2}{12}f^{(4)}(\xi)$.

> **[Calculus]** The method of undetermined coefficients is extremely versatile. By choosing which points to include and which derivative to approximate, you can derive any finite difference stencil. The system of equations is always linear in the unknown coefficients, making it straightforward to solve.

---

<br>

## 9. Differentiable Programming (Brief Overview)

**Differentiable programming** is a modern programming paradigm where everything in the program is written so that it can be **differentiated automatically**.

Key points:

- You can compute how small changes in input affect the output
- Enables **gradient-based optimization** (e.g., training neural networks)
- Represents a **new hybrid approach** combining traditional numerical methods and data-driven methods
- Provides **exact gradients** using automatic differentiation (AD), which is more accurate than finite differences and more efficient than symbolic differentiation

**Popular tools** for differentiable programming:
- **PyTorch** (Python)
- **TensorFlow** (Python)
- **JAX** (Python, by Google)
- **Zygote** (Julia)

> **[Calculus]** Automatic differentiation is fundamentally different from numerical differentiation. While numerical differentiation approximates derivatives using finite differences (and suffers from truncation and round-off errors), AD computes derivatives to machine precision by systematically applying the chain rule through every operation in the program.

---

<br>

## 10. Python Implementation

### Forward, Backward, and Central Differences

```python
import numpy as np

def forward_diff(f, x, h):
    """First derivative using forward difference -- O(h)"""
    return (f(x + h) - f(x)) / h

def backward_diff(f, x, h):
    """First derivative using backward difference -- O(h)"""
    return (f(x) - f(x - h)) / h

def central_diff(f, x, h):
    """First derivative using central difference -- O(h^2)"""
    return (f(x + h) - f(x - h)) / (2 * h)

def central_diff_2nd(f, x, h):
    """Second derivative using central difference -- O(h^2)"""
    return (f(x + h) - 2 * f(x) + f(x - h)) / h**2

def forward_diff_2nd_order(f, x, h):
    """First derivative using 2nd-order forward difference -- O(h^2)"""
    return (-f(x + 2*h) + 4*f(x + h) - 3*f(x)) / (2 * h)
```

### Richardson Extrapolation for Differentiation

```python
def richardson_diff(f, x, h, levels=3):
    """
    Richardson extrapolation for central difference.
    Returns a triangular table of increasingly accurate estimates.
    """
    # Build step sizes: h, 2h, 4h, ...
    D = np.zeros((levels, levels))

    # First column: central differences with step sizes h, 2h, 4h, ...
    for i in range(levels):
        hi = h * (2 ** i)
        D[i, 0] = (f(x + hi) - f(x - hi)) / (2 * hi)

    # Richardson extrapolation (note: entries go from fine to coarse)
    # We rearrange so D[0,0] uses smallest h (most accurate base)
    # Rebuild with D[i,0] using step size h * 2^(levels-1-i)
    R = np.zeros((levels, levels))
    for i in range(levels):
        hi = h * (2 ** (levels - 1 - i))
        R[i, 0] = (f(x + hi) - f(x - hi)) / (2 * hi)

    for j in range(1, levels):
        for i in range(j, levels):
            R[i, j] = R[i, j-1] + (R[i, j-1] - R[i-1, j-1]) / (4**j - 1)

    return R

# --- Example ---
f = lambda x: np.exp(x)   # f(x) = e^x, true derivative = e^x
x0 = 1.0
h = 0.1
true_val = np.exp(x0)

print(f"True f'({x0}) = {true_val:.10f}")
print(f"Forward diff:           {forward_diff(f, x0, h):.10f}  "
      f"Error: {abs(forward_diff(f, x0, h) - true_val):.2e}")
print(f"Central diff:           {central_diff(f, x0, h):.10f}  "
      f"Error: {abs(central_diff(f, x0, h) - true_val):.2e}")
print(f"2nd-order forward diff: {forward_diff_2nd_order(f, x0, h):.10f}  "
      f"Error: {abs(forward_diff_2nd_order(f, x0, h) - true_val):.2e}")

R = richardson_diff(f, x0, h, levels=3)
print(f"Richardson (level 1):   {R[1,1]:.10f}  "
      f"Error: {abs(R[1,1] - true_val):.2e}")
print(f"Richardson (level 2):   {R[2,2]:.10f}  "
      f"Error: {abs(R[2,2] - true_val):.2e}")
```

### Derivative of Unequally Spaced Data (Lagrange)

```python
def lagrange_deriv(x_pts, y_pts, x_eval):
    """
    Compute the derivative at x_eval using the derivative
    of the Lagrange interpolating polynomial.

    Parameters:
        x_pts: array of x data points
        y_pts: array of y data points (same length)
        x_eval: point at which to evaluate the derivative

    Returns:
        Approximate value of f'(x_eval)
    """
    n = len(x_pts)
    deriv = 0.0

    for i in range(n):
        # Compute N_i'(x_eval)
        Ni_prime = 0.0
        for j in range(n):
            if j == i:
                continue
            # Product of all (x_eval - x_k) / (x_i - x_k) for k != i, k != j
            prod = 1.0
            for k in range(n):
                if k == i or k == j:
                    continue
                prod *= (x_eval - x_pts[k]) / (x_pts[i] - x_pts[k])
            Ni_prime += prod / (x_pts[i] - x_pts[j])
        deriv += y_pts[i] * Ni_prime

    return deriv

# --- Example: unequally spaced data ---
x_data = np.array([1.0, 2.5, 4.0])
y_data = np.exp(x_data)  # f(x) = e^x
x_eval = 2.5

approx = lagrange_deriv(x_data, y_data, x_eval)
true = np.exp(x_eval)
print(f"\nUnequally spaced data derivative at x={x_eval}:")
print(f"  Lagrange approx: {approx:.8f}")
print(f"  True value:      {true:.8f}")
print(f"  Error:           {abs(approx - true):.2e}")
```

### Numerical Partial Derivatives

```python
def partial_x_central(f, x, y, dx):
    """Central difference for df/dx"""
    return (f(x + dx, y) - f(x - dx, y)) / (2 * dx)

def partial_y_central(f, x, y, dy):
    """Central difference for df/dy"""
    return (f(x, y + dy) - f(x, y - dy)) / (2 * dy)

def mixed_partial_xy(f, x, y, dx, dy):
    """Central difference for d^2f/(dx dy)"""
    return (f(x+dx, y+dy) - f(x+dx, y-dy)
          - f(x-dx, y+dy) + f(x-dx, y-dy)) / (4 * dx * dy)

# --- Example ---
g = lambda x, y: x**2 * y + np.sin(x * y)  # test function
x0, y0 = 1.0, 2.0
h = 1e-5

print(f"\nPartial derivatives of f(x,y) = x^2*y + sin(x*y) at ({x0},{y0}):")
print(f"  df/dx (central):      {partial_x_central(g, x0, y0, h):.8f}")
print(f"  df/dy (central):      {partial_y_central(g, x0, y0, h):.8f}")
print(f"  d^2f/(dx dy) (mixed): {mixed_partial_xy(g, x0, y0, h, h):.8f}")
```

---

<br>

## Summary

| Method | Formula | Accuracy | Points Used |
|:-------|:--------|:---------|:------------|
| **Forward Difference** ($f'$) | $\dfrac{f_{i+1} - f_i}{h}$ | $O(h)$ | $x_i, x_{i+1}$ |
| **Backward Difference** ($f'$) | $\dfrac{f_i - f_{i-1}}{h}$ | $O(h)$ | $x_{i-1}, x_i$ |
| **Central Difference** ($f'$) | $\dfrac{f_{i+1} - f_{i-1}}{2h}$ | $O(h^2)$ | $x_{i-1}, x_{i+1}$ |
| **2nd-Order Forward** ($f'$) | $\dfrac{-f_{i+2} + 4f_{i+1} - 3f_i}{2h}$ | $O(h^2)$ | $x_i, x_{i+1}, x_{i+2}$ |
| **Central Difference** ($f''$) | $\dfrac{f_{i+1} - 2f_i + f_{i-1}}{h^2}$ | $O(h^2)$ | $x_{i-1}, x_i, x_{i+1}$ |
| **Forward Difference** ($f''$) | $\dfrac{f_{i+2} - 2f_{i+1} + f_i}{h^2}$ | $O(h)$ | $x_i, x_{i+1}, x_{i+2}$ |
| **Richardson Extrap. (1 level)** | $\dfrac{4}{3}D(h) - \dfrac{1}{3}D(2h)$ | $O(h^4)$ | From $D(h)$ |
| **Richardson Extrap. (2 levels)** | $\dfrac{16}{15}\bar{D}(h) - \dfrac{1}{15}\bar{D}(2h)$ | $O(h^6)$ | From $\bar{D}(h)$ |
| **Mixed Partial** ($f_{xy}$) | $\dfrac{f(x\!+\!\Delta x, y\!+\!\Delta y) - f(x\!+\!\Delta x, y\!-\!\Delta y) - f(x\!-\!\Delta x, y\!+\!\Delta y) + f(x\!-\!\Delta x, y\!-\!\Delta y)}{4\Delta x \Delta y}$ | $O(\Delta x^2, \Delta y^2)$ | 4 diagonal pts |

**Key Takeaways**:
- All finite difference formulas are derived from **Taylor series expansions**
- **Central differences** are more accurate than forward/backward differences at the same step size due to error cancellation from symmetry
- **Richardson extrapolation** systematically eliminates leading error terms by combining estimates at different step sizes, increasing accuracy from $O(h^2) \to O(h^4) \to O(h^6) \to \cdots$
- For **unequally spaced data**, differentiate the Lagrange interpolating polynomial
- The **method of undetermined coefficients** provides a systematic way to derive any finite difference stencil
- **Differentiable programming** (PyTorch, JAX, TensorFlow) provides exact gradients via automatic differentiation, representing a modern alternative to numerical differentiation
