# Chapter 22 Lecture — Initial-Value Problems for ODEs

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 22

> **Prerequisites**: [Calculus] ODEs, derivatives (Ch 1-21).
>
> **Learning Objectives**:
> 1. Apply Euler's method for initial value problems
> 2. Implement Runge-Kutta methods (RK2, RK4)
> 3. Analyze stability and accuracy of ODE solvers

---

<br>

## Table of Contents

- [1. Overview of ODEs and Initial-Value Problems](#1-overview-of-odes-and-initial-value-problems)
  - [1.1 Ordinary Differential Equations](#11-ordinary-differential-equations)
  - [1.2 The Initial-Value Problem (IVP)](#12-the-initial-value-problem-ivp)
  - [1.3 Systems of ODEs](#13-systems-of-odes)
  - [1.4 Integral Formulation](#14-integral-formulation)
- [2. Euler Method](#2-euler-method)
  - [2.1 Derivation](#21-derivation)
  - [2.2 Geometric Interpretation](#22-geometric-interpretation)
- [3. Local and Global Truncation Errors](#3-local-and-global-truncation-errors)
  - [3.1 Global Error](#31-global-error)
  - [3.2 Big-O Notation for Errors](#32-big-o-notation-for-errors)
  - [3.3 Taylor Series Analysis of Euler's Method](#33-taylor-series-analysis-of-eulers-method)
  - [3.4 Approximate Local Truncation Error](#34-approximate-local-truncation-error)
  - [3.5 From Local to Global Error](#35-from-local-to-global-error)
- [4. Runge-Kutta (RK) Methods](#4-runge-kutta-rk-methods)
  - [4.1 Core Idea](#41-core-idea)
  - [4.2 The Increment Function](#42-the-increment-function)
  - [4.3 Explicit RK Methods — General Framework](#43-explicit-rk-methods--general-framework)
  - [4.4 Butcher Tableau](#44-butcher-tableau)
- [5. 2nd-Order Runge-Kutta Methods](#5-2nd-order-runge-kutta-methods)
  - [5.1 General 2-Stage RK Form](#51-general-2-stage-rk-form)
  - [5.2 Deriving the Order Conditions](#52-deriving-the-order-conditions)
  - [5.3 Order Conditions for 2nd-Order Accuracy](#53-order-conditions-for-2nd-order-accuracy)
  - [5.4 Heun's Method (b2 = 1/2)](#54-heuns-method-b2--12)
  - [5.5 Modified Euler / Midpoint Method (b2 = 1)](#55-modified-euler--midpoint-method-b2--1)
- [6. Classical 4th-Order Runge-Kutta Method (RK4)](#6-classical-4th-order-runge-kutta-method-rk4)
  - [6.1 The Formula](#61-the-formula)
  - [6.2 Butcher Tableau for RK4](#62-butcher-tableau-for-rk4)
  - [6.3 Python Implementation](#63-python-implementation)
- [7. Systems of ODEs](#7-systems-of-odes)
  - [7.1 Vector Form](#71-vector-form)
  - [7.2 Higher-Order ODEs as Systems](#72-higher-order-odes-as-systems)
  - [7.3 RK4 for Systems — Python Implementation](#73-rk4-for-systems--python-implementation)
- [8. Numerical Stability](#8-numerical-stability)
  - [8.1 The Model Problem](#81-the-model-problem)
  - [8.2 Stability of the Exact Solution](#82-stability-of-the-exact-solution)
  - [8.3 Stability of Euler's Method (Forward Euler)](#83-stability-of-eulers-method-forward-euler)
  - [8.4 Stability Region of Forward Euler](#84-stability-region-of-forward-euler)
  - [8.5 Stability of Backward Euler Method](#85-stability-of-backward-euler-method)
  - [8.6 Stability of the Trapezoidal Method](#86-stability-of-the-trapezoidal-method)
  - [8.7 A-Stability and L-Stability](#87-a-stability-and-l-stability)
- [Summary](#summary)

---

<br>

## 1. Overview of ODEs and Initial-Value Problems

### 1.1 Ordinary Differential Equations

Let $y = y(t)$ be an unknown function of the independent variable $t$. An **ordinary differential equation (ODE)** relates $y$ and its derivatives:

$$\frac{dy}{dt} = f(t, y)$$

The function $f(t, y)$ is known and gives the **slope** (tendency) of the solution at any point $(t, y)$ in the $t$-$y$ plane. The solution $y(t)$ traces a **solution curve** through this slope field.

> **[Calculus]** Unlike algebraic equations where the unknown is a number, the unknown in a differential equation is a *function*. The ODE $\frac{dy}{dt} = f(t,y)$ states that the rate of change of $y$ at any time $t$ is determined by the current values of $t$ and $y$.

### 1.2 The Initial-Value Problem (IVP)

Given an initial condition $y(t_0) = y_0$, solve for $y(t)$ for $t > t_0$. This is the **initial-value problem**:

$$\begin{cases} \dfrac{dy}{dt} = f(t, y), & t > t_0 \\[6pt] y(t_0) = g \end{cases}$$

The initial condition selects a unique solution curve from the family of all solutions to the ODE.

### 1.3 Systems of ODEs

Consider a pair of coupled ODEs:

$$\frac{du}{dt} = p(t, u, v), \quad u(t_0) = u_0$$

$$\frac{dv}{dt} = q(t, u, v), \quad v(t_0) = v_0$$

In **vector form**, define $\mathbf{y} := (u, v)^T$ and $\mathbf{f} := (p, q)^T$:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}), \quad \mathbf{y}(t_0) = (u_0, v_0)^T$$

This is a **system of ODEs**. All single-equation methods (Euler, Runge-Kutta, etc.) extend naturally to systems by replacing scalar operations with vector operations.

### 1.4 Integral Formulation

Integrating the ODE from $t_0$ to $t_1$:

$$\int_{t_0}^{t_1} \frac{dy}{dt}\,dt = y(t_1) - y(t_0) = \int_{t_0}^{t_1} f(t, y(t))\,dt$$

Therefore:

$$y(t_1) = y(t_0) + \int_{t_0}^{t_1} f(t, y(t))\,dt$$

The key question is: **how can we numerically approximate this integral** when we do not know $y(t)$ in closed form?

---

<br>

## 2. Euler Method

### 2.1 Derivation

Decompose the time interval $[t_0, T]$ into $n$ subintervals with uniform step size:

$$h = \frac{T - t_0}{n}, \quad t_{n+1} = t_n + h$$

From the integral formulation over one step:

$$y(t_{n+1}) = y(t_n) + \int_{t_n}^{t_{n+1}} f(t, y(t))\,dt$$

The simplest approximation is to treat $f$ as constant over $[t_n, t_{n+1}]$, evaluated at the left endpoint:

$$\int_{t_n}^{t_{n+1}} f(t, y(t))\,dt \approx f(t_n, y_n) \cdot \int_{t_n}^{t_{n+1}} dt = f(t_n, y_n) \cdot h$$

This gives **Euler's method** (Forward Euler):

$$\boxed{y_{n+1} = y_n + h\,f_n, \quad f_n := f(t_n, y_n)}$$

### 2.2 Geometric Interpretation

At each time step, Euler's method advances the solution along the tangent line at the current point $(t_n, y_n)$. The slope of this tangent is $f(t_n, y_n)$, so moving a horizontal distance $h$ produces a vertical change of $h \cdot f(t_n, y_n)$.

> **[Calculus]** Euler's method is the simplest member of a vast family of numerical ODE solvers. It corresponds to approximating the integral $\int f\,dt$ by a left-endpoint rectangle rule. Better quadrature rules lead to better ODE methods (this is the Runge-Kutta idea).

---

<br>

## 3. Local and Global Truncation Errors

### 3.1 Global Error

The **global error** at step $n$ is the difference between the true solution and the numerical approximation:

$$e_n = y(t_n) - y_n$$

Note that $e_0 = y(t_0) - y_0 = 0$ by the initial condition. The global error is the **accumulation** of local errors over all steps from $t_0$ to $t_n$.

### 3.2 Big-O Notation for Errors

A quantity $z$ is $O(h^p)$ if there exist constants $h_0, C > 0$ such that:

$$|z| \le C h^p \quad \forall\; 0 < h < h_0$$

This means $z$ converges to zero as $h \to 0$, and the **rate of convergence** is $p$.

### 3.3 Taylor Series Analysis of Euler's Method

Recall the Taylor expansion of $y$ about $(t_n, y_n)$ with $h = t_{n+1} - t_n$:

$$y(t_n + h) = y(t_n) + y'(t_n)\,h + y''(t_n)\frac{h^2}{2} + \cdots + y^{(n)}(t_n)\frac{h^n}{n!} + R_n$$

where $R_n = y^{(n+1)}(\xi)\frac{h^{n+1}}{(n+1)!}$ for some $\xi \in [t_n, t_n + h]$.

Since $y'(t_n) = f(t_n, y_n)$, we can rewrite the Taylor expansion using derivatives of $f$:

$$y(t_n + h) = \underbrace{y_n + f(t_n, y_n)\,h}_{\text{Euler approx}} + \underbrace{f'(t_n, y_n)\frac{h^2}{2} + \cdots + f^{(n)}(t_n, y_n)\frac{h^n}{n!} + O(h^{n+1})}_{\text{truncated terms (local truncation error)}}$$

The **true local truncation error** of Euler's method is $O(h^2)$ — the leading term of what is discarded.

### 3.4 Approximate Local Truncation Error

The leading term of the local truncation error for Euler's method is:

$$E_a = f'(t_n, y_n)\frac{h^2}{2} = O(h^2)$$

This is the **approximate local truncation error** at each step.

### 3.5 From Local to Global Error

With a uniform step size $h$ over $[0, T]$, we have $n = \frac{T}{h}$ steps. Summing the local errors:

$$\left|\sum_{i=1}^{n} E_{a,i}\right| = \left|(f'_0 + f'_1 + \cdots + f'_{n-1})\frac{h^2}{2}\right| \le \tilde{f}' \cdot \frac{h^2}{2} \cdot n$$

where $\tilde{f}' = \max_{0 \le i \le n-1} |f'_i|$.

Substituting $n = \frac{T}{h}$:

$$\left|\sum E_{a,i}\right| \le \tilde{f}' \cdot \frac{h^2}{2} \cdot \frac{T}{h} = \frac{\tilde{f}'\,T}{2}\,h = O(h)$$

**Key result: Euler's method has local truncation error $O(h^2)$ and global truncation error $O(h)$.** Euler's method is a **first-order method**.

> **[Calculus]** The local error loses one order of $h$ when accumulated into global error because the number of steps $n = T/h$ grows as $h$ decreases. This is a universal pattern: for a method with local error $O(h^{p+1})$, the global error is $O(h^p)$.

---

<br>

## 4. Runge-Kutta (RK) Methods

### 4.1 Core Idea

Runge-Kutta methods **improve accuracy by combining multiple slope estimates** from different points within the interval $[t_n, t_{n+1}]$, without requiring explicit computation of higher derivatives of $f$.

Starting from the integral form and performing a change of variable $t = (t_{n+1} - t_n)z + t_n$ where $z \in [0, 1]$:

$$y_{n+1} = y_n + \int_{t_n}^{t_{n+1}} f(t, y)\,dt = y_n + h \int_0^1 f(z, y(z))\,dz$$

### 4.2 The Increment Function

Approximate the integral $\int_0^1 f(z, y(z))\,dz$ using a quadrature rule:

$$\int_0^1 f(z, y(z))\,dz \approx \phi$$

where $\phi$ is the **increment function** (an averaged slope):

$$y_{n+1} = y_n + \phi \cdot h$$

The increment function is a weighted sum of slope evaluations:

$$\phi = \sum_{i=1}^{s} b_i k_i = \sum_{i=1}^{s} b_i\,f(t_n + c_i h,\; y(t_n + c_i h))$$

This is essentially a **quadrature** of $f$ over the step, where $b_i$ are weights and $c_i$ are nodes.

### 4.3 Explicit RK Methods — General Framework

Since we do not know $y(t_n + c_i h)$ in general, we approximate it by an **internal stage** value $Y_i$:

**Internal stages** (computed sequentially):

$$Y_1 = y_n \quad \Rightarrow \quad k_1 = f(t_n, Y_1)$$

$$Y_2 = y_n + h\,a_{21}\,k_1 \quad \Rightarrow \quad k_2 = f(t_n + c_2 h, Y_2)$$

$$Y_3 = y_n + h(a_{31}\,k_1 + a_{32}\,k_2) \quad \Rightarrow \quad k_3 = f(t_n + c_3 h, Y_3)$$

$$\vdots$$

$$Y_s = y_n + h \sum_{j=1}^{s-1} a_{sj}\,k_j \quad \Rightarrow \quad k_s = f(t_n + c_s h, Y_s)$$

**Step completion:**

$$y_{n+1} = y_n + h \sum_{i=1}^{s} b_i\,k_i$$

where $s$ is the **number of stages**.

> **[Calculus]** The key insight of RK methods is replacing the intractable problem of evaluating $f$ at unknown intermediate solution values with a bootstrapping procedure: use Euler-like predictions to estimate intermediate values, evaluate $f$ there, then combine those slopes.

### 4.4 Butcher Tableau

An explicit RK method is fully specified by the matrix $A$ (strictly lower triangular), the weight vector $\mathbf{b} = (b_1, b_2, \ldots, b_s)^T$, and the node vector $\mathbf{c} = (c_1, c_2, \ldots, c_s)^T$.

In matrix form for the internal stages:

$$\begin{pmatrix} Y_1 \\ Y_2 \\ Y_3 \\ \vdots \\ Y_s \end{pmatrix} = \begin{pmatrix} y_n \\ y_n \\ y_n \\ \vdots \\ y_n \end{pmatrix} + h \begin{pmatrix} 0 & 0 & 0 & \cdots & 0 \\ a_{21} & 0 & 0 & \cdots & 0 \\ a_{31} & a_{32} & 0 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ a_{s1} & a_{s2} & a_{s3} & \cdots & 0 \end{pmatrix} \begin{pmatrix} k_1 \\ k_2 \\ k_3 \\ \vdots \\ k_s \end{pmatrix}$$

The matrix $A$ being strictly lower triangular is what makes the method **explicit** — each stage can be computed from previously computed stages.

The **Butcher tableau** is the standard compact notation:

$$\begin{array}{c|c} \mathbf{c} & A \\ \hline & \mathbf{b}^T \end{array}$$

---

<br>

## 5. 2nd-Order Runge-Kutta Methods

### 5.1 General 2-Stage RK Form

A 2-stage ($s = 2$) explicit RK method has the form:

$$y_{n+1} = y_n + h(b_1 k_1 + b_2 k_2)$$

where:

$$k_1 = f(t_n, y_n)$$

$$k_2 = f(t_n + ch,\; y_n + h\,a\,k_1)$$

The parameters to determine are $b_1, b_2, c, a$ (4 unknowns).

### 5.2 Deriving the Order Conditions

**Taylor expansion of the true solution** (using total derivative with respect to $t$):

$$y(t_n + h) = y_n + h\,f + \frac{h^2}{2}(f_t + f_y f) + O(h^3)$$

where all $f, f_t, f_y$ are evaluated at $(t_n, y_n)$.

> **[Calculus]** The second derivative $y'' = \frac{d}{dt}f(t, y) = f_t + f_y \cdot \frac{dy}{dt} = f_t + f_y f$ uses the **chain rule** for the total derivative. This is because $f$ depends on $t$ both directly and indirectly through $y(t)$.

**Taylor expansion of the RK2 formula:**

Expanding $k_2 = f(t_n + ch, y_n + a\,f(t_n,y_n)\,h)$ using a multivariable Taylor expansion with $\Delta t = ch$ and $\Delta y = a\,f(t_n,y_n)\,h$:

$$f(t_n + \Delta t, y_n + \Delta y) = f + f_t \cdot ch + f_y \cdot a\,f\,h + O(h^2)$$

Substituting back into the RK2 formula:

$$y_{n+1} = y_n + h(b_1 + b_2)f + \frac{h^2}{2} \cdot 2b_2(c\,f_t + a\,f_y f) + O(h^3)$$

### 5.3 Order Conditions for 2nd-Order Accuracy

Matching with the Taylor expansion up to $O(h^2)$ terms gives **3 equations in 4 unknowns** (one degree of freedom):

$$b_1 + b_2 = 1$$

$$b_2 c = \frac{1}{2}$$

$$b_2 a = \frac{1}{2}$$

Note: the last two conditions imply $c = a$.

### 5.4 Heun's Method (b2 = 1/2)

Choose $b_2 = 1/2$:

$$b_1 = \frac{1}{2}, \quad a = 1, \quad c = 1$$

$$Y_1 = y_n \quad \Rightarrow \quad k_1 = f(t_n, y_n)$$

$$Y_2 = y_n + h\,k_1 \quad \Rightarrow \quad k_2 = f(t_n + h, Y_2)$$

$$\boxed{y_{n+1} = y_n + h\left(\frac{k_1 + k_2}{2}\right)}$$

This is **Heun's method** (also called the improved Euler or trapezoidal predictor-corrector). The Butcher tableau is:

$$\begin{array}{c|cc} 0 & & \\ 1 & 1 & \\ \hline & 1/2 & 1/2 \end{array}$$

> **[Calculus]** Heun's method first predicts $Y_2$ using a full Euler step to $t_{n+1}$, evaluates the slope there, then averages the slopes at both endpoints. This is analogous to the trapezoidal rule for integration.

### 5.5 Modified Euler / Midpoint Method (b2 = 1)

Choose $b_2 = 1$:

$$b_1 = 0, \quad a = \frac{1}{2}, \quad c = \frac{1}{2}$$

$$k_1 = f(t_n, y_n)$$

$$k_2 = f\!\left(t_n + \frac{h}{2},\; y_n + \frac{h}{2}\,k_1\right)$$

$$\boxed{y_{n+1} = y_n + h\,k_2}$$

This is the **modified Euler** (midpoint) method. The Butcher tableau is:

$$\begin{array}{c|cc} 0 & & \\ 1/2 & 1/2 & \\ \hline & 0 & 1 \end{array}$$

> **[Calculus]** The midpoint method uses an Euler half-step to estimate the slope at the interval's midpoint, then uses that midpoint slope to advance the full step. This is analogous to the midpoint rule for integration.

---

<br>

## 6. Classical 4th-Order Runge-Kutta Method (RK4)

### 6.1 The Formula

The classical RK4 is a 4-stage method with local truncation error $O(h^5)$ and global error $O(h^4)$:

$$k_1 = f(t_n,\; y_n)$$

$$k_2 = f\!\left(t_n + \frac{h}{2},\; y_n + \frac{h}{2}\,k_1\right)$$

$$k_3 = f\!\left(t_n + \frac{h}{2},\; y_n + \frac{h}{2}\,k_2\right)$$

$$k_4 = f(t_n + h,\; y_n + h\,k_3)$$

$$\boxed{y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)}$$

The weights $\frac{1}{6}, \frac{2}{6}, \frac{2}{6}, \frac{1}{6}$ correspond to **Simpson's 1/3 rule** applied to the four slope estimates.

> **[Calculus]** RK4 is the workhorse of ODE solving. It achieves 4th-order accuracy using only 4 function evaluations per step, without requiring any derivatives of $f$. For many practical problems, it provides an excellent balance between accuracy and computational cost.

### 6.2 Butcher Tableau for RK4

$$\begin{array}{c|cccc} 0 & & & & \\ 1/2 & 1/2 & & & \\ 1/2 & 0 & 1/2 & & \\ 1 & 0 & 0 & 1 & \\ \hline & 1/6 & 1/3 & 1/3 & 1/6 \end{array}$$

### 6.3 Python Implementation

```python
import numpy as np

def rk4(f, t0, y0, T, h):
    """
    Classical 4th-order Runge-Kutta method.

    Parameters
    ----------
    f  : callable, f(t, y) -> dy/dt
    t0 : float, initial time
    y0 : float, initial value
    T  : float, final time
    h  : float, step size

    Returns
    -------
    t_arr : ndarray, time values
    y_arr : ndarray, solution values
    """
    n = int((T - t0) / h)
    t_arr = np.linspace(t0, T, n + 1)
    y_arr = np.zeros(n + 1)
    y_arr[0] = y0

    for i in range(n):
        t = t_arr[i]
        y = y_arr[i]
        k1 = f(t, y)
        k2 = f(t + h/2, y + h/2 * k1)
        k3 = f(t + h/2, y + h/2 * k2)
        k4 = f(t + h, y + h * k3)
        y_arr[i+1] = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    return t_arr, y_arr
```

---

<br>

## 7. Systems of ODEs

### 7.1 Vector Form

A system of $m$ first-order ODEs is written in vector form as:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}), \quad \mathbf{y}(t_0) = \mathbf{y}_0$$

where $\mathbf{y} = (y_1, y_2, \ldots, y_m)^T$ and $\mathbf{f} = (f_1, f_2, \ldots, f_m)^T$.

All Euler and Runge-Kutta formulas apply **component-by-component** — simply replace scalar $y$, $f$, and $k$ with vectors $\mathbf{y}$, $\mathbf{f}$, and $\mathbf{k}$.

### 7.2 Higher-Order ODEs as Systems

An $m$-th order ODE can be converted to a system of $m$ first-order ODEs. For example, the second-order ODE:

$$y'' = g(t, y, y')$$

is rewritten by introducing $u = y$ and $v = y'$:

$$\frac{du}{dt} = v, \quad \frac{dv}{dt} = g(t, u, v)$$

with initial conditions $u(t_0) = y(t_0)$ and $v(t_0) = y'(t_0)$.

> **[Calculus]** This reduction technique is fundamental: *any* numerical method for first-order systems can solve higher-order ODEs by first converting them. For instance, Newton's second law $m\,x'' = F(t, x, x')$ becomes a first-order system in position $x$ and velocity $v = x'$.

### 7.3 RK4 for Systems — Python Implementation

```python
import numpy as np

def rk4_system(f, t0, y0, T, h):
    """
    RK4 for a system of ODEs.

    Parameters
    ----------
    f  : callable, f(t, y) -> ndarray of shape (m,)
    t0 : float, initial time
    y0 : ndarray of shape (m,), initial conditions
    T  : float, final time
    h  : float, step size

    Returns
    -------
    t_arr : ndarray of shape (n+1,)
    y_arr : ndarray of shape (n+1, m)
    """
    n = int((T - t0) / h)
    m = len(y0)
    t_arr = np.linspace(t0, T, n + 1)
    y_arr = np.zeros((n + 1, m))
    y_arr[0] = y0

    for i in range(n):
        t = t_arr[i]
        y = y_arr[i]
        k1 = f(t, y)
        k2 = f(t + h/2, y + h/2 * k1)
        k3 = f(t + h/2, y + h/2 * k2)
        k4 = f(t + h, y + h * k3)
        y_arr[i+1] = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    return t_arr, y_arr


# Example: Solve y'' + y = 0  (simple harmonic oscillator)
# Convert to system: u' = v, v' = -u
# with u(0) = 1, v(0) = 0  => exact solution: u(t) = cos(t)

def harmonic(t, y):
    u, v = y
    return np.array([v, -u])

t, sol = rk4_system(harmonic, 0, np.array([1.0, 0.0]), 10, 0.01)
# sol[:, 0] is u(t) ≈ cos(t)
# sol[:, 1] is v(t) ≈ -sin(t)
```

---

<br>

## 8. Numerical Stability

### 8.1 The Model Problem

A numerical solution is said to be **unstable** if errors grow exponentially. To analyze stability, consider the **model (test) problem**:

$$\frac{dy}{dt} = \lambda y, \quad y(0) = y_0$$

The exact solution is:

$$y(t) = y_0 e^{\lambda t}$$

### 8.2 Stability of the Exact Solution

The behavior of the exact solution depends on $\lambda$:

| Value of $\lambda$ | Behavior of $y(t) = y_0 e^{\lambda t}$ |
|:---|:---|
| $\lambda > 0$ | Exponential growth (unstable) |
| $\lambda = 0$ | Constant (neutrally stable) |
| $\lambda < 0$ | Exponential decay (stable) |

The solution asymptotically approaches zero when $\lambda < 0$. A **stable** numerical method should reproduce this decaying behavior for $\lambda < 0$.

### 8.3 Stability of Euler's Method (Forward Euler)

Applying Euler's method to $y' = \lambda y$:

$$y_{n+1} = y_n + h \lambda y_n = y_n(1 + \lambda h)$$

The factor $(1 + \lambda h)$ is the **amplification factor** (or **stability function** $\phi(\lambda h)$). By induction:

$$y_1 = y_0(1 + \lambda h)$$

$$y_2 = y_0(1 + \lambda h)^2$$

$$\vdots$$

$$y_n = y_0(1 + \lambda h)^n$$

For stability (bounded solution), we need:

$$|1 + \lambda h| \le 1$$

### 8.4 Stability Region of Forward Euler

Let $z = \lambda h \in \mathbb{C}$. The stability condition is:

$$|1 + z| \le 1$$

This is a **disk of radius 1 centered at $(-1, 0)$** in the complex $z$-plane.

On the **real axis**, the condition reduces to:

$$-1 \le 1 + z \le 1 \quad \Rightarrow \quad -2 \le z \le 0$$

Since $z = \lambda h$ with $h > 0$, this requires:

$$\lambda \le 0 \quad \text{and} \quad h \le \frac{-2}{\lambda} = \frac{2}{|\lambda|}$$

**The forward Euler method is conditionally stable** — the step size $h$ must be chosen sufficiently small to ensure stability. For strongly negative $\lambda$ (stiff problems), $h$ must be very small.

### 8.5 Stability of Backward Euler Method

The **backward Euler method** (implicit) evaluates the slope at the *next* time step:

$$x_{n+1} = x_n + h\,f(x_{n+1})$$

Applied to $\frac{dx}{dt} = \lambda x$:

$$x_{n+1} = x_n + \lambda h\,x_{n+1}$$

$$(1 - \lambda h)\,x_{n+1} = x_n$$

$$x_{n+1} = \frac{1}{1 - \lambda h}\,x_n$$

By induction:

$$x_n = (1 - \lambda h)^{-n}\,x_0$$

The stability function is $\phi(z) = \frac{1}{1 - z}$. For stability:

$$\frac{1}{|1 - z|} < 1 \quad \Rightarrow \quad |1 - z| > 1$$

The stability region is the **exterior** of the disk of radius 1 centered at $(1, 0)$ in the complex $z$-plane. This region contains the **entire left half-plane** $\text{Re}(z) < 0$.

For any $\lambda < 0$, we can choose the step size $h$ **arbitrarily large** and the method remains stable. This means the backward Euler method does not require step size restrictions for stability — only accuracy considerations limit $h$.

> **[Calculus]** The price of this unconditional stability is that backward Euler is **implicit**: at each step we must solve an equation (possibly nonlinear) for $x_{n+1}$. For the linear test problem, this reduces to a simple algebraic solve, but for general nonlinear $f$, Newton's method or similar techniques are required.

### 8.6 Stability of the Trapezoidal Method

The **trapezoidal method** uses the averaged slope at both endpoints:

$$x_{n+1} = x_n + \frac{h}{2}(f_n + f_{n+1})$$

Applied to $\frac{dx}{dt} = \lambda x$:

$$x_{n+1} = x_n + \frac{\lambda h}{2}(x_n + x_{n+1})$$

$$\left(1 - \frac{\lambda h}{2}\right)x_{n+1} = \left(1 + \frac{\lambda h}{2}\right)x_n$$

By induction:

$$x_n = \left(\frac{1 + \frac{\lambda h}{2}}{1 - \frac{\lambda h}{2}}\right)^n x_0$$

The stability function is $\phi(z) = \frac{1 + z/2}{1 - z/2}$. For stability:

$$\left|\frac{1 + z/2}{1 - z/2}\right| < 1 \quad \Rightarrow \quad |1 + z/2| < |1 - z/2| \quad \Rightarrow \quad |2 + z| < |2 - z|$$

This condition is satisfied for the **entire left half-plane** $\text{Re}(z) < 0$, meaning the trapezoidal method is also stable for any $\lambda < 0$ regardless of step size.

### 8.7 A-Stability and L-Stability

**A-stability:** A numerical method is **A-stable** if its stability region includes the **entire left half-plane** $\{z \in \mathbb{C} : \text{Re}(z) \le 0\}$.

- Forward Euler: **NOT** A-stable (finite stability region)
- Backward Euler: **A-stable**
- Trapezoidal method: **A-stable**

**L-stability:** A numerical method is **L-stable** if:
1. The method is A-stable, **and**
2. The stability function $\phi(z) \to 0$ as $z \to \infty$

L-stability is **more restrictive** than A-stability. It concerns the **asymptotic behavior** — an L-stable method correctly captures the rapid decay of very stiff components.

- Backward Euler: $\phi(z) = \frac{1}{1-z} \to 0$ as $z \to -\infty$ — **L-stable**
- Trapezoidal method: $\phi(z) = \frac{1 + z/2}{1 - z/2} \to -1$ as $z \to -\infty$ — A-stable but **NOT** L-stable

> **[Calculus]** For stiff problems (where $|\lambda|$ is very large), L-stability is desirable because it ensures the numerical solution of rapidly decaying components also decays rapidly. The trapezoidal method, while A-stable, may exhibit spurious oscillations for very stiff problems because its stability function approaches $-1$ rather than $0$ for large $|z|$.

---

<br>

## Summary

| Topic | Key Formula / Result |
|:------|:---------------------|
| **IVP** | $\frac{dy}{dt} = f(t,y)$, $y(t_0) = y_0$ |
| **Euler method** | $y_{n+1} = y_n + h\,f(t_n, y_n)$ |
| **Euler local error** | $O(h^2)$ |
| **Euler global error** | $O(h)$ — 1st-order method |
| **RK general form** | $y_{n+1} = y_n + h \sum b_i k_i$ |
| **RK2 order conditions** | $b_1 + b_2 = 1$, $b_2 c = 1/2$, $b_2 a = 1/2$ |
| **Heun's method** | $y_{n+1} = y_n + \frac{h}{2}(k_1 + k_2)$, $k_2 = f(t_n+h, y_n+hk_1)$ |
| **Midpoint method** | $y_{n+1} = y_n + h\,k_2$, $k_2 = f(t_n + h/2, y_n + \frac{h}{2}k_1)$ |
| **RK4** | $y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$ — 4th order |
| **Systems of ODEs** | Replace scalars with vectors; same formulas apply |
| **Test equation** | $y' = \lambda y$ with $y(t) = y_0 e^{\lambda t}$ |
| **Forward Euler stability** | $\|1 + \lambda h\| \le 1$ — **conditionally stable** |
| **Backward Euler stability** | $\|1 - \lambda h\|^{-1} < 1$ — **unconditionally stable** (A-stable, L-stable) |
| **Trapezoidal stability** | $\left\|\frac{1+z/2}{1-z/2}\right\| < 1$ — **A-stable but NOT L-stable** |
| **A-stable** | Stability region contains the entire left half-plane |
| **L-stable** | A-stable and $\phi(z) \to 0$ as $z \to \infty$ |
