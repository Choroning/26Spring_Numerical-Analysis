# Chapter 24 — Boundary-Value Problems

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 24

> **Learning Objectives:**
> 1. Formulate boundary value problems for ODEs
> 2. Obtain the analytical solution of a linear BVP as the sum of the homogeneous and particular solutions
> 3. Apply the shooting method and the finite difference method
> 4. Analyze the accuracy and convergence of BVP solvers

---

<br>

## Table of Contents

- [1. Introduction to Boundary-Value Problems (BVP)](#1-introduction-to-boundary-value-problems-bvp)
  - [1.1 Motivating Example: Heated Rod](#11-motivating-example-heated-rod)
  - [1.2 Governing Equation](#12-governing-equation)
  - [1.3 BVP vs. IVP](#13-bvp-vs-ivp)
- [2. Analytical Solution](#2-analytical-solution)
  - [2.1 Rearranging the Governing Equation](#21-rearranging-the-governing-equation)
  - [2.2 Homogeneous Solution](#22-homogeneous-solution)
  - [2.3 Particular Solution](#23-particular-solution)
  - [2.4 Total (General) Solution](#24-total-general-solution)
  - [2.5 Determining the Constants from Boundary Conditions](#25-determining-the-constants-from-boundary-conditions)
- [3. The Shooting Method (Section 24.2)](#3-the-shooting-method-section-242)
  - [3.1 Core Idea: Convert BVP into IVP](#31-core-idea-convert-bvp-into-ivp)
  - [3.2 Converting a 2nd-Order ODE to Two 1st-Order ODEs](#32-converting-a-2nd-order-ode-to-two-1st-order-odes)
  - [3.3 The Guessing Procedure](#33-the-guessing-procedure)
  - [3.4 Algorithm Summary](#34-algorithm-summary)
  - [3.5 Python Implementation](#35-python-implementation)
- [4. Finite Difference Methods (Section 24.3)](#4-finite-difference-methods-section-243)
  - [4.1 Core Idea: Discretize the Domain](#41-core-idea-discretize-the-domain)
  - [4.2 Second Derivative Approximation](#42-second-derivative-approximation)
  - [4.3 Generating the Finite Difference Equations](#43-generating-the-finite-difference-equations)
  - [4.4 Counting Unknowns and Equations](#44-counting-unknowns-and-equations)
  - [4.5 Assembling the Tridiagonal System](#45-assembling-the-tridiagonal-system)
  - [4.6 Matrix Formulation](#46-matrix-formulation)
  - [4.7 Python Implementation](#47-python-implementation)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Introduction to Boundary-Value Problems (BVP)

### 1.1 Motivating Example: Heated Rod

Consider a **steady-state temperature distribution** for a long, thin rod positioned between two constant-temperature walls:

- The left wall is held at temperature $`T_a`$.
- The right wall is held at temperature $`T_b`$.
- The surrounding air temperature is $`T_\infty`$.
- The rod has length $`L`$.

```text
                 Ambient air T∞
 ┃   ↕    ↕    ↕    ↕    ↕    ↕    ↕    ↕   ┃   ↕ Convection
 ┃██████████████████████████████████████████┃   ◀─▶ Conduction (along the rod)
T_a                                        T_b
 ├──────────────────────────────────────────┤──▶ x
 0                                          L
```

*Figure 1. Conduction and convection in a heated rod between two walls*

Two mechanisms of heat transfer act on the rod:

| Mechanism | Direction | Description |
|:----------|:----------|:------------|
| **Conduction** | Along the rod (horizontal) | Heat flows through the rod material from one end to the other |
| **Convection** | Perpendicular to the rod (vertical) | Heat is exchanged between the rod surface and the surrounding air at $`T_\infty`$ |

*Table 1. Two heat transfer mechanisms acting on the heated rod*

The fundamental question is: **What is the temperature distribution $`T(x)`$ along the rod?**

The temperature $`T(x)`$ varies from $`T_a`$ at $`x = 0`$ to $`T_b`$ at $`x = L`$, and the profile between these endpoints depends on the balance of conduction and convection.

> **Note:** **Physics background.** Conduction is governed by Fourier's law and depends on the thermal conductivity of the rod material. Convection depends on the heat transfer coefficient between the rod surface and the surrounding fluid. The parameter $`h'`$ in the governing equation captures the relative strength of convective and conductive effects.

### 1.2 Governing Equation

The governing equation that combines heat conduction along the rod and convection to the surrounding air is:

```math
0 = \frac{d^2 T}{dx^2} + h'(T_\infty - T) \quad \cdots (**)
```

with **boundary conditions**:

```math
T(x = 0) = T_a, \quad T(x = L) = T_b
```

Here $`h'`$ is a parameter related to the convective heat transfer coefficient, the rod geometry, and the thermal conductivity.

Given two boundary values at the endpoints of the domain, we can solve for the temperature distribution $`T(x)`$.

### 1.3 BVP vs. IVP

| Feature | Initial-Value Problem (IVP) | Boundary-Value Problem (BVP) |
|:--------|:---------------------------|:-----------------------------|
| **Conditions given at** | A single point ($`x = x_0`$) | Two or more distinct points |
| **Typical conditions** | $`y(x_0) = y_0`$, $`y'(x_0) = y'_0`$ | $`y(a) = \alpha`$, $`y(b) = \beta`$ |
| **Solution approach** | March forward from $`x_0`$ | Cannot simply march from one end; needs iterative guessing (shooting) or a global solve (finite differences) |
| **Order** | Needs $`n`$ conditions at one point for $`n`$-th order ODE | Conditions split across boundary |

*Table 2. Comparison of initial-value and boundary-value problems*

> **[Calculus]** A second-order ODE generally requires two conditions to determine a unique solution. In an IVP, both conditions are specified at the same point (e.g., position and velocity at $`t = 0`$). In a BVP, the conditions are specified at different spatial locations (e.g., temperatures at both ends of a rod). A linear IVP with continuous coefficients has exactly one solution, but a BVP may have no solution or infinitely many even with two conditions. For example, $`y'' + y = 0`$, $`y(0) = 0`$, $`y(\pi) = 0`$ has the solution $`y = C \sin x`$ for every constant $`C`$.

---

<br>

## 2. Analytical Solution

### 2.1 Rearranging the Governing Equation

Express equation $`(**)`$ as:

```math
\frac{d^2 T}{dx^2} - h' T = -h' T_\infty
```

This is a **second-order linear ODE with constant coefficients** and a non-homogeneous term $`-h' T_\infty`$ on the right-hand side.

The total solution is the sum of:
1. The **homogeneous** solution
2. A **particular** solution

### 2.2 Homogeneous Solution

First, solve the homogeneous equation:

```math
\frac{d^2 T}{dx^2} - h' T = 0 \quad \cdots (***)
```

Assume a trial solution of the form $`T = C e^{\lambda x}`$.

Substituting into $`(***)`$:

```math
\lambda^2 C e^{\lambda x} - h' C e^{\lambda x} = (\lambda^2 - h') C e^{\lambda x} = 0
```

Since we seek a non-trivial solution ($`C \neq 0`$) and $`e^{\lambda x}`$ is never zero, the characteristic equation must be satisfied:

```math
\lambda^2 - h' = 0
```

```math
\lambda = \pm \sqrt{h'}
```

Thus the **homogeneous solution** is:

```math
T_h = C_1 e^{\lambda x} + C_2 e^{-\lambda x}, \quad \lambda = \sqrt{h'}
```

> **[Calculus]** The characteristic equation $`\lambda^2 - h' = 0`$ has two distinct real roots when $`h' > 0`$ (which is always the case for physical heat transfer). This gives two linearly independent solutions $`e^{\lambda x}`$ and $`e^{-\lambda x}`$, and their linear combination forms the complete homogeneous solution.

### 2.3 Particular Solution

A particular solution is any one solution that satisfies the full non-homogeneous equation. Since the right-hand side is constant, a constant function is a natural trial, and if $`T = T_\infty`$:

```math
\frac{d^2 T_\infty}{dx^2} - h' T_\infty = 0 - h' T_\infty = -h' T_\infty \quad \checkmark
```

Therefore the particular solution is:

```math
T_p = T_\infty
```

### 2.4 Total (General) Solution

The total solution is the sum of the homogeneous and particular solutions:

```math
T(x) = T_\infty + C_1 e^{\lambda x} + C_2 e^{-\lambda x}, \quad \lambda = \sqrt{h'}
```

### 2.5 Determining the Constants from Boundary Conditions

To determine $`C_1`$ and $`C_2`$, apply the two boundary conditions:

**At $`x = 0`$:**

```math
T_a = T(0) = T_\infty + C_1 + C_2
```

**At $`x = L`$:**

```math
T_b = T(L) = T_\infty + C_1 e^{\lambda L} + C_2 e^{-\lambda L}
```

This is a system of two equations in two unknowns ($`C_1`$, $`C_2`$).

**Solving for $`C_2`$:** Multiply the first equation by $`e^{\lambda L}`$ and subtract the second:

```math
e^{\lambda L} T_a - T_b = T_\infty (e^{\lambda L} - 1) + C_2 (e^{\lambda L} - e^{-\lambda L})
```

```math
\therefore \; C_2 = \frac{e^{\lambda L}(T_a - T_\infty) + (T_\infty - T_b)}{e^{\lambda L} - e^{-\lambda L}}
```

**Solving for $`C_1`$:** Multiply the first equation by $`e^{-\lambda L}`$ and subtract the second:

```math
e^{-\lambda L} T_a - T_b = e^{-\lambda L} T_\infty - T_\infty + C_1 (e^{-\lambda L} - e^{\lambda L})
```

```math
\therefore \; C_1 = \frac{e^{-\lambda L}(T_a - T_\infty) + (T_\infty - T_b)}{e^{-\lambda L} - e^{\lambda L}}
```

> **[Calculus]** The denominator $`e^{\lambda L} - e^{-\lambda L} = 2\sinh(\lambda L)`$ is always non-zero for $`\lambda > 0`$ and $`L > 0`$, guaranteeing a unique solution to the BVP. This confirms the physical expectation that the temperature distribution is uniquely determined by the boundary temperatures.

---

<br>

## 3. The Shooting Method (Section 24.2)

### 3.1 Core Idea: Convert BVP into IVP

The shooting method transforms a boundary-value problem into an initial-value problem by:

1. **Using** the known boundary condition at one end as the initial condition.
2. **Guessing** the missing initial condition (the derivative at that end).
3. **Integrating** the ODE from one boundary to the other using any IVP solver (e.g., Euler, RK4).
4. **Comparing** the computed result at the far boundary with the known boundary value.
5. **Adjusting** the guess and repeating until the far boundary condition is satisfied.

The name "shooting method" comes from the analogy of firing a projectile: you adjust the launch angle (the initial slope guess) until the projectile hits the target (the far boundary condition).

### 3.2 Converting a 2nd-Order ODE to Two 1st-Order ODEs

Starting from the governing equation:

```math
0 = \frac{d^2 T}{dx^2} + h'(T_\infty - T)
```

Introduce a new variable $`z`$:

```math
z := \frac{dT}{dx}
```

Then the single second-order ODE becomes a **system of two first-order ODEs**:

```math
\begin{cases} \dfrac{dT}{dx} = z \\[10pt] \dfrac{dz}{dx} = h'(T - T_\infty) \end{cases}
```

> **[Calculus]** This is a general technique: any $`n`$-th order ODE can be converted into a system of $`n`$ first-order ODEs by introducing auxiliary variables for each derivative up to order $`n-1`$. This is essential because most numerical ODE solvers (Euler, Runge-Kutta, etc.) are designed for first-order systems.

### 3.3 The Guessing Procedure

If we know both $`z(x = 0)`$ and $`T(x = 0)`$, we can integrate the system of ODEs from $`x = 0`$ to $`x = L`$.

However, in a BVP:
- $`T(x = 0) = T_a`$ is **known** (boundary condition).
- $`z(x = 0) = \frac{dT}{dx}\big|_{x=0}`$ is **unknown**.

**Procedure:**

1. **Guess** the initial slope as $`z(x = 0) = z_{a,1}`$ (first guess).
2. **Integrate** the pair of first-order ODEs from $`x = 0`$ to $`x = L`$ using a numerical method (e.g., RK4).
3. Obtain the computed temperature at the right boundary, $`T(x = L) = T_{b,1}`$.
4. **Compare** $`T_{b,1}`$ with the known $`T_b`$.
5. In general, $`T_{b,1} \neq T_b`$.
6. **Adjust** the guess to $`z(x = 0) = z_{a,2}`$ (new guess).
7. **Repeat** the integration and comparison until $`T(L) \approx T_b`$ within a desired tolerance.

The adjustment of the guess can be done systematically using root-finding methods:

- **Bisection method**: Bracket the correct slope between two guesses and keep narrowing the bracket.
- **Secant method / linear interpolation**: Use two guesses and their results to interpolate the next guess.

```math
z_{a,\text{new}} = z_{a,1} + \frac{(T_b - T_{b,1})(z_{a,2} - z_{a,1})}{T_{b,2} - T_{b,1}}
```

> **Note:** When the ODE is linear, as in this example, $`T(L)`$ is a linear function of $`z(0)`$. The $`z_{a,\text{new}}`$ obtained by linear interpolation of two trial shots is therefore the exact initial slope apart from the numerical error of the integration, and it satisfies the boundary condition without further iteration. With the example values of Section 4.7, the implementation in Section 3.5 converges at the first interpolation. For a nonlinear ODE this relation does not hold, so the interpolation must be repeated.

### 3.4 Algorithm Summary

```text
1. Set T(0) = Ta  (known boundary condition)
2. Guess z(0) = z_a1
3. Integrate the system of ODEs from x = 0 to x = L
4. Compute T(L) = T_b1
5. If |T_b1 - Tb| < tolerance: DONE
6. Else: adjust z(0) using root-finding, go to step 3
```

### 3.5 Python Implementation

```python
import numpy as np

def shooting_method(h_prime, T_inf, Ta, Tb, L, n_steps=100, tol=1e-6, max_iter=50):
    """
    Solve the BVP: d^2T/dx^2 + h'(T_inf - T) = 0
    with T(0) = Ta, T(L) = Tb using the shooting method.

    Parameters
    ----------
    h_prime : float  - heat transfer parameter
    T_inf   : float  - ambient temperature
    Ta      : float  - temperature at x = 0
    Tb      : float  - temperature at x = L
    L       : float  - length of the rod
    n_steps : int    - number of integration steps
    tol     : float  - convergence tolerance
    max_iter: int    - maximum number of shooting iterations
    """
    dx = L / n_steps

    def integrate(z0):
        """Integrate the system using RK4 from x=0 to x=L."""
        T = Ta
        z = z0
        T_vals = [T]
        x_vals = [0.0]

        for i in range(n_steps):
            x = i * dx

            # RK4 for the system dT/dx = z, dz/dx = h'(T - T_inf)
            k1_T = z
            k1_z = h_prime * (T - T_inf)

            k2_T = z + 0.5 * dx * k1_z
            k2_z = h_prime * ((T + 0.5 * dx * k1_T) - T_inf)

            k3_T = z + 0.5 * dx * k2_z
            k3_z = h_prime * ((T + 0.5 * dx * k2_T) - T_inf)

            k4_T = z + dx * k3_z
            k4_z = h_prime * ((T + dx * k3_T) - T_inf)

            T = T + (dx / 6) * (k1_T + 2*k2_T + 2*k3_T + k4_T)
            z = z + (dx / 6) * (k1_z + 2*k2_z + 2*k3_z + k4_z)

            T_vals.append(T)
            x_vals.append((i + 1) * dx)

        return np.array(x_vals), np.array(T_vals), T  # T at x = L

    # Two initial guesses for z(0)
    z1 = 0.0
    _, _, Tb1 = integrate(z1)

    z2 = 1.0
    _, _, Tb2 = integrate(z2)

    # Secant method to find the correct z(0)
    for iteration in range(max_iter):
        if abs(Tb2 - Tb1) < 1e-15:
            break

        # Linear interpolation for next guess
        z_new = z1 + (Tb - Tb1) * (z2 - z1) / (Tb2 - Tb1)
        x_vals, T_vals, Tb_new = integrate(z_new)

        if abs(Tb_new - Tb) < tol:
            print(f"Converged in {iteration + 1} iterations, z(0) = {z_new:.6f}")
            return x_vals, T_vals

        # Update for next iteration
        z1, Tb1 = z2, Tb2
        z2, Tb2 = z_new, Tb_new

    print(f"Warning: did not converge after {max_iter} iterations")
    return x_vals, T_vals
```

---

<br>

## 4. Finite Difference Methods (Section 24.3)

### 4.1 Core Idea: Discretize the Domain

Instead of converting the BVP into an IVP, the finite difference (FD) method **directly discretizes** the governing equation at a set of grid points, producing a system of algebraic equations that can be solved simultaneously.

Starting from:

```math
\frac{d^2 T}{dx^2} - h' T = -h' T_\infty
```

Note that the right-hand side $`-h' T_\infty`$ is a **known value**.

### 4.2 Second Derivative Approximation

Divide the domain $`[0, L]`$ into $`n`$ equal subintervals with spacing $`\Delta x = L / n`$, creating nodes $`x_0, x_1, \ldots, x_n`$.

```text
T_0 = T_a                                         T_n = T_b
  ●───────○───────○── ... ──○───────○───────○── ... ──●
 x_0     x_1     x_2     x_{i-1}   x_i   x_{i+1}     x_n
  ├── Δx ─┤
 ● Boundary node (known value)   ○ Interior node (unknown)
```

*Figure 2. Boundary and interior nodes of the finite difference grid*

The **central difference approximation** for the second derivative at node $`x_i`$ is:

```math
\left. \frac{d^2 T}{dx^2} \right|_{x_i} \approx \frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2}
```

> **Note:** This approximation is the central difference formula for the second derivative derived in Section 3.4 of Chapter 4 and Section 3.1 of Chapter 21, obtained by adding the forward and backward Taylor expansions around $`x_i`$. The leading error term is $`O((\Delta x)^2)`$, meaning the approximation has second-order accuracy.

### 4.3 Generating the Finite Difference Equations

Substituting the central difference approximation into the governing equation at node $`x_i`$:

```math
\frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2} - h' T_i = -h' T_\infty
```

This equation holds for each **interior** node $`i = 1, 2, \ldots, n-1`$.

Combined with the boundary conditions:
- $`T_0 = T_a`$ (left boundary)
- $`T_n = T_b`$ (right boundary)

This gives a total of $`(n+1)`$ equations for $`(n+1)`$ unknowns $`T_0, T_1, \ldots, T_n`$.

### 4.4 Counting Unknowns and Equations

| Item | Count |
|:-----|:------|
| Total unknowns: $`T_0, T_1, \ldots, T_n`$ | $`n + 1`$ |
| Boundary conditions: $`T_0 = T_a`$, $`T_n = T_b`$ | $`2`$ |
| Interior FD equations (at $`x_1, x_2, \ldots, x_{n-1}`$) | $`n - 1`$ |
| **Total equations** | $`n + 1`$ |

*Table 3. Number of unknowns and equations*

Since $`T_0`$ and $`T_n`$ are known from the boundary conditions, the system reduces to $`(n - 1)`$ equations in $`(n - 1)`$ unknowns $`T_1, T_2, \ldots, T_{n-1}`$.

### 4.5 Assembling the Tridiagonal System

Define the following constants for convenience:

```math
G = 2 + (\Delta x)^2 h', \quad F = h' T_\infty (\Delta x)^2
```

Multiply the FD equation at each node by $`(\Delta x)^2`$ and rearrange:

**i) First interior node ($`i = 1`$):**

```math
T_0 - 2T_1 + T_2 - (\Delta x)^2 h' T_1 = -h' T_\infty (\Delta x)^2
```

Since $`T_0`$ is known, move it to the right-hand side:

```math
(-2 - (\Delta x)^2 h') T_1 + T_2 = -h' T_\infty (\Delta x)^2 - T_0
```

```math
-G \cdot T_1 + T_2 = -F - T_0
```

**ii) General interior nodes ($`i = 2, 3, \ldots, n-2`$):**

```math
T_{i-1} - 2T_i + T_{i+1} - (\Delta x)^2 h' T_i = -h' T_\infty (\Delta x)^2
```

```math
T_{i-1} - (2 + (\Delta x)^2 h') T_i + T_{i+1} = -h' T_\infty (\Delta x)^2
```

```math
T_{i-1} - G \cdot T_i + T_{i+1} = -F
```

**iii) Last interior node ($`i = n-1`$):**

```math
T_{n-2} - 2T_{n-1} + T_n - (\Delta x)^2 h' T_{n-1} = -h' T_\infty (\Delta x)^2
```

Since $`T_n`$ is known, move it to the right-hand side:

```math
T_{n-2} - (2 + (\Delta x)^2 h') T_{n-1} = -h' T_\infty (\Delta x)^2 - T_n
```

```math
T_{n-2} - G \cdot T_{n-1} = -F - T_n
```

If $`n = 2`$, there is only one interior node $`T_1`$, so i) and iii) become the same equation, $`-G \cdot T_1 = -F - T_0 - T_2`$. The general interior equation ii) appears only when $`n \geq 4`$.

### 4.6 Matrix Formulation

Combining all the interior equations into matrix form:

```math
\begin{bmatrix} -G & 1 & & & \\ 1 & -G & 1 & & \\ & 1 & -G & 1 & \\ & & \ddots & \ddots & \ddots \\ & & & 1 & -G \end{bmatrix} \begin{bmatrix} T_1 \\ T_2 \\ T_3 \\ \vdots \\ T_{n-1} \end{bmatrix} = \begin{bmatrix} -F \\ -F \\ -F \\ \vdots \\ -F \end{bmatrix} - \begin{bmatrix} T_0 \\ 0 \\ 0 \\ \vdots \\ T_n \end{bmatrix}
```

This is a **tridiagonal system** of $`(n-1)`$ equations that can be solved using:
- **Thomas algorithm** (tridiagonal matrix algorithm), which solves it efficiently in $`O(n)`$ time.
- **NumPy's `linalg.solve`** or any standard linear algebra solver, which is simple to use but treats the matrix as dense and takes $`O(n^3)`$ time.

> **Note:** The coefficient matrix is tridiagonal, symmetric, and strictly diagonally dominant. For $`h' > 0`$, the magnitude of each diagonal entry, $`G = 2 + (\Delta x)^2 h'`$, exceeds 2, while each row has at most two off-diagonal entries, each equal to 1. A strictly diagonally dominant matrix is nonsingular, so the system has a unique solution, and by the convergence condition of Section 3 of Chapter 12, iterative methods like Gauss-Seidel converge.

### 4.7 Python Implementation

```python
import numpy as np

def finite_difference_bvp(h_prime, T_inf, Ta, Tb, L, n=10):
    """
    Solve the BVP: d^2T/dx^2 + h'(T_inf - T) = 0
    with T(0) = Ta, T(L) = Tb using finite differences.

    Parameters
    ----------
    h_prime : float  - heat transfer parameter
    T_inf   : float  - ambient temperature
    Ta      : float  - temperature at x = 0 (T_0)
    Tb      : float  - temperature at x = L (T_n)
    L       : float  - length of the rod
    n       : int    - number of subintervals
    """
    dx = L / n
    G = 2 + dx**2 * h_prime       # diagonal coefficient
    F = h_prime * T_inf * dx**2   # forcing term

    # Build the (n-1) x (n-1) tridiagonal coefficient matrix
    size = n - 1
    A = np.zeros((size, size))
    b = np.full(size, -F)

    for i in range(size):
        A[i, i] = -G                          # main diagonal
        if i > 0:
            A[i, i - 1] = 1.0                 # lower diagonal
        if i < size - 1:
            A[i, i + 1] = 1.0                 # upper diagonal

    # Adjust RHS for boundary conditions
    b[0] -= Ta        # T_0 = Ta is known
    b[-1] -= Tb       # T_n = Tb is known

    # Solve the tridiagonal system
    T_interior = np.linalg.solve(A, b)

    # Assemble full solution including boundary nodes
    T_full = np.concatenate(([Ta], T_interior, [Tb]))
    x_full = np.linspace(0, L, n + 1)

    return x_full, T_full
```

**Example usage:** comparing both methods with the analytical solution.

```python
import matplotlib.pyplot as plt

# Problem parameters
h_prime = 0.01    # heat transfer parameter
T_inf = 20.0      # ambient temperature
Ta = 40.0         # left boundary temperature
Tb = 200.0        # right boundary temperature
L = 10.0          # rod length

# Analytical solution
lam = np.sqrt(h_prime)
C1 = (np.exp(-lam*L) * (Ta - T_inf) + (T_inf - Tb)) / (np.exp(-lam*L) - np.exp(lam*L))
C2 = (np.exp(lam*L) * (Ta - T_inf) + (T_inf - Tb)) / (np.exp(lam*L) - np.exp(-lam*L))

x_exact = np.linspace(0, L, 200)
T_exact = T_inf + C1 * np.exp(lam * x_exact) + C2 * np.exp(-lam * x_exact)

# Shooting method
x_shoot, T_shoot = shooting_method(h_prime, T_inf, Ta, Tb, L, n_steps=100)

# Finite difference method
x_fd, T_fd = finite_difference_bvp(h_prime, T_inf, Ta, Tb, L, n=10)

# Plot comparison
plt.figure(figsize=(10, 6))
plt.plot(x_exact, T_exact, 'k-', label='Analytical', linewidth=2)
plt.plot(x_shoot, T_shoot, 'b--o', label='Shooting method', markersize=4)
plt.plot(x_fd, T_fd, 'r--s', label='Finite difference', markersize=6)
plt.xlabel('x')
plt.ylabel('T(x)')
plt.title('BVP Solution: Heated Rod')
plt.legend()
plt.grid(True)
plt.show()
```

Running this example gives the analytical constants $`C_1 \approx 73.4523`$ and $`C_2 \approx -53.4523`$. The shooting method converges at the first interpolation with $`z(0) \approx 12.6905`$, and the finite difference solution with $`n = 10`$ gives $`108.6904`$ at $`x = 5`$ (analytical $`108.6819`$), with an error below $`0.009`$ at every node.

---

<br>

## Summary

| Concept | Key Summary |
|:------|:-----------|
| Boundary-Value Problem | An ODE problem with conditions specified at two or more distinct points, unlike an IVP where all conditions are given at one point. |
| Heated Rod BVP | The governing equation is $`\dfrac{d^2 T}{dx^2} + h'(T_\infty - T) = 0`$, with boundary conditions $`T(0) = T_a`$ and $`T(L) = T_b`$. |
| Analytical Solution | The solution is $`T(x) = T_\infty + C_1 e^{\lambda x} + C_2 e^{-\lambda x}`$ with $`\lambda = \sqrt{h'}`$, and the constants $`C_1`$, $`C_2`$ are determined from the boundary conditions. |
| Shooting Method | It converts the BVP to an IVP by guessing the unknown initial slope $`z(0) = dT/dx\vert_{x=0}`$, integrates forward, and adjusts the guess by root-finding until $`T(L) = T_b`$. |
| ODE Conversion | A single $`n`$-th order ODE can be converted to a system of $`n`$ first-order ODEs by introducing auxiliary variables. |
| Finite Difference Method | It discretizes the domain into a grid, approximates $`d^2T/dx^2 \approx \frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2}`$, and solves the resulting tridiagonal linear system. |
| FD Matrix Structure | The matrix is tridiagonal with diagonal entries $`-G = -(2 + (\Delta x)^2 h')`$ and off-diagonal entries $`1`$, and the Thomas algorithm solves it in $`O(n)`$. |
| Key Difference | The shooting method iterates by guess and check, while the FD method builds an algebraic system directly and solves a linear BVP in one step. |

---

<br>

## Self-Check Questions

1. **BVP and IVP:** Explain why a BVP cannot be integrated forward from one end like an IVP, and describe how the shooting method resolves this.

   <details>
   <summary>Answer</summary>

   Integrating a second-order ODE forward requires both $`T`$ and $`dT/dx`$ at the starting point. In a BVP, $`T`$ is given at $`x = 0`$ and at $`x = L`$, so the starting slope $`z(0)`$ is unknown. The shooting method guesses $`z(0)`$, integrates the resulting IVP, and adjusts the guess with a root-finding method until the difference between the computed $`T(L)`$ and the given $`T_b`$ is within the tolerance.

   </details>

2. **Constants of the analytical solution:** For $`h' = 0.01`$, $`L = 10`$, $`T_\infty = 20`$, $`T_a = 40`$, and $`T_b = 200`$, find $`\lambda`$, $`C_1`$, and $`C_2`$, and compute $`T(5)`$.

   <details>
   <summary>Answer</summary>

   $`\lambda = \sqrt{0.01} = 0.1`$ and $`\lambda L = 1`$. Substituting $`e \approx 2.718282`$ and $`e^{-1} \approx 0.367879`$ gives $`C_1 = (0.367879 \times 20 - 180)/(0.367879 - 2.718282) = (-172.6424)/(-2.350402) \approx 73.4523`$ and $`C_2 = (2.718282 \times 20 - 180)/(2.718282 - 0.367879) = (-125.6344)/2.350402 \approx -53.4523`$. As a check, $`C_1 + C_2 = 20 = T_a - T_\infty`$. Therefore $`T(5) = 20 + 73.4523 e^{0.5} - 53.4523 e^{-0.5} \approx 20 + 121.1023 - 32.4204 = 108.6819`$.

   </details>

3. **Linear interpolation in the shooting method:** Under the conditions of Question 2, integrating with $`z(0) = 10`$ gives $`T(10) \approx 168.3817`$, and integrating with $`z(0) = 20`$ gives $`T(10) \approx 285.9019`$. Find the next guess $`z_{a,\text{new}}`$, and explain why it matches the value of $`dT/dx`$ at $`x = 0`$ from the analytical solution of Question 2.

   <details>
   <summary>Answer</summary>

   $`z_{a,\text{new}} = 10 + (200 - 168.3817)(20 - 10)/(285.9019 - 168.3817) = 10 + 316.183/117.5202 \approx 12.6905`$. From the analytical solution, $`dT/dx\big|_{x=0} = \lambda(C_1 - C_2) = 0.1 \times (73.4523 + 53.4523) \approx 12.6905`$, so the two values agree. Because this ODE is linear, $`T(L) = T_\infty + (T_a - T_\infty)\cosh(\lambda L) + (z(0)/\lambda)\sinh(\lambda L)`$ is a linear function of $`z(0)`$. Interpolating along the straight line through the two trial points therefore gives the exact initial slope.

   </details>

4. **Finite difference calculation:** Under the conditions of Question 2, apply the finite difference method with $`n = 2`$ ($`\Delta x = 5`$) to find $`T_1 \approx T(5)`$, and compare it with the analytical solution and with the result for $`n = 10`$ ($`108.6904`$).

   <details>
   <summary>Answer</summary>

   $`G = 2 + 25 \times 0.01 = 2.25`$ and $`F = 0.01 \times 20 \times 25 = 5`$. With only one interior node, $`-2.25 T_1 = -5 - 40 - 200 = -245`$, so $`T_1 = 245/2.25 \approx 108.8889`$. The error against the analytical value $`108.6819`$ is about $`0.2070`$, while the error for $`n = 10`$ is about $`0.0085`$. When $`\Delta x`$ shrinks from 5 to 1 (a factor of 1/5), the error drops by about 1/24, which nearly matches the factor 1/25 expected from second-order accuracy, where the error is proportional to $`(\Delta x)^2`$.

   </details>

5. **Properties of the coefficient matrix:** For $`n = 10`$, state the size of the linear system that the finite difference method solves, and explain why the system of Section 4.6 has a unique solution and why Gauss-Seidel iteration converges.

   <details>
   <summary>Answer</summary>

   $`T_0`$ and $`T_{10}`$ are known from the boundary conditions, so the unknowns are $`T_1, \ldots, T_9`$ and the system is $`9 \times 9`$. In each row, the magnitude of the diagonal entry is $`G = 2 + (\Delta x)^2 h'`$, which exceeds 2 when $`h' > 0`$, while the magnitudes of the off-diagonal entries sum to at most 2. The matrix is therefore strictly diagonally dominant, and a strictly diagonally dominant matrix is nonsingular, so the solution is unique and Gauss-Seidel iteration converges.

   </details>
