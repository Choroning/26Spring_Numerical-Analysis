# Chapter 1 Lecture — Mathematical Modeling and Engineering Problem Solving

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. What is a Mathematical Model?](#1-what-is-a-mathematical-model)
  - [1.1 Definition and General Form](#11-definition-and-general-form)
  - [1.2 Analytical vs. Numerical Solutions](#12-analytical-vs-numerical-solutions)
- [2. Newton's 2nd Law — Bungee Jumper Example](#2-newtons-2nd-law--bungee-jumper-example)
  - [2.1 System Setup and Force Balance](#21-system-setup-and-force-balance)
  - [2.2 Deriving the ODE](#22-deriving-the-ode)
  - [2.3 Analytical Solution and Terminal Velocity](#23-analytical-solution-and-terminal-velocity)
- [3. Euler's Method — Numerical Solution](#3-eulers-method--numerical-solution)
  - [3.1 Core Idea: Finite Difference Approximation](#31-core-idea-finite-difference-approximation)
  - [3.2 Step-by-Step Calculation](#32-step-by-step-calculation)
  - [3.3 Effect of Step Size](#33-effect-of-step-size)
- [4. Types of Engineering Problems Solved Numerically](#4-types-of-engineering-problems-solved-numerically)
- [Summary](#summary)

---

<br>

## 1. What is a Mathematical Model?

### 1.1 Definition and General Form

A **mathematical model** is a formulation that expresses the features of a physical system using mathematical equations. Rather than building a physical prototype or relying solely on experiments, we translate a real-world phenomenon into the language of mathematics so that we can analyze, predict, and optimize its behavior.

The general form of a mathematical model is:

$$\text{Dependent variable} = f(\text{independent variables},\ \text{parameters},\ \text{forcing functions})$$

- **Dependent variable**: the quantity we want to predict or understand (e.g., velocity, temperature, concentration)
- **Independent variables**: dimensions over which the system evolves (e.g., time $t$, spatial coordinates $x, y, z$)
- **Parameters**: properties of the system that are typically constant for a given scenario (e.g., mass $m$, drag coefficient $c$)
- **Forcing functions**: external influences acting on the system (e.g., gravity, applied voltage, wind)

### 1.2 Analytical vs. Numerical Solutions

Many real-world engineering and science problems cannot be solved analytically — that is, we cannot always find a closed-form expression for the solution. This is where **numerical methods** become essential: they provide approximate solutions using computers.

| Aspect | Analytical Solution | Numerical Solution |
|:-------|:-------------------|:-------------------|
| Form | Exact closed-form expression | Approximate discrete values |
| Applicability | Simple, idealized problems | Complex, realistic problems |
| Tools | Pen, paper, symbolic math | Computers, algorithms |
| Example | $v(t) = \frac{gm}{c}(1 - e^{-(c/m)t})$ | Table of $v$ values at discrete $t_i$ |

> **Note:** Numerical methods do not replace analytical methods — they complement them. When an analytical solution exists, it serves as a benchmark to verify numerical results. When it does not exist, numerical methods are often the only practical approach.

---

<br>

## 2. Newton's 2nd Law — Bungee Jumper Example

### 2.1 System Setup and Force Balance

Consider a **bungee jumper** falling through the air. We want to model how the jumper's velocity changes over time during free fall (before the cord engages).

Two forces act on the jumper:

1. **Gravity (downward)**: $F_D = mg$, where $m$ is the mass of the jumper and $g$ is gravitational acceleration ($9.81\ \text{m/s}^2$)
2. **Air resistance / Drag (upward)**: opposes the motion
   - **Linear drag model**: $F_U = cv$, where $c$ is the drag coefficient (kg/s) and $v$ is velocity
   - **Quadratic drag model**: $F_U = c_d\, v|v|$, where $c_d$ is the drag coefficient for quadratic resistance

The **net force** on the jumper is:

$$F_{\text{net}} = F_{\text{down}} - F_{\text{up}} = mg - cv$$

(using the linear drag model)

### 2.2 Deriving the ODE

Applying **Newton's 2nd law** ($F = ma$):

$$ma = mg - cv$$

Since acceleration is the time derivative of velocity ($a = dv/dt$):

$$m\frac{dv}{dt} = mg - cv$$

Dividing both sides by $m$:

$$\frac{dv}{dt} = g - \frac{c}{m}v$$

This is a **first-order linear ordinary differential equation (ODE)** that describes the velocity of the jumper as a function of time.

> **[Physics]** Newton's 2nd law: $F = ma$. The net force on the jumper is the difference between gravitational pull (downward) and air resistance/drag (upward). At terminal velocity, these forces balance and acceleration becomes zero.

### 2.3 Analytical Solution and Terminal Velocity

The ODE $\frac{dv}{dt} = g - \frac{c}{m}v$ can be solved analytically with the initial condition $v(0) = 0$ (starting from rest):

$$v(t) = \frac{gm}{c}\left(1 - e^{-(c/m)t}\right)$$

**Terminal velocity** is found by considering the long-time behavior. As $t \to \infty$, the exponential term $e^{-(c/m)t} \to 0$, so:

$$v_{\text{terminal}} = \frac{gm}{c}$$

Equivalently, terminal velocity occurs when acceleration equals zero:

$$\frac{dv}{dt} = 0 \implies g - \frac{c}{m}v = 0 \implies v = \frac{gm}{c}$$

> **[Calculus]** This is a first-order linear ODE. The analytical solution is obtained by separation of variables and integration. The exponential decay term $e^{-(c/m)t}$ shows how the velocity approaches the terminal value asymptotically.

---

<br>

## 3. Euler's Method — Numerical Solution

### 3.1 Core Idea: Finite Difference Approximation

**Euler's method** is the simplest numerical technique for solving ordinary differential equations. The core idea is to approximate the continuous derivative with a **finite difference**:

$$\frac{dv}{dt}\bigg|_{t_i} \approx \frac{v(t_{i+1}) - v(t_i)}{h}$$

where $h = \Delta t$ is the **step size** (time increment between successive points).

Rearranging, we obtain the **Euler formula**:

$$v(t_{i+1}) = v(t_i) + \frac{dv}{dt}\bigg|_{t_i} \cdot h$$

This says: the new value equals the old value plus the slope at the current point times the step size. Geometrically, we are following the tangent line at each point to estimate the next value.

> **[Calculus]** Euler's method replaces the continuous derivative $dv/dt$ with a discrete approximation: $\frac{\Delta v}{\Delta t} = \frac{v(t_{i+1}) - v(t_i)}{h}$. This is the simplest numerical ODE solver, connecting to the formal definition of a derivative as $h \to 0$.

### 3.2 Step-by-Step Calculation

Using the bungee jumper ODE $\frac{dv}{dt} = g - \frac{c}{m}v$ with:

- $m = 68.1$ kg
- $c = 12.5$ kg/s
- $g = 9.81$ m/s$^2$
- Initial condition: $v(0) = 0$ m/s
- Step size: $h = 2$ s

**Step 1** — From $t = 0$ to $t = 2$:

$$\frac{dv}{dt}\bigg|_{t=0} = 9.81 - \frac{12.5}{68.1}(0) = 9.81\ \text{m/s}^2$$

$$v(2) = 0 + 9.81 \times 2 = 19.62\ \text{m/s}$$

**Step 2** — From $t = 2$ to $t = 4$:

$$\frac{dv}{dt}\bigg|_{t=2} = 9.81 - \frac{12.5}{68.1}(19.62) = 9.81 - 3.60 = 6.21\ \text{m/s}^2$$

$$v(4) = 19.62 + 6.21 \times 2 = 32.04\ \text{m/s}$$

**Step 3** — From $t = 4$ to $t = 6$:

$$\frac{dv}{dt}\bigg|_{t=4} = 9.81 - \frac{12.5}{68.1}(32.04) = 9.81 - 5.88 = 3.93\ \text{m/s}^2$$

$$v(6) = 32.04 + 3.93 \times 2 = 39.90\ \text{m/s}$$

Continuing this process:

| $t$ (s) | $v_{\text{numerical}}$ (m/s) | $v_{\text{analytical}}$ (m/s) |
|:--------|:----------------------------|:-----------------------------|
| 0 | 0.00 | 0.00 |
| 2 | 19.62 | 16.40 |
| 4 | 32.04 | 27.77 |
| 6 | 39.90 | 35.64 |
| 8 | 44.87 | 41.10 |
| 10 | 48.02 | 44.87 |
| 12 | 50.01 | 47.49 |

At $t = 12$ s, the numerical result ($50.01$ m/s) is reasonably close to the analytical result, but there is a noticeable discrepancy due to the relatively large step size.

### 3.3 Effect of Step Size

The accuracy of Euler's method depends heavily on the **step size** $h$:

| Step Size ($h$) | $v(12)$ Numerical | Error vs. Analytical |
|:----------------|:------------------|:---------------------|
| 2.0 s | ~50.01 m/s | Larger error |
| 1.0 s | Closer to analytical | Moderate error |
| 0.5 s | Even closer | Smaller error |

**Key insight**: Reducing the step size $h$ improves accuracy because the finite difference approximation becomes closer to the true derivative. However, smaller $h$ means more computation steps — this is the fundamental **trade-off between accuracy and computational cost** in numerical methods.

As $h \to 0$, the numerical solution converges to the exact analytical solution.

---

<br>

## 4. Types of Engineering Problems Solved Numerically

Numerical methods provide tools for solving a wide range of engineering and scientific problems. The major categories include:

1. **Roots of Equations**: Finding $x$ such that $f(x) = 0$
   - Example: determining the drag coefficient that produces a desired terminal velocity

2. **Systems of Linear Algebraic Equations**: Solving $[A]\{x\} = \{b\}$
   - Example: structural analysis, circuit analysis, mass balance in chemical systems

3. **Curve Fitting / Regression**: Finding a function that best fits observed data
   - Example: fitting experimental measurements to a theoretical model

4. **Integration**: Computing definite integrals $\int_a^b f(x)\,dx$
   - Example: calculating work, area, or total accumulated quantity

5. **Ordinary Differential Equations (ODEs)**: Solving equations involving derivatives of one independent variable
   - Example: the bungee jumper problem, population dynamics, radioactive decay

6. **Partial Differential Equations (PDEs)**: Solving equations involving partial derivatives of multiple independent variables
   - Example: heat conduction, fluid flow, electromagnetic wave propagation

7. **Optimization**: Finding parameter values that minimize or maximize an objective function
   - Example: minimizing cost, maximizing efficiency, optimal design

> **Note:** This course progresses through these topics systematically, building from simple root-finding techniques to more advanced PDE and optimization methods. Each topic relies on the foundations established in earlier chapters.

---

<br>

## Summary

| Topic | Key Point |
|:------|:----------|
| Mathematical Model | Dependent variable = $f$(independent variables, parameters, forcing functions) |
| Analytical vs. Numerical | Analytical = exact closed-form; Numerical = approximate via computation |
| Newton's 2nd Law (Bungee Jumper) | ODE: $dv/dt = g - (c/m)v$; Terminal velocity: $v_t = gm/c$ |
| Euler's Method | $v(t_{i+1}) = v(t_i) + (dv/dt) \cdot h$; simplest numerical ODE solver |
| Step Size Trade-off | Smaller $h$ = better accuracy but more computation |
| Problem Categories | Roots, linear systems, curve fitting, integration, ODEs, PDEs, optimization |

---
