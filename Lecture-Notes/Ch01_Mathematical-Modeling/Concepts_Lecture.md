# Chapter 1 — Mathematical Modeling and Engineering Problem Solving

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 1

> **Learning Objectives:**
> 1. Explain the role of mathematical modeling in engineering
> 2. Identify sources of error in numerical computations
> 3. Compute an approximate solution of an ordinary differential equation step by step with Euler's method
> 4. Describe the trade-off between accuracy and computational cost that depends on the step size

---

## Table of Contents

- [1. What is a Mathematical Model?](#1-what-is-a-mathematical-model)
  - [1.1 Definition and General Form](#11-definition-and-general-form)
  - [1.2 Analytical vs. Numerical Solutions](#12-analytical-vs-numerical-solutions)
- [2. Newton's 2nd Law: Bungee Jumper Example](#2-newtons-2nd-law-bungee-jumper-example)
  - [2.1 System Setup and Force Balance](#21-system-setup-and-force-balance)
  - [2.2 Deriving the ODE](#22-deriving-the-ode)
  - [2.3 Analytical Solution and Terminal Velocity](#23-analytical-solution-and-terminal-velocity)
- [3. Euler's Method: Numerical Solution](#3-eulers-method-numerical-solution)
  - [3.1 Core Idea: Finite Difference Approximation](#31-core-idea-finite-difference-approximation)
  - [3.2 Step-by-Step Calculation](#32-step-by-step-calculation)
  - [3.3 Effect of Step Size](#33-effect-of-step-size)
- [4. Types of Engineering Problems Solved Numerically](#4-types-of-engineering-problems-solved-numerically)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. What is a Mathematical Model?

### 1.1 Definition and General Form

A **mathematical model** is a formulation that expresses the features of a physical system using mathematical equations. Rather than building a physical prototype or relying solely on experiments, we translate a real-world phenomenon into the language of mathematics so that we can analyze, predict, and optimize its behavior.

The general form of a mathematical model is:

```math
\text{Dependent variable} = f(\text{independent variables},\ \text{parameters},\ \text{forcing functions})
```

- **Dependent variable**: the quantity we want to predict or understand (e.g., velocity, temperature, concentration).
- **Independent variables**: the dimensions over which the system evolves (e.g., time $`t`$, spatial coordinates $`x, y, z`$).
- **Parameters**: properties of the system that are typically constant for a given scenario (e.g., mass $`m`$, drag coefficient $`c`$).
- **Forcing functions**: external influences acting on the system (e.g., gravity, applied voltage, wind).

### 1.2 Analytical vs. Numerical Solutions

Many real-world engineering and science problems cannot be solved analytically. That is, we cannot always find a closed-form expression for the solution. This is where **numerical methods**, which provide approximate solutions using computers, become essential.

| Aspect | Analytical Solution | Numerical Solution |
|:-------|:-------------------|:-------------------|
| Form | Exact closed-form expression | Approximate discrete values |
| Applicability | Simple, idealized problems | Complex, realistic problems |
| Tools | Pen, paper, symbolic math | Computers, algorithms |
| Example | $`v(t) = \frac{gm}{c}(1 - e^{-(c/m)t})`$ | Table of $`v`$ values at discrete $`t_i`$ |

*Table 1. Comparison of analytical and numerical solutions*

> **Note:** Numerical methods do not replace analytical methods; they complement them. When an analytical solution exists, it serves as a benchmark to verify numerical results. When it does not exist, numerical methods are often the only practical approach.

---

<br>

## 2. Newton's 2nd Law: Bungee Jumper Example

### 2.1 System Setup and Force Balance

Consider a **bungee jumper** falling through the air. We want to model how the jumper's velocity changes over time during free fall (before the cord becomes taut).

Two forces act on the jumper:

1. **Gravity (downward)**: $`F_D = mg`$, where $`m`$ is the mass of the jumper and $`g`$ is gravitational acceleration ($`9.81\ \text{m/s}^2`$).
2. **Air resistance or drag (upward)**: it opposes the motion.
   - **Linear drag model**: $`F_U = cv`$, where $`c`$ is the drag coefficient (kg/s) and $`v`$ is velocity.
   - **Quadratic drag model**: $`F_U = c_d\, v|v|`$, where $`c_d`$ is the drag coefficient for quadratic resistance (kg/m).

```text
              ▲  F_U = cv   (air resistance, upward)
              │
           ┌──┴──┐
           │  m  │          jumper
           └──┬──┘
              │
              ▼  F_D = mg   (gravity, downward)
```

*Figure 1. Gravity and air resistance acting on the bungee jumper*

The **net force** on the jumper is the downward force $`F_D`$ minus the upward force $`F_U`$. With the linear drag model, it is:

```math
F_{\text{net}} = F_D - F_U = mg - cv
```

The subscripts D and U stand for the downward and upward directions.

### 2.2 Deriving the ODE

Applying **Newton's 2nd law** ($`F = ma`$) gives:

```math
ma = mg - cv
```

Since acceleration is the time derivative of velocity ($`a = dv/dt`$), we can write:

```math
m\frac{dv}{dt} = mg - cv
```

Dividing both sides by $`m`$ gives:

```math
\frac{dv}{dt} = g - \frac{c}{m}v
```

This is a **first-order linear ordinary differential equation (ODE)** that describes the velocity of the jumper as a function of time.

> **Note:** **Physics background.** Newton's 2nd law states $`F = ma`$. The net force on the jumper is the difference between gravity (downward) and air resistance (upward). At terminal velocity, these forces balance and the acceleration becomes zero.

### 2.3 Analytical Solution and Terminal Velocity

Solving the ODE $`\frac{dv}{dt} = g - \frac{c}{m}v`$ analytically with the initial condition $`v(0) = 0`$ (starting from rest) gives:

```math
v(t) = \frac{gm}{c}\left(1 - e^{-(c/m)t}\right)
```

**Terminal velocity** is found by considering the long-time behavior. As $`t \to \infty`$, the exponential term $`e^{-(c/m)t} \to 0`$, so:

```math
v_{\text{terminal}} = \frac{gm}{c}
```

The same result follows from the condition that the acceleration is zero:

```math
\frac{dv}{dt} = 0 \implies g - \frac{c}{m}v = 0 \implies v = \frac{gm}{c}
```

> **[Calculus]** This is a first-order linear ODE. The analytical solution is obtained by separation of variables and integration. The exponential decay term $`e^{-(c/m)t}`$ shows how the velocity approaches the terminal value asymptotically.

---

<br>

## 3. Euler's Method: Numerical Solution

### 3.1 Core Idea: Finite Difference Approximation

**Euler's method** is the simplest numerical technique for solving ordinary differential equations. The core idea is to approximate the continuous derivative with a **finite difference**:

```math
\frac{dv}{dt}\bigg|_{t_i} \approx \frac{v(t_{i+1}) - v(t_i)}{h}
```

where $`h = \Delta t`$ is the **step size**, the time increment between two successive points.

Rearranging, we obtain the **Euler formula**:

```math
v(t_{i+1}) = v(t_i) + \frac{dv}{dt}\bigg|_{t_i} \cdot h
```

This says that the new value equals the old value plus the slope at the current point times the step size. Geometrically, we follow the tangent line at each point to estimate the next value.

> **[Calculus]** Euler's method replaces the continuous derivative $`dv/dt`$ with a discrete approximation: $`\frac{\Delta v}{\Delta t} = \frac{v(t_{i+1}) - v(t_i)}{h}`$. This is the simplest numerical ODE solver, and as $`h \to 0`$ this difference quotient leads to the formal definition of a derivative.

### 3.2 Step-by-Step Calculation

Using the bungee jumper ODE $`\frac{dv}{dt} = g - \frac{c}{m}v`$ with:

- $`m = 68.1`$ kg
- $`c = 12.5`$ kg/s
- $`g = 9.81`$ m/s$`^2`$
- Initial condition: $`v(0) = 0`$ m/s
- Step size: $`h = 2`$ s

**Step 1:** From $`t = 0`$ to $`t = 2`$.

```math
\frac{dv}{dt}\bigg|_{t=0} = 9.81 - \frac{12.5}{68.1}(0) = 9.81\ \text{m/s}^2
```

```math
v(2) = 0 + 9.81 \times 2 = 19.62\ \text{m/s}
```

**Step 2:** From $`t = 2`$ to $`t = 4`$.

```math
\frac{dv}{dt}\bigg|_{t=2} = 9.81 - \frac{12.5}{68.1}(19.62) = 9.81 - 3.60 = 6.21\ \text{m/s}^2
```

```math
v(4) = 19.62 + 6.21 \times 2 = 32.04\ \text{m/s}
```

**Step 3:** From $`t = 4`$ to $`t = 6`$.

```math
\frac{dv}{dt}\bigg|_{t=4} = 9.81 - \frac{12.5}{68.1}(32.04) = 9.81 - 5.88 = 3.93\ \text{m/s}^2
```

```math
v(6) = 32.04 + 3.93 \times 2 = 39.90\ \text{m/s}
```

Continuing this process gives the following table. The analytical column evaluates the formula of Section 2.3 with the same parameters ($`g = 9.81`$ m/s$`^2`$).

| $`t`$ (s) | $`v_{\text{numerical}}`$ (m/s) | $`v_{\text{analytical}}`$ (m/s) |
|:--------|:----------------------------|:-----------------------------|
| 0 | 0.00 | 0.00 |
| 2 | 19.62 | 16.42 |
| 4 | 32.04 | 27.80 |
| 6 | 39.90 | 35.68 |
| 8 | 44.87 | 41.14 |
| 10 | 48.02 | 44.92 |
| 12 | 50.01 | 47.54 |

*Table 2. Euler's method and the analytical solution for h = 2 s*

At $`t = 12`$ s, the numerical result ($`50.01`$ m/s) is reasonably close to the analytical result ($`47.54`$ m/s), but there is a noticeable discrepancy of about $`2.47`$ m/s due to the relatively large step size. The numerical solution is larger because the slope decreases as the velocity grows, yet Euler's method applies the slope from the start of each interval, the largest one in that interval, to the whole interval.

### 3.3 Effect of Step Size

The accuracy of Euler's method depends heavily on the **step size** $`h`$:

| Step Size ($`h`$) | $`v(12)`$ Numerical | Error vs. Analytical |
|:----------------|:------------------|:---------------------|
| 2.0 s | 50.01 m/s | 2.47 m/s (larger error) |
| 1.0 s | 48.76 m/s (closer to analytical) | 1.22 m/s (moderate error) |
| 0.5 s | 48.14 m/s (even closer) | 0.60 m/s (smaller error) |

*Table 3. Numerical solution and error at t = 12 s for different step sizes*

> **Key Point:** Reducing the step size $`h`$ improves accuracy because the finite difference approximation becomes closer to the true derivative. However, smaller $`h`$ means more computation steps. This is the fundamental **trade-off between accuracy and computational cost** in numerical methods.

In Table 3, each halving of $`h`$ roughly halves the error, while the number of steps doubles from 6 to 12 to 24. This is because the global error of Euler's method is proportional to $`h`$ (first-order accuracy).

Ignoring round-off error, the numerical solution converges to the exact analytical solution as $`h \to 0`$. On a real computer, an excessively small $`h`$ lets round-off error accumulate, which is covered in Chapter 3.

---

<br>

## 4. Types of Engineering Problems Solved Numerically

Numerical methods provide tools for solving a wide range of engineering and scientific problems. The major categories include:

1. **Roots of Equations**: Finding $`x`$ such that $`f(x) = 0`$.
   - Example: determining the drag coefficient $`c`$ for which the velocity reaches 40 m/s after 10 s of free fall.

2. **Systems of Linear Algebraic Equations**: Solving $`[A]\{x\} = \{b\}`$.
   - Example: structural analysis, circuit analysis, mass balance in chemical systems

3. **Curve Fitting / Regression**: Finding a function that best fits observed data.
   - Example: fitting experimental measurements to a theoretical model

4. **Integration**: Computing definite integrals $`\int_a^b f(x)\,dx`$.
   - Example: calculating work, area, or total accumulated quantity

5. **Ordinary Differential Equations (ODEs)**: Solving equations involving derivatives with respect to one independent variable.
   - Example: the bungee jumper problem, population dynamics, radioactive decay

6. **Partial Differential Equations (PDEs)**: Solving equations involving partial derivatives with respect to two or more independent variables.
   - Example: heat conduction, fluid flow, electromagnetic wave propagation

7. **Optimization**: Finding parameter values that minimize or maximize an objective function.
   - Example: minimizing cost, maximizing efficiency, optimal design

> **Note:** Of these categories, roots of equations are covered in Chapters 5-6, optimization in Chapter 7, systems of linear algebraic equations in Chapters 8-11, curve fitting in Chapters 14-18, integration in Chapters 19-20, and ODEs in Chapters 22-24. Each topic builds on the foundations established in earlier chapters.

---

<br>

## Summary

| Concept | Key Summary |
|:------|:----------|
| Mathematical Model | A model expresses the dependent variable as a function $`f`$ of independent variables, parameters, and forcing functions. |
| Analytical vs. Numerical | An analytical solution is an exact closed-form expression, while a numerical solution is an approximation obtained by computation. |
| Newton's 2nd Law (Bungee Jumper) | With linear drag the ODE is $`dv/dt = g - (c/m)v`$, and the terminal velocity is $`v_t = gm/c`$. |
| Euler's Method | It is the simplest numerical ODE solver, updating values by $`v(t_{i+1}) = v(t_i) + (dv/dt) \cdot h`$. |
| Step Size Trade-off | A smaller $`h`$ gives better accuracy but more computation, and the error of Euler's method is proportional to $`h`$. |
| Problem Categories | Numerical methods solve roots, linear systems, curve fitting, integration, ODEs, PDEs, and optimization problems. |

---

<br>

## Self-Check Questions

1. **Model components:** In the bungee jumper model $`dv/dt = g - (c/m)v`$, classify $`v`$, $`t`$, $`m`$, $`c`$, and $`g`$ as the dependent variable, independent variable, parameters, or forcing function.

   <details>
   <summary>Answer</summary>

   $`v`$ is the dependent variable, $`t`$ is the independent variable, and $`m`$ and $`c`$ are parameters. The gravitational acceleration $`g`$ is the forcing function, an external influence acting on the system.

   </details>

2. **Terminal velocity:** Find the terminal velocity of the linear drag model for $`m = 68.1`$ kg, $`c = 12.5`$ kg/s, and $`g = 9.81`$ m/s$`^2`$.

   <details>
   <summary>Answer</summary>

   $`v_{\text{terminal}} = gm/c = (9.81 \times 68.1)/12.5 = 668.061/12.5 \approx 53.44`$ m/s. The analytical value of 47.54 m/s at $`t = 12`$ s in Table 2 is still only about 89% of the terminal velocity.

   </details>

3. **One Euler step:** Under the conditions of Section 3.2, start from $`v(6) = 39.90`$ m/s and compute $`v(8)`$ with Euler's method.

   <details>
   <summary>Answer</summary>

   The slope at $`t = 6`$ is $`9.81 - (12.5/68.1)(39.90) = 9.81 - 7.324 = 2.486`$ m/s$`^2`$. Hence $`v(8) = 39.90 + 2.486 \times 2 \approx 44.87`$ m/s, which matches Table 2.

   </details>

4. **Step size and error:** Under the conditions of Section 3.2, how do the error at $`t = 12`$ s and the number of steps change when the step size is reduced from $`h = 0.5`$ s to $`h = 0.25`$ s?

   <details>
   <summary>Answer</summary>

   The error of Euler's method is proportional to $`h`$, so the error drops to about half of 0.60 m/s. The actual computation gives $`v(12) \approx 47.84`$ m/s with an error of about 0.30 m/s. The number of steps doubles from 24 to 48.

   </details>

5. **Problem category:** Finding the drag coefficient $`c`$ for which a 68.1 kg jumper reaches 40 m/s after 10 s of free fall belongs to which category in Section 4? Write the equation to be solved.

   <details>
   <summary>Answer</summary>

   It is a roots-of-equations problem. We must find $`c`$ such that $`f(c) = \frac{gm}{c}\left(1 - e^{-(c/m) \cdot 10}\right) - 40 = 0`$. Since this equation cannot be solved for $`c`$ in closed form, the numerical root-finding methods of Chapters 5-6 are used; with $`g = 9.81`$ m/s$`^2`$ the root is about 14.80 kg/s.

   </details>
