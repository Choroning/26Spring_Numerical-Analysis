# Chapter 13 Lecture — Eigenvalue Methods: Symmetric Matrices

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 13

> **Prerequisites**: [Linear Algebra] Eigenvalues basics (Ch 12).
>
> **Learning Objectives**:
> 1. Apply Jacobi method for symmetric eigenvalue problems
> 2. Implement QR algorithm for eigenvalue computation
> 3. Analyze special properties of symmetric matrices

---

<br>

## Table of Contents

- [1. Eigenvalues and Eigenvectors](#1-eigenvalues-and-eigenvectors)
  - [1.1 Definition](#11-definition)
  - [1.2 Geometric Interpretation](#12-geometric-interpretation)
  - [1.3 Example: Projection Matrix](#13-example-projection-matrix)
  - [1.4 Example: Reflection Matrix](#14-example-reflection-matrix)
- [2. Computing Eigenvalues and Eigenvectors](#2-computing-eigenvalues-and-eigenvectors)
  - [2.1 The Characteristic Equation](#21-the-characteristic-equation)
  - [2.2 The 2x2 Case](#22-the-2x2-case)
  - [2.3 Example: Finding Eigenvalues and Eigenvectors](#23-example-finding-eigenvalues-and-eigenvectors)
- [3. Applications: Differential Equations](#3-applications-differential-equations)
  - [3.1 First-Order ODE — Scalar Case](#31-first-order-ode--scalar-case)
  - [3.2 Converting 2nd-Order to 1st-Order System](#32-converting-2nd-order-to-1st-order-system)
  - [3.3 Matrix Exponential and Stability](#33-matrix-exponential-and-stability)
  - [3.4 Eigenvalue Location and System Behavior](#34-eigenvalue-location-and-system-behavior)
- [4. Example 13.2 — Stability of a First-Order System](#4-example-132--stability-of-a-first-order-system)
- [5. Eigenvalues and Ordinary Differential Equations](#5-eigenvalues-and-ordinary-differential-equations)
  - [5.1 Pure Oscillations — Single 2nd-Order ODE](#51-pure-oscillations--single-2nd-order-ode)
  - [5.2 Eigenvectors of the Oscillation System](#52-eigenvectors-of-the-oscillation-system)
- [6. Coupled 2nd-Order ODE Systems](#6-coupled-2nd-order-ode-systems)
  - [6.1 General Form](#61-general-form)
  - [6.2 Proposed Periodic Solution](#62-proposed-periodic-solution)
  - [6.3 Characteristic Equation for 2nd-Order Systems](#63-characteristic-equation-for-2nd-order-systems)
- [7. Example 13.3 — Behavior of a Coupled System](#7-example-133--behavior-of-a-coupled-system)
  - [7.1 Eigenvalue Computation](#71-eigenvalue-computation)
  - [7.2 Eigenvector Computation](#72-eigenvector-computation)
  - [7.3 Frequencies and General Solution](#73-frequencies-and-general-solution)
  - [7.4 Applying Initial Conditions](#74-applying-initial-conditions)
- [8. Physical Settings: Mass-Spring Systems](#8-physical-settings-mass-spring-systems)
  - [8.1 Two-Mass, Three-Spring System](#81-two-mass-three-spring-system)
  - [8.2 Eigenvalue Formula for the Mass-Spring System](#82-eigenvalue-formula-for-the-mass-spring-system)
- [9. Example 13.4 — Mass-Spring System](#9-example-134--mass-spring-system)
- [10. The Power Method](#10-the-power-method)
  - [10.1 Three-Mass, Four-Spring System (Ex 13.5)](#101-three-mass-four-spring-system-ex-135)
  - [10.2 Algorithm](#102-algorithm)
  - [10.3 Iteration Walkthrough](#103-iteration-walkthrough)
  - [10.4 Convergence and Finding the Smallest Eigenvalue](#104-convergence-and-finding-the-smallest-eigenvalue)
- [11. Python Implementation](#11-python-implementation)
- [Summary](#summary)

---

<br>

## 1. Eigenvalues and Eigenvectors

### 1.1 Definition

We work on **square matrices**. Given an $n \times n$ matrix $[A]$, consider the equation:

$$[A]\{x\} = \lambda \{x\}$$

where $\{x\}$ is an $n \times 1$ vector and $\lambda \in \mathbb{R}$ or $\lambda \in \mathbb{C}$.

For any value of $\lambda$, a **trivial solution** always exists: $\{x\} = \{0\}$.

If there is a **non-trivial solution** $\{x\} \neq \{0\}$ such that $[A]\{x\} = \lambda \{x\}$, then:

- $\{x\}$ is called an **eigenvector** of $[A]$
- $\lambda$ is called the **eigenvalue** associated with $\{x\}$

> **[Linear Algebra]** The word "eigen" comes from German, meaning "own" or "characteristic." An eigenvector is a special direction in which the matrix acts simply as a scalar multiplier. An $n \times n$ matrix has at most $n$ eigenvalues (counting multiplicities), which are the roots of the characteristic polynomial.

### 1.2 Geometric Interpretation

The eigenvector $\{x\}$ is **parallel** to $[A]\{x\}$. That is, multiplying $\{x\}$ by the matrix $[A]$ does not change its direction — it only scales it by the factor $\lambda$.

$$\lambda \{x\} \quad \longrightarrow \quad \text{(same direction as } [A]\{x\}\text{)}$$

### 1.3 Example: Projection Matrix

Suppose $P$ is the projection matrix onto a plane.

**Case 1:** For any vector $\{x\}$ that lies **in the plane**:

$$[P]\{x\} = \{x\} = 1 \cdot \{x\}$$

Here $\lambda = 1$ is the eigenvalue, and $\{x\}$ (any vector in the plane) is the eigenvector.

**Case 2:** If a vector $\{x\}$ is **perpendicular** to the plane:

$$[P]\{x\} = \{0\} = 0 \cdot \{x\}$$

Here $\lambda = 0$ is the eigenvalue, and $\{x\}$ (normal to the plane) is the eigenvector.

### 1.4 Example: Reflection Matrix

Consider the reflection matrix:

$$A = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}$$

This is a reflection map. It has two eigenvectors: $\begin{Bmatrix} 1 \\ 1 \end{Bmatrix}$ and $\begin{Bmatrix} 1 \\ -1 \end{Bmatrix}$.

The two eigenvectors **span the whole space**, meaning any vector in $\mathbb{R}^2$ can be written as a linear combination of these eigenvectors.

> **[Linear Algebra]** When the eigenvectors of an $n \times n$ matrix span the entire $\mathbb{R}^n$, the matrix is said to be **diagonalizable**. Symmetric matrices are always diagonalizable and have real eigenvalues — a key result known as the Spectral Theorem.

---

<br>

## 2. Computing Eigenvalues and Eigenvectors

### 2.1 The Characteristic Equation

Starting from the eigenvalue equation:

$$[A]\{x\} = \lambda \{x\}$$

Rearranging:

$$[A]\{x\} - \lambda [I]\{x\} = \{0\}$$

$$[A - \lambda I]\{x\} = \{0\}$$

For a **non-trivial solution** $\{x\} \neq \{0\}$ to exist, the matrix $[A - \lambda I]$ must be **singular** (non-invertible). This requires:

$$\det(A - \lambda I) = |A - \lambda I| = 0$$

This is called the **characteristic equation** (or characteristic polynomial). For an $n \times n$ matrix, this yields a polynomial of degree $n$ in $\lambda$.

### 2.2 The 2x2 Case

For a $2 \times 2$ matrix:

$$\begin{vmatrix} a_{11} - \lambda & a_{12} \\ a_{21} & a_{22} - \lambda \end{vmatrix} = 0$$

Expanding:

$$\lambda^2 - (a_{11} + a_{22})\lambda + (a_{11}a_{22} - a_{12}a_{21}) = 0$$

Using the quadratic formula:

$$\lambda = \frac{(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}$$

> **[Linear Algebra]** The quantity $a_{11} + a_{22}$ is the **trace** of $A$, and $a_{11}a_{22} - a_{12}a_{21}$ is the **determinant** of $A$. For any $2 \times 2$ matrix, the sum of eigenvalues equals the trace, and the product of eigenvalues equals the determinant.

### 2.3 Example: Finding Eigenvalues and Eigenvectors

Consider:

$$A = \begin{bmatrix} 3 & 1 \\ 1 & 3 \end{bmatrix}$$

**Step 1: Characteristic equation**

$$|A - \lambda I| = \begin{vmatrix} 3 - \lambda & 1 \\ 1 & 3 - \lambda \end{vmatrix} = (3 - \lambda)^2 - 1 = \lambda^2 - 6\lambda + 8 = (\lambda - 4)(\lambda - 2) = 0$$

Therefore: $\lambda_1 = 4$ and $\lambda_2 = 2$.

**Step 2: Eigenvector for $\lambda_1 = 4$**

$$\begin{bmatrix} 3 - 4 & 1 \\ 1 & 3 - 4 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}$$

$$-x_1 + x_2 = 0 \quad \text{and} \quad x_1 - x_2 = 0 \quad \Rightarrow \quad \text{same equation}$$

Let $x_1 = 1$, then $x_2 = 1$:

$$\{x\}_1 = \begin{Bmatrix} 1 \\ 1 \end{Bmatrix}$$

**Step 3: Eigenvector for $\lambda_2 = 2$**

$$\begin{bmatrix} 3 - 2 & 1 \\ 1 & 3 - 2 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}$$

$$x_1 + x_2 = 0$$

Let $x_1 = 1$, then $x_2 = -1$:

$$\{x\}_2 = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}$$

---

<br>

## 3. Applications: Differential Equations

### 3.1 First-Order ODE — Scalar Case

Consider the scalar ODE:

$$\frac{dx}{dt} = ax, \quad a \in \mathbb{R}$$

The solution is:

$$x(t) = e^{at} C, \quad C \text{ is a constant}$$

The behavior depends on $a$:

| Value of $a$ | Solution behavior |
|---|---|
| $a > 0$ (real, positive) | Exponential growth (unstable) |
| $a < 0$ (real, negative) | Exponential decay (stable) |
| $a = 0$ | Constant: $x(t) = C$ |
| $a = i\beta$ (imaginary) | Oscillation: $x(t) = e^{it}C = (\cos t + i\sin t)C$ |

### 3.2 Converting 2nd-Order to 1st-Order System

A second-order ODE with constant coefficients:

$$\frac{d^2 y}{dt^2} + a\frac{dy}{dt} + by = 0$$

can be converted to a system of **first-order** ODEs by the substitution:

$$x_1 = y, \quad x_2 = \frac{dx_1}{dt}$$

This gives:

$$\frac{dx_1}{dt} = x_2$$

$$\frac{dx_2}{dt} = -bx_1 - ax_2$$

In matrix form:

$$\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 0 & 1 \\ -b & -a \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}$$

We have converted one 2nd-order ODE into **two first-order** ODEs.

### 3.3 Matrix Exponential and Stability

For the vector ODE:

$$\frac{d\mathbf{x}}{dt} = A\mathbf{x}, \quad \mathbf{x} \in \mathbb{R}^n, \quad A \in \mathbb{R}^{n \times n}$$

The solution is:

$$\mathbf{x}(t) = e^{At}\mathbf{C}, \quad \mathbf{C} \in \mathbb{R}^n$$

Here $e^{At}$ is the **matrix exponential**. If $A$ can be diagonalized as $A = U D U^T$ (where $D$ is the diagonal matrix of eigenvalues and $U$ contains the eigenvectors), then:

$$e^{At} = U \, e^{Dt} \, U^T$$

This means:

- **Negative real eigenvalues** produce stable (decaying) solutions
- **Positive real eigenvalues** produce unstable (growing) solutions
- **Imaginary eigenvalues** produce oscillatory solutions

### 3.4 Eigenvalue Location and System Behavior

The complex plane provides an intuitive map for stability:

| Region | Eigenvalue type | System behavior |
|---|---|---|
| Left half-plane ($\text{Re}(\lambda) < 0$) | Negative real part | **Stable** — solution decays |
| Right half-plane ($\text{Re}(\lambda) > 0$) | Positive real part | **Unstable** — solution grows |
| Imaginary axis ($\text{Re}(\lambda) = 0$) | Pure imaginary | **Oscillation** — bounded periodic motion |
| Far left | Large negative real | Faster convergence |
| Far right | Large positive real | Faster divergence |

---

<br>

## 4. Example 13.2 — Stability of a First-Order System

Consider the system:

$$\frac{dx_1}{dt} = -3x_1 + x_2$$

$$\frac{dx_2}{dt} = x_1 - 3x_2$$

In matrix form:

$$\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -3 & 1 \\ 1 & -3 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}$$

**Characteristic equation:**

$$\begin{vmatrix} -3 - \lambda & 1 \\ 1 & -3 - \lambda \end{vmatrix} = (3 + \lambda)^2 - 1 = \lambda^2 + 6\lambda + 8 = (\lambda + 4)(\lambda + 2) = 0$$

$$\therefore \lambda = -2 \text{ or } -4$$

**Both eigenvalues are negative**, so the solutions will be **stable** and die out to zero **without oscillation**.

---

<br>

## 5. Eigenvalues and Ordinary Differential Equations

### 5.1 Pure Oscillations — Single 2nd-Order ODE

Consider the undamped oscillator:

$$\ddot{y} = -ay \quad \Longleftrightarrow \quad \frac{d^2 y}{dt^2} = -ay$$

Cast into first-order ODEs with $x_1 = y$, $x_2 = \frac{dx_1}{dt}$:

$$\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 0 & 1 \\ -a & 0 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}$$

**Finding eigenvalues of $A$:**

$$|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ -a & -\lambda \end{vmatrix} = \lambda^2 + a = 0$$

$$\therefore \lambda = \pm \sqrt{-a} = \pm \sqrt{a}\, i$$

The eigenvalues are **pure imaginary** (lie on the imaginary axis), so the solution will be **oscillatory** with frequency $\sqrt{a}/2\pi$.

Recall that:

$$e^{i\sqrt{a}\,t} = \cos(\sqrt{a}\,t) + i\sin(\sqrt{a}\,t)$$

$$\text{Period} = \frac{2\pi}{\sqrt{a}}, \qquad \text{Frequency} = \frac{1}{\text{Period}} = \frac{\sqrt{a}}{2\pi}$$

> **[Linear Algebra]** Eigenvalues are directly related to the **frequency of oscillations**. This is why eigenvalue problems appear throughout physics and engineering whenever we analyze vibrating systems.

### 5.2 Eigenvectors of the Oscillation System

Solve $[A - \lambda I]\{x\} = \{0\}$:

$$\begin{bmatrix} -\lambda & 1 \\ -a & -\lambda \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}$$

**Case i)** $\lambda = \sqrt{a}\, i$:

$$-\sqrt{a}\,i \cdot x_1 + x_2 = 0$$

Take $x_1 = 1$, then $x_2 = \sqrt{a}\,i$:

$$\{x\} = \begin{Bmatrix} 1 \\ \sqrt{a}\,i \end{Bmatrix}$$

**Case ii)** $\lambda = -\sqrt{a}\, i$:

$$\sqrt{a}\,i \cdot x_1 + x_2 = 0$$

Take $x_1 = 1$, then $x_2 = -\sqrt{a}\,i$:

$$\{x\} = \begin{Bmatrix} 1 \\ -\sqrt{a}\,i \end{Bmatrix}$$

---

<br>

## 6. Coupled 2nd-Order ODE Systems

### 6.1 General Form

A system of two coupled second-order ODEs:

$$\frac{d^2 y_1}{dt^2} = -a_{11}y_1 + a_{12}y_2$$

$$\frac{d^2 y_2}{dt^2} = a_{21}y_1 - a_{22}y_2$$

In matrix form:

$$\frac{d^2}{dt^2}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix} \begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} \quad \cdots \quad (*)$$

### 6.2 Proposed Periodic Solution

Based on the periodic solution form, we propose:

$$\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t}$$

where $\omega$ is the frequency. Substituting into $(*)$:

$$-\omega^2 \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t} = \begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t}$$

Cancel $e^{i\omega t}$:

$$\underbrace{-\omega^2}_{\lambda} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} = \underbrace{\begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix}}_{A} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix}$$

**This is an eigenvalue problem!** The eigenvalue is $\lambda = -\omega^2$, so $\omega^2 = -\lambda$.

### 6.3 Characteristic Equation for 2nd-Order Systems

From $|A - \lambda I| = 0$:

$$\begin{vmatrix} -a_{11} - \lambda & a_{12} \\ a_{21} & -a_{22} - \lambda \end{vmatrix} = 0$$

$$(a_{11} + \lambda)(a_{22} + \lambda) - a_{12}a_{21} = 0$$

$$\lambda^2 + (a_{11} + a_{22})\lambda + a_{11}a_{22} - a_{12}a_{21} = 0$$

$$\therefore \lambda = \frac{-(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}$$

Equivalently, in terms of $\omega^2 = -\lambda$:

$$|A + \omega^2 I| = \begin{vmatrix} -a_{11} + \omega^2 & a_{12} \\ a_{21} & -a_{22} + \omega^2 \end{vmatrix} = 0$$

$$(\omega^2 - a_{11})(\omega^2 - a_{22}) - a_{12}a_{21} = 0$$

$$\omega^4 - (a_{11} + a_{22})\omega^2 + a_{11}a_{22} - a_{12}a_{21} = 0$$

$$\therefore \omega^2 = \frac{(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}$$

---

<br>

## 7. Example 13.3 — Behavior of a Coupled System

### 7.1 Eigenvalue Computation

Given:

$$\frac{d^2 y_1}{dt^2} = -5y_1 + 2y_2, \qquad \frac{d^2 y_2}{dt^2} = 2y_1 - 2y_2$$

Matrix form:

$$\frac{d^2}{dt^2}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \underbrace{\begin{bmatrix} -5 & 2 \\ 2 & -2 \end{bmatrix}}_{A}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix}$$

Characteristic equation:

$$\begin{vmatrix} -5 - \lambda & 2 \\ 2 & -2 - \lambda \end{vmatrix} = (5 + \lambda)(2 + \lambda) - 4 = \lambda^2 + 7\lambda + 6 = (\lambda + 6)(\lambda + 1) = 0$$

$$\therefore \lambda_1 = -1, \quad \lambda_2 = -6$$

### 7.2 Eigenvector Computation

**For $\lambda_1 = -1$:**

$$\begin{bmatrix} -5 + 1 & 2 \\ 2 & -2 + 1 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix} \quad \Rightarrow \quad -4x_1 + 2x_2 = 0 \quad \Rightarrow \quad 2x_1 - x_2 = 0$$

Take $x_1 = 1$, $x_2 = 2$:

$$\{x\}_1 = \begin{Bmatrix} 1 \\ 2 \end{Bmatrix}$$

**For $\lambda_2 = -6$:**

$$\begin{bmatrix} -5 + 6 & 2 \\ 2 & -2 + 6 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix} \quad \Rightarrow \quad x_1 + 2x_2 = 0$$

Take $x_2 = -1$, $x_1 = 2$:

$$\{x\}_2 = \begin{Bmatrix} 2 \\ -1 \end{Bmatrix}$$

### 7.3 Frequencies and General Solution

From the eigenvalues, $\omega^2 = -\lambda$:

$$\omega_1^2 = 1 \quad \Rightarrow \quad \omega_1 = \pm 1$$

$$\omega_2^2 = 6 \quad \Rightarrow \quad \omega_2 = \pm \sqrt{6}$$

The solution behaves as:

$$\{x\}_1 e^{\pm it} = \{x\}_1 (\cos t \pm i\sin t)$$

$$\{x\}_2 e^{\pm i\sqrt{6}t} = \{x\}_2 (\cos\sqrt{6}t \pm i\sin\sqrt{6}t)$$

General solution:

$$\mathbf{y} = C_1 \{x\}_1 e^{+it} + C_2 \{x\}_2 e^{+i\sqrt{6}t} + C_3 \{x\}_1 e^{-it} + C_4 \{x\}_2 e^{-i\sqrt{6}t}$$

Its time derivative:

$$\dot{\mathbf{y}} = iC_1 \{x\}_1 e^{+it} + i\sqrt{6}C_2 \{x\}_2 e^{+i\sqrt{6}t} - iC_3 \{x\}_1 e^{-it} - i\sqrt{6}C_4 \{x\}_2 e^{-i\sqrt{6}t}$$

### 7.4 Applying Initial Conditions

Given initial conditions:

$$\mathbf{y}(t=0) = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}, \quad \dot{\mathbf{y}}(t=0) = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}$$

From the zero-velocity condition, we require $C_1 = C_3$ and $C_2 = C_4$. Let $\tilde{C}_1 = C_1 = C_3$ and $\tilde{C}_2 = C_2 = C_4$.

Substituting $t = 0$ into the position equations:

From $y_1$: $1 = 2\tilde{C}_1 + 4\tilde{C}_2$

From $y_2$: $-1 = 4\tilde{C}_1 - 2\tilde{C}_2$

Solving: multiply the second equation by 2 and add:

$$-2 = 8\tilde{C}_1 - 4\tilde{C}_2$$

$$1 = 2\tilde{C}_1 + 4\tilde{C}_2$$

Adding: $-1 = 10\tilde{C}_1 \Rightarrow \tilde{C}_1 = -\frac{1}{10}$

From $1 = 2(-\frac{1}{10}) + 4\tilde{C}_2 \Rightarrow \tilde{C}_2 = \frac{1}{4}(1 + \frac{1}{5}) = \frac{3}{10}$

**Final solution:**

$$\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = -\frac{1}{5}\begin{Bmatrix} 1 \\ 2 \end{Bmatrix}\cos(t) + \frac{3}{5}\begin{Bmatrix} 2 \\ -1 \end{Bmatrix}\cos(\sqrt{6}\,t)$$

---

<br>

## 8. Physical Settings: Mass-Spring Systems

### 8.1 Two-Mass, Three-Spring System

Consider a two-mass, three-spring system with frictionless rollers, vibrating between two walls.

Two masses $m_1$ and $m_2$ are connected by springs of stiffness $k$. Newton's second law for each mass:

$$m_1 \frac{d^2 x_1}{dt^2} = -kx_1 + k(x_2 - x_1)$$

$$m_2 \frac{d^2 x_2}{dt^2} = -k(x_2 - x_1) - kx_2$$

In matrix form:

$$\frac{d^2}{dt^2}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -2k/m_1 & k/m_1 \\ k/m_2 & -2k/m_2 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}$$

### 8.2 Eigenvalue Formula for the Mass-Spring System

From $|A - \lambda I| = 0$:

$$\begin{vmatrix} -2k/m_1 - \lambda & k/m_1 \\ k/m_2 & -2k/m_2 - \lambda \end{vmatrix} = 0$$

$$(\lambda + 2k/m_1)(\lambda + 2k/m_2) - k^2/(m_1 m_2) = 0$$

$$\lambda^2 + \lambda(2k/m_1 + 2k/m_2) + 4k^2/(m_1 m_2) - k^2/(m_1 m_2) = 0$$

Solving via the quadratic formula:

$$-\omega^2 = \lambda = -\frac{k}{m_1 m_2}(m_1 + m_2) \pm \frac{k}{m_1 m_2}\sqrt{m_2^2 + m_1^2 - m_1 m_2}$$

---

<br>

## 9. Example 13.4 — Mass-Spring System

Given: $m_1 = m_2 = 40$ kg, $k = 200$ N/m.

**Eigenvalue computation:**

$$\lambda = -\frac{k}{m_1 m_2}(m_1 + m_2) \pm \frac{k}{m_1 m_2}\sqrt{m_2^2 + m_1^2 - m_1 m_2}$$

$$= -\frac{200}{1600}(80) \pm \frac{1}{8}\sqrt{1600}$$

$$= -10 \pm 5$$

$$\therefore \lambda_1 = -5, \quad \lambda_2 = -15$$

**Frequencies:**

$$\omega^2 = -\lambda \quad \Rightarrow \quad \omega_1 = \sqrt{5}, \quad f_1 = \frac{\sqrt{5}}{2\pi}$$

$$\omega_2 = \sqrt{15}, \quad f_2 = \frac{\sqrt{15}}{2\pi}$$

(Recall: $\omega = 2\pi f$)

**Eigenvectors:**

For $\lambda_1 = -5$: With $A = \begin{bmatrix} -10 & 5 \\ 5 & -10 \end{bmatrix}$:

$$\begin{bmatrix} -10 - (-5) & 5 \\ 5 & -10 - (-5) \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -5 & 5 \\ 5 & -5 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}$$

$$-x_1 + x_2 = 0 \quad \Rightarrow \quad \{x\}_1 = \begin{Bmatrix} 1 \\ 1 \end{Bmatrix}$$

Both masses move **in the same direction** (in-phase mode).

For $\lambda_2 = -15$:

$$\begin{bmatrix} -10 + 15 & 5 \\ 5 & -10 + 15 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 5 & 5 \\ 5 & 5 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}$$

$$x_1 + x_2 = 0 \quad \Rightarrow \quad \{x\}_2 = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}$$

The masses move in **opposite directions** (out-of-phase mode).

**General solution:**

$$\mathbf{y} = C_1 \{x\}_1 \cos\sqrt{5}\,t + C_2 \{x\}_1 \sin\sqrt{5}\,t + C_3 \{x\}_2 \cos\sqrt{15}\,t + C_4 \{x\}_2 \sin\sqrt{15}\,t$$

The coefficients $C_1, C_2, C_3, C_4$ are determined by the initial conditions.

Expanding for each component:

$$y_1 = C_1 \cdot 1 \cdot \cos\sqrt{5}\,t + C_2 \cdot 1 \cdot \sin\sqrt{5}\,t + C_3 \cdot 1 \cdot \cos\sqrt{15}\,t + C_4 \cdot 1 \cdot \sin\sqrt{15}\,t$$

$$y_2 = C_1 \cdot 1 \cdot \cos\sqrt{5}\,t + C_2 \cdot 1 \cdot \sin\sqrt{5}\,t + C_3 \cdot (-1) \cdot \cos\sqrt{15}\,t + C_4 \cdot (-1) \cdot \sin\sqrt{15}\,t$$

---

<br>

## 10. The Power Method

The power method **iteratively** determines the **largest** (dominant) eigenvalue and its corresponding eigenvector.

### 10.1 Three-Mass, Four-Spring System (Ex 13.5)

Consider a three-mass, four-spring system between two walls. With $m = m_1 = m_2 = m_3 = 1$ kg and $k = 20$ N/m:

Newton's second law for each mass:

$$m_1 \frac{d^2 x_1}{dt^2} = -kx_1 + k(x_2 - x_1)$$

$$m_2 \frac{d^2 x_2}{dt^2} = -k(x_2 - x_1) + k(x_3 - x_2)$$

$$m_3 \frac{d^2 x_3}{dt^2} = -k(x_3 - x_2) - kx_3$$

Assuming a periodic solution $\{x_i\} = \{X_i\}e^{i\omega t}$ and substituting:

$$\begin{bmatrix} 40 - \omega^2 & -20 & 0 \\ -20 & 40 - \omega^2 & -20 \\ 0 & -20 & 40 - \omega^2 \end{bmatrix} \begin{Bmatrix} X_1 \\ X_2 \\ X_3 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \\ 0 \end{Bmatrix}$$

This is equivalent to the eigenvalue problem $A\{X\} = \lambda\{X\}$ with:

$$A = \begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix}, \quad \lambda = \omega^2$$

### 10.2 Algorithm

1. **Start** with an initial guess vector $\{x\}^{(0)}$
2. **Multiply**: $\{y\}^{(k)} = [A]\{x\}^{(k-1)}$
3. **Normalize**: extract the largest element (in absolute value) as the eigenvalue estimate $\lambda^{(k)}$, and divide $\{y\}^{(k)}$ by this element to get the new normalized vector $\{x\}^{(k)}$
4. **Repeat** until the eigenvalue estimate converges

The normalization factor at each step is our **estimate of the dominant eigenvalue**.

### 10.3 Iteration Walkthrough

**Iteration 1:** Take $X_1 = X_2 = X_3 = 1$ as initial guess:

$$\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ 1 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 20 \\ 0 \\ 20 \end{Bmatrix} = 20 \begin{Bmatrix} 1 \\ 0 \\ 1 \end{Bmatrix}$$

Eigenvalue estimate: $\lambda \approx 20$.

**Iteration 2:** Use $\{1, 0, 1\}^T$:

$$\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ 0 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 40 \\ -40 \\ 40 \end{Bmatrix} = 40 \begin{Bmatrix} 1 \\ -1 \\ 1 \end{Bmatrix}$$

Eigenvalue estimate: $\lambda \approx 40$. Approximate error: $\varepsilon_a = \left|\frac{40 - 20}{40}\right| \times 100\% = 50\%$.

**Iteration 3:** Use $\{1, -1, 1\}^T$:

$$\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ -1 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 60 \\ -80 \\ 60 \end{Bmatrix} = -80 \begin{Bmatrix} -3/4 \\ 1 \\ -3/4 \end{Bmatrix}$$

Eigenvalue estimate: $\lambda \approx -80$. Error: $\varepsilon_a = \left|\frac{-80 - 40}{-80}\right| \times 100\% = 150\%$.

**Iteration 4:** Use $\{-3/4, 1, -3/4\}^T$:

$$A\begin{Bmatrix} -3/4 \\ 1 \\ -3/4 \end{Bmatrix} = \begin{Bmatrix} -50 \\ 70 \\ -50 \end{Bmatrix} = 70\begin{Bmatrix} -5/7 \\ 1 \\ -5/7 \end{Bmatrix}$$

Eigenvalue estimate: $\lambda \approx 70$. Error: $\varepsilon_a \approx 21.4\%$.

**Iteration 5:** Use $\{-5/7, 1, -5/7\}^T$:

$$A\begin{Bmatrix} -5/7 \\ 1 \\ -5/7 \end{Bmatrix} = \begin{Bmatrix} -48.52 \\ 68.52 \\ -48.52 \end{Bmatrix} = 68.52\begin{Bmatrix} -0.71 \\ 1 \\ -0.71 \end{Bmatrix}$$

Eigenvalue estimate: $\lambda \approx 68.52$. Error: $\varepsilon_a = \left|\frac{68.52 - 70}{68.52}\right| \approx 2.08\%$.

### 10.4 Convergence and Finding the Smallest Eigenvalue

The eigenvalue is converging. After sufficient iterations:

$$\lambda = 68.28427..., \quad \{x\} = \begin{Bmatrix} -0.70711 \\ 1 \\ -0.70711 \end{Bmatrix}$$

> **[Linear Algebra]** The power method converges to the eigenvalue with the **largest absolute value**. The rate of convergence depends on the ratio $|\lambda_2/\lambda_1|$, where $\lambda_1$ is the dominant eigenvalue and $\lambda_2$ is the second-largest. The closer this ratio is to 1, the slower the convergence.

**For finding the smallest eigenvalue**, one can apply the power method to the **inverse** of $A$:

$$[A]^{-1}\{x\} = \frac{1}{\lambda}\{x\}$$

Since the largest eigenvalue of $A^{-1}$ is $1/\lambda_{\min}$, the power method applied to $A^{-1}$ converges to $1/\lambda_{\min}$. This is called the **Inverse Power Method**.

---

<br>

## 11. Python Implementation

### Eigenvalue computation using NumPy

```python
import numpy as np

# Example: 2x2 matrix from Section 2.3
A = np.array([[3, 1],
              [1, 3]])

eigenvalues, eigenvectors = np.linalg.eig(A)
print("Eigenvalues:", eigenvalues)       # [4. 2.]
print("Eigenvectors:\n", eigenvectors)   # columns are eigenvectors
```

### Power Method implementation

```python
import numpy as np

def power_method(A, x0, tol=1e-6, max_iter=100):
    """
    Power method to find the dominant eigenvalue and eigenvector.

    Parameters
    ----------
    A : ndarray
        Square matrix (n x n).
    x0 : ndarray
        Initial guess vector (n x 1).
    tol : float
        Convergence tolerance.
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    lam : float
        Dominant eigenvalue estimate.
    x : ndarray
        Corresponding eigenvector (normalized).
    """
    x = x0.copy().astype(float)
    lam_old = 0.0

    for k in range(max_iter):
        y = A @ x
        # Find the element with the largest absolute value
        idx = np.argmax(np.abs(y))
        lam = y[idx]
        x = y / lam

        # Check convergence
        if abs(lam - lam_old) / abs(lam) < tol:
            print(f"Converged in {k + 1} iterations")
            return lam, x
        lam_old = lam

    print("Did not converge within max_iter iterations")
    return lam, x


# Example 13.5: Three-mass four-spring system
A = np.array([[ 40, -20,   0],
              [-20,  40, -20],
              [  0, -20,  40]])

x0 = np.array([1, 1, 1])
lam, x = power_method(A, x0)
print(f"Dominant eigenvalue: {lam:.5f}")
print(f"Eigenvector: {x}")
```

### Mass-Spring System (Ex 13.4)

```python
import numpy as np

# Parameters
m1, m2 = 40.0, 40.0  # kg
k = 200.0             # N/m

# Coefficient matrix
A = np.array([[-2*k/m1,    k/m1],
              [   k/m2, -2*k/m2]])

eigenvalues, eigenvectors = np.linalg.eig(A)
print("Eigenvalues (lambda):", eigenvalues)  # [-5. -15.]

# Natural frequencies
omega = np.sqrt(-eigenvalues)
freq = omega / (2 * np.pi)
print("Angular frequencies (omega):", omega)
print("Frequencies (Hz):", freq)
print("Eigenvectors:\n", eigenvectors)
```

### Inverse Power Method for smallest eigenvalue

```python
import numpy as np

def inverse_power_method(A, x0, tol=1e-6, max_iter=100):
    """
    Inverse power method to find the smallest eigenvalue.
    """
    x = x0.copy().astype(float)
    lam_old = 0.0

    for k in range(max_iter):
        # Solve A @ y = x instead of computing A_inv @ x
        y = np.linalg.solve(A, x)
        idx = np.argmax(np.abs(y))
        lam_inv = y[idx]
        x = y / lam_inv

        lam = 1.0 / lam_inv
        if abs(lam - lam_old) / abs(lam) < tol:
            print(f"Converged in {k + 1} iterations")
            return lam, x
        lam_old = lam

    print("Did not converge within max_iter iterations")
    return lam, x


# Example: find smallest eigenvalue of the 3x3 spring system
A = np.array([[ 40, -20,   0],
              [-20,  40, -20],
              [  0, -20,  40]])

x0 = np.array([1, 1, 1])
lam_min, x_min = inverse_power_method(A, x0)
print(f"Smallest eigenvalue: {lam_min:.5f}")
print(f"Eigenvector: {x_min}")
```

---

<br>

## Summary

| Topic | Key Formula / Concept |
|---|---|
| **Eigenvalue equation** | $[A]\{x\} = \lambda\{x\}$ |
| **Characteristic equation** | $\det(A - \lambda I) = 0$ |
| **2x2 eigenvalues** | $\lambda = \frac{(a_{11}+a_{22}) \pm \sqrt{(a_{11}+a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}$ |
| **Non-trivial solution condition** | $[A - \lambda I]$ must be singular |
| **Eigenvector property** | $\{x\}$ is parallel to $[A]\{x\}$ |
| **Scalar ODE** $dx/dt = ax$ | Solution: $x(t) = Ce^{at}$ |
| **Matrix ODE** $d\mathbf{x}/dt = A\mathbf{x}$ | Solution: $\mathbf{x}(t) = e^{At}\mathbf{C}$; $e^{At} = Ue^{Dt}U^T$ |
| **Stability** | $\text{Re}(\lambda) < 0$: stable; $\text{Re}(\lambda) > 0$: unstable; $\text{Re}(\lambda) = 0$: oscillatory |
| **2nd-order ODE to eigenvalue** | $\ddot{y} = Ay \Rightarrow \lambda = -\omega^2$, so $\omega = \sqrt{-\lambda}$ |
| **Frequency relation** | $\text{Period} = 2\pi/\omega$; $f = \omega/(2\pi)$ |
| **Mass-spring system** | $m\ddot{x} = -kx + k(x_{\text{neighbor}} - x) \Rightarrow$ eigenvalue problem |
| **Power Method** | Iteratively computes the **dominant** (largest) eigenvalue |
| **Inverse Power Method** | Apply power method to $A^{-1}$ to find the **smallest** eigenvalue |
| **Convergence rate** | Depends on $\|\lambda_2/\lambda_1\|$ — closer to 1 means slower convergence |
