# Chapter 13 — Eigenvalue Methods: Symmetric Matrices

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 13

> **Learning Objectives:**
> 1. Compute eigenvalues and eigenvectors from the characteristic equation
> 2. Determine the stability and oscillation of a linear ODE system from the location of its eigenvalues
> 3. Formulate coupled oscillators and mass-spring systems as eigenvalue problems and find their natural frequencies and modes
> 4. Find the largest and smallest eigenvalues with the power method and the inverse power method

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
  - [3.1 First-Order ODE: Scalar Case](#31-first-order-ode-scalar-case)
  - [3.2 Converting 2nd-Order to 1st-Order System](#32-converting-2nd-order-to-1st-order-system)
  - [3.3 Matrix Exponential and Stability](#33-matrix-exponential-and-stability)
  - [3.4 Eigenvalue Location and System Behavior](#34-eigenvalue-location-and-system-behavior)
- [4. Example 13.2: Stability of a First-Order System](#4-example-132-stability-of-a-first-order-system)
- [5. Eigenvalues and Ordinary Differential Equations](#5-eigenvalues-and-ordinary-differential-equations)
  - [5.1 Pure Oscillations: Single 2nd-Order ODE](#51-pure-oscillations-single-2nd-order-ode)
  - [5.2 Eigenvectors of the Oscillation System](#52-eigenvectors-of-the-oscillation-system)
- [6. Coupled 2nd-Order ODE Systems](#6-coupled-2nd-order-ode-systems)
  - [6.1 General Form](#61-general-form)
  - [6.2 Proposed Periodic Solution](#62-proposed-periodic-solution)
  - [6.3 Characteristic Equation for 2nd-Order Systems](#63-characteristic-equation-for-2nd-order-systems)
- [7. Example 13.3: Behavior of a Coupled System](#7-example-133-behavior-of-a-coupled-system)
  - [7.1 Eigenvalue Computation](#71-eigenvalue-computation)
  - [7.2 Eigenvector Computation](#72-eigenvector-computation)
  - [7.3 Frequencies and General Solution](#73-frequencies-and-general-solution)
  - [7.4 Applying Initial Conditions](#74-applying-initial-conditions)
- [8. Physical Settings: Mass-Spring Systems](#8-physical-settings-mass-spring-systems)
  - [8.1 Two-Mass, Three-Spring System](#81-two-mass-three-spring-system)
  - [8.2 Eigenvalue Formula for the Mass-Spring System](#82-eigenvalue-formula-for-the-mass-spring-system)
- [9. Example 13.4: Mass-Spring System](#9-example-134-mass-spring-system)
- [10. The Power Method](#10-the-power-method)
  - [10.1 Three-Mass, Four-Spring System (Ex 13.5)](#101-three-mass-four-spring-system-ex-135)
  - [10.2 Algorithm](#102-algorithm)
  - [10.3 Iteration Walkthrough](#103-iteration-walkthrough)
  - [10.4 Convergence and Finding the Smallest Eigenvalue](#104-convergence-and-finding-the-smallest-eigenvalue)
- [11. Python Implementation](#11-python-implementation)
  - [11.1 Eigenvalue computation using NumPy](#111-eigenvalue-computation-using-numpy)
  - [11.2 Power Method implementation](#112-power-method-implementation)
  - [11.3 Mass-Spring System (Ex 13.4)](#113-mass-spring-system-ex-134)
  - [11.4 Inverse Power Method for smallest eigenvalue](#114-inverse-power-method-for-smallest-eigenvalue)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Eigenvalues and Eigenvectors

### 1.1 Definition

We work on **square matrices**. Given an $`n \times n`$ matrix $`[A]`$, consider the equation:

```math
[A]\{x\} = \lambda \{x\}
```

where $`\{x\}`$ is an $`n \times 1`$ vector and $`\lambda \in \mathbb{R}`$ or $`\lambda \in \mathbb{C}`$.

For any value of $`\lambda`$, a **trivial solution** $`\{x\} = \{0\}`$ always exists.

If there is a **non-trivial solution** $`\{x\} \neq \{0\}`$ such that $`[A]\{x\} = \lambda \{x\}`$, then:

- $`\{x\}`$ is called an **eigenvector** of $`[A]`$.
- $`\lambda`$ is called the **eigenvalue** associated with $`\{x\}`$.

> **[Linear Algebra]** The word "eigen" comes from German, meaning "own" or "characteristic." An eigenvector is a special direction in which the matrix acts simply as a scalar multiplier. An $`n \times n`$ matrix has exactly $`n`$ (complex) eigenvalues counting multiplicities, which are the roots of the characteristic polynomial, so it has at most $`n`$ distinct eigenvalues.

### 1.2 Geometric Interpretation

The eigenvector $`\{x\}`$ is **parallel** to $`[A]\{x\}`$. That is, multiplying $`\{x\}`$ by the matrix $`[A]`$ keeps it on the same line through the origin and only scales it by the factor $`\lambda`$: it points in the same direction if $`\lambda > 0`$, in the opposite direction if $`\lambda < 0`$, and becomes the zero vector if $`\lambda = 0`$.

```math
\lambda \{x\} \quad \longrightarrow \quad \text{(same direction as } [A]\{x\}\text{)}
```

### 1.3 Example: Projection Matrix

Suppose $`P`$ is the projection matrix onto a plane.

**Case 1:** For any vector $`\{x\}`$ that lies **in the plane**:

```math
[P]\{x\} = \{x\} = 1 \cdot \{x\}
```

Here $`\lambda = 1`$ is the eigenvalue, and $`\{x\}`$ (any nonzero vector in the plane) is the eigenvector.

**Case 2:** If a vector $`\{x\}`$ is **perpendicular** to the plane:

```math
[P]\{x\} = \{0\} = 0 \cdot \{x\}
```

Here $`\lambda = 0`$ is the eigenvalue, and $`\{x\}`$ (a vector normal to the plane) is the eigenvector.

### 1.4 Example: Reflection Matrix

Consider the reflection matrix:

```math
A = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}
```

This is the reflection across the line $`x_1 = x_2`$. Its eigenvectors are $`\begin{Bmatrix} 1 \\ 1 \end{Bmatrix}`$ (eigenvalue 1) and $`\begin{Bmatrix} 1 \\ -1 \end{Bmatrix}`$ (eigenvalue $`-1`$).

The two eigenvectors **span the whole space**, meaning any vector in $`\mathbb{R}^2`$ can be written as a linear combination of these eigenvectors.

> **[Linear Algebra]** When the eigenvectors of an $`n \times n`$ matrix span the entire $`\mathbb{R}^n`$, the matrix is said to be **diagonalizable**. Real symmetric matrices are always diagonalizable by an orthogonal matrix and have real eigenvalues, a key result known as the **spectral theorem**.

---

<br>

## 2. Computing Eigenvalues and Eigenvectors

### 2.1 The Characteristic Equation

Starting from the eigenvalue equation:

```math
[A]\{x\} = \lambda \{x\}
```

Rearranging:

```math
[A]\{x\} - \lambda [I]\{x\} = \{0\}
```

```math
[A - \lambda I]\{x\} = \{0\}
```

For a **non-trivial solution** $`\{x\} \neq \{0\}`$ to exist, the matrix $`[A - \lambda I]`$ must be non-invertible, that is, **singular**. This requires:

```math
\det(A - \lambda I) = |A - \lambda I| = 0
```

This is called the **characteristic equation**, and the determinant on the left, viewed as a polynomial in $`\lambda`$, is the characteristic polynomial. For an $`n \times n`$ matrix, the characteristic polynomial has degree $`n`$ in $`\lambda`$.

### 2.2 The 2x2 Case

For a $`2 \times 2`$ matrix:

```math
\begin{vmatrix} a_{11} - \lambda & a_{12} \\ a_{21} & a_{22} - \lambda \end{vmatrix} = 0
```

Expanding:

```math
\lambda^2 - (a_{11} + a_{22})\lambda + (a_{11}a_{22} - a_{12}a_{21}) = 0
```

Using the quadratic formula:

```math
\lambda = \frac{(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}
```

> **[Linear Algebra]** The quantity $`a_{11} + a_{22}`$ is the **trace** of $`A`$, and $`a_{11}a_{22} - a_{12}a_{21}`$ is the **determinant** of $`A`$. For any $`2 \times 2`$ matrix, the sum of eigenvalues equals the trace and the product of eigenvalues equals the determinant, and these relations hold for general $`n \times n`$ matrices as well (Chapter 12 Lab, Section 1.1).

### 2.3 Example: Finding Eigenvalues and Eigenvectors

Consider:

```math
A = \begin{bmatrix} 3 & 1 \\ 1 & 3 \end{bmatrix}
```

**Step 1: Characteristic equation**

```math
|A - \lambda I| = \begin{vmatrix} 3 - \lambda & 1 \\ 1 & 3 - \lambda \end{vmatrix} = (3 - \lambda)^2 - 1 = \lambda^2 - 6\lambda + 8 = (\lambda - 4)(\lambda - 2) = 0
```

Therefore: $`\lambda_1 = 4`$ and $`\lambda_2 = 2`$.

**Step 2: Eigenvector for $`\lambda_1 = 4`$**

```math
\begin{bmatrix} 3 - 4 & 1 \\ 1 & 3 - 4 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
-x_1 + x_2 = 0 \quad \text{and} \quad x_1 - x_2 = 0 \quad \Rightarrow \quad \text{same equation}
```

An eigenvector is determined only up to a scalar multiple, so let $`x_1 = 1`$; then $`x_2 = 1`$ and the eigenvector is:

```math
\{x\}_1 = \begin{Bmatrix} 1 \\ 1 \end{Bmatrix}
```

**Step 3: Eigenvector for $`\lambda_2 = 2`$**

```math
\begin{bmatrix} 3 - 2 & 1 \\ 1 & 3 - 2 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
x_1 + x_2 = 0
```

Let $`x_1 = 1`$; then $`x_2 = -1`$ and the eigenvector is:

```math
\{x\}_2 = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}
```

---

<br>

## 3. Applications: Differential Equations

### 3.1 First-Order ODE: Scalar Case

Consider the scalar ODE:

```math
\frac{dx}{dt} = ax, \quad a \in \mathbb{R}
```

The solution is:

```math
x(t) = e^{at} C, \quad C \text{ is a constant}
```

The behavior depends on $`a`$ as shown below. The last row extends $`a`$ to complex values.

| Value of $`a`$ | Solution behavior |
|---|---|
| $`a > 0`$ (real, positive) | Exponential growth (unstable) |
| $`a < 0`$ (real, negative) | Exponential decay (stable) |
| $`a = 0`$ | Constant: $`x(t) = C`$ |
| $`a = i\beta`$ (imaginary) | Oscillation: $`x(t) = e^{i\beta t}C = (\cos \beta t + i\sin \beta t)C`$ |

*Table 1. Behavior of the scalar ODE solution by the value of a*

### 3.2 Converting 2nd-Order to 1st-Order System

Consider the following second-order ODE with constant coefficients:

```math
\frac{d^2 y}{dt^2} + a\frac{dy}{dt} + by = 0
```

It can be converted to a system of **first-order** ODEs by the substitution:

```math
x_1 = y, \quad x_2 = \frac{dx_1}{dt}
```

This gives:

```math
\frac{dx_1}{dt} = x_2
```

```math
\frac{dx_2}{dt} = -bx_1 - ax_2
```

In matrix form:

```math
\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 0 & 1 \\ -b & -a \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

We have converted one 2nd-order ODE into **two first-order** ODEs.

### 3.3 Matrix Exponential and Stability

For the vector ODE:

```math
\frac{d\mathbf{x}}{dt} = A\mathbf{x}, \quad \mathbf{x} \in \mathbb{R}^n, \quad A \in \mathbb{R}^{n \times n}
```

The solution is:

```math
\mathbf{x}(t) = e^{At}\mathbf{C}, \quad \mathbf{C} \in \mathbb{R}^n
```

Here $`e^{At}`$ is the **matrix exponential**. If $`A`$ is symmetric, it can be diagonalized as $`A = U D U^T`$, where $`D`$ is the diagonal matrix of eigenvalues and $`U`$ has orthonormal eigenvectors as its columns, and then the following holds (for a general diagonalizable matrix, $`U^{-1}`$ replaces $`U^T`$):

```math
e^{At} = U \, e^{Dt} \, U^T
```

This means:

- **Negative real eigenvalues** produce stable (decaying) solutions.
- **Positive real eigenvalues** produce unstable (growing) solutions.
- **Imaginary eigenvalues** produce oscillatory solutions.

### 3.4 Eigenvalue Location and System Behavior

The complex plane provides an intuitive map for stability:

| Region | Eigenvalue type | System behavior |
|---|---|---|
| Left half-plane ($`\text{Re}(\lambda) < 0`$) | Negative real part | **Stable** (solution decays) |
| Right half-plane ($`\text{Re}(\lambda) > 0`$) | Positive real part | **Unstable** (solution grows) |
| Imaginary axis ($`\text{Re}(\lambda) = 0`$) | Pure imaginary | **Oscillation** (bounded periodic motion) |
| Far left | Large negative real | Faster convergence |
| Far right | Large positive real | Faster divergence |

*Table 2. Eigenvalue location in the complex plane and system behavior*

---

<br>

## 4. Example 13.2: Stability of a First-Order System

Consider the system:

```math
\frac{dx_1}{dt} = -3x_1 + x_2
```

```math
\frac{dx_2}{dt} = x_1 - 3x_2
```

In matrix form:

```math
\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -3 & 1 \\ 1 & -3 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

**Characteristic equation:**

```math
\begin{vmatrix} -3 - \lambda & 1 \\ 1 & -3 - \lambda \end{vmatrix} = (3 + \lambda)^2 - 1 = \lambda^2 + 6\lambda + 8 = (\lambda + 4)(\lambda + 2) = 0
```

```math
\therefore \lambda = -2 \text{ or } -4
```

**Both eigenvalues are negative real numbers**, so the solutions will be **stable**, and since they have no imaginary part, the solutions die out to zero **without oscillation**.

---

<br>

## 5. Eigenvalues and Ordinary Differential Equations

### 5.1 Pure Oscillations: Single 2nd-Order ODE

Consider the undamped oscillator ($`a > 0`$):

```math
\ddot{y} = -ay \quad \Longleftrightarrow \quad \frac{d^2 y}{dt^2} = -ay
```

Cast into first-order ODEs with $`x_1 = y`$, $`x_2 = \frac{dx_1}{dt}`$:

```math
\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 0 & 1 \\ -a & 0 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

**Finding eigenvalues of $`A`$:**

```math
|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ -a & -\lambda \end{vmatrix} = \lambda^2 + a = 0
```

```math
\therefore \lambda = \pm \sqrt{-a} = \pm \sqrt{a}\, i
```

The eigenvalues are **pure imaginary** (lie on the imaginary axis), so the solution will be **oscillatory** with frequency $`\sqrt{a}/2\pi`$.

Euler's formula and the relations for the period and frequency are:

```math
e^{i\sqrt{a}\,t} = \cos(\sqrt{a}\,t) + i\sin(\sqrt{a}\,t)
```

```math
\text{Period} = \frac{2\pi}{\sqrt{a}}, \qquad \text{Frequency} = \frac{1}{\text{Period}} = \frac{\sqrt{a}}{2\pi}
```

> **Note:** **Eigenvalues and frequency.** Eigenvalues are directly related to the **frequency of oscillations**. This is why eigenvalue problems appear throughout physics and engineering whenever we analyze vibrating systems.

### 5.2 Eigenvectors of the Oscillation System

Solve $`[A - \lambda I]\{x\} = \{0\}`$:

```math
\begin{bmatrix} -\lambda & 1 \\ -a & -\lambda \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

**Case i)** $`\lambda = \sqrt{a}\, i`$:

```math
-\sqrt{a}\,i \cdot x_1 + x_2 = 0
```

Take $`x_1 = 1`$; then $`x_2 = \sqrt{a}\,i`$:

```math
\{x\} = \begin{Bmatrix} 1 \\ \sqrt{a}\,i \end{Bmatrix}
```

**Case ii)** $`\lambda = -\sqrt{a}\, i`$:

```math
\sqrt{a}\,i \cdot x_1 + x_2 = 0
```

Take $`x_1 = 1`$; then $`x_2 = -\sqrt{a}\,i`$:

```math
\{x\} = \begin{Bmatrix} 1 \\ -\sqrt{a}\,i \end{Bmatrix}
```

---

<br>

## 6. Coupled 2nd-Order ODE Systems

### 6.1 General Form

A system of two coupled second-order ODEs:

```math
\frac{d^2 y_1}{dt^2} = -a_{11}y_1 + a_{12}y_2
```

```math
\frac{d^2 y_2}{dt^2} = a_{21}y_1 - a_{22}y_2
```

In matrix form:

```math
\frac{d^2}{dt^2}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix} \begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} \quad \cdots \quad (*)
```

### 6.2 Proposed Periodic Solution

Based on the periodic solution form, we propose:

```math
\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t}
```

where $`\omega`$ is the angular frequency. Substituting into $`(*)`$:

```math
-\omega^2 \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t} = \begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t}
```

Cancel $`e^{i\omega t}`$:

```math
\underbrace{-\omega^2}_{\lambda} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} = \underbrace{\begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix}}_{A} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix}
```

**This is an eigenvalue problem.** The eigenvalue is $`\lambda = -\omega^2`$, so $`\omega^2 = -\lambda`$.

### 6.3 Characteristic Equation for 2nd-Order Systems

From $`|A - \lambda I| = 0`$:

```math
\begin{vmatrix} -a_{11} - \lambda & a_{12} \\ a_{21} & -a_{22} - \lambda \end{vmatrix} = 0
```

```math
(a_{11} + \lambda)(a_{22} + \lambda) - a_{12}a_{21} = 0
```

```math
\lambda^2 + (a_{11} + a_{22})\lambda + a_{11}a_{22} - a_{12}a_{21} = 0
```

```math
\therefore \lambda = \frac{-(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}
```

Equivalently, in terms of $`\omega^2 = -\lambda`$:

```math
|A + \omega^2 I| = \begin{vmatrix} -a_{11} + \omega^2 & a_{12} \\ a_{21} & -a_{22} + \omega^2 \end{vmatrix} = 0
```

```math
(\omega^2 - a_{11})(\omega^2 - a_{22}) - a_{12}a_{21} = 0
```

```math
\omega^4 - (a_{11} + a_{22})\omega^2 + a_{11}a_{22} - a_{12}a_{21} = 0
```

```math
\therefore \omega^2 = \frac{(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}
```

---

<br>

## 7. Example 13.3: Behavior of a Coupled System

### 7.1 Eigenvalue Computation

Given:

```math
\frac{d^2 y_1}{dt^2} = -5y_1 + 2y_2, \qquad \frac{d^2 y_2}{dt^2} = 2y_1 - 2y_2
```

Matrix form:

```math
\frac{d^2}{dt^2}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \underbrace{\begin{bmatrix} -5 & 2 \\ 2 & -2 \end{bmatrix}}_{A}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix}
```

Characteristic equation:

```math
\begin{vmatrix} -5 - \lambda & 2 \\ 2 & -2 - \lambda \end{vmatrix} = (5 + \lambda)(2 + \lambda) - 4 = \lambda^2 + 7\lambda + 6 = (\lambda + 6)(\lambda + 1) = 0
```

```math
\therefore \lambda_1 = -1, \quad \lambda_2 = -6
```

### 7.2 Eigenvector Computation

**For $`\lambda_1 = -1`$:**

```math
\begin{bmatrix} -5 + 1 & 2 \\ 2 & -2 + 1 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix} \quad \Rightarrow \quad -4x_1 + 2x_2 = 0 \quad \Rightarrow \quad 2x_1 - x_2 = 0
```

Take $`x_1 = 1`$, $`x_2 = 2`$; the eigenvector is:

```math
\{x\}_1 = \begin{Bmatrix} 1 \\ 2 \end{Bmatrix}
```

**For $`\lambda_2 = -6`$:**

```math
\begin{bmatrix} -5 + 6 & 2 \\ 2 & -2 + 6 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix} \quad \Rightarrow \quad x_1 + 2x_2 = 0
```

Take $`x_2 = -1`$, $`x_1 = 2`$; the eigenvector is:

```math
\{x\}_2 = \begin{Bmatrix} 2 \\ -1 \end{Bmatrix}
```

### 7.3 Frequencies and General Solution

From the eigenvalues, $`\omega^2 = -\lambda`$:

```math
\omega_1^2 = 1 \quad \Rightarrow \quad \omega_1 = \pm 1
```

```math
\omega_2^2 = 6 \quad \Rightarrow \quad \omega_2 = \pm \sqrt{6}
```

The solution behaves as:

```math
\{x\}_1 e^{\pm it} = \{x\}_1 (\cos t \pm i\sin t)
```

```math
\{x\}_2 e^{\pm i\sqrt{6}t} = \{x\}_2 (\cos\sqrt{6}t \pm i\sin\sqrt{6}t)
```

General solution:

```math
\mathbf{y} = C_1 \{x\}_1 e^{+it} + C_2 \{x\}_2 e^{+i\sqrt{6}t} + C_3 \{x\}_1 e^{-it} + C_4 \{x\}_2 e^{-i\sqrt{6}t}
```

Its time derivative:

```math
\dot{\mathbf{y}} = iC_1 \{x\}_1 e^{+it} + i\sqrt{6}C_2 \{x\}_2 e^{+i\sqrt{6}t} - iC_3 \{x\}_1 e^{-it} - i\sqrt{6}C_4 \{x\}_2 e^{-i\sqrt{6}t}
```

### 7.4 Applying Initial Conditions

Given initial conditions:

```math
\mathbf{y}(t=0) = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}, \quad \dot{\mathbf{y}}(t=0) = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

At $`t = 0`$, the velocity is $`i(C_1 - C_3)\{x\}_1 + i\sqrt{6}(C_2 - C_4)\{x\}_2 = \{0\}`$, and the two eigenvectors are linearly independent, so we require $`C_1 = C_3`$ and $`C_2 = C_4`$. Let $`\tilde{C}_1 = C_1 = C_3`$ and $`\tilde{C}_2 = C_2 = C_4`$.

Substituting $`t = 0`$ into the position equations gives $`\mathbf{y}(0) = 2\tilde{C}_1\{x\}_1 + 2\tilde{C}_2\{x\}_2`$.

From $`y_1`$: $`1 = 2\tilde{C}_1 + 4\tilde{C}_2`$

From $`y_2`$: $`-1 = 4\tilde{C}_1 - 2\tilde{C}_2`$

Multiplying the second equation by 2 and writing it alongside the first:

```math
-2 = 8\tilde{C}_1 - 4\tilde{C}_2
```

```math
1 = 2\tilde{C}_1 + 4\tilde{C}_2
```

Adding the two equations: $`-1 = 10\tilde{C}_1 \Rightarrow \tilde{C}_1 = -\frac{1}{10}`$

Substituting into $`1 = 2(-\frac{1}{10}) + 4\tilde{C}_2`$ gives $`\tilde{C}_2 = \frac{1}{4}(1 + \frac{1}{5}) = \frac{3}{10}`$.

**Final solution:** Using $`e^{it} + e^{-it} = 2\cos t`$:

```math
\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = -\frac{1}{5}\begin{Bmatrix} 1 \\ 2 \end{Bmatrix}\cos(t) + \frac{3}{5}\begin{Bmatrix} 2 \\ -1 \end{Bmatrix}\cos(\sqrt{6}\,t)
```

---

<br>

## 8. Physical Settings: Mass-Spring Systems

### 8.1 Two-Mass, Three-Spring System

Consider a two-mass, three-spring system with frictionless rollers, vibrating between two walls.

```text
Wall ─[k]─ m₁ ─[k]─ m₂ ─[k]─ Wall
```

*Figure 1. Two-mass, three-spring system between two walls*

Two masses $`m_1`$ and $`m_2`$ are connected by springs of stiffness $`k`$. Newton's second law for each mass:

```math
m_1 \frac{d^2 x_1}{dt^2} = -kx_1 + k(x_2 - x_1)
```

```math
m_2 \frac{d^2 x_2}{dt^2} = -k(x_2 - x_1) - kx_2
```

In matrix form:

```math
\frac{d^2}{dt^2}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -2k/m_1 & k/m_1 \\ k/m_2 & -2k/m_2 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

### 8.2 Eigenvalue Formula for the Mass-Spring System

From $`|A - \lambda I| = 0`$:

```math
\begin{vmatrix} -2k/m_1 - \lambda & k/m_1 \\ k/m_2 & -2k/m_2 - \lambda \end{vmatrix} = 0
```

```math
(\lambda + 2k/m_1)(\lambda + 2k/m_2) - k^2/(m_1 m_2) = 0
```

```math
\lambda^2 + \lambda(2k/m_1 + 2k/m_2) + 4k^2/(m_1 m_2) - k^2/(m_1 m_2) = 0
```

Solving via the quadratic formula:

```math
-\omega^2 = \lambda = -\frac{k}{m_1 m_2}(m_1 + m_2) \pm \frac{k}{m_1 m_2}\sqrt{m_2^2 + m_1^2 - m_1 m_2}
```

---

<br>

## 9. Example 13.4: Mass-Spring System

Given: $`m_1 = m_2 = 40`$ kg, $`k = 200`$ N/m.

**Eigenvalue computation:**

```math
\lambda = -\frac{k}{m_1 m_2}(m_1 + m_2) \pm \frac{k}{m_1 m_2}\sqrt{m_2^2 + m_1^2 - m_1 m_2}
```

```math
= -\frac{200}{1600}(80) \pm \frac{1}{8}\sqrt{1600}
```

```math
= -10 \pm 5
```

```math
\therefore \lambda_1 = -5, \quad \lambda_2 = -15
```

**Frequencies:**

```math
\omega^2 = -\lambda \quad \Rightarrow \quad \omega_1 = \sqrt{5}, \quad f_1 = \frac{\sqrt{5}}{2\pi}
```

```math
\omega_2 = \sqrt{15}, \quad f_2 = \frac{\sqrt{15}}{2\pi}
```

Here $`\omega = 2\pi f`$.

**Eigenvectors:**

For $`\lambda_1 = -5`$, with $`A = \begin{bmatrix} -10 & 5 \\ 5 & -10 \end{bmatrix}`$:

```math
\begin{bmatrix} -10 - (-5) & 5 \\ 5 & -10 - (-5) \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -5 & 5 \\ 5 & -5 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
-x_1 + x_2 = 0 \quad \Rightarrow \quad \{x\}_1 = \begin{Bmatrix} 1 \\ 1 \end{Bmatrix}
```

Both masses move **in the same direction** (in-phase mode).

For $`\lambda_2 = -15`$:

```math
\begin{bmatrix} -10 + 15 & 5 \\ 5 & -10 + 15 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 5 & 5 \\ 5 & 5 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
x_1 + x_2 = 0 \quad \Rightarrow \quad \{x\}_2 = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}
```

The masses move in **opposite directions** (out-of-phase mode).

**General solution:**

```math
\mathbf{y} = C_1 \{x\}_1 \cos\sqrt{5}\,t + C_2 \{x\}_1 \sin\sqrt{5}\,t + C_3 \{x\}_2 \cos\sqrt{15}\,t + C_4 \{x\}_2 \sin\sqrt{15}\,t
```

The coefficients $`C_1, C_2, C_3, C_4`$ are determined by the initial conditions.

Expanding for each component:

```math
y_1 = C_1 \cdot 1 \cdot \cos\sqrt{5}\,t + C_2 \cdot 1 \cdot \sin\sqrt{5}\,t + C_3 \cdot 1 \cdot \cos\sqrt{15}\,t + C_4 \cdot 1 \cdot \sin\sqrt{15}\,t
```

```math
y_2 = C_1 \cdot 1 \cdot \cos\sqrt{5}\,t + C_2 \cdot 1 \cdot \sin\sqrt{5}\,t + C_3 \cdot (-1) \cdot \cos\sqrt{15}\,t + C_4 \cdot (-1) \cdot \sin\sqrt{15}\,t
```

---

<br>

## 10. The Power Method

The power method **iteratively** determines the dominant eigenvalue, the one **largest** in absolute value, and its corresponding eigenvector.

### 10.1 Three-Mass, Four-Spring System (Ex 13.5)

Consider a three-mass, four-spring system between two walls, with $`m = m_1 = m_2 = m_3 = 1`$ kg and $`k = 20`$ N/m.

```text
Wall ─[k]─ m₁ ─[k]─ m₂ ─[k]─ m₃ ─[k]─ Wall
```

*Figure 2. Three-mass, four-spring system between two walls*

Newton's second law for each mass:

```math
m_1 \frac{d^2 x_1}{dt^2} = -kx_1 + k(x_2 - x_1)
```

```math
m_2 \frac{d^2 x_2}{dt^2} = -k(x_2 - x_1) + k(x_3 - x_2)
```

```math
m_3 \frac{d^2 x_3}{dt^2} = -k(x_3 - x_2) - kx_3
```

Assuming a periodic solution $`\{x_i\} = \{X_i\}e^{i\omega t}`$ and substituting:

```math
\begin{bmatrix} 40 - \omega^2 & -20 & 0 \\ -20 & 40 - \omega^2 & -20 \\ 0 & -20 & 40 - \omega^2 \end{bmatrix} \begin{Bmatrix} X_1 \\ X_2 \\ X_3 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \\ 0 \end{Bmatrix}
```

This is equivalent to the eigenvalue problem $`A\{X\} = \lambda\{X\}`$ with the following $`A`$. In this section, the sign is flipped so that $`\lambda = \omega^2`$.

```math
A = \begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix}, \quad \lambda = \omega^2
```

### 10.2 Algorithm

1. **Start** with an initial guess vector $`\{x\}^{(0)}`$.
2. **Multiply:** compute $`\{y\}^{(k)} = [A]\{x\}^{(k-1)}`$.
3. **Normalize:** take the element largest in absolute value as the eigenvalue estimate $`\lambda^{(k)}`$, and divide $`\{y\}^{(k)}`$ by this element to get the new normalized vector $`\{x\}^{(k)}`$.
4. **Repeat** until the eigenvalue estimate converges.

The normalization factor at each step is our **estimate of the dominant eigenvalue**.

### 10.3 Iteration Walkthrough

**Iteration 1:** Take $`X_1 = X_2 = X_3 = 1`$ as the initial guess.

```math
\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ 1 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 20 \\ 0 \\ 20 \end{Bmatrix} = 20 \begin{Bmatrix} 1 \\ 0 \\ 1 \end{Bmatrix}
```

The eigenvalue estimate is $`\lambda \approx 20`$.

**Iteration 2:** Use $`\{1, 0, 1\}^T`$.

```math
\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ 0 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 40 \\ -40 \\ 40 \end{Bmatrix} = 40 \begin{Bmatrix} 1 \\ -1 \\ 1 \end{Bmatrix}
```

The eigenvalue estimate is $`\lambda \approx 40`$, and the approximate error is $`\varepsilon_a = \left|\frac{40 - 20}{40}\right| \times 100\% = 50\%`$.

**Iteration 3:** Use $`\{1, -1, 1\}^T`$.

```math
\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ -1 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 60 \\ -80 \\ 60 \end{Bmatrix} = -80 \begin{Bmatrix} -3/4 \\ 1 \\ -3/4 \end{Bmatrix}
```

The eigenvalue estimate is $`\lambda \approx -80`$, and the error is $`\varepsilon_a = \left|\frac{-80 - 40}{-80}\right| \times 100\% = 150\%`$.

**Iteration 4:** Use $`\{-3/4, 1, -3/4\}^T`$.

```math
A\begin{Bmatrix} -3/4 \\ 1 \\ -3/4 \end{Bmatrix} = \begin{Bmatrix} -50 \\ 70 \\ -50 \end{Bmatrix} = 70\begin{Bmatrix} -5/7 \\ 1 \\ -5/7 \end{Bmatrix}
```

The eigenvalue estimate is $`\lambda \approx 70`$, and the error is $`\varepsilon_a = \left|\frac{70 - (-80)}{70}\right| \times 100\% \approx 214\%`$.

**Iteration 5:** Use $`\{-5/7, 1, -5/7\}^T`$.

```math
A\begin{Bmatrix} -5/7 \\ 1 \\ -5/7 \end{Bmatrix} = \begin{Bmatrix} -48.57 \\ 68.57 \\ -48.57 \end{Bmatrix} = 68.57\begin{Bmatrix} -0.71 \\ 1 \\ -0.71 \end{Bmatrix}
```

The eigenvalue estimate is $`\lambda \approx 68.57`$, and the error is $`\varepsilon_a = \left|\frac{68.57 - 70}{68.57}\right| \times 100\% \approx 2.08\%`$.

### 10.4 Convergence and Finding the Smallest Eigenvalue

The eigenvalue is converging. After sufficient iterations:

```math
\lambda = 68.28427..., \quad \{x\} = \begin{Bmatrix} -0.70711 \\ 1 \\ -0.70711 \end{Bmatrix}
```

This value equals $`40 + 20\sqrt{2}`$, and the angular frequency of the highest mode is $`\omega = \sqrt{68.284} \approx 8.263`$ rad/s.

> **Note:** **Convergence rate of the power method.** The power method converges to the eigenvalue with the **largest absolute value**. The rate of convergence depends on the ratio $`|\lambda_2/\lambda_1|`$, where $`\lambda_1`$ is the dominant eigenvalue and $`\lambda_2`$ is the second-largest in absolute value. The closer this ratio is to 1, the slower the convergence; in this example it is $`40/68.284 \approx 0.586`$. The starting vector must also have a component along the eigenvector of $`\lambda_1`$.

**For finding the smallest eigenvalue**, one can apply the power method to the **inverse** of $`A`$:

```math
[A]^{-1}\{x\} = \frac{1}{\lambda}\{x\}
```

Since the eigenvalue of $`A^{-1}`$ largest in absolute value is $`1/\lambda_{\min}`$, the power method applied to $`A^{-1}`$ converges to $`1/\lambda_{\min}`$, where $`\lambda_{\min}`$ is the eigenvalue smallest in absolute value. This is called the **inverse power method**.

---

<br>

## 11. Python Implementation

### 11.1 Eigenvalue computation using NumPy

Applying `np.linalg.eig` to the matrix of Section 2.3 gives the eigenvalues 4 and 2 together with eigenvectors normalized to unit length in each column.

```python
import numpy as np

# Example: 2x2 matrix from Section 2.3
A = np.array([[3, 1],
              [1, 3]])

eigenvalues, eigenvectors = np.linalg.eig(A)
print("Eigenvalues:", eigenvalues)       # [4. 2.]
print("Eigenvectors:\n", eigenvectors)   # columns are eigenvectors
```

### 11.2 Power Method implementation

The function implements the algorithm of Section 10.2 directly, normalizing by the element largest in absolute value. Applied to the matrix of Example 13.5, it gives the dominant eigenvalue 68.28428 and the eigenvector $`[-0.70710681, 1, -0.70710681]`$ after 11 iterations.

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

### 11.3 Mass-Spring System (Ex 13.4)

The coefficient matrix of Section 9 gives the eigenvalues $`-5`$ and $`-15`$, and $`\omega = \sqrt{-\lambda}`$ gives the angular frequencies 2.23607 and 3.87298 rad/s and the frequencies 0.35588 and 0.61640 Hz.

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

### 11.4 Inverse Power Method for smallest eigenvalue

Instead of forming the inverse, each iteration solves $`A\mathbf{y} = \mathbf{x}`$ with `np.linalg.solve`. Applied to the matrix of Example 13.5, it gives the smallest eigenvalue 11.71573 ($`= 40 - 20\sqrt{2}`$) and the eigenvector $`[0.70710681, 1, 0.70710681]`$ after 9 iterations. Because `np.linalg.solve` factors the matrix again on every call, factoring once by LU as in the Chapter 12 Lab is more efficient for large matrices.

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

| Concept | Key Summary |
|---|---|
| Eigenvalue equation | A nonzero $`\{x\}`$ satisfying $`[A]\{x\} = \lambda\{x\}`$ is an eigenvector, and $`\lambda`$ is its eigenvalue. |
| Characteristic equation | The eigenvalues are the roots of $`\det(A - \lambda I) = 0`$. |
| 2x2 eigenvalues | The eigenvalues of a $`2 \times 2`$ matrix are $`\lambda = \frac{(a_{11}+a_{22}) \pm \sqrt{(a_{11}+a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}`$. |
| Non-trivial solution condition | A non-trivial solution exists only if $`[A - \lambda I]`$ is singular. |
| Eigenvector property | An eigenvector $`\{x\}`$ is parallel to $`[A]\{x\}`$. |
| Scalar ODE $`dx/dt = ax`$ | The solution is $`x(t) = Ce^{at}`$. |
| Matrix ODE $`d\mathbf{x}/dt = A\mathbf{x}`$ | The solution is $`\mathbf{x}(t) = e^{At}\mathbf{C}`$, and $`e^{At} = Ue^{Dt}U^T`$ when $`A`$ is symmetric. |
| Stability | The system is stable if $`\text{Re}(\lambda) < 0`$, unstable if $`\text{Re}(\lambda) > 0`$, and oscillatory if the eigenvalues are pure imaginary ($`\text{Re}(\lambda) = 0`$, $`\lambda \neq 0`$). |
| 2nd-order ODE to eigenvalue | For $`\ddot{y} = Ay`$, $`\lambda = -\omega^2`$, so $`\omega = \sqrt{-\lambda}`$. |
| Frequency relation | The period is $`2\pi/\omega`$ and the frequency is $`f = \omega/(2\pi)`$. |
| Mass-spring system | Equations of motion of the form $`m\ddot{x} = -kx + k(x_{\text{neighbor}} - x)`$ lead to an eigenvalue problem. |
| Power Method | Repeated multiplication and normalization compute the dominant eigenvalue, the one largest in absolute value. |
| Inverse Power Method | Applying the power method to $`A^{-1}`$ finds the eigenvalue smallest in absolute value. |
| Convergence rate | Convergence of the power method depends on $`\vert\lambda_2/\lambda_1\vert`$ and slows as this ratio approaches 1. |

---

<br>

## Self-Check Questions

1. **Eigenvalues and eigenvectors:** Find the eigenvalues and eigenvectors of $`A = \begin{bmatrix} 4 & 1 \\ 2 & 3 \end{bmatrix}`$, and check them with the trace and determinant.

   <details>
   <summary>Answer</summary>

   The characteristic equation is $`(4-\lambda)(3-\lambda) - 2 = \lambda^2 - 7\lambda + 10 = (\lambda - 5)(\lambda - 2) = 0`$, so $`\lambda_1 = 5`$ and $`\lambda_2 = 2`$. For $`\lambda_1 = 5`$, $`-x_1 + x_2 = 0`$ gives $`\{x\}_1 = \{1, 1\}^T`$, and for $`\lambda_2 = 2`$, $`2x_1 + x_2 = 0`$ gives $`\{x\}_2 = \{1, -2\}^T`$. The sum of the eigenvalues, 7, equals the trace $`4 + 3`$, and their product, 10, equals the determinant $`4 \cdot 3 - 1 \cdot 2`$.

   </details>

2. **Stability:** Find the eigenvalues of $`\frac{d\mathbf{x}}{dt} = \begin{bmatrix} -1 & 2 \\ -2 & -1 \end{bmatrix}\mathbf{x}`$ and describe the behavior of the solution.

   <details>
   <summary>Answer</summary>

   The characteristic equation is $`(-1-\lambda)^2 + 4 = 0`$, so $`\lambda = -1 \pm 2i`$. The real part $`-1 < 0`$ makes the solution stable with decay proportional to $`e^{-t}`$, and the imaginary part makes it oscillate with angular frequency 2 (period $`\pi`$). The solution therefore oscillates while converging to 0.

   </details>

3. **Mass-spring system:** For the two-mass, three-spring system of Section 8 with $`m_1 = m_2 = 10`$ kg and $`k = 100`$ N/m, find the eigenvalues, angular frequencies, and mode shapes.

   <details>
   <summary>Answer</summary>

   Since $`k/(m_1 m_2) = 100/100 = 1`$, the formula of Section 8.2 gives $`\lambda = -1 \cdot 20 \pm 1 \cdot \sqrt{100 + 100 - 100} = -20 \pm 10`$, that is, $`\lambda_1 = -10`$ and $`\lambda_2 = -30`$. The angular frequencies are $`\omega_1 = \sqrt{10} \approx 3.162`$ rad/s and $`\omega_2 = \sqrt{30} \approx 5.477`$ rad/s. The coefficient matrix is $`\begin{bmatrix} -20 & 10 \\ 10 & -20 \end{bmatrix}`$, so as in Example 13.4 the mode of $`\lambda_1`$ is $`\{1, 1\}^T`$ (in phase) and the mode of $`\lambda_2`$ is $`\{1, -1\}^T`$ (out of phase).

   </details>

4. **Inverse power method and the lowest mode:** To which eigenvalue does the inverse power method converge for the matrix of Example 13.5? Find the corresponding angular frequency and compare it with the value found by the power method.

   <details>
   <summary>Answer</summary>

   The eigenvalues are $`40 - 20\sqrt{2} \approx 11.716`$, $`40`$, and $`40 + 20\sqrt{2} \approx 68.284`$, so the inverse power method converges to 11.716, the one smallest in absolute value. In this example $`\lambda = \omega^2`$, so the angular frequency of the lowest mode is $`\omega = \sqrt{11.716} \approx 3.423`$ rad/s, which is smaller than the angular frequency of the highest mode found by the power method, $`\sqrt{68.284} \approx 8.263`$ rad/s.

   </details>
