# Chapter 7 Lecture — Optimization

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. Introduction to Optimization](#1-introduction-to-optimization)
  - [1.1 What Is Optimization?](#11-what-is-optimization)
  - [1.2 Motivating Example: Aluminum Can Design](#12-motivating-example-aluminum-can-design)
  - [1.3 Mathematical Description](#13-mathematical-description)
  - [1.4 Minimization vs. Maximization](#14-minimization-vs-maximization)
  - [1.5 Feasible Region and Constraints](#15-feasible-region-and-constraints)
- [2. Key Topological Concepts](#2-key-topological-concepts)
  - [2.1 Interior Points](#21-interior-points)
  - [2.2 Boundary Points](#22-boundary-points)
  - [2.3 Feasible Directions](#23-feasible-directions)
- [3. Local and Global Optimality](#3-local-and-global-optimality)
  - [3.1 Global Minimum](#31-global-minimum)
  - [3.2 Infimum and Supremum](#32-infimum-and-supremum)
  - [3.3 Local Minimum](#33-local-minimum)
  - [3.4 Strict Local Minimum](#34-strict-local-minimum)
- [4. Existence of an Optimal Solution](#4-existence-of-an-optimal-solution)
  - [4.1 Weierstrass Theorem](#41-weierstrass-theorem)
  - [4.2 Coercive Functions](#42-coercive-functions)
  - [4.3 Proof Sketch for Coercive Functions](#43-proof-sketch-for-coercive-functions)
- [5. One-Dimensional Optimization](#5-one-dimensional-optimization)
- [6. Golden-Section Search](#6-golden-section-search)
  - [6.1 The Golden Ratio](#61-the-golden-ratio)
  - [6.2 Algorithm: Interval Reduction](#62-algorithm-interval-reduction)
  - [6.3 Reuse of Function Evaluations](#63-reuse-of-function-evaluations)
  - [6.4 Error Analysis](#64-error-analysis)
- [7. Parabolic Interpolation](#7-parabolic-interpolation)
  - [7.1 Constructing the Parabola](#71-constructing-the-parabola)
  - [7.2 Optimal Point Formula](#72-optimal-point-formula)
  - [7.3 Updating Points for the Next Iteration](#73-updating-points-for-the-next-iteration)
- [8. Level Sets and Gradient](#8-level-sets-and-gradient)
  - [8.1 Level Set Definition](#81-level-set-definition)
  - [8.2 Gradient Is Orthogonal to Level Sets](#82-gradient-is-orthogonal-to-level-sets)
- [9. Convex Problems](#9-convex-problems)
  - [9.1 Convex Combination](#91-convex-combination)
  - [9.2 Convex Set](#92-convex-set)
  - [9.3 Hyperplane and Polyhedral Set](#93-hyperplane-and-polyhedral-set)
  - [9.4 Extreme Point](#94-extreme-point)
  - [9.5 Convex Function](#95-convex-function)
- [10. Steepest Gradient (Descent) Method](#10-steepest-gradient-descent-method)
  - [10.1 Directional Derivative and Cauchy-Schwarz](#101-directional-derivative-and-cauchy-schwarz)
  - [10.2 The Steepest Descent Update Rule](#102-the-steepest-descent-update-rule)
  - [10.3 Convergence and Orthogonality Property](#103-convergence-and-orthogonality-property)
  - [10.4 Example](#104-example)
- [Summary](#summary)

---

<br>

## 1. Introduction to Optimization

### 1.1 What Is Optimization?

**Optimization** is a methodology to find the **best** among all **feasible solutions**. An **optimal solution** is the feasible solution that yields the best objective function value.

A function $f(x)$ can have several types of extreme points:

- **Global maximum**: the largest value of $f$ over the entire domain
- **Global minimum**: the smallest value of $f$ over the entire domain
- **Local maximum**: a point where $f$ is larger than at all nearby points
- **Local minimum**: a point where $f$ is smaller than at all nearby points

### 1.2 Motivating Example: Aluminum Can Design

A company needs to design an aluminum can in the shape of a cylinder with diameter $d$ and height $h$.

**Formulas:**

$$\text{Volume} = \pi \left(\frac{d}{2}\right)^2 h = \frac{\pi d^2 h}{4}$$

$$\text{Surface Area} = \frac{\pi d^2}{2} + \pi d h$$

Both volume and surface area are functions of the diameter $d$ and height $h$.

**Goal:** Design a can that requires the **minimum possible amount of aluminum** (minimize surface area), while satisfying the following conditions:

1. The height must be at least 50% greater than the diameter: $h \geq \frac{3}{2}d \iff \frac{3}{2}d - h \leq 0$
2. The height cannot be more than twice the diameter: $h \leq 2d \iff -2d + h \leq 0$
3. The volume must be at least 330 mL: $\frac{\pi d^2 h}{4} \geq 330$
4. $h > 0, \; d > 0$ (already guaranteed by conditions 1 and 2)

**Summary of the optimization problem:**

$$\text{minimize} \quad \frac{\pi d^2}{2} + \pi d h$$

$$\text{subject to} \quad \frac{3}{2}d - h \leq 0$$

$$-2d + h \leq 0$$

$$\frac{\pi d^2 h}{4} \geq 330$$

### 1.3 Mathematical Description

Let $\mathbb{X} \subseteq \mathbb{R}^n$ be a set of $n$-dimensional real vectors, so that $\mathbb{X} \ni x = [x_1, x_2, \ldots, x_n]^T$.

Let $f : \mathbb{X} \to \mathbb{R}$.

An optimization problem has the general form:

$$\min_{x \in \mathbb{X}} f(x)$$

equivalently written as:

$$\text{minimize} \quad f(x) \quad \text{subject to} \quad x \in \mathbb{X}$$

where:

- $\mathbb{X}$ is the **feasible (admissible) region**
- $f$ is the **objective function**
- $x$ is a **feasible solution** or a **decision variable**

### 1.4 Minimization vs. Maximization

A minimization problem $\min_{x \in \mathbb{X}} f(x)$ can be converted into an equivalent maximization problem:

$$\max_{x \in \mathbb{X}} (-f(x))$$

and vice versa. Therefore, without loss of generality, we can always consider only minimization problems.

### 1.5 Feasible Region and Constraints

The feasible region $\mathbb{X}$ is described by a set of **constraints**:

- **Equality constraints:** $h_i(x) = 0$, where $h_i : \mathbb{R}^n \to \mathbb{R}$, $i \in \mathcal{E}_i$ (index set)
- **Inequality constraints:** $g_j(x) \leq 0$, where $g_j : \mathbb{R}^n \to \mathbb{R}$, $j \in \mathcal{E}_j$ (index set)

**Example:** Consider a feasible set $\mathbb{X}$ defined by two inequality constraints:

$$\mathbb{X} = \{x \in \mathbb{R}^2 : x_1^2 + x_2^2 \leq 1, \; (x_1 - 1)^2 + x_2^2 \leq 1\}$$

Here $g_1(x_1, x_2) = x_1^2 + x_2^2 - 1 \leq 0$ and $g_2(x_1, x_2) = (x_1 - 1)^2 + x_2^2 - 1 \leq 0$. The feasible region is the intersection (overlap) of the two disks.

---

<br>

## 2. Key Topological Concepts

### 2.1 Interior Points

**Definition.** A point $\hat{x} \in \mathbb{X}$ is an **interior point** of $\mathbb{X}$ if there exists $\varepsilon > 0$ such that $\hat{x}$ is included in $\mathbb{X}$ together with the $\varepsilon$-ball centered at $\hat{x}$:

$$B(\hat{x}, \varepsilon) \subset \mathbb{X}$$

**Definition.** $\text{int}(\mathbb{X})$ is the set of all interior points of $\mathbb{X}$.

### 2.2 Boundary Points

**Definition.** A point $\bar{x}$ is a **boundary point** of $\mathbb{X}$ if:

$$\bar{x} \in \mathbb{X} \quad \text{and} \quad \bar{x} \notin \text{int}(\mathbb{X})$$

In other words, $\bar{x}$ belongs to the set, but no $\varepsilon$-ball around it is entirely contained within the set. Geometrically, boundary points lie on the "edge" of the set, while interior points lie "inside" with room to spare in every direction.

### 2.3 Feasible Directions

**Definition.** A vector $d \in \mathbb{R}^n$ is a **feasible direction** for the set $\mathbb{X} \subseteq \mathbb{R}^n$ at $\bar{x} \in \mathbb{X}$ if there exists $\delta > 0$ such that:

$$\bar{x} + \alpha d \in \mathbb{X} \quad \text{for any } \alpha \leq \delta$$

> **[Geometry]** A feasible direction is one in which we can move a small positive step from the current point and remain within the feasible set. At interior points, every direction is feasible. At boundary points, only directions pointing "inward" (or along the boundary) are feasible.

---

<br>

## 3. Local and Global Optimality

Consider a minimization problem:

$$\text{minimize} \quad f(x) \quad \text{subject to} \quad x \in \mathbb{X}$$

where $f : \mathbb{R}^n \to \mathbb{R}$.

### 3.1 Global Minimum

**Definition (Global Minimum).** The point $x^* \in \mathbb{X}$ is a point of **global minimum** (global minimizer) for the problem $\min_{x \in \mathbb{X}} f(x)$, i.e., $x^* = \arg\min_{x \in \mathbb{X}} f(x)$, if:

$$f(x^*) \leq f(x) \quad \forall \; x \in \mathbb{X}$$

A global minimizer $x^*$ is a **strict global minimizer** if:

$$f(x^*) < f(x) \quad \forall \; x \in \mathbb{X} \setminus \{x^*\}$$

### 3.2 Infimum and Supremum

Not every function has a global minimizer. Consider $f(x) = e^x$ with $x \in \mathbb{R}$:

- $f(x) \geq 0$ (bounded from below)
- $\lim_{x \to -\infty} f(x) \to 0$ (greatest lower bound), but $f(x) = 0$ is never achieved

**Definition (Infimum, Supremum).** For a function $f : \mathbb{X} \to \mathbb{R}$:

- Its **greatest lower bound** is called the **infimum**: $\inf_{x \in \mathbb{X}} f(x)$
- Its **least upper bound** is called the **supremum**: $\sup_{x \in \mathbb{X}} f(x)$

**Examples:**

- $\sup_{x \in \mathbb{R}} e^x = +\infty$
- $\inf_{x \in (0, +\infty)} \ln(x) = -\infty$

> **[Analysis]** The minimum of a function is the infimum that is actually attained: $\min f = \inf f$ when there exists some $x^*$ with $f(x^*) = \inf f$. When the infimum is not attained, the minimum does not exist even though the infimum is well-defined.

### 3.3 Local Minimum

**Definition (Local Minimum).** $x^* \in \mathbb{X}$ is a point of **local minimum** (local minimizer) for the problem $\min_{x \in \mathbb{X}} f(x)$ if there exists $\varepsilon > 0$ such that:

$$f(x) \geq f(x^*) \quad \forall \; x \in \mathbb{X} \text{ with } \|x - x^*\| \leq \varepsilon$$

### 3.4 Strict Local Minimum

A local minimizer $x^*$ is a **strict local minimizer** if there exists $\varepsilon > 0$ such that:

$$f(x) > f(x^*) \quad \forall \; x \in \mathbb{X} \setminus \{x^*\} \text{ with } \|x - x^*\| \leq \varepsilon$$

> **[Key Distinction]** A global minimizer has the smallest function value over the **entire** feasible set, while a local minimizer only needs to be optimal within some neighborhood. Every global minimizer is a local minimizer, but the converse is not necessarily true.

---

<br>

## 4. Existence of an Optimal Solution

Checking if a global minimizer exists is an extremely difficult problem in general. However, in some situations, we can guarantee the existence of a global minimizer.

### 4.1 Weierstrass Theorem

**Theorem (Weierstrass).** The problems $\min_{x \in \mathbb{X}} f(x)$ and $\max_{x \in \mathbb{X}} f(x)$, where $\mathbb{X} \subset \mathbb{R}^n$ is a **compact set** (closed and bounded) and $f : \mathbb{X} \to \mathbb{R}$ is a **continuous function**, have **global optimal solutions**.

> **[Topology]** A compact set in $\mathbb{R}^n$ is one that is both closed (contains all its boundary points) and bounded (fits inside a ball of finite radius). The Weierstrass theorem guarantees that continuous functions attain their extreme values on compact sets.

### 4.2 Coercive Functions

**Definition (Coercive Function).** A function $f : \mathbb{R}^n \to \mathbb{R}$ is called **coercive** if:

$$\lim_{\|x\| \to \infty} f(x) = +\infty$$

As $x$ moves arbitrarily far away from the origin, the function value $f(x)$ grows to infinity. Examples include $f(x) = x^2$ (coercive) vs. $f(x) = \sin(x)$ (not coercive).

**Theorem.** Any continuous coercive function $f : \mathbb{R}^n \to \mathbb{R}$ has a **global minimizer** in $\mathbb{R}^n$.

### 4.3 Proof Sketch for Coercive Functions

1. Consider $\hat{x} \in \mathbb{R}^n$. Since $f$ is coercive, there exists $C > 0$ such that $\|x\| > C$ implies $f(x) > f(\hat{x})$.
2. Represent $\mathbb{R}^n = \mathbb{X}_1 \cup \mathbb{X}_2$ as a union of disjoint sets: $\mathbb{X}_1 = \{x : \|x\| \leq C\}$ and $\mathbb{X}_2 = \{x : \|x\| > C\}$.
3. Take $\hat{x} \in \mathbb{X}_1$ and any $x \in \mathbb{X}_2$. Then $f(\hat{x}) < f(x)$.
4. Since $\mathbb{X}_1$ is a compact set, by the Weierstrass theorem, there exists $x^* \in \mathbb{X}_1$ such that $f(x^*) = \min_{x \in \mathbb{X}_1} f(x)$.
5. Thus $f(\hat{x}) \geq f(x^*)$, which yields $f(x) \geq f(x^*)$ for all $x \in \mathbb{R}^n$.
6. This implies $x^*$ is a global minimizer.

---

<br>

## 5. One-Dimensional Optimization

In one-dimensional optimization, the goal is to find the minimum (or maximum) of a function $f(x)$ where $x \in \mathbb{R}$.

**Objectives:**

- Understand one-dimensional and multi-dimensional optimization
- Distinguish between global and local optima
- Cast maximization problems into minimization problems
- Learn **golden-section search** and **parabolic interpolation** methods

---

<br>

## 6. Golden-Section Search

The golden-section search is an efficient bracketing method for finding the minimum of a unimodal function on a closed interval. It systematically narrows the interval that contains the minimum.

### 6.1 The Golden Ratio

**Definition (Golden Ratio).** Divide a line into two segments $l_1$ and $l_2$ ($l_1 > l_2 > 0$) such that:

$$\frac{l_1 + l_2}{l_1} = \frac{l_1}{l_2} =: \phi \quad \text{(golden ratio)}$$

The ratio of the **whole line** to the **larger segment** equals the ratio of the **larger segment** to the **smaller segment**.

**Derivation:** From the defining relation:

$$l_1 l_2 + l_2^2 = l_1^2$$

Dividing by $l_2^2$:

$$\frac{l_1}{l_2} + 1 = \left(\frac{l_1}{l_2}\right)^2$$

$$\phi^2 - \phi - 1 = 0$$

$$\phi = \frac{1 \pm \sqrt{5}}{2}$$

Taking $\phi > 0$:

$$\phi = \frac{1 + \sqrt{5}}{2} \approx 1.61803$$

A useful reciprocal relation follows from dividing $(\ast)$ by $\phi$:

$$\frac{1}{\phi} = \phi - 1 = \frac{\sqrt{5} - 1}{2} \approx 0.61803$$

### 6.2 Algorithm: Interval Reduction

Given an interval $[x_l, x_u]$ that brackets a minimum, define:

$$d = (\phi - 1)(x_u - x_l)$$

Choose two interior points:

$$x_1^{(1)} = x_l + d = x_l(2 - \phi) + x_u(\phi - 1)$$

$$x_2^{(1)} = x_u - d = x_u(2 - \phi) + x_l(\phi - 1)$$

Note that $x_2 < x_1$.

**Decision rule:**

- If $f(x_2) > f(x_1)$: the minimum is in $[x_2, x_u]$, so set $x_l \leftarrow x_2^{(1)}$
- If $f(x_2) < f(x_1)$: the minimum is in $[x_l, x_1]$, so set $x_u \leftarrow x_1^{(1)}$

### 6.3 Reuse of Function Evaluations

The beauty of this method is that after each iteration, one of the two interior points from the previous step coincides with one of the new interior points. Specifically, when $f(x_2) > f(x_1)$ and we set $x_l \leftarrow x_2^{(1)}$:

$$x_2^{(2)} = x_u - d_{\text{new}} = x_1^{(1)}$$

This can be verified algebraically using the property $\phi^2 = \phi + 1$. Therefore, we **do not need to recalculate** the function evaluation at $x_2^{(2)}$ because it is the same as $f(x_1^{(1)})$.

> **[Efficiency]** Each iteration of the golden-section search requires only **one new function evaluation** (instead of two), making it very efficient. The interval shrinks by a constant factor of $\frac{1}{\phi} \approx 0.61803$ at each step.

### 6.4 Error Analysis

After one iteration, the estimate $x^*$ lies between $x_2$ and $x_1$.

**Case 1:** If $f(x_2) > f(x_1)$, then $x_2 \leq x^* \leq x_1$, and the interval is $[x_2, x_1, x_u]$. The maximum distance from the estimate to $x^*$:

$$\Delta x_a = x_1 - x_2 = (2\phi - 3)(x_u - x_l) \approx 0.2361 \cdot (x_u - x_l)$$

**Case 2:** If $f(x_2) < f(x_1)$, then $x_1 \leq x^* \leq x_u$, and the interval is $[x_l, x_2, x_1]$. The maximum distance:

$$\Delta x_b = x_u - x_1 = (2 - \phi)(x_u - x_l) \approx 0.3820 \cdot (x_u - x_l)$$

The **relative error** is defined as:

$$\varepsilon_a = (2 - \phi) \left| \frac{x_u - x_l}{x_{\text{opt}}} \right|$$

where $x_{\text{opt}}$ is the current best estimate of the optimum.

---

<br>

## 7. Parabolic Interpolation

Parabolic interpolation uses three points to fit a quadratic (parabola) and then locates the optimum of the parabola as the next estimate. This approach can converge faster than golden-section search for smooth functions.

### 7.1 Constructing the Parabola

Given three points $(x_1, y_1)$, $(x_2, y_2)$, $(x_3, y_3)$ where $y_i = f(x_i)$, construct the parabolic equation:

$$y = ax^2 + bx + c$$

This leads to the linear system:

$$\begin{bmatrix} y_1 \\ y_2 \\ y_3 \end{bmatrix} = \begin{bmatrix} x_1^2 & x_1 & 1 \\ x_2^2 & x_2 & 1 \\ x_3^2 & x_3 & 1 \end{bmatrix} \begin{bmatrix} a \\ b \\ c \end{bmatrix}$$

Using **Cramer's rule**, let $D = \begin{vmatrix} x_1^2 & x_1 & 1 \\ x_2^2 & x_2 & 1 \\ x_3^2 & x_3 & 1 \end{vmatrix}$, then:

$$a = \frac{1}{D}\begin{vmatrix} y_1 & x_1 & 1 \\ y_2 & x_2 & 1 \\ y_3 & x_3 & 1 \end{vmatrix}, \quad b = \frac{1}{D}\begin{vmatrix} x_1^2 & y_1 & 1 \\ x_2^2 & y_2 & 1 \\ x_3^2 & y_3 & 1 \end{vmatrix}, \quad c = \frac{1}{D}\begin{vmatrix} x_1^2 & x_1 & y_1 \\ x_2^2 & x_2 & y_2 \\ x_3^2 & x_3 & y_3 \end{vmatrix}$$

### 7.2 Optimal Point Formula

The optimal point of the parabola (vertex) is found by setting $\frac{dy}{dx} = 0$:

$$2ax_4 + b = 0 \implies x_4 = -\frac{b}{2a}$$

Substituting the Cramer's rule expressions:

$$x_4 = x_2 - \frac{1}{2} \cdot \frac{(x_2 - x_1)^2(y_2 - y_3) - (x_2 - x_3)^2(y_2 - y_1)}{(x_2 - x_1)(y_2 - y_3) - (x_2 - x_3)(y_2 - y_1)}$$

> **[Note]** For a **maximization** problem, the parabola should open downward ($a < 0$), and the vertex gives the approximate maximum. For a **minimization** problem, the parabola should open upward ($a > 0$), and the vertex gives the approximate minimum.

### 7.3 Updating Points for the Next Iteration

Once we compute $x_4$, we need to determine the new triplet $(x_1, x_2, x_3)$ for the next iteration. The decision depends on where $x_4$ falls relative to $x_2$:

**If $x_4$ is between $x_2$ and $x_3$ (i.e., $x_2 \leq x_4 \leq x_3$):**

| Condition | Action |
|-----------|--------|
| $f(x_4) < f(x_2)$ (Case 4) | $x_2 \Rightarrow x_1$, $x_4 \Rightarrow x_2$, $x_3$ stays |
| $f(x_4) \geq f(x_2)$ (Case 3) | $x_1$ stays, $x_2$ stays, $x_4 \Rightarrow x_3$ |

**If $x_4$ is between $x_1$ and $x_2$ (i.e., $x_1 \leq x_4 \leq x_2$):**

| Condition | Action |
|-----------|--------|
| $f(x_4) < f(x_2)$ (Case 2) | $x_2 \Rightarrow x_3$, $x_4 \Rightarrow x_2$, $x_1$ stays |
| $f(x_4) \geq f(x_2)$ (Case 1) | $x_4 \Rightarrow x_1$, $x_2$ stays, $x_3$ stays |

> **[Convergence]** Parabolic interpolation has **superlinear convergence** (approximately order 1.324), which is faster than the linear convergence of golden-section search. However, it can fail if the three points do not bracket a minimum well or if the fitted parabola degenerates.

---

<br>

## 8. Level Sets and Gradient

In optimization, level set functions play a crucial role in analyzing and visualizing objective functions.

### 8.1 Level Set Definition

**Definition (Level Set).** For a function $f : \mathbb{R}^n \to \mathbb{R}$ and a constant $c$:

- The set $\{x \in \mathbb{R}^n : f(x) = c\}$ is called the **level set** of $f$ at the level $c$
- The set $\{x \in \mathbb{R}^n : f(x) \leq c\}$ is called the **lower level set** of $f$ at $c$
- The set $\{x \in \mathbb{R}^n : f(x) \geq c\}$ is called the **upper level set** of $f$ at $c$

**Example:** For $f(x_1, x_2) = x_1^2 + x_2^2$, the level set at $c = 1$ is the unit circle $x_1^2 + x_2^2 = 1$. The lower level set $f(x,y) \leq 1$ is the unit disk, and the upper level set $f(x,y) \geq 1$ is the exterior of the unit disk.

### 8.2 Gradient Is Orthogonal to Level Sets

Consider a parametric curve $\gamma = \{x(t) : t \in (a, b)\} \subset S$ lying on a level set $S = \{x : f(x) = c\}$, where $x(t) : (a, b) \to S$ is a continuous function.

Since $f(x(t)) = c$ for $t \in (a, b)$, differentiating with respect to $t$:

$$\frac{df}{dt} = \nabla f^T \frac{dx}{dt} = 0$$

where:

$$\nabla f = \begin{bmatrix} \frac{\partial f}{\partial x_1} \\ \frac{\partial f}{\partial x_2} \end{bmatrix}, \quad \frac{dx}{dt} = \begin{bmatrix} \frac{dx_1}{dt} \\ \frac{dx_2}{dt} \end{bmatrix}$$

Since $\frac{dx}{dt}$ is the tangent vector to the curve $x(t)$, and $\nabla f^T \frac{dx}{dt} = 0$, this means:

$$\nabla f \perp \frac{dx}{dt}$$

**The gradient $\nabla f$ is orthogonal to the level set.**

**Example:** For $f(x_1, x_2) = x_1^2 + x_2^2$:

$$\nabla f = \begin{bmatrix} 2x_1 \\ 2x_2 \end{bmatrix} = 2\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}$$

The gradient points radially outward, perpendicular to the circular level sets, confirming the orthogonality.

> **[Geometry]** The gradient always points in the direction of steepest ascent. It is perpendicular to the level set at every point, pointing "uphill" toward higher values of $f$. This is the foundation for gradient-based optimization methods.

---

<br>

## 9. Convex Problems

### 9.1 Convex Combination

**Definition (Convex Combination).** Given points $x_1, x_2, \ldots, x_m$ in $\mathbb{R}^n$ and real numbers $\alpha_i \geq 0$ with $\sum_{i=1}^m \alpha_i = 1$, the point:

$$\sum_{i=1}^m \alpha_i x_i$$

is called a **convex combination** of these points.

### 9.2 Convex Set

**Definition (Convex Set).** A set $\mathbb{X} \subseteq \mathbb{R}^n$ is said to be **convex** if:

$$\alpha x + (1 - \alpha)y \in \mathbb{X}$$

holds for any $x, y \in \mathbb{X}$ and any $\alpha \in (0, 1)$.

Geometrically, a set is convex if the line segment connecting any two points in the set lies entirely within the set.

**Example:** Show that $\mathbb{X} = \{x \in \mathbb{R}^n : A_1 x = b_1, \; A_2 x \leq b_2\}$ is convex.

For any $x, y \in \mathbb{X}$ and any $\alpha \in (0, 1)$, let $z = \alpha x + (1 - \alpha)y$:

$$A_1 z = \alpha A_1 x + (1 - \alpha) A_1 y = \alpha b_1 + (1 - \alpha) b_1 = b_1$$

$$A_2 z = \alpha A_2 x + (1 - \alpha) A_2 y \leq \alpha b_2 + (1 - \alpha) b_2 = b_2$$

Therefore $z \in \mathbb{X}$, confirming that $\mathbb{X}$ is a convex set.

### 9.3 Hyperplane and Polyhedral Set

**Definition (Hyperplane).** A hyperplane $\mathbb{H} \subset \mathbb{R}^n$ is a set of the form:

$$\mathbb{H} = \{x \in \mathbb{R}^n : c^T x = b\}$$

where $c \in \mathbb{R}^n \setminus \{0\}$ and $b \in \mathbb{R}$.

A hyperplane is the intersection of two **half-spaces**:

$$\mathbb{H}_+ = \{x \in \mathbb{R}^n : c^T x \geq b\}, \quad \mathbb{H}_- = \{x \in \mathbb{R}^n : c^T x \leq b\}$$

**Proof that $\mathbb{H}$ is convex:** Take $x, y \in \mathbb{H}$ so that $c^T x = b$ and $c^T y = b$. For $z = \alpha x + (1 - \alpha) y$:

$$c^T z = \alpha c^T x + (1 - \alpha) c^T y = \alpha b + (1 - \alpha)b = b$$

Therefore $z \in \mathbb{H}$, which implies $\mathbb{H}$ is convex.

**Definition (Polyhedral Set).** A set defined by linear equations and/or inequalities is called a **polyhedral set** (or polyhedron).

### 9.4 Extreme Point

**Definition (Extreme Point).** Given a convex set $\mathbb{X}$, a point $x \in \mathbb{X}$ is called an **extreme point** of $\mathbb{X}$ if it cannot be represented as a convex combination of two **distinct** points in $\mathbb{X}$. That is, there do not exist distinct points $x', x'' \in \mathbb{X}$ and $\alpha \in (0, 1)$ such that $x = \alpha x' + (1 - \alpha) x''$.

**Examples:**

- For a circle (disk), every point on the boundary is an extreme point
- For a rectangle, only the four **vertices (corners)** are extreme points. Points on the edges (but not vertices) are NOT extreme points because they can be written as convex combinations of the two endpoints

For a convex feasible region with a linear objective function, the optimal solution (if it exists) occurs at an extreme point. Therefore, we only need to check extreme points for optimality.

### 9.5 Convex Function

**Definition (Convex Function).** Given a set $\mathbb{X} \subseteq \mathbb{R}^n$, a function $f : \mathbb{X} \to \mathbb{R}$ is **convex** if:

$$f(\alpha x + (1 - \alpha)y) \leq \alpha f(x) + (1 - \alpha) f(y) \quad \forall \; x, y \in \mathbb{X}, \; \alpha \in (0, 1)$$

We say $f$ is **strictly convex** if:

$$f(\alpha x + (1 - \alpha)y) < \alpha f(x) + (1 - \alpha) f(y) \quad \forall \; x, y \in \mathbb{X}, \; x \neq y, \; \alpha \in (0, 1)$$

Geometrically, a function is convex if the line segment connecting any two points on its graph lies above (or on) the graph.

**Concave functions:** A function $f$ is **concave** if $-f$ is convex.

**Examples:**

- $f(x) = \frac{1}{x}$ for $x > 0$ is **convex**
- $f(x) = \ln(x)$ for $x > 0$ is **concave**

**Proof that $f(x) = \frac{1}{x}$ is convex on $(0, \infty)$:** Need to show $\frac{1}{z} \leq \frac{\alpha}{x_1} + \frac{(1 - \alpha)}{x_2}$ where $z = \alpha x_1 + (1 - \alpha)x_2$.

This is equivalent to showing $x_1 x_2 \leq z(\alpha x_2 + (1 - \alpha)x_1)$, which reduces to $0 \leq (1 - \alpha)\alpha(x_1^2 - 2x_1 x_2 + x_2^2) = (1 - \alpha)\alpha(x_1 - x_2)^2$, which is always true.

**Properties of convex functions:**

- A convex function is **continuous** in the interior $\text{int}(\mathbb{X})$
- A convex function is **not necessarily differentiable** in $\text{int}(\mathbb{X})$
  - Example: $\mathbb{X} = \mathbb{R}$, $f(x) = |x|$ is convex but not differentiable at $x = 0$

> **[Optimization]** The importance of convexity in optimization cannot be overstated: for a convex function on a convex set, every local minimizer is also a global minimizer. This makes convex optimization problems fundamentally easier to solve.

---

<br>

## 10. Steepest Gradient (Descent) Method

### 10.1 Directional Derivative and Cauchy-Schwarz

Consider the unconstrained minimization problem:

$$\min_{x \in \mathbb{R}^n} f(x)$$

where $f$ is a continuously differentiable function.

Given $x^{(0)} \in \mathbb{R}^n$ and a direction $d \in \mathbb{R}^n$, the **directional derivative** of $f$ at $x^{(0)}$ in direction $d$ is:

$$\nabla f\big|_{x=x^{(0)}}^T d$$

If $\|d\| = 1$, this is the rate of increase of $f$ at $x^{(0)}$ in the direction $d$.

By the **Cauchy-Schwarz inequality**:

$$\nabla f^T d \leq \|\nabla f\| \cdot \|d\|$$

Choosing $d = \alpha \nabla f$ for some $\alpha > 0$:

$$\alpha \nabla f^T \nabla f = \alpha \|\nabla f\| \cdot \|\nabla f\|$$

This achieves equality in Cauchy-Schwarz, so $\alpha \nabla f$ is the direction of **maximum rate of increase** for $f$.

Similarly, for the direction of **maximum rate of decrease**, choose $d = -\alpha \nabla f$:

$$-\nabla f^T d \leq \|\nabla f\| \cdot \|d\|$$

$$\nabla f^T d \geq -\|\nabla f\| \cdot \|d\|$$

Choosing $d = -\alpha \nabla f$: $-\alpha \nabla f^T \nabla f = -\alpha \|\nabla f\| \cdot \|\nabla f\|$, which achieves the lower bound.

Therefore, $-\nabla f$ is the **best direction to take in a minimization method**.

### 10.2 The Steepest Descent Update Rule

The idea is to iteratively update:

$$x^{(k+1)} = x^{(k)} - \alpha_k \nabla f\big|_{x=x^{(k)}}$$

where $\alpha_k \geq 0$ is chosen such that:

$$f(x^{(k+1)}) < f(x^{(k)}), \quad k \geq 0$$

The step size $\alpha_k$ is typically found by **line search** (minimizing $\phi_k(\alpha) = f(x^{(k)} - \alpha \nabla f_k)$ with respect to $\alpha$).

```python
import numpy as np

def steepest_descent(f, grad_f, x0, tol=1e-8, max_iter=10000):
    """
    Steepest descent method with exact line search (golden-section).

    Parameters:
        f      : objective function
        grad_f : gradient function
        x0     : initial point (numpy array)
        tol    : tolerance for convergence
        max_iter: maximum number of iterations

    Returns:
        x      : approximate minimizer
        k      : number of iterations
    """
    x = np.array(x0, dtype=float)

    for k in range(max_iter):
        g = grad_f(x)
        if np.linalg.norm(g) < tol:
            break

        # Line search: minimize f(x - alpha * g)
        # Using golden-section search on alpha in [0, 1]
        phi = lambda a: f(x - a * g)
        alpha = golden_section_search(phi, 0, 1, tol=1e-10)

        x = x - alpha * g

    return x, k
```

### 10.3 Convergence and Orthogonality Property

**Theorem.** If $x^{(k)} \to x^*$ where $\{x^{(k)} : k \geq 0\}$ is the sequence generated by the steepest descent method, then:

$$\nabla f(x^*) = 0$$

**Proof sketch:** Let $\phi_k(\alpha) = f(x^{(k)} - \alpha \nabla f_k)$ and $\nabla f_k := \nabla f\big|_{x=x^{(k)}}$. We seek $\alpha$ that minimizes $\phi_k(\alpha)$:

$$\frac{d\phi_k}{d\alpha} = \frac{df(x^{(k+1)})}{dx^{(k+1)}} \cdot \frac{dx^{(k+1)}}{d\alpha} = -\nabla f_{k+1}^T \nabla f_k = 0$$

This yields the **orthogonality property**:

$$\nabla f_{k+1}^T \nabla f_k = 0, \quad k \geq 0$$

The gradient of $f$ at two **consecutive** points generated by the steepest descent method are **orthogonal to each other**.

Taking the limit as $k \to \infty$: $\nabla f\big|_{x^*}^T \nabla f\big|_{x^*} = \|\nabla f(x^*)\|^2 = 0$.

> **[Zigzag Behavior]** The orthogonality property explains the characteristic "zigzag" pattern of steepest descent: each step is perpendicular to the previous one. This can lead to slow convergence, especially near elongated (ill-conditioned) level sets, as the method repeatedly overshoots and corrects.

### 10.4 Example

$$f(x, y) = 3x^2 - 10x - 4xy + 2y^2 - 5y + 18$$

Starting from $(x^{(0)}, y^{(0)}) = (-2, 3)$ with tolerance $= 10^{-8}$.

The gradient is:

$$\nabla f = \begin{bmatrix} 6x - 10 - 4y \\ -4x + 4y - 5 \end{bmatrix}$$

The steepest descent method generates a sequence of iterates that zigzag toward the minimum, with consecutive gradients orthogonal to each other.

---

<br>

## Summary

| Topic | Key Concept |
|-------|-------------|
| **Optimization** | Find the best (min or max) among feasible solutions |
| **Feasible Region** | Defined by equality ($h_i(x)=0$) and inequality ($g_j(x) \leq 0$) constraints |
| **Min vs. Max** | $\min f(x) \iff \max(-f(x))$; any max problem can be recast as min |
| **Interior / Boundary** | Interior: $\varepsilon$-ball fits inside set; Boundary: in set but not interior |
| **Global Minimum** | $f(x^*) \leq f(x)$ for all $x \in \mathbb{X}$ |
| **Local Minimum** | $f(x^*) \leq f(x)$ for all $x$ in a neighborhood |
| **Weierstrass Theorem** | Continuous $f$ on compact $\mathbb{X}$ has global optimizer |
| **Coercive Function** | $f(x) \to \infty$ as $\|x\| \to \infty$; guarantees global min on $\mathbb{R}^n$ |
| **Golden-Section Search** | Bracket-based 1D method; reduces interval by factor $\frac{1}{\phi} \approx 0.618$ per step; 1 new function eval per iteration |
| **Parabolic Interpolation** | Fits quadratic through 3 points; vertex $x_4 = -b/(2a)$; superlinear convergence |
| **Level Set** | $\{x : f(x) = c\}$; gradient is orthogonal to level sets |
| **Convex Set** | Line segment between any two points in the set stays inside |
| **Convex Function** | $f(\alpha x + (1-\alpha)y) \leq \alpha f(x) + (1-\alpha)f(y)$; local min = global min |
| **Hyperplane** | $\{x : c^T x = b\}$; always convex |
| **Extreme Point** | Cannot be expressed as convex combination of two distinct points in the set |
| **Steepest Descent** | Update: $x^{(k+1)} = x^{(k)} - \alpha_k \nabla f_k$; $-\nabla f$ is direction of steepest decrease |
| **Orthogonality Property** | Consecutive steepest descent gradients are perpendicular: $\nabla f_{k+1}^T \nabla f_k = 0$ |
