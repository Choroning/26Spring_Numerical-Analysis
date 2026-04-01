# Chapter 8 Lecture — Solving Linear Algebraic Equations

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 8

> **Prerequisites**: [Linear Algebra] Matrix operations. [Calculus] Linear systems (Ch 1-7).
>
> **Learning Objectives**:
> 1. Formulate engineering problems as linear algebraic equations
> 2. Identify solution existence and uniqueness conditions
> 3. Apply graphical interpretation of linear systems

---

<br>

## Table of Contents

- [1. Motivation](#1-motivation)
- [2. Small Systems (2x2): Graphical Method](#2-small-systems-2x2-graphical-method)
  - [2.1 Two Equations as Two Lines](#21-two-equations-as-two-lines)
  - [2.2 Three Possibilities](#22-three-possibilities)
- [3. Determinant and Solvability](#3-determinant-and-solvability)
  - [3.1 Determinant of a 2x2 System](#31-determinant-of-a-2x2-system)
  - [3.2 Singular and Ill-Conditioned Systems](#32-singular-and-ill-conditioned-systems)
- [4. Cramer's Rule](#4-cramers-rule)
  - [4.1 Derivation for 2x2](#41-derivation-for-2x2)
  - [4.2 General Form](#42-general-form)
  - [4.3 Limitation](#43-limitation)
- [5. Naive Gauss Elimination](#5-naive-gauss-elimination)
  - [5.1 Forward Elimination](#51-forward-elimination)
  - [5.2 Back Substitution](#52-back-substitution)
- [6. Computational Cost](#6-computational-cost)
- [7. Example: 3x3 System (Complete Step-by-Step)](#7-example-3x3-system-complete-step-by-step)
  - [7.1 Problem Setup](#71-problem-setup)
  - [7.2 Forward Elimination — Step 1](#72-forward-elimination--step-1)
  - [7.3 Forward Elimination — Step 2](#73-forward-elimination--step-2)
  - [7.4 Back Substitution](#74-back-substitution)
- [8. Pitfalls of Naive Gauss Elimination](#8-pitfalls-of-naive-gauss-elimination)
- [Summary](#summary)

---

<br>

## 1. Motivation

Many engineering and scientific problems ultimately reduce to solving a system of simultaneous linear algebraic equations of the form:

$$[A]\{x\} = \{b\}$$

where $[A]$ is an $n \times n$ coefficient matrix, $\{x\}$ is the vector of $n$ unknowns, and $\{b\}$ is the right-hand side vector of known constants.

These systems arise naturally across a wide range of disciplines:

- **Mass-spring systems**: Hooke's law ($F = kx$) applied to each spring, combined with force equilibrium at each mass, produces a system of linear equations. The coefficient matrix is the **stiffness matrix** $[K]$ encoding spring constants and connectivity.
- **Structural truss analysis**: Equilibrium of forces at each joint of a truss structure yields linear equations relating member forces to applied loads.
- **Circuit analysis**: Kirchhoff's current and voltage laws applied to electrical networks produce linear systems relating currents and voltages through resistances.

> **[Physics]** In a mass-spring system, Hooke's law ($F = kx$) for each spring combined with force equilibrium at each mass yields a system of linear equations. The stiffness matrix $[K]$ encodes the spring constants and connectivity.

The size of these systems can range from a handful of equations (simple textbook problems) to millions of equations (finite element models of aircraft structures or computational fluid dynamics simulations). Efficient and reliable methods for solving such systems are therefore essential tools in the engineer's toolkit.

---

<br>

## 2. Small Systems (2x2): Graphical Method

### 2.1 Two Equations as Two Lines

For a system of two equations in two unknowns:

$$a_{11}x_1 + a_{12}x_2 = b_1$$
$$a_{21}x_1 + a_{22}x_2 = b_2$$

each equation represents a straight line in the $x_1$-$x_2$ plane. Solving the system means finding the point(s) where both lines are simultaneously satisfied.

### 2.2 Three Possibilities

There are exactly three possible geometric configurations for two lines in 2D:

1. **Unique solution** — The two lines intersect at exactly one point. The system has a single solution $(x_1, x_2)$. This is the typical case when the equations are independent and the determinant is nonzero.

2. **No solution** — The two lines are parallel (same slope, different intercepts). The system is **inconsistent**: no point satisfies both equations simultaneously. This occurs when the rows of $[A]$ are proportional but $\{b\}$ is not similarly proportional.

3. **Infinite solutions** — The two lines are identical (one equation is a scalar multiple of the other). Every point on the line is a solution. The system is **dependent**: the equations carry redundant information.

These three cases correspond precisely to the rank conditions of the augmented matrix $[A|b]$. In the unique-solution case, $\text{rank}([A]) = \text{rank}([A|b]) = n$. When the system is inconsistent, $\text{rank}([A]) < \text{rank}([A|b])$. When the system is dependent, $\text{rank}([A]) = \text{rank}([A|b]) < n$.

---

<br>

## 3. Determinant and Solvability

### 3.1 Determinant of a 2x2 System

For the $2 \times 2$ coefficient matrix:

$$[A] = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix}$$

the **determinant** is:

$$D = a_{11}a_{22} - a_{12}a_{21}$$

The determinant provides a single scalar that characterizes the solvability of the system.

### 3.2 Singular and Ill-Conditioned Systems

- **$D = 0$**: The system is **singular** — no unique solution exists. The rows of $[A]$ are linearly dependent (one row is a scalar multiple of the other), and the system is either inconsistent or dependent.

- **$D$ very small relative to the coefficients**: The system is **ill-conditioned**. A unique solution technically exists, but it is extremely sensitive to small perturbations in the coefficients or the right-hand side. Tiny round-off errors in the input data or during computation can produce wildly inaccurate results.

> **[Linear Algebra]** The determinant measures the "volume" of the parallelogram/parallelepiped formed by the row (or column) vectors. Zero determinant means the vectors are linearly dependent — the system of equations is redundant.

Geometrically, ill-conditioning for a $2 \times 2$ system means the two lines intersect at a very shallow angle. In such cases, small changes in the coefficients (shifting a line slightly) can move the intersection point dramatically. The ratio of the largest to smallest singular value of $[A]$ — the **condition number** — quantifies this sensitivity and will be discussed further in Chapter 11.

---

<br>

## 4. Cramer's Rule

### 4.1 Derivation for 2x2

Consider the $2 \times 2$ system:

$$a_{11}x_1 + a_{12}x_2 = b_1$$
$$a_{21}x_1 + a_{22}x_2 = b_2$$

To solve for $x_1$, we eliminate $x_2$ algebraically:

1. Multiply the first equation by $a_{22}$:
   $$a_{11}a_{22}x_1 + a_{12}a_{22}x_2 = a_{22}b_1$$

2. Multiply the second equation by $a_{12}$:
   $$a_{21}a_{12}x_1 + a_{12}a_{22}x_2 = a_{12}b_2$$

3. Subtract the second from the first to eliminate $x_2$:
   $$(a_{11}a_{22} - a_{12}a_{21})x_1 = a_{22}b_1 - a_{12}b_2$$

4. Solve for $x_1$ (provided $D \neq 0$):
   $$x_1 = \frac{a_{22}b_1 - a_{12}b_2}{a_{11}a_{22} - a_{12}a_{21}}$$

Similarly, eliminating $x_1$ yields:

$$x_2 = \frac{a_{11}b_2 - a_{21}b_1}{a_{11}a_{22} - a_{12}a_{21}}$$

### 4.2 General Form

In compact notation using determinants:

$$x_1 = \frac{a_{22}b_1 - a_{12}b_2}{a_{11}a_{22} - a_{12}a_{21}}, \quad x_2 = \frac{a_{11}b_2 - a_{21}b_1}{a_{11}a_{22} - a_{12}a_{21}}$$

More generally, for an $n \times n$ system, Cramer's rule states:

$$x_i = \frac{\det(A_i)}{\det(A)}$$

where $A_i$ is the matrix formed by replacing the $i$-th column of $A$ with the right-hand side vector $b$. The numerator determinant $\det(A_i)$ is computed by substituting $b$ into the appropriate column position, and the denominator is always $\det(A)$.

For example, in a $3 \times 3$ system:

$$x_1 = \frac{\det\begin{bmatrix} b_1 & a_{12} & a_{13} \\ b_2 & a_{22} & a_{23} \\ b_3 & a_{32} & a_{33} \end{bmatrix}}{\det(A)}, \quad x_2 = \frac{\det\begin{bmatrix} a_{11} & b_1 & a_{13} \\ a_{21} & b_2 & a_{23} \\ a_{31} & b_3 & a_{33} \end{bmatrix}}{\det(A)}, \quad x_3 = \frac{\det\begin{bmatrix} a_{11} & a_{12} & b_1 \\ a_{21} & a_{22} & b_2 \\ a_{31} & a_{32} & b_3 \end{bmatrix}}{\det(A)}$$

### 4.3 Limitation

Cramer's rule requires evaluating $n + 1$ determinants, each of which involves $O(n!)$ operations when computed by cofactor expansion. The total cost is therefore $O((n+1) \cdot n!) \approx O(n!)$, which grows astronomically fast:

| $n$ | $n!$ operations |
|:----|:----------------|
| 3 | 6 |
| 5 | 120 |
| 10 | 3,628,800 |
| 20 | $\approx 2.4 \times 10^{18}$ |

For $n > 3$, Cramer's rule is computationally impractical. A computer performing $10^9$ operations per second would need approximately **77 years** to solve a $20 \times 20$ system using Cramer's rule. This motivates the development of more efficient algorithms.

---

<br>

## 5. Naive Gauss Elimination

Gauss elimination is the workhorse algorithm for solving linear systems. The method consists of two phases: **forward elimination** (reduce the system to upper triangular form) and **back substitution** (solve the triangular system from bottom to top).

### 5.1 Forward Elimination

**Goal:** Transform the augmented matrix $[A|b]$ into upper triangular form, so that all entries below the main diagonal are zero.

The procedure works column by column from left to right:

**Step 1 — Pivot on row 1:** Use $a_{11}$ as the **pivot element**. For each row $j = 2, 3, \ldots, n$ below the pivot:

$$\text{factor} = \frac{a_{j1}}{a_{11}}$$

$$\text{Row}_j \leftarrow \text{Row}_j - \text{factor} \times \text{Row}_1$$

This zeros out all entries in column 1 below the diagonal.

**Step 2 — Pivot on row 2:** Use $a_{22}^{(1)}$ (the updated entry after Step 1) as the pivot element. Eliminate entries in column 2 from rows $3, 4, \ldots, n$.

**Continue** through column $n-1$: At each step $i$, use $a_{ii}^{(i-1)}$ as the pivot element and eliminate all entries below it.

**Pseudocode:**

```
for i = 1 to n-1:          // pivot column
    for j = i+1 to n:       // rows below pivot
        factor = a[j][i] / a[i][i]
        for k = i to n:
            a[j][k] = a[j][k] - factor * a[i][k]
        b[j] = b[j] - factor * b[i]
```

After forward elimination, the system has the form:

$$\begin{bmatrix} a_{11} & a_{12} & a_{13} & \cdots & a_{1n} \\ 0 & a_{22}^{(1)} & a_{23}^{(1)} & \cdots & a_{2n}^{(1)} \\ 0 & 0 & a_{33}^{(2)} & \cdots & a_{3n}^{(2)} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & a_{nn}^{(n-1)} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ \vdots \\ x_n \end{bmatrix} = \begin{bmatrix} b_1 \\ b_2^{(1)} \\ b_3^{(2)} \\ \vdots \\ b_n^{(n-1)} \end{bmatrix}$$

where the superscripts indicate modified values after each elimination step.

### 5.2 Back Substitution

Once the system is in upper triangular form, we solve from the bottom up.

The last equation involves only $x_n$:

$$x_n = \frac{b_n^{(n-1)}}{a_{nn}^{(n-1)}}$$

For each remaining unknown, working upward from $i = n-1$ to $i = 1$:

$$x_i = \frac{b_i^{(i-1)} - \sum_{j=i+1}^{n} a_{ij}^{(i-1)} x_j}{a_{ii}^{(i-1)}}, \quad i = n-1, n-2, \ldots, 1$$

The idea is straightforward: since all the $x_j$ values for $j > i$ have already been computed, we substitute them in and solve for $x_i$.

**Pseudocode:**

```
x[n] = b[n] / a[n][n]
for i = n-1 down to 1:
    sum = 0
    for j = i+1 to n:
        sum = sum + a[i][j] * x[j]
    x[i] = (b[i] - sum) / a[i][i]
```

---

<br>

## 6. Computational Cost

The efficiency of Gauss elimination can be quantified by counting the number of floating-point multiplications and divisions (which dominate the cost):

**Forward elimination:**

The number of multiplications/divisions in the forward elimination phase is:

$$\sum_{k=1}^{n-1} \sum_{i=k+1}^{n} (n - k + 1) = \frac{2n^3 + 3n^2 - 5n}{6} \approx \frac{2}{3}n^3 + O(n^2)$$

**Back substitution:**

$$\sum_{i=1}^{n-1} (n - i) + n = \frac{n^2 + n}{2} \approx \frac{1}{2}n^2 + O(n)$$

**Total cost:**

$$\text{Total} = \frac{2}{3}n^3 + O(n^2) + \frac{1}{2}n^2 + O(n)$$

For large $n$, the total cost is **dominated by the $(2/3)n^3$ term** from forward elimination. Back substitution contributes only an $O(n^2)$ term, which is negligible for large systems.

| Phase | Dominant Term | Growth |
|:------|:-------------|:-------|
| Forward Elimination | $(2/3)n^3$ | Cubic |
| Back Substitution | $(1/2)n^2$ | Quadratic |
| **Total** | **$(2/3)n^3$** | **Cubic** |

This is a dramatic improvement over Cramer's rule ($O(n!)$). For $n = 20$:
- Cramer's rule: $\approx 2.4 \times 10^{18}$ operations
- Gauss elimination: $\approx (2/3)(20)^3 = 5333$ operations

---

<br>

## 7. Example: 3x3 System (Complete Step-by-Step)

### 7.1 Problem Setup

Solve the following system:

$$3x_1 - 0.1x_2 - 0.2x_3 = 7.85$$
$$0.1x_1 + 7x_2 - 0.3x_3 = -19.3$$
$$0.3x_1 - 0.2x_2 + 10x_3 = 71.4$$

The augmented matrix is:

$$[A|b] = \left[\begin{array}{ccc|c} 3 & -0.1 & -0.2 & 7.85 \\ 0.1 & 7 & -0.3 & -19.3 \\ 0.3 & -0.2 & 10 & 71.4 \end{array}\right]$$

### 7.2 Forward Elimination — Step 1

**Pivot element:** $a_{11} = 3$

**Eliminate $x_1$ from Row 2:**

$$\text{factor}_{21} = \frac{a_{21}}{a_{11}} = \frac{0.1}{3} = 0.0333\overline{3}$$

$$\text{Row}_2 \leftarrow \text{Row}_2 - 0.0333\overline{3} \times \text{Row}_1$$

- $a_{22}^{(1)} = 7 - 0.0333\overline{3}(-0.1) = 7.00333\overline{3}$
- $a_{23}^{(1)} = -0.3 - 0.0333\overline{3}(-0.2) = -0.29333\overline{3}$
- $b_2^{(1)} = -19.3 - 0.0333\overline{3}(7.85) = -19.5617$

**Eliminate $x_1$ from Row 3:**

$$\text{factor}_{31} = \frac{a_{31}}{a_{11}} = \frac{0.3}{3} = 0.1$$

$$\text{Row}_3 \leftarrow \text{Row}_3 - 0.1 \times \text{Row}_1$$

- $a_{32}^{(1)} = -0.2 - 0.1(-0.1) = -0.19$
- $a_{33}^{(1)} = 10 - 0.1(-0.2) = 10.02$
- $b_3^{(1)} = 71.4 - 0.1(7.85) = 70.615$

After Step 1:

$$\left[\begin{array}{ccc|c} 3 & -0.1 & -0.2 & 7.85 \\ 0 & 7.00333 & -0.29333 & -19.5617 \\ 0 & -0.19 & 10.02 & 70.615 \end{array}\right]$$

### 7.3 Forward Elimination — Step 2

**Pivot element:** $a_{22}^{(1)} = 7.00333\overline{3}$

**Eliminate $x_2$ from Row 3:**

$$\text{factor}_{32} = \frac{a_{32}^{(1)}}{a_{22}^{(1)}} = \frac{-0.19}{7.00333} = -0.027130$$

$$\text{Row}_3 \leftarrow \text{Row}_3 - (-0.027130) \times \text{Row}_2$$

- $a_{33}^{(2)} = 10.02 - (-0.027130)(-0.29333) = 10.02 - 0.007960 = 10.0120$
- $b_3^{(2)} = 70.615 - (-0.027130)(-19.5617) = 70.615 - 0.530786 = 70.0842$

After Step 2 (upper triangular form):

$$\left[\begin{array}{ccc|c} 3 & -0.1 & -0.2 & 7.85 \\ 0 & 7.00333 & -0.29333 & -19.5617 \\ 0 & 0 & 10.0120 & 70.0842 \end{array}\right]$$

### 7.4 Back Substitution

**Solve for $x_3$:**

$$x_3 = \frac{b_3^{(2)}}{a_{33}^{(2)}} = \frac{70.0842}{10.0120} = 7.0003$$

**Solve for $x_2$:**

$$x_2 = \frac{b_2^{(1)} - a_{23}^{(1)}x_3}{a_{22}^{(1)}} = \frac{-19.5617 - (-0.29333)(7.0003)}{7.00333} = \frac{-19.5617 + 2.05343}{7.00333} = \frac{-17.5083}{7.00333} = -2.5000$$

**Solve for $x_1$:**

$$x_1 = \frac{b_1 - a_{12}x_2 - a_{13}x_3}{a_{11}} = \frac{7.85 - (-0.1)(-2.5) - (-0.2)(7.0003)}{3} = \frac{7.85 - 0.25 + 1.40006}{3} = \frac{9.00006}{3} = 3.0000$$

**Final solution:**

$$x_1 = 3, \quad x_2 = -2.5, \quad x_3 = 7.0003$$

> **Note:** The slight deviation of $x_3$ from $7.0$ is due to round-off error accumulated during the elimination process. The exact solution (to more decimal places) is $x_3 \approx 7.00003$, which is very close to the exact value of $7$.

---

<br>

## 8. Pitfalls of Naive Gauss Elimination

The term "naive" refers to the fact that this basic algorithm has no safeguards against numerical difficulties. Three major pitfalls exist:

**1. Division by zero if the pivot element is zero:**

If $a_{ii}^{(i-1)} = 0$ at any step, the algorithm fails because it attempts to divide by zero when computing the factor. This can happen even when the system has a perfectly valid unique solution — it just means the equations are ordered in an unfavorable way.

**Example:** Consider the system where $a_{11} = 0$:

$$0 \cdot x_1 + 2x_2 = 4$$
$$3x_1 + x_2 = 5$$

The solution is $x_1 = 1, x_2 = 2$, but naive Gauss elimination fails at the first step because the pivot $a_{11} = 0$.

**2. Near-zero pivots amplify round-off error:**

When a pivot element is very small but nonzero, the computed factor $a_{ji}/a_{ii}$ becomes very large. Multiplying row $i$ by this large factor and subtracting magnifies any round-off errors present in row $i$, potentially corrupting the result severely.

**3. Large factor values magnify errors:**

Related to the previous point, when the factor is much larger than 1, the subtraction $\text{Row}_j - \text{factor} \times \text{Row}_i$ involves subtracting a large multiple of one row from another. This can lead to catastrophic cancellation and significant loss of significant digits.

These three pitfalls all stem from the same root cause: **the pivot element is too small relative to other entries in the column**. The remedy is **partial pivoting** — at each step, swap rows so that the largest available entry in the current column becomes the pivot. This strategy, covered in Chapter 9, effectively eliminates all three problems and is used in all production-quality implementations of Gauss elimination.

---

<br>

## Summary

| Topic | Key Point |
|:------|:----------|
| System form | $[A]\{x\} = \{b\}$; arises in springs, trusses, circuits |
| Graphical (2x2) | Unique, no solution (parallel), or infinite solutions (same line) |
| Determinant | $D = a_{11}a_{22} - a_{12}a_{21}$; $D = 0$ means singular |
| Cramer's Rule | $x_i = \det(A_i)/\det(A)$; $O(n!)$ — impractical for $n > 3$ |
| Forward Elimination | Zero out below-diagonal entries column by column |
| Back Substitution | Solve from $x_n$ upward using upper triangular form |
| Computational Cost | $(2/3)n^3$ total (dominated by forward elimination) |
| Pitfalls | Zero/small pivots → division by zero or error amplification |
| Remedy | Partial pivoting (Chapter 9) |

---
