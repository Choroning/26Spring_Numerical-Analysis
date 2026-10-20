# Chapter 8 — Solving Linear Algebraic Equations

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 8, 9

> **Learning Objectives:**
> 1. Formulate engineering problems as linear algebraic equations
> 2. Identify solution existence and uniqueness conditions
> 3. Apply graphical interpretation of linear systems
> 4. Solve linear systems with naive Gauss elimination and analyze its computational cost

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
  - [7.2 Forward Elimination: Step 1](#72-forward-elimination-step-1)
  - [7.3 Forward Elimination: Step 2](#73-forward-elimination-step-2)
  - [7.4 Back Substitution](#74-back-substitution)
- [8. Pitfalls of Naive Gauss Elimination](#8-pitfalls-of-naive-gauss-elimination)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Motivation

Many engineering and scientific problems ultimately reduce to solving a system of simultaneous linear algebraic equations of the form:

```math
[A]\{x\} = \{b\}
```

where $`[A]`$ is an $`n \times n`$ coefficient matrix, $`\{x\}`$ is the vector of $`n`$ unknowns, and $`\{b\}`$ is the right-hand side vector of known constants.

These systems arise naturally across a wide range of disciplines:

- **Mass-spring systems**: Hooke's law ($`F = kx`$) applied to each spring, combined with force equilibrium at each mass, produces a system of linear equations. The coefficient matrix is the **stiffness matrix** $`[K]`$ encoding spring constants and connectivity.
- **Structural truss analysis**: Equilibrium of forces at each joint of a truss structure yields linear equations relating member forces to applied loads.
- **Circuit analysis**: Kirchhoff's current and voltage laws applied to electrical networks produce linear systems relating currents and voltages through resistances.

The size of these systems can range from a handful of equations (simple textbook problems) to millions of equations (finite element models of aircraft structures or computational fluid dynamics simulations). Efficient and reliable methods for solving such systems are therefore essential tools in the engineer's toolkit.

---

<br>

## 2. Small Systems (2x2): Graphical Method

### 2.1 Two Equations as Two Lines

Consider a system of two equations in two unknowns:

```math
a_{11}x_1 + a_{12}x_2 = b_1
```
```math
a_{21}x_1 + a_{22}x_2 = b_2
```

Each equation represents a straight line in the $`x_1`$-$`x_2`$ plane. Solving the system means finding the point(s) that satisfy both line equations simultaneously.

### 2.2 Three Possibilities

There are exactly three possible geometric configurations for two lines in 2D:

1. **Unique solution**: The two lines intersect at exactly one point. The system has a single solution $`(x_1, x_2)`$. This is the typical case when the equations are independent and the determinant is nonzero.

2. **No solution**: The two lines are parallel (same slope, different intercepts). The system is **inconsistent**: no point satisfies both equations simultaneously. This occurs when the rows of $`[A]`$ are proportional but $`\{b\}`$ is not similarly proportional.

3. **Infinite solutions**: The two lines are identical (one equation is a scalar multiple of the other). Every point on the line is a solution. The system is **dependent**: the equations carry redundant information.

These three cases correspond precisely to the rank conditions of the augmented matrix $`[A|b]`$. In the unique-solution case, $`\text{rank}([A]) = \text{rank}([A|b]) = n`$. When the system is inconsistent, $`\text{rank}([A]) < \text{rank}([A|b])`$. When the system is dependent, $`\text{rank}([A]) = \text{rank}([A|b]) < n`$.

---

<br>

## 3. Determinant and Solvability

### 3.1 Determinant of a 2x2 System

For the $`2 \times 2`$ coefficient matrix:

```math
[A] = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix}
```

the **determinant** is:

```math
D = a_{11}a_{22} - a_{12}a_{21}
```

The determinant provides a single scalar that characterizes the solvability of the system.

### 3.2 Singular and Ill-Conditioned Systems

- **$`D = 0`$**: The system is **singular**, and no unique solution exists. The rows of $`[A]`$ are linearly dependent (one row is a scalar multiple of the other), and the system is either inconsistent or dependent.

- **$`D`$ very small relative to the coefficients**: The system is **ill-conditioned**. A unique solution technically exists, but it is extremely sensitive to small perturbations in the coefficients or the right-hand side. Tiny round-off errors in the input data or during computation can produce wildly inaccurate results. Because the magnitude of the determinant changes when an equation is merely multiplied by a constant, the comparison is meaningful only after each equation is normalized.

> **[Linear Algebra]** The absolute value of the determinant is the area of the parallelogram (2D) or the volume of the parallelepiped (3D) formed by the row (or column) vectors. A zero determinant means the vectors are linearly dependent and the area or volume collapses to zero; the redundant rows of the coefficient matrix in the singular systems above are exactly this case.

Geometrically, ill-conditioning for a $`2 \times 2`$ system means the two lines intersect at a very shallow angle. In such cases, small changes in the coefficients (shifting a line slightly) can move the intersection point dramatically. This sensitivity is quantified by the **condition number**, which in the 2-norm is the ratio of the largest to the smallest singular value of $`[A]`$. The condition number is discussed further in Chapter 11.

---

<br>

## 4. Cramer's Rule

### 4.1 Derivation for 2x2

Consider the $`2 \times 2`$ system:

```math
a_{11}x_1 + a_{12}x_2 = b_1
```
```math
a_{21}x_1 + a_{22}x_2 = b_2
```

To solve for $`x_1`$, we eliminate $`x_2`$ algebraically:

1. Multiply the first equation by $`a_{22}`$:
   ```math
   a_{11}a_{22}x_1 + a_{12}a_{22}x_2 = a_{22}b_1
   ```

2. Multiply the second equation by $`a_{12}`$:
   ```math
   a_{21}a_{12}x_1 + a_{12}a_{22}x_2 = a_{12}b_2
   ```

3. Subtract the second from the first to eliminate $`x_2`$:
   ```math
   (a_{11}a_{22} - a_{12}a_{21})x_1 = a_{22}b_1 - a_{12}b_2
   ```

4. Solve for $`x_1`$ (provided $`D \neq 0`$):
   ```math
   x_1 = \frac{a_{22}b_1 - a_{12}b_2}{a_{11}a_{22} - a_{12}a_{21}}
   ```

Similarly, eliminating $`x_1`$ yields:

```math
x_2 = \frac{a_{11}b_2 - a_{21}b_1}{a_{11}a_{22} - a_{12}a_{21}}
```

### 4.2 General Form

Written together, the two solutions are as follows, where each numerator and denominator is a $`2 \times 2`$ determinant:

```math
x_1 = \frac{a_{22}b_1 - a_{12}b_2}{a_{11}a_{22} - a_{12}a_{21}}, \quad x_2 = \frac{a_{11}b_2 - a_{21}b_1}{a_{11}a_{22} - a_{12}a_{21}}
```

More generally, for an $`n \times n`$ system, Cramer's rule states:

```math
x_i = \frac{\det(A_i)}{\det(A)}
```

where $`A_i`$ is the matrix formed by replacing the $`i`$-th column of $`A`$ with the right-hand side vector $`b`$. The numerator determinant $`\det(A_i)`$ is computed by substituting $`b`$ into the appropriate column position, and the denominator is always $`\det(A)`$.

For example, in a $`3 \times 3`$ system:

```math
x_1 = \frac{\det\begin{bmatrix} b_1 & a_{12} & a_{13} \\ b_2 & a_{22} & a_{23} \\ b_3 & a_{32} & a_{33} \end{bmatrix}}{\det(A)}, \quad x_2 = \frac{\det\begin{bmatrix} a_{11} & b_1 & a_{13} \\ a_{21} & b_2 & a_{23} \\ a_{31} & b_3 & a_{33} \end{bmatrix}}{\det(A)}, \quad x_3 = \frac{\det\begin{bmatrix} a_{11} & a_{12} & b_1 \\ a_{21} & a_{22} & b_2 \\ a_{31} & a_{32} & b_3 \end{bmatrix}}{\det(A)}
```

### 4.3 Limitation

Cramer's rule requires evaluating $`n + 1`$ determinants, and computing each one by cofactor expansion takes about $`n!`$ operations. The total cost is therefore about $`(n+1) \cdot n! = (n+1)!`$ operations, which grows factorially, that is, astronomically fast. The table below lists the operation count for one determinant:

| $`n`$ | $`n!`$ operations |
|:----|:----------------|
| 3 | 6 |
| 5 | 120 |
| 10 | 3,628,800 |
| 20 | $`\approx 2.4 \times 10^{18}`$ |

*Table 1. Operation count for one determinant by cofactor expansion*

For $`n > 3`$, Cramer's rule is computationally impractical. Even a computer performing $`10^9`$ operations per second would need approximately **77 years** to evaluate a single $`20 \times 20`$ determinant ($`20! \approx 2.4 \times 10^{18}`$ operations), and Cramer's rule needs 21 such determinants. This motivates the development of more efficient algorithms.

---

<br>

## 5. Naive Gauss Elimination

Gauss elimination is the workhorse algorithm for solving linear systems. The method consists of two phases. **Forward elimination** reduces the system to upper triangular form, and **back substitution** solves the triangular system from bottom to top.

### 5.1 Forward Elimination

**Goal:** Transform the augmented matrix $`[A|b]`$ into upper triangular form, so that all entries below the main diagonal are zero.

The procedure works column by column from left to right:

**Step 1 (pivot on row 1):** Use $`a_{11}`$ as the **pivot element**. For each row $`j = 2, 3, \ldots, n`$ below the pivot:

```math
\text{factor} = \frac{a_{j1}}{a_{11}}
```

```math
\text{Row}_j \leftarrow \text{Row}_j - \text{factor} \times \text{Row}_1
```

This zeros out all entries in column 1 below the diagonal.

**Step 2 (pivot on row 2):** Use $`a_{22}^{(1)}`$ (the updated entry after Step 1) as the pivot element. Eliminate entries in column 2 from rows $`3, 4, \ldots, n`$.

**Continue through column $`n-1`$:** At each step $`i`$, use $`a_{ii}^{(i-1)}`$ as the pivot element and eliminate all entries below it. For a $`3 \times 3`$ augmented matrix, the process proceeds as in Figure 1.

```text
[ * * * | * ]       [ * * * | * ]       [ * * * | * ]
[ * * * | * ]  -->  [ 0 * * | * ]  -->  [ 0 * * | * ]
[ * * * | * ]       [ 0 * * | * ]       [ 0 0 * | * ]
```

*Figure 1. A 3x3 augmented matrix reaching upper triangular form through the two steps of forward elimination*

**Pseudocode:**

```text
for i = 1 to n-1:          // pivot column
    for j = i+1 to n:       // rows below pivot
        factor = a[j][i] / a[i][i]
        for k = i to n:
            a[j][k] = a[j][k] - factor * a[i][k]
        b[j] = b[j] - factor * b[i]
```

After forward elimination, the system has the form:

```math
\begin{bmatrix} a_{11} & a_{12} & a_{13} & \cdots & a_{1n} \\ 0 & a_{22}^{(1)} & a_{23}^{(1)} & \cdots & a_{2n}^{(1)} \\ 0 & 0 & a_{33}^{(2)} & \cdots & a_{3n}^{(2)} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & a_{nn}^{(n-1)} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ \vdots \\ x_n \end{bmatrix} = \begin{bmatrix} b_1 \\ b_2^{(1)} \\ b_3^{(2)} \\ \vdots \\ b_n^{(n-1)} \end{bmatrix}
```

where the superscripts indicate modified values after each elimination step.

### 5.2 Back Substitution

Once the system is in upper triangular form, we solve from the bottom up.

The last equation involves only $`x_n`$:

```math
x_n = \frac{b_n^{(n-1)}}{a_{nn}^{(n-1)}}
```

For each remaining unknown, working upward from $`i = n-1`$ to $`i = 1`$:

```math
x_i = \frac{b_i^{(i-1)} - \sum_{j=i+1}^{n} a_{ij}^{(i-1)} x_j}{a_{ii}^{(i-1)}}, \quad i = n-1, n-2, \ldots, 1
```

Since all the $`x_j`$ values for $`j > i`$ have already been computed, we substitute them in and solve for $`x_i`$.

**Pseudocode:**

```text
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

The efficiency of Gauss elimination can be quantified by counting the number of floating-point operations (multiplications and divisions, additions and subtractions).

**Forward elimination:**

At step $`k`$, each row below the pivot uses one division for the factor, $`n - k`$ multiplications to update the coefficient matrix, and one multiplication to update the right-hand side, so the number of multiplications/divisions is:

```math
\sum_{k=1}^{n-1} \sum_{i=k+1}^{n} (n - k + 2) = \frac{2n^3 + 3n^2 - 5n}{6} \approx \frac{1}{3}n^3 + O(n^2)
```

Each row also needs $`n - k + 1`$ additions/subtractions, so their number is:

```math
\sum_{k=1}^{n-1} \sum_{i=k+1}^{n} (n - k + 1) = \frac{n^3 - n}{3} \approx \frac{1}{3}n^3 + O(n^2)
```

Forward elimination therefore needs about $`\frac{2}{3}n^3 + O(n^2)`$ operations in total.

**Back substitution:**

The number of multiplications/divisions is given below, and there are $`\sum_{i=1}^{n-1} (n - i) = \frac{n^2 - n}{2}`$ additions/subtractions, so back substitution takes $`n^2`$ operations in total:

```math
\sum_{i=1}^{n-1} (n - i) + n = \frac{n^2 + n}{2} \approx \frac{1}{2}n^2 + O(n)
```

**Total cost:**

```math
\text{Total} = \frac{2}{3}n^3 + O(n^2) + n^2 + O(n) \approx \frac{2}{3}n^3
```

For large $`n`$, the total cost is **dominated by the $`(2/3)n^3`$ term** from forward elimination. Back substitution contributes only an $`O(n^2)`$ term, which is negligible for large systems.

> **[Algorithms]** Asymptotic notation $`O(\cdot)`$ keeps only the growth rate as the input size $`n`$ increases and ignores constant factors and lower-order terms. In the operation count $`\frac{2}{3}n^3 + O(n^2)`$ above, $`O(n^2)`$ collects the terms that grow in proportion to $`n^2`$, so for large $`n`$ the cubic term determines the cost, and multiplying $`n`$ by 10 multiplies the cost by about 1000.

| Phase | Dominant Term (total operations) | Growth |
|:------|:-------------|:-------|
| Forward Elimination | $`(2/3)n^3`$ | Cubic |
| Back Substitution | $`n^2`$ | Quadratic |
| **Total** | **$`(2/3)n^3`$** | **Cubic** |

*Table 2. Dominant terms of the operation count in each phase of Gauss elimination*

This is a dramatic improvement over Cramer's rule. For $`n = 20`$:

- Cramer's rule: about $`2.4 \times 10^{18}`$ operations for just one determinant.
- Gauss elimination: about $`(2/3)(20)^3 = 5333`$ operations.

---

<br>

## 7. Example: 3x3 System (Complete Step-by-Step)

### 7.1 Problem Setup

Solve the following system:

```math
3x_1 - 0.1x_2 - 0.2x_3 = 7.85
```
```math
0.1x_1 + 7x_2 - 0.3x_3 = -19.3
```
```math
0.3x_1 - 0.2x_2 + 10x_3 = 71.4
```

The augmented matrix is:

```math
[A|b] = \left[\begin{array}{ccc|c} 3 & -0.1 & -0.2 & 7.85 \\ 0.1 & 7 & -0.3 & -19.3 \\ 0.3 & -0.2 & 10 & 71.4 \end{array}\right]
```

### 7.2 Forward Elimination: Step 1

**Pivot element:** $`a_{11} = 3`$

**Eliminate $`x_1`$ from Row 2:**

```math
\text{factor}_{21} = \frac{a_{21}}{a_{11}} = \frac{0.1}{3} = 0.0333\overline{3}
```

```math
\text{Row}_2 \leftarrow \text{Row}_2 - 0.0333\overline{3} \times \text{Row}_1
```

- $`a_{22}^{(1)} = 7 - 0.0333\overline{3}(-0.1) = 7.00333\overline{3}`$
- $`a_{23}^{(1)} = -0.3 - 0.0333\overline{3}(-0.2) = -0.29333\overline{3}`$
- $`b_2^{(1)} = -19.3 - 0.0333\overline{3}(7.85) = -19.5617`$

**Eliminate $`x_1`$ from Row 3:**

```math
\text{factor}_{31} = \frac{a_{31}}{a_{11}} = \frac{0.3}{3} = 0.1
```

```math
\text{Row}_3 \leftarrow \text{Row}_3 - 0.1 \times \text{Row}_1
```

- $`a_{32}^{(1)} = -0.2 - 0.1(-0.1) = -0.19`$
- $`a_{33}^{(1)} = 10 - 0.1(-0.2) = 10.02`$
- $`b_3^{(1)} = 71.4 - 0.1(7.85) = 70.615`$

After Step 1, the augmented matrix is:

```math
\left[\begin{array}{ccc|c} 3 & -0.1 & -0.2 & 7.85 \\ 0 & 7.00333 & -0.29333 & -19.5617 \\ 0 & -0.19 & 10.02 & 70.615 \end{array}\right]
```

### 7.3 Forward Elimination: Step 2

**Pivot element:** $`a_{22}^{(1)} = 7.00333\overline{3}`$

**Eliminate $`x_2`$ from Row 3:**

```math
\text{factor}_{32} = \frac{a_{32}^{(1)}}{a_{22}^{(1)}} = \frac{-0.19}{7.00333} = -0.027130
```

```math
\text{Row}_3 \leftarrow \text{Row}_3 - (-0.027130) \times \text{Row}_2
```

- $`a_{33}^{(2)} = 10.02 - (-0.027130)(-0.29333) = 10.02 - 0.007958 = 10.0120`$
- $`b_3^{(2)} = 70.615 - (-0.027130)(-19.5617) = 70.615 - 0.530709 = 70.0843`$

After Step 2 (upper triangular form), the augmented matrix is:

```math
\left[\begin{array}{ccc|c} 3 & -0.1 & -0.2 & 7.85 \\ 0 & 7.00333 & -0.29333 & -19.5617 \\ 0 & 0 & 10.0120 & 70.0843 \end{array}\right]
```

### 7.4 Back Substitution

**Solve for $`x_3`$:**

```math
x_3 = \frac{b_3^{(2)}}{a_{33}^{(2)}} = \frac{70.0843}{10.0120} = 7.00003
```

**Solve for $`x_2`$:**

```math
x_2 = \frac{b_2^{(1)} - a_{23}^{(1)}x_3}{a_{22}^{(1)}} = \frac{-19.5617 - (-0.293333)(7.00003)}{7.00333} = -2.50000
```

**Solve for $`x_1`$:**

```math
x_1 = \frac{b_1 - a_{12}x_2 - a_{13}x_3}{a_{11}} = \frac{7.85 - (-0.1)(-2.50000) - (-0.2)(7.00003)}{3} = \frac{7.85 - 0.25 + 1.40001}{3} = \frac{9.00001}{3} = 3.00000
```

**Final solution:**

```math
x_1 = 3.00000, \quad x_2 = -2.50000, \quad x_3 = 7.00003
```

> **Note:** $`x_3`$ comes out as 7.00003 instead of the exact value 7 because of round-off error accumulated by rounding to six significant digits at each step. The exact solution of this system is $`x_1 = 3`$, $`x_2 = -2.5`$, $`x_3 = 7`$, and the computed solution is very close to it.

---

<br>

## 8. Pitfalls of Naive Gauss Elimination

The term "naive" refers to the fact that this basic algorithm has no safeguards against numerical difficulties. Three major pitfalls exist:

**1. Division by zero if the pivot element is zero:**

If $`a_{ii}^{(i-1)} = 0`$ at any step, the algorithm fails because it attempts to divide by zero when computing the factor. This can happen even when the system has a perfectly valid unique solution; it just means the equations are ordered in an unfavorable way.

**Example:** Consider the following system where $`a_{11} = 0`$:

```math
0 \cdot x_1 + 2x_2 = 4
```
```math
3x_1 + x_2 = 5
```

The solution is $`x_1 = 1, x_2 = 2`$, but naive Gauss elimination fails at the first step because the pivot $`a_{11} = 0`$.

**2. Near-zero pivots amplify round-off error:**

When a pivot element is very small but nonzero, the computed factor $`a_{ji}/a_{ii}`$ becomes very large. Multiplying row $`i`$ by this large factor and subtracting magnifies any round-off errors present in row $`i`$, potentially corrupting the result severely.

**3. Large factor values magnify errors:**

Related to the previous point, when the factor is much larger than 1, the operation $`\text{Row}_j - \text{factor} \times \text{Row}_i`$ subtracts a large multiple of one row from another. This can lead to catastrophic cancellation and significant loss of significant digits.

These three pitfalls all stem from the same root cause: **the pivot element is too small relative to other entries in the column**. The remedy is **partial pivoting**. At each step, swap rows so that the available entry with the largest absolute value in the current column becomes the pivot. This strategy, covered in Chapter 9, effectively mitigates all three problems and is used in virtually all production-quality implementations of Gauss elimination.

---

<br>

## Summary

| Concept | Key Summary |
|:------|:----------|
| System form | It has the form $`[A]\{x\} = \{b\}`$ and arises in spring, truss, and circuit analysis. |
| Graphical (2x2) | Two lines meeting at one point give a unique solution, parallel lines give no solution, and identical lines give infinitely many solutions. |
| Determinant | It is $`D = a_{11}a_{22} - a_{12}a_{21}`$, and $`D = 0`$ means the system is singular. |
| Cramer's Rule | It computes $`x_i = \det(A_i)/\det(A)`$, but its cost grows as $`O((n+1)!)`$, making it impractical for $`n > 3`$. |
| Forward Elimination | It zeros out the below-diagonal entries column by column to reach upper triangular form. |
| Back Substitution | It solves from $`x_n`$ upward using the upper triangular form. |
| Computational Cost | The total operation count is about $`(2/3)n^3`$, dominated by forward elimination. |
| Pitfalls | A zero or small pivot causes division by zero or error amplification. |
| Remedy | Partial pivoting swaps rows to avoid small pivots and is covered in Chapter 9. |

---

<br>

## Self-Check Questions

1. **Cramer's rule:** Solve $`2x_1 + x_2 = 5`$, $`x_1 - x_2 = 1`$ by Cramer's rule.

   <details>
   <summary>Answer</summary>

   $`D = (2)(-1) - (1)(1) = -3`$. Then $`x_1 = \frac{a_{22}b_1 - a_{12}b_2}{D} = \frac{(-1)(5) - (1)(1)}{-3} = 2`$ and $`x_2 = \frac{a_{11}b_2 - a_{21}b_1}{D} = \frac{(2)(1) - (1)(5)}{-3} = 1`$. The two lines meet at the point $`(2, 1)`$.

   </details>

2. **Solvability:** Find the determinant of $`x_1 + 2x_2 = 3`$, $`2x_1 + 4x_2 = 7`$ and decide how many solutions the system has. What happens if the right-hand side of the second equation is 6?

   <details>
   <summary>Answer</summary>

   $`D = (1)(4) - (2)(2) = 0`$, so the system is singular. The coefficient rows are proportional (by a factor of 2), but $`7 \neq 2 \times 3`$, so the lines are parallel and there is no solution (inconsistent). If the right-hand side is 6, the second equation is twice the first, the lines coincide, and there are infinitely many solutions (dependent).

   </details>

3. **Naive Gauss elimination:** Solve the following system by naive Gauss elimination: $`2x_1 + x_2 - x_3 = 8`$, $`-3x_1 - x_2 + 2x_3 = -11`$, $`-2x_1 + x_2 + 2x_3 = -3`$.

   <details>
   <summary>Answer</summary>

   Eliminating column 1 with factors $`-3/2 = -1.5`$ and $`-2/2 = -1`$ turns row 2 into $`0.5x_2 + 0.5x_3 = 1`$ and row 3 into $`2x_2 + x_3 = 5`$. Eliminating column 2 with the factor $`2/0.5 = 4`$ turns row 3 into $`-x_3 = 1`$. Back substitution gives $`x_3 = -1`$, $`x_2 = (1 - 0.5(-1))/0.5 = 3`$, and $`x_1 = (8 - x_2 + x_3)/2 = (8 - 3 - 1)/2 = 2`$.

   </details>

4. **Growth of the cost:** If $`n`$ grows tenfold from 100 to 1000, by roughly what factor does the operation count of Gauss elimination grow? Also give the operation count in each case.

   <details>
   <summary>Answer</summary>

   The operation count is about $`(2/3)n^3`$, so it is about $`6.7 \times 10^5`$ for $`n = 100`$ and about $`6.7 \times 10^8`$ for $`n = 1000`$. A tenfold increase in $`n`$ multiplies the cost by $`10^3 = 1000`$. The $`n^2`$ term of back substitution, $`10^4`$ and $`10^6`$ respectively, is small compared with forward elimination.

   </details>
