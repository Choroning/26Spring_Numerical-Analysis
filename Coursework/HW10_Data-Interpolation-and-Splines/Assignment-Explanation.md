# HW10 Assignment Explanation — Data Interpolation and Splines

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 17, 18

> **Learning Objectives:**
> 1. Implement the product formula of the Lagrange basis functions with array operations
> 2. Verify the implementation with the properties that each basis function is 1 or 0 at the nodes and that the basis functions sum to 1
> 3. Explain why Lagrange and Newton interpolation give the same value
> 4. Interpret the basis-function values as weights of the data values and compute interpolated values by hand

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem: Lagrange Basis Functions](#3-problem-lagrange-basis-functions)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
  - [3.4 Results Used by Later Cells](#34-results-used-by-later-cells)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the single blank of the lab notebook `Lab11-Spline.ipynb` distributed by the instructor and run it. The notebook has seven parts: why interpolation is needed, Newton's divided differences, Lagrange interpolation, Runge's phenomenon, quadratic splines, cubic splines, and `CubicSpline`, and the other cells were distributed with their code filled in. The function completed in the blank is used in the cell of Section 3 that plots the basis functions and compares them with Newton interpolation, and in Exercise 3.1. The problem statement is in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab11-Spline.ipynb](Original/Lab11-Spline.ipynb), and the notebook with the blank filled in and run is [Solution/Lab11-Spline_Completed.ipynb](Solution/Lab11-Spline_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem | 3.1 Theory of Lagrange interpolation | The `Li` update in the function `lagrange_basis` |

*Table 1. The blank that makes up the assignment*

The Newton and Lagrange interpolating polynomials are covered in Sections 2 and 3 of the [Chapter 17 lecture note](../../Lecture-Notes/Ch17_Interpolation/Concepts_Lecture.md), and the splines treated in the later part of the notebook are covered in the [Chapter 18 lecture note](../../Lecture-Notes/Ch18_Splines/Concepts_Lecture.md). This document repeats only what the solution needs.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of the cell stays as distributed, and only the blank is filled in.

1. `lagrange_basis(x, i, xx)` must return $`L_i(xx) = \prod_{j \neq i} (xx - x_j)/(x_i - x_j)`$. The blank runs once for every iteration with $`j \neq i`$, so it must multiply in one factor at a time.
2. `np.atleast_1d` turns `xx` into an array, so the computation must be elementwise over every entry of `xx`. The cell that plots the basis functions as curves passes an array of 200 points.
3. `lagrange_interp` computes $`\sum_i y_i L_i(xx)`$, so the value at 350 K must equal that of Newton interpolation.

---

<br>

## 3. Problem: Lagrange Basis Functions

### 3.1 Approach

The product $`\prod_{j \neq i}`$ is computed as a running product in the loop. `Li` starts from `np.ones(xx.shape)`; 1 is the multiplicative identity, so it is the value of the empty product before any factor is multiplied in. Then, for every node with $`j \neq i`$, the factor $`(xx - x_j)/(x_i - x_j)`$ is multiplied in.

$$
L_i \leftarrow L_i \cdot \frac{xx - x_j}{x_i - x_j}
$$

`xx` and `Li` are arrays of the same shape and `x[i]` and `x[j]` are scalars, so NumPy's elementwise operations evaluate all points at once. The denominator $`x_i - x_j`$ is nonzero as long as the nodes are distinct. The Lagrange interpolating polynomial is covered in Section 3 of the [Chapter 17 lecture note](../../Lecture-Notes/Ch17_Interpolation/Concepts_Lecture.md).

### 3.2 Code Written

```python
            Li = Li * (xx - x[j]) / (x[i] - x[j])   # Lagrange basis
```

### 3.3 Results and Verification

```text
Lagrange at 350 K = 0.56762500
Newton   at 350 K = 0.56762500
Difference        = 0.00e+00  (should be ~machine precision)
```

The output matches the output in the distributed notebook. With nodes at 300, 400, and 500 K, the basis-function values at 350 K computed by hand are as follows.

$$
L_0 = \frac{(350-400)(350-500)}{(300-400)(300-500)} = 0.375, \quad L_1 = \frac{(350-300)(350-500)}{(400-300)(400-500)} = 0.75, \quad L_2 = \frac{(350-300)(350-400)}{(500-300)(500-400)} = -0.125
$$

The interpolated value is therefore $`0.616 \times 0.375 + 0.525 \times 0.75 + 0.457 \times (-0.125) = 0.567625`$, matching the output. Three further checks confirm the result.

- **Values at the nodes:** Evaluating each basis function at the three nodes with the same function gives 1 at its own node and 0 at the others. In the basis-function plot drawn by the cell, each curve passing through 1 at its own node and 0 at the remaining nodes shows this property.
- **Sum equal to 1:** At 350 K, $`0.375 + 0.75 - 0.125 = 1`$, and over all 200 evaluation points of the plot the sum of the three basis functions differs from 1 by at most $`3.3 \times 10^{-16}`$. Interpolating the constant function $`y = 1`$ gives the constant 1 itself as the unique interpolating polynomial, so $`\sum_i L_i(x) = 1`$ holds for every $`x`$.
- **Agreement with Newton interpolation:** Only one polynomial of degree at most 2 passes through three points, so the Lagrange and Newton forms are the same polynomial. In this case even the floating-point results coincide and the difference prints as exactly 0, but in general the two computations use different orders of operations and can differ at the level of machine epsilon. The comment in the cell expects that level.

### 3.4 Results Used by Later Cells

Exercise 3.1 adds the point at 600 K, builds a cubic Lagrange polynomial through the four points, and interpolates at 350, 450, and 550 K.

```text
ρ(350 K) ≈ 0.56712 kg/m³
ρ(450 K) ≈ 0.48863 kg/m³
ρ(550 K) ≈ 0.42912 kg/m³
```

Repeating the computation gives the following basis-function values at each evaluation point.

| Evaluation point | $`L_0`$ | $`L_1`$ | $`L_2`$ | $`L_3`$ | Interpolated value |
|:-----|:-----|:-----|:-----|:-----|:-----|
| 350 K | 0.3125 | 0.9375 | -0.3125 | 0.0625 | 0.567125 |
| 450 K | -0.0625 | 0.5625 | 0.5625 | -0.0625 | 0.488625 |
| 550 K | 0.0625 | -0.3125 | 0.9375 | 0.3125 | 0.429125 |

*Table 2. Basis-function values and interpolated values for four nodes*

Each row sums to 1. At 450 K, the middle of four equally spaced nodes, the weights are the symmetric $`(-1, 9, 9, -1)/16`$, and the outer nodes receive negative weights. The value at 350 K changes by 0.0005, from 0.567625 with three points to 0.567125, because adding a node redefines every basis function in the Lagrange form. This is what the comparison in Table 1 means by saying that Lagrange interpolation must recompute everything when a point is added.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Lagrange basis function | $`L_i(x)`$ is the product of the factors $`(x - x_j)/(x_i - x_j)`$ over the nodes with $`j \neq i`$, computed as a running product starting from 1. |
| Property at the nodes | Since $`L_i(x_j) = \delta_{ij}`$, the interpolating polynomial $`\sum_i y_i L_i(x)`$ passes exactly through every node. |
| Sum equal to 1 | The interpolating polynomial of a constant function is that constant itself, so the basis functions sum to 1 at every point, which can be used to verify the implementation. |
| Uniqueness | Only one polynomial of degree at most $`n`$ passes through $`n+1`$ points, so Lagrange and Newton interpolation give the same value. |
| Weight interpretation | The interpolated value is a weighted sum of the data values, the weights can be negative, and adding a node changes every weight. |

---

<br>

## Self-Check Questions

1. **Starting from 1:** In `lagrange_basis`, why must `Li` start from `np.ones` rather than `np.zeros`?

   <details>
   <summary>Answer</summary>

   `Li` is a running product built by multiplying in factors one after another, so it must start from 1, the multiplicative identity. Starting from 0 would make the result 0 whatever factors are multiplied in, so every basis function would be 0.

   </details>

2. **Why the sum is 1:** Explain why $`\sum_i L_i(x) = 1`$ holds for every $`x`$, whatever the nodes are.

   <details>
   <summary>Answer</summary>

   Interpolating data whose value is 1 at every node gives $`\sum_i 1 \cdot L_i(x)`$, and the constant function 1 is also a polynomial of degree 0 through those nodes. The polynomial of degree at most $`n`$ through $`n+1`$ points is unique, so the two polynomials are equal, and the basis functions therefore always sum to 1.

   </details>

3. **Hand calculation:** With nodes at 300, 400, and 500 K, compute $`L_1(350)`$.

   <details>
   <summary>Answer</summary>

   $`L_1(350) = \frac{(350 - 300)(350 - 500)}{(400 - 300)(400 - 500)} = \frac{50 \times (-150)}{100 \times (-100)} = 0.75`$.

   </details>

4. **Negative weights:** In Exercise 3.1, the data at 300 K and 600 K receive negative weights when the value at 450 K is interpolated. What does this mean?

   <details>
   <summary>Answer</summary>

   The interpolated value is not a convex combination of the data values, so there is no guarantee that it lies between the smallest and largest data values. The weights $`(-1, 9, 9, -1)/16`$ add to the average of the two near points, $`(y_1 + y_2)/2`$, a correction $`(y_1 + y_2 - y_0 - y_3)/16`$ built from the differences with the two outer points, which accounts for curvature. Because of the negative weights, the values of a high-degree polynomial can leave the range of the data between the nodes, and Runge's phenomenon is an extreme example.

   </details>
