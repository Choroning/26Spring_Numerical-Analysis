# Chapter 8 Lab — Linear Systems and Gaussian Elimination

[Repository Guide](../../README.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **Last Updated:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 8

> **Learning Objectives:**
> 1. Formulate engineering problems as linear algebraic equations
> 2. Identify solution existence and uniqueness conditions
> 3. Apply graphical interpretation of linear systems
> 4. Implement naive Gauss elimination and confirm its computational cost experimentally

---

<br>

## Table of Contents

- [1. Setting Up $`Ax = b`$](#1-setting-up-ax--b)
  - [1.1 Matrix Form](#11-matrix-form)
  - [1.2 Geometric Interpretation (2x2)](#12-geometric-interpretation-2x2)
- [2. Spring-Mass System Application](#2-spring-mass-system-application)
  - [2.1 Problem Setup](#21-problem-setup)
  - [2.2 Solving and Verification](#22-solving-and-verification)
- [3. Naive Gaussian Elimination Implementation](#3-naive-gaussian-elimination-implementation)
  - [3.1 Algorithm Implementation](#31-algorithm-implementation)
  - [3.2 Testing with the Textbook Example](#32-testing-with-the-textbook-example)
- [4. Computational Cost Analysis](#4-computational-cost-analysis)
  - [4.1 Empirical Timing](#41-empirical-timing)
  - [4.2 Verifying $`O(n^3)`$ Scaling](#42-verifying-on3-scaling)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Setting Up $`Ax = b`$

### 1.1 Matrix Form

A system of $`n`$ linear equations in $`n`$ unknowns can be written compactly in matrix-vector form:

```math
[A]\{x\} = \{b\}
```

where:
- $`[A]`$ is the $`n \times n`$ **coefficient matrix** containing the known coefficients.
- $`\{x\}`$ is the $`n \times 1`$ **unknown vector** containing the variables we wish to find.
- $`\{b\}`$ is the $`n \times 1`$ **right-hand side (RHS) vector** containing the known constants.

For example, the system:

```math
2x_1 + x_2 = 5
```
```math
x_1 - x_2 = 1
```

is expressed as:

```math
\begin{bmatrix} 2 & 1 \\ 1 & -1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 5 \\ 1 \end{bmatrix}
```

### 1.2 Geometric Interpretation (2x2)

For a $`2 \times 2`$ system, each equation defines a line in the plane. The three possible outcomes, namely a unique solution, no solution, and infinitely many solutions, can therefore be visualized directly:

```python
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Case 1: Unique solution (2x + y = 5, x - y = 1)
x = np.linspace(-1, 4, 100)
axes[0].plot(x, 5 - 2*x, label='2x + y = 5')
axes[0].plot(x, x - 1, label='x - y = 1')
axes[0].plot(2, 1, 'ko', markersize=8)
axes[0].set_title('Unique Solution')
axes[0].legend()
axes[0].grid(True)

# Case 2: No solution (parallel: y = 2x + 1, y = 2x + 3)
axes[1].plot(x, 2*x + 1, label='y = 2x + 1')
axes[1].plot(x, 2*x + 3, label='y = 2x + 3')
axes[1].set_title('No Solution (Parallel)')
axes[1].legend()
axes[1].grid(True)

# Case 3: Infinite solutions (y = 2x + 1, 2y = 4x + 2)
axes[2].plot(x, 2*x + 1, label='y = 2x + 1', linewidth=3)
axes[2].plot(x, 2*x + 1, '--', label='2y = 4x + 2', linewidth=1)
axes[2].set_title('Infinite Solutions')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.show()
```

In the unique-solution case, the intersection point $`(2, 1)`$ is marked with a black dot. The parallel-lines case shows two lines with identical slope but different intercepts, which never meet. The infinite-solutions case shows two equations that describe the exact same line (one is simply double the other).

---

<br>

## 2. Spring-Mass System Application

### 2.1 Problem Setup

Consider a system of 3 masses connected by 4 springs between two fixed walls (Figure 1). The springs have constants $`k_1 = 100`$, $`k_2 = 200`$, $`k_3 = 150`$, $`k_4 = 250`$ N/m, and the applied forces are $`F_1 = 50`$, $`F_2 = 0`$, $`F_3 = -30`$ N. The displacements $`x_i`$ and forces $`F_i`$ are positive in the same direction.

```text
||--k1--[m1]--k2--[m2]--k3--[m3]--k4--||
```

*Figure 1. Three masses connected by four springs between two fixed walls*

Applying Hooke's law ($`F = kx`$) and force equilibrium at each mass yields the stiffness matrix equation $`[K]\{x\} = \{F\}`$:

```math
\begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} F_1 \\ F_2 \\ F_3 \end{bmatrix}
```

Each diagonal entry $`(k_i + k_{i+1})`$ represents the total stiffness acting on mass $`i`$. The off-diagonal entries $`(-k_j)`$ represent the coupling between adjacent masses through the shared spring.

### 2.2 Solving and Verification

```python
import numpy as np

# 3 masses, 4 springs between walls
# k1=100, k2=200, k3=150, k4=250 (N/m)
# F1=50, F2=0, F3=-30 (N)

# Stiffness matrix
K = np.array([
    [100+200, -200, 0],
    [-200, 200+150, -150],
    [0, -150, 150+250]
], dtype=float)

F = np.array([50.0, 0.0, -30.0])

x = np.linalg.solve(K, F)
print(f'Displacements: x = {x}')

# Verify: residual should be near zero
residual = K @ x - F
print(f'Residual norm: {np.linalg.norm(residual):.2e}')
```

The residual $`\|[K]\{x\} - \{F\}\|`$ should be near machine epsilon, confirming that the computed solution satisfies the original equations to full floating-point precision. For this system, the displacements are $`x_1 \approx 0.2584`$, $`x_2 \approx 0.1377`$, and $`x_3 \approx -0.0234`$ m, and the residual norm is about $`4 \times 10^{-16}`$.

> **Note:** **Properties of the stiffness matrix.** When all spring constants are positive and the system is anchored to the walls, the stiffness matrix K is symmetric positive definite (SPD). Each diagonal entry is the sum of spring constants connected to that mass. Off-diagonal entries represent coupling between masses through shared springs.

The SPD property guarantees that the system always has a unique solution, the displacements are finite for finite forces, and specialized solvers (such as Cholesky decomposition) can exploit the symmetry for improved efficiency.

---

<br>

## 3. Naive Gaussian Elimination Implementation

### 3.1 Algorithm Implementation

The following function implements naive Gaussian elimination with forward elimination and back substitution. It operates on the **augmented matrix** $`[A|b]`$ to avoid managing $`A`$ and $`b`$ separately:

```python
import numpy as np

def gauss_naive(A, b):
    """Solve Ax = b using naive Gaussian elimination.

    Args:
        A: n×n coefficient matrix (not modified)
        b: n×1 right-hand side vector (not modified)
    Returns:
        x: solution vector
    """
    A = A.astype(float).copy()
    b = b.astype(float).copy()
    n = len(b)

    # Augmented matrix [A|b]
    Aug = np.hstack([A, b.reshape(-1, 1)])
    nb = n + 1

    # Forward Elimination
    for k in range(n - 1):
        for i in range(k + 1, n):
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]

    # Back Substitution
    x = np.zeros(n)
    x[n-1] = Aug[n-1, nb-1] / Aug[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, nb-1] - Aug[i, i+1:n] @ x[i+1:n]) / Aug[i, i]

    return x
```

Key implementation details:

- **Copying inputs**: `A.astype(float).copy()` ensures the original arrays are not modified and that integer arrays are promoted to float for correct division.
- **Augmented matrix**: Combining $`[A|b]`$ into a single array ensures that row operations are applied consistently to both the coefficient matrix and the right-hand side.
- **Vectorized inner loop**: The line `Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]` uses NumPy slicing to update an entire row at once, avoiding an explicit inner loop over columns.

### 3.2 Testing with the Textbook Example

We test the implementation on the $`3 \times 3`$ system solved by hand in Section 7 of the Chapter 8 lecture note:

```python
# Test with the 3x3 system from the lecture
A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
b = np.array([7.85, -19.3, 71.4])
x = gauss_naive(A, b)
print(f'Solution: {x}')  # [ 3.  -2.5  7. ]

# Compare with NumPy's built-in solver
x_np = np.linalg.solve(A, b)
print(f'NumPy solution: {x_np}')
print(f'Difference: {np.abs(x - x_np)}')
```

In double precision, the result is the exact solution $`(3, -2.5, 7)`$ to within machine epsilon, rather than the 7.00003 of the hand calculation. Both methods produce essentially the same result, and the small differences (on the order of machine epsilon) arise from different internal implementations. `np.linalg.solve` uses LAPACK's optimized routines with partial pivoting, while this naive implementation does not pivot.

---

<br>

## 4. Computational Cost Analysis

### 4.1 Empirical Timing

We can empirically verify that `np.linalg.solve` scales as $`O(n^3)`$ by timing it on matrices of increasing size:

```python
import numpy as np
import time

sizes = [50, 100, 200, 400, 800]
times = []

for n in sizes:
    A = np.random.rand(n, n)
    b = np.random.rand(n)

    start = time.time()
    for _ in range(10):
        np.linalg.solve(A, b)
    elapsed = (time.time() - start) / 10
    times.append(elapsed)
    print(f'n={n:4d}: {elapsed:.6f}s')
```

### 4.2 Verifying $`O(n^3)`$ Scaling

On a log-log plot, an $`O(n^3)`$ algorithm appears as a straight line with slope approximately 3:

```python
import matplotlib.pyplot as plt

plt.loglog(sizes, times, 'bo-')
plt.xlabel('Matrix size n')
plt.ylabel('Time (s)')
plt.title('np.linalg.solve: O(n³) scaling')
plt.grid(True)
plt.show()
```

The slope of the log-log plot can be estimated from any two data points. If the time at size $`n_1`$ is $`t_1`$ and at size $`n_2`$ is $`t_2`$, then:

```math
\text{slope} = \frac{\log(t_2/t_1)}{\log(n_2/n_1)} \approx 3
```

For small matrices, overhead costs (function call, memory allocation) may dominate, causing the slope to appear less than 3. For very large matrices, the slope approaches 3 as the cubic term dominates all lower-order terms.

---

<br>

## Summary

| Concept | Key Summary |
|:-------|:----------|
| Matrix form | A system of $`n`$ linear equations is written compactly as $`[A]\{x\} = \{b\}`$. |
| Geometry of 2x2 systems | Lines meeting at one point give a unique solution, parallel lines give none, and identical lines give infinitely many. |
| Stiffness matrix | The stiffness matrix of a spring-mass system is symmetric positive definite, so the system has a unique solution. |
| Residual check | A residual near machine epsilon confirms that the computed solution satisfies the original equations. |
| Cramer's Rule | Its cost is $`O((n+1)!)`$, pivoting does not apply, and it is used only for tiny systems. |
| Naive Gauss | Its cost is $`O(n^3)`$, and because it does not pivot it may fail with a zero pivot. |
| `np.linalg.solve` | It costs $`O(n^3)`$, uses an optimized LAPACK routine with partial pivoting, and is preferred in practice. |
| Scaling check | A slope of about 3 for run time on a log-log plot confirms $`O(n^3)`$ scaling. |

---

<br>

## Self-Check Questions

1. **Building a stiffness matrix:** Two masses are connected by three springs ($`k_1 = k_2 = k_3 = 100`$ N/m) between two fixed walls, with $`F_1 = 10`$ N and $`F_2 = 0`$ N. Write the stiffness matrix and find the displacements.

   <details>
   <summary>Answer</summary>

   $`[K] = \begin{bmatrix} 200 & -100 \\ -100 & 200 \end{bmatrix}`$. Since $`\det [K] = 200^2 - 100^2 = 30000`$, Cramer's rule gives $`x_1 = (200 \times 10 - (-100) \times 0)/30000 \approx 0.0667`$ m and $`x_2 = (200 \times 0 - (-100) \times 10)/30000 \approx 0.0333`$ m. Mass 1, where the force is applied, moves more, and mass 2 is pulled along through spring 2.

   </details>

2. **Log-log slope:** `np.linalg.solve` took 0.010 s at $`n = 400`$ and 0.075 s at $`n = 800`$. Estimate the scaling slope.

   <details>
   <summary>Answer</summary>

   $`\text{slope} = \log(0.075/0.010)/\log(800/400) = \log 7.5/\log 2 \approx 2.015/0.693 \approx 2.91`$. This is close to 3, consistent with $`O(n^3)`$ scaling.

   </details>

3. **Zero pivot:** What happens when $`A = \begin{bmatrix} 0 & 1 \\ 1 & 1 \end{bmatrix}`$ and $`b = (1, 2)`$ are passed to `gauss_naive`? What is the solution of this system?

   <details>
   <summary>Answer</summary>

   The first pivot `Aug[0, 0]` is zero, so computing the factor divides by zero. NumPy issues a warning and produces `inf` and `nan`, and the result is `[nan nan]`. The system itself is nonsingular, however, with solution $`x_2 = 1`$ and $`x_1 = 2 - 1 = 1`$. Swapping the two rows makes it solvable, and the partial pivoting of Chapter 9 systematizes this idea.

   </details>

4. **Limits of the residual:** If the residual is as small as machine epsilon, is the error of the solution always small?

   <details>
   <summary>Answer</summary>

   No. The residual only shows how well the computed solution satisfies the equations. For an ill-conditioned system, the relative error of the solution can be as large as the condition number times the relative residual, so the error can be large even when the residual is small. The condition number is covered in the Chapter 9 lab and in Chapter 11.

   </details>
