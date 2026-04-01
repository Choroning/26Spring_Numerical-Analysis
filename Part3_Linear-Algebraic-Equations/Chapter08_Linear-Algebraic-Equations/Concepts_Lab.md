# Chapter 8 Lab — Linear Systems and Gaussian Elimination

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 8

> **Prerequisites**: [Programming Language] MATLAB/Python. [Linear Algebra] Matrix operations. [Calculus] Linear systems (Ch 1-7).
>
> **Learning Objectives**:
> 1. Formulate engineering problems as linear algebraic equations
> 2. Identify solution existence and uniqueness conditions
> 3. Apply graphical interpretation of linear systems

---

<br>

## Table of Contents

- [1. Setting Up $Ax = b$](#1-setting-up-ax--b)
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
  - [4.2 Verifying O(n^3) Scaling](#42-verifying-on3-scaling)
- [Summary](#summary)

---

<br>

## 1. Setting Up $Ax = b$

### 1.1 Matrix Form

A system of $n$ linear equations in $n$ unknowns can be written compactly in matrix-vector form:

$$[A]\{x\} = \{b\}$$

where:
- $[A]$ is the $n \times n$ **coefficient matrix** containing the known coefficients
- $\{x\}$ is the $n \times 1$ **unknown vector** containing the variables we wish to find
- $\{b\}$ is the $n \times 1$ **right-hand side (RHS) vector** containing the known constants

For example, the system:

$$2x_1 + x_2 = 5$$
$$x_1 - x_2 = 1$$

is expressed as:

$$\begin{bmatrix} 2 & 1 \\ 1 & -1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 5 \\ 1 \end{bmatrix}$$

### 1.2 Geometric Interpretation (2x2)

For a $2 \times 2$ system, each equation defines a line in the plane. The three possible outcomes — unique solution, no solution, and infinite solutions — can be visualized directly:

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

In the unique-solution case, the intersection point $(2, 1)$ is marked with a black dot. The parallel-lines case shows two lines with identical slope but different intercepts — they never meet. The infinite-solutions case shows two equations that describe the exact same line (one is simply double the other).

---

<br>

## 2. Spring-Mass System Application

### 2.1 Problem Setup

Consider a system of 3 masses connected by 4 springs between two fixed walls. The springs have constants $k_1 = 100$, $k_2 = 200$, $k_3 = 150$, $k_4 = 250$ N/m, and the applied forces are $F_1 = 50$, $F_2 = 0$, $F_3 = -30$ N.

Applying Hooke's law ($F = kx$) and force equilibrium at each mass yields the stiffness matrix equation $[K]\{x\} = \{F\}$:

$$\begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} F_1 \\ F_2 \\ F_3 \end{bmatrix}$$

Each diagonal entry $(k_i + k_{i+1})$ represents the total stiffness acting on mass $i$. The off-diagonal entries $(-k_j)$ represent the coupling between adjacent masses through the shared spring.

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

The residual $\|[K]\{x\} - \{F\}\|$ should be near machine epsilon, confirming that the computed solution satisfies the original equations to full floating-point precision.

> **[Physics]** The stiffness matrix K is symmetric positive definite (SPD) for physically valid spring systems. Each diagonal entry is the sum of spring constants connected to that mass. Off-diagonal entries represent coupling between masses through shared springs.

The SPD property guarantees that the system always has a unique solution, the displacements are finite for finite forces, and specialized solvers (such as Cholesky decomposition) can exploit the symmetry for improved efficiency.

---

<br>

## 3. Naive Gaussian Elimination Implementation

### 3.1 Algorithm Implementation

The following function implements naive Gaussian elimination with forward elimination and back substitution. It operates on the **augmented matrix** $[A|b]$ to avoid managing $A$ and $b$ separately:

```python
import numpy as np

def gauss_naive(A, b):
    """Solve Ax = b using naive Gaussian elimination.

    Args:
        A: n×n coefficient matrix (will be modified)
        b: n×1 right-hand side vector (will be modified)
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
- **Augmented matrix**: Combining $[A|b]$ into a single array ensures that row operations are applied consistently to both the coefficient matrix and the right-hand side.
- **Vectorized inner loop**: The line `Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]` uses NumPy slicing to update an entire row at once, avoiding an explicit inner loop over columns.

### 3.2 Testing with the Textbook Example

```python
# Test with the 3x3 system from the lecture
A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
b = np.array([7.85, -19.3, 71.4])
x = gauss_naive(A, b)
print(f'Solution: {x}')  # [3.0, -2.5, 7.0003...]

# Compare with NumPy's built-in solver
x_np = np.linalg.solve(A, b)
print(f'NumPy solution: {x_np}')
print(f'Difference: {np.abs(x - x_np)}')
```

Both methods should produce essentially the same result. The small numerical differences (on the order of machine epsilon) arise from different internal implementations — `np.linalg.solve` uses LAPACK's optimized routines with partial pivoting, while our naive implementation does not.

---

<br>

## 4. Computational Cost Analysis

### 4.1 Empirical Timing

We can empirically verify that `np.linalg.solve` scales as $O(n^3)$ by timing it on matrices of increasing size:

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

### 4.2 Verifying O(n^3) Scaling

On a log-log plot, an $O(n^3)$ algorithm appears as a straight line with slope approximately 3:

```python
import matplotlib.pyplot as plt

plt.loglog(sizes, times, 'bo-')
plt.xlabel('Matrix size n')
plt.ylabel('Time (s)')
plt.title('np.linalg.solve: O(n³) scaling')
plt.grid(True)
plt.show()
```

The slope of the log-log plot can be estimated from any two data points. If the time at size $n_1$ is $t_1$ and at size $n_2$ is $t_2$, then:

$$\text{slope} = \frac{\log(t_2/t_1)}{\log(n_2/n_1)} \approx 3$$

For small matrices, overhead costs (function call, memory allocation) may dominate, causing the slope to appear less than 3. For very large matrices, the slope approaches 3 as the cubic term dominates all lower-order terms.

---

<br>

## Summary

| Method | Complexity | Pivoting | Notes |
|:-------|:----------|:---------|:------|
| Cramer's Rule | $O(n!)$ | N/A | Only for tiny systems |
| Naive Gauss | $O(n^3)$ | No | May fail with zero pivot |
| `np.linalg.solve` | $O(n^3)$ | Yes | Optimized LAPACK, preferred |

---
