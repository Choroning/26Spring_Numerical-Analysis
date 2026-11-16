# HW5 Assignment Explanation — Eigenvalues and Eigenvectors

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 13

> **Learning Objectives:**
> 1. Compute eigenpairs with `numpy.linalg.eig` and verify them with the trace, the determinant, and the residuals
> 2. Complete the power method with infinity-norm normalization and explain its convergence rate by the ratio of eigenvalues
> 3. Complete the shifted power method that reuses one LU factorization and explain why the iteration count depends on the shift
> 4. Compute natural frequencies and periods from the result of a generalized eigenvalue problem
> 5. Find the modal constants of a linear ODE system and interpret its stability and time constants from the eigenvalues

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Computing Eigenvalues and Eigenvectors](#3-problem-1-computing-eigenvalues-and-eigenvectors)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Power Method](#4-problem-2-power-method)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: Inverse and Shifted Power Methods](#5-problem-3-inverse-and-shifted-power-methods)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Problem 4: Natural Frequencies of a Vibrating System](#6-problem-4-natural-frequencies-of-a-vibrating-system)
  - [6.1 Approach](#61-approach)
  - [6.2 Code Written](#62-code-written)
  - [6.3 Results and Verification](#63-results-and-verification)
- [7. Problem 5: Modal Constants of an ODE System](#7-problem-5-modal-constants-of-an-ode-system)
  - [7.1 Approach](#71-approach)
  - [7.2 Code Written](#72-code-written)
  - [7.3 Results and Verification](#73-results-and-verification)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the five blanks of the lab notebook `Lab06-Eigen.ipynb` distributed by the instructor and run it. The notebook has six parts: eigenvalues and eigenvectors, the power method, the inverse and shifted power methods, symmetric matrices and spectral decomposition, an application to vibration analysis, and an application to ODE systems and stability, and each part except the one on symmetric matrices contains one blank. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab06-Eigen.ipynb](Original/Lab06-Eigen.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab06-Eigen_Completed.ipynb](Solution/Lab06-Eigen_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 1.3 Using `numpy.linalg.eig` | `eigvals`, `eigvecs` |
| Problem 2 | 2.1 Power method | `lam` and `x` in the function `power_method` |
| Problem 3 | 3.1, 3.2 Inverse and shifted power methods | `lu_obj`, `y`, and `lam` in the function `shifted_power_method` |
| Problem 4 | 5.1 Vibration analysis | `omega_n`, `freq_hz` |
| Problem 5 | 6.1, 6.2 ODE systems and stability | `c_ode` |

*Table 1. The five blanks that make up the assignment*

The power, inverse power, and shifted power methods are covered in the [Chapter 12 lab note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lab.md), and symmetric eigenvalue problems, vibration analysis, and the stability of ODE systems are covered in the [Chapter 13 lab note](../../Lecture-Notes/Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lab.md) and the [Chapter 13 lecture note](../../Lecture-Notes/Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lecture.md). This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell (output format, loop structure, plots) stays as distributed, and only the blanks are filled in.

1. **Problem 1:** `eigvals` must be the array of eigenvalues of the matrix `A`, and `eigvecs` must be the matrix whose columns are the eigenvectors. The rest of the cell treats `eigvecs[:, i]` as the eigenvector of the $`i`$-th eigenvalue when computing the residuals.
2. **Problem 2:** In each iteration, the component of `y` with the largest magnitude must become the eigenvalue estimate `lam`, and `y` divided by that value must become the `x` of the next iteration.
3. **Problem 3:** `B` must be LU-factored only once before the loop, each iteration must solve $`B\mathbf{y} = \mathbf{x}`$ with that factorization, and the eigenvalue of the original matrix must be recovered from the normalization factor `mu`.
4. **Problem 4:** `omega_n` must be the array of angular frequencies (rad/s) and `freq_hz` the array of frequencies (Hz). The output cell computes the period as `1/freq_hz`.
5. **Problem 5:** `c_ode` must be the modal constants satisfying $`V\mathbf{c} = \mathbf{y}_0`$, and the rest of the cell builds the time solution from these values.

---

<br>

## 3. Problem 1: Computing Eigenvalues and Eigenvectors

### 3.1 Approach

`np.linalg.eig(A)` returns two values, the array of eigenvalues and the matrix of eigenvectors, where the $`i`$-th column of the matrix is the eigenvector of the $`i`$-th eigenvalue and each column has 2-norm 1. This order matches the verification code `A @ eigvecs[:, i] - eigvals[i] * eigvecs[:, i]` in the cell, so the two values are unpacked as they are. Eigenvalues and the characteristic equation are covered in the [Chapter 13 lecture note](../../Lecture-Notes/Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lecture.md).

### 3.2 Code Written

```python
eigvals, eigvecs = np.linalg.eig(A)
```

### 3.3 Results and Verification

```text
Eigenvalues:
  lambda_1 = 21.72915618
  lambda_2 = 13.18232626
  lambda_3 = 6.08851756

Eigenvectors (columns):
[[ 0.89819   0.427012  0.104481]
 [ 0.30972  -0.44601  -0.839732]
 [ 0.311976 -0.786598  0.532855]]

Verification  max|A v - lambda v|:
  lambda_1: 3.55e-15
  lambda_2: 7.11e-15
  lambda_3: 8.88e-16
```

The eigenvalues and eigenvectors match the output of the distributed notebook. The last residual differs slightly from the distributed `1.33e-15`, but both are at the level of roundoff, so the difference can be attributed to the order of floating-point operations in the two environments. Three checks confirm the result.

- **Residuals:** For all three eigenpairs, $`\max|A\mathbf{v} - \lambda\mathbf{v}|`$ is on the order of $`10^{-15}`$, about machine epsilon times the size of the matrix entries.
- **Trace and determinant:** In the output of the following cell, $`\text{tr}(A) = 20 + 9 + 12 = 41`$ equals the sum of the eigenvalues $`21.729 + 13.182 + 6.089 = 41.000`$, and $`\det(A) = 1744`$ equals their product.
- **Orthogonality:** Since $`A`$ is symmetric, its eigenvectors are mutually orthogonal, and the $`V^T V`$ printed by the following cell is the identity matrix.

---

<br>

## 4. Problem 2: Power Method

### 4.1 Approach

The blanks transcribe steps 2 and 3 of the algorithm. The index $`i^* = \arg\max_i |y_i|`$ is `np.argmax(np.abs(y))`, and the eigenvalue estimate is the component at that index itself, `y[i*]`. Using the signed component rather than its absolute value matters. When the dominant eigenvalue is negative, the sign of $`\mathbf{y} = A\mathbf{x}`$ flips at every iteration, and only division by the signed component keeps the sign of $`\mathbf{x}`$ fixed and lets the estimate converge to the negative eigenvalue. Dividing with `x = y / lam` makes the largest component 1, which keeps $`\|\mathbf{x}\|_\infty = 1`$ and prevents overflow however many iterations run. The power method is covered in Section 2 of the [Chapter 12 lab note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lab.md).

### 4.2 Code Written

```python
        lam = y[np.argmax(np.abs(y))]
        x   = y / lam
```

### 4.3 Results and Verification

```text
Converged in 51 iterations
Dominant eigenvalue (Power Method): 21.729156182199
True value (numpy)                : 21.729156182169
Absolute error                    : 3.02e-11

Theoretical convergence rate |λ₂/λ₁| = 0.606665
  -> error halves every 1.4 iterations
```

The output matches the output of the distributed notebook. The convergence rate is $`|\lambda_2/\lambda_1| = 13.1823/21.7292 = 0.6067`$, so the error shrinks by about 0.61 per iteration and halves every $`\ln 2 / (-\ln 0.6067) = 1.39`$ iterations.

The final error of $`3.02 \times 10^{-11}`$ is also consistent with the stopping criterion. If the error shrinks as $`e_k \approx C r^k`$ $`(r = 0.6067)`$, the change between successive estimates is $`e_{k-1} - e_k = e_k (1/r - 1)`$, so the error is $`r/(1-r) = 1.54`$ times the change. The loop stops when the relative change falls below $`10^{-12}`$, that is, when the change falls below $`21.73 \times 10^{-12}`$, so the error must be smaller than about $`1.54 \times 2.17 \times 10^{-11} = 3.35 \times 10^{-11}`$, and the actual value $`3.02 \times 10^{-11}`$ lies within this bound. The slower the convergence (the closer $`r`$ is to 1), the more a stopping criterion based on the change underestimates the actual error.

---

<br>

## 5. Problem 3: Inverse and Shifted Power Methods

### 5.1 Approach

The power method is applied to the inverse of $`B = A - \sigma I`$ without forming the inverse. Each iteration solves $`B\mathbf{y} = \mathbf{x}`$ instead of computing $`\mathbf{y} = B^{-1}\mathbf{x}`$, and since $`B`$ stays the same throughout, it is factored once with `sc.linalg.lu_factor(B)` before the loop, and each iteration performs only forward and back substitution with `sc.linalg.lu_solve(lu_obj, x)`. The factorization costs $`O(n^3)`$ and the substitution $`O(n^2)`$.

The normalization factor `mu` converges to the dominant eigenvalue $`1/(\lambda - \sigma)`$ of $`B^{-1}`$, so the original eigenvalue is recovered as follows.

$$
\mu = \frac{1}{\lambda - \sigma} \quad\Longrightarrow\quad \lambda = \sigma + \frac{1}{\mu}
$$

The inverse and shifted power methods are covered in Section 3 of the [Chapter 12 lab note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lab.md).

### 5.2 Code Written

```python
    lu_obj = sc.linalg.lu_factor(B)
```

```python
        y   = sc.linalg.lu_solve(lu_obj, x)
        mu  = y[np.argmax(np.abs(y))]                # largest element of y  (= 1 / (lambda - sigma))
        x   = y / mu
        lam = sigma + 1.0 / mu
```

### 5.3 Results and Verification

```text
All eigenvalues (numpy, ascending): [ 6.088518 13.182326 21.729156]

Method                                Found lambda   True lambda       Error   Iters
------------------------------------------------------------------------------------
Power Method   (dominant)              21.72915618   21.72915618    3.02e-11      51
sigma=0        (smallest)               6.08851756    6.08851756    2.68e-12      38
sigma=10       (middle)                13.18232626   13.18232626    5.06e-12     127
sigma=24       (largest)               21.72915618   21.72915618    1.34e-12      18
```

The output matches the output of the distributed notebook, and each shift found the eigenvalue closest to $`\sigma`$. The convergence rate in $`B^{-1}`$ explains the difference in iteration counts. The eigenvalues of $`B^{-1}`$ are $`1/(\lambda_i - \sigma)`$, so the convergence rate is the distance from $`\sigma`$ to the closest eigenvalue divided by the distance to the second closest one.

| Shift | Distance to the closest eigenvalue | Distance to the second closest eigenvalue | Convergence rate | Iterations |
|:-----|:-----|:-----|:-----|:-----|
| $`\sigma = 0`$ | 6.0885 | 13.1823 | 0.462 | 38 |
| $`\sigma = 10`$ | 3.1823 | 3.9115 | 0.814 | 127 |
| $`\sigma = 24`$ | 2.2708 | 10.8177 | 0.210 | 18 |

*Table 2. Convergence rate and iteration count for each shift*

The shift $`\sigma = 10`$ lies between the eigenvalues 6.09 and 13.18 at similar distances, so convergence is slow, while $`\sigma = 24`$ is close to 21.73 and far from 13.18, which makes it faster even than the power method (rate 0.607, 51 iterations). The closer the shift is to the target eigenvalue, the faster the convergence.

---

<br>

## 6. Problem 4: Natural Frequencies of a Vibrating System

### 6.1 Approach

The `omega_sq` returned by `sc.linalg.eigh(K_vib, M_vib)` is $`\lambda = \omega^2`$, so the angular frequency is its square root, and the frequency is the angular frequency divided by the angle of one cycle, $`2\pi`$.

$$
\omega_n = \sqrt{\lambda}, \qquad f = \frac{\omega_n}{2\pi}, \qquad T = \frac{1}{f}
$$

`omega_sq` is an array, so `np.sqrt` computes all three modes at once. Generalized eigenvalue problems and natural frequencies are covered in Section 2 of the [Chapter 13 lab note](../../Lecture-Notes/Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lab.md).

### 6.2 Code Written

```python
omega_n = np.sqrt(omega_sq)
freq_hz = omega_n / (2 * np.pi)
```

### 6.3 Results and Verification

```text
=== Natural Frequencies ===
 Mode         omega^2   omega (rad/s)      f (Hz)       T (s)
------------------------------------------------------------
    1          9.0470          3.0078      0.4787      2.0889
    2         32.3636          5.6889      0.9054      1.1045
    3         56.9227          7.5447      1.2008      0.8328

Mode shapes (columns = eigenvectors, M-normalized):
[[-0.329371  0.484607 -0.395817]
 [-0.569877  0.085174  0.578492]
 [-0.543958 -0.720715 -0.429743]]

Verification  max|K v - omega^2 M v|:
  Mode 1: 4.44e-15
  Mode 2: 7.11e-15
  Mode 3: 1.42e-14
```

The output matches the output of the distributed notebook. Checking the first mode, $`\sqrt{9.0470} = 3.0078`$ rad/s, $`3.0078 / (2\pi) = 0.4787`$ Hz, and $`1/0.4787 = 2.0889`$ s. The residuals are on the order of $`10^{-14}`$, and the mode shapes returned by `eigh` are normalized with respect to $`M`$, satisfying $`V^T M V = I`$. The eigenvalues are positive, so the square roots are real, because $`K`$ and $`M`$ are symmetric positive definite.

---

<br>

## 7. Problem 5: Modal Constants of an ODE System

### 7.1 Approach

At $`t = 0`$ the solution is $`\mathbf{y}(0) = \sum_i c_i \mathbf{v}_i = V\mathbf{c}`$, so the modal constants are the solution of the linear system $`V\mathbf{c} = \mathbf{y}_0`$. Solving it directly with `np.linalg.solve` instead of forming $`V^{-1}`$ and multiplying costs less and introduces less roundoff error. This is the same golden rule, solve instead of inverting, covered in Section 1.2 of the [Chapter 11 lab note](../../Lecture-Notes/Ch11_Matrix-Inverse-and-Condition/Concepts_Lab.md). The stability of ODE systems is covered in Section 3 of the [Chapter 13 lab note](../../Lecture-Notes/Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lab.md) and in the [Chapter 13 lecture note](../../Lecture-Notes/Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lecture.md).

### 7.2 Code Written

```python
c_ode = np.linalg.solve(V_ode, y0)
```

### 7.3 Results and Verification

```text
Eigenvalues:
  lambda_1 = -3.989933   -> mode decays
  lambda_2 = -302.010067   -> mode decays
System is STABLE

Modal constants: c1 = 53.6413,  c2 = 82.8879

Time constants:  tau_1 = 1/|lambda_1| = 0.2506 s  (slow mode)
                 tau_2 = 1/|lambda_2| = 0.003311 s  (fast mode)
Stiffness ratio  |lambda_2 / lambda_1| = 75.7
```

The output matches the output of the distributed notebook. The results are checked as follows.

- **Eigenvalues:** From the trace $`-306`$ and the determinant $`(-5)(-301) - 3 \cdot 100 = 1205`$, the characteristic equation is $`\lambda^2 + 306\lambda + 1205 = 0`$, giving $`\lambda = (-306 \pm \sqrt{88816})/2 = -3.9899, -302.0101`$. Both real parts are negative, so by Table 1 the system is asymptotically stable.
- **Modal constants:** With the eigenvectors $`\mathbf{v}_1 = (0.9477, 0.3191)`$ and $`\mathbf{v}_2 = (-0.0101, 0.99995)`$ returned by `np.linalg.eig`, computing $`V\mathbf{c}`$ gives $`0.9477 \times 53.6413 - 0.0101 \times 82.8879 = 50.0`$ and $`0.3191 \times 53.6413 + 0.99995 \times 82.8879 = 100.0`$, which recovers $`\mathbf{y}_0`$.
- **Time constants and stiffness ratio:** $`1/3.9899 = 0.2506`$ s and $`1/302.01 = 0.003311`$ s, and the ratio 75.7 between the two time constants shows a stiff system whose fast and slow modes have very different time scales.

> **Note:** An eigenvector is free in its length and sign, so the values of the modal constants $`c_i`$ depend on how the eigenvectors are normalized. `np.linalg.eig` returns eigenvectors of 2-norm 1, and another normalization changes $`c_i`$ but leaves the products $`c_i \mathbf{v}_i`$ and the solution $`\mathbf{y}(t)`$ unchanged.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Computing eigenpairs | `np.linalg.eig` returns the eigenvalue array and the matrix whose columns are the eigenvectors, and the residuals and the trace and determinant identities verify the result. |
| Power method | The component with the largest magnitude, with its sign, serves as the eigenvalue estimate and as the normalization divisor, and the error shrinks by about $`\vert\lambda_2/\lambda_1\vert`$ per iteration. |
| Stopping criterion and actual error | When the error shrinks at rate $`r`$, the actual error is about $`r/(1-r)`$ times the change between successive estimates. |
| Shifted power method | It LU-factors $`A - \sigma I`$ once and performs only substitution in each iteration, recovering the eigenvalue closest to $`\sigma`$ as $`\lambda = \sigma + 1/\mu`$. |
| Shift and convergence rate | The convergence rate is the ratio of the distances from $`\sigma`$ to the closest and the second closest eigenvalues, so a shift closer to the target converges faster. |
| Natural frequencies | From the generalized eigenvalue $`\lambda = \omega^2`$ one obtains $`\omega_n = \sqrt{\lambda}`$, $`f = \omega_n/(2\pi)`$, and $`T = 1/f`$. |
| Modal constants and stability | The modal constants come from solving $`V\mathbf{c} = \mathbf{y}_0`$, and the system is asymptotically stable when all eigenvalues have negative real parts. |

---

<br>

## Self-Check Questions

1. **Signed eigenvalue estimate:** In Problem 2, what happens for a matrix with a negative dominant eigenvalue if the code uses `lam = np.max(np.abs(y))`?

   <details>
   <summary>Answer</summary>

   The estimate is always positive, so the negative eigenvalue can never be obtained. Moreover, the sign of $`\mathbf{y} = A\mathbf{x}`$ flips at every iteration while the division is always by a positive number, so $`\mathbf{x}`$ also alternates in sign and does not converge. Only division by the signed component `y[np.argmax(np.abs(y))]` keeps $`\mathbf{x}`$ fixed and lets the estimate converge to the negative eigenvalue.

   </details>

2. **Placement of the LU factorization:** In Problem 3, calling `lu_factor` inside the loop gives the same result. Why place it outside the loop anyway?

   <details>
   <summary>Answer</summary>

   $`B = A - \sigma I`$ does not change during the iterations, so one factorization suffices. The factorization costs $`O(n^3)`$ and the substitution using it costs $`O(n^2)`$, so factoring once outside the loop costs only $`O(n^2)`$ per iteration, while factoring inside the loop costs $`O(n^3)`$ per iteration.

   </details>

3. **Choosing the shift:** For the matrix of Problem 3, how should the shift be chosen to find the middle eigenvalue 13.18 faster than with $`\sigma = 10`$? Compute the convergence rate for $`\sigma = 13`$.

   <details>
   <summary>Answer</summary>

   Place the shift closer to 13.18. With $`\sigma = 13`$, the distance to the closest eigenvalue 13.18 is 0.18 and the distance to the second closest, 6.09, is 6.91, so the convergence rate is $`0.18/6.91 \approx 0.026`$, far smaller than the 0.814 of $`\sigma = 10`$, and the method converges in a few iterations.

   </details>

4. **Computing the period:** For the second mode of Problem 4 with $`\omega^2 = 32.3636`$, compute the frequency and the period by hand.

   <details>
   <summary>Answer</summary>

   $`\omega = \sqrt{32.3636} = 5.6889`$ rad/s, $`f = 5.6889/(2\pi) = 0.9054`$ Hz, and $`T = 1/0.9054 = 1.1045`$ s.

   </details>

5. **Stability check:** If $`a_{22} = -301`$ in the matrix of Problem 5 is changed to $`+301`$, is the system stable?

   <details>
   <summary>Answer</summary>

   The determinant becomes $`(-5)(301) - 3 \cdot 100 = -1805 < 0`$, so the product of the two eigenvalues is negative. One eigenvalue is therefore positive, and its mode grows exponentially, so the system is unstable.

   </details>
