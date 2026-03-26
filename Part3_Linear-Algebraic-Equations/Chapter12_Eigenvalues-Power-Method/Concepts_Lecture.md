# Chapter 12 Lecture — Iterative Methods for Linear Systems

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. Introduction to Iterative Methods](#1-introduction-to-iterative-methods)
- [2. Gauss-Seidel Method](#2-gauss-seidel-method)
  - [2.1 Problem Setup](#21-problem-setup)
  - [2.2 Derivation for a 3x3 System](#22-derivation-for-a-3x3-system)
  - [2.3 Iteration Scheme with Superscripts](#23-iteration-scheme-with-superscripts)
  - [2.4 Convergence Criterion](#24-convergence-criterion)
  - [2.5 Example 12.1: Gauss-Seidel Iteration](#25-example-121-gauss-seidel-iteration)
- [3. Convergence and Diagonal Dominance](#3-convergence-and-diagonal-dominance)
  - [3.1 Diagonal Dominance Condition](#31-diagonal-dominance-condition)
  - [3.2 Example: 3x3 Diagonal Dominance Check](#32-example-3x3-diagonal-dominance-check)
  - [3.3 Sufficient vs. Necessary Condition](#33-sufficient-vs-necessary-condition)
- [4. Gauss-Seidel Method in Matrix Form](#4-gauss-seidel-method-in-matrix-form)
  - [4.1 Matrix Splitting: L, D, U Decomposition](#41-matrix-splitting-l-d-u-decomposition)
  - [4.2 Matrix Form of the Iteration](#42-matrix-form-of-the-iteration)
- [5. Jacobi Iteration](#5-jacobi-iteration)
  - [5.1 Comparison with Gauss-Seidel](#51-comparison-with-gauss-seidel)
  - [5.2 Component Form](#52-component-form)
  - [5.3 Matrix Form](#53-matrix-form)
- [6. Relaxation Method](#6-relaxation-method)
  - [6.1 Relaxation Formula](#61-relaxation-formula)
  - [6.2 Types of Relaxation](#62-types-of-relaxation)
  - [6.3 Relaxation Applied to 3x3 System](#63-relaxation-applied-to-3x3-system)
  - [6.4 Matrix Form of Relaxation (SOR)](#64-matrix-form-of-relaxation-sor)
  - [6.5 Example 12.2: Relaxation with lambda = 1.2](#65-example-122-relaxation-with-lambda--12)
- [7. Preconditioning](#7-preconditioning)
  - [7.1 General Framework](#71-general-framework)
  - [7.2 Preconditioners for Each Method](#72-preconditioners-for-each-method)
  - [7.3 Error Analysis and Convergence](#73-error-analysis-and-convergence)
  - [7.4 Spectral Radius](#74-spectral-radius)
  - [7.5 Eigenvalue Relationship for I - A](#75-eigenvalue-relationship-for-i---a)
- [8. Python Implementations](#8-python-implementations)
- [Summary](#summary)

---

<br>

## 1. Introduction to Iterative Methods

Iterative methods solve a linear system $[A]\{x\} = \{b\}$ by generating a sequence of approximations:

$$\mathbf{x}^{k+1} := \Psi(\mathbf{x}^k), \quad k \geq 0$$

Starting from an initial guess $\mathbf{x}^0$, the method repeatedly applies a mapping $\Psi$ until convergence to the true solution. Unlike direct methods (Gauss elimination, LU decomposition) that produce an exact answer in a finite number of operations, iterative methods approach the solution gradually and can be stopped when the approximation is sufficiently accurate.

### Learning Objectives

1. Understand the **Gauss-Seidel** and **Jacobi** methods
2. Know how to assess **diagonal dominance**
3. Understand how **relaxation** improves the convergence of iterative methods

> **[Linear Algebra]** Direct methods cost $O(n^3)$ and work well for small-to-moderate systems. For very large, sparse systems (common in engineering and scientific computing), iterative methods are often far more efficient because they only require matrix-vector products and can exploit sparsity.

---

<br>

## 2. Gauss-Seidel Method

### 2.1 Problem Setup

Suppose we are given $n$ equations in the standard form:

$$[A]_{n \times n} \{x\}_{n \times 1} = \{b\}_{n \times 1}$$

where $[A]$ is the coefficient matrix, $\{x\}$ is the unknown vector, and $\{b\}$ is the right-hand side vector.

### 2.2 Derivation for a 3x3 System

For simplicity, consider a $3 \times 3$ system:

$$\begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \\ x_3 \end{Bmatrix} = \begin{Bmatrix} b_1 \\ b_2 \\ b_3 \end{Bmatrix}$$

Expanding the rows:

$$a_{11} x_1 = b_1 - a_{12} x_2 - a_{13} x_3$$

$$a_{22} x_2 = b_2 - a_{21} x_1 - a_{23} x_3$$

$$a_{33} x_3 = b_3 - a_{31} x_1 - a_{32} x_2$$

Solving each equation for the diagonal unknown:

$$x_1 = \frac{b_1 - a_{12} x_2 - a_{13} x_3}{a_{11}}$$

$$x_2 = \frac{b_2 - a_{21} x_1 - a_{23} x_3}{a_{22}}$$

$$x_3 = \frac{b_3 - a_{31} x_1 - a_{32} x_2}{a_{33}}$$

### 2.3 Iteration Scheme with Superscripts

The Gauss-Seidel method introduces iteration indices. At iteration $j$, we compute:

$$x_1^{j} = \frac{b_1 - a_{12} x_2^{j-1} - a_{13} x_3^{j-1}}{a_{11}}$$

$$x_2^{j} = \frac{b_2 - a_{21} x_1^{j} - a_{23} x_3^{j-1}}{a_{22}}$$

$$x_3^{j} = \frac{b_3 - a_{31} x_1^{j} - a_{32} x_2^{j}}{a_{33}}$$

The key feature of Gauss-Seidel is that it **uses the most recently computed values** as soon as they are available. When computing $x_2^{j}$, we already use $x_1^{j}$ (just computed in the current iteration), not $x_1^{j-1}$. Similarly, when computing $x_3^{j}$, we use both $x_1^{j}$ and $x_2^{j}$.

> **[Linear Algebra]** This "sequential update" property is what distinguishes Gauss-Seidel from the Jacobi method. It typically results in faster convergence because new information is incorporated immediately, but it also means the computation is inherently sequential -- $x_2^j$ cannot be computed until $x_1^j$ is known.

### 2.4 Convergence Criterion

We repeat the procedure until our solution converges to the true solution. The **approximate relative error** for each variable is:

$$\varepsilon_{a,i} = \left| \frac{x_i^{j} - x_i^{j-1}}{x_i^{j}} \right| \times 100\%$$

The iteration stops when $\varepsilon_{a,i} \leq \varepsilon_s$ for **all** variables, where $\varepsilon_s$ is the user-specified tolerance.

### 2.5 Example 12.1: Gauss-Seidel Iteration

**Given system:**

$$3x_1 - 0.1x_2 - 0.2x_3 = 7.85$$

$$0.1x_1 + 7x_2 - 0.3x_3 = -19.3$$

$$0.3x_1 - 0.2x_2 + 10x_3 = 71.4$$

**Step 1 -- Rearrange** each equation for its unknown on the diagonal:

$$x_1 = \frac{7.85 + 0.1x_2 + 0.2x_3}{3}$$

$$x_2 = \frac{-19.3 - 0.1x_1 + 0.3x_3}{7}$$

$$x_3 = \frac{71.4 - 0.3x_1 + 0.2x_2}{10}$$

**Step 2 -- First iteration** ($j = 1$) with initial guess $\{x\}^0 = \{0, 0, 0\}$:

$$x_1^1 = \frac{7.85 + 0.1(0) + 0.2(0)}{3} = 2.616667$$

$$x_2^1 = \frac{-19.3 - 0.1(2.616667) + 0.3(0)}{7} = -2.794524$$

$$x_3^1 = \frac{71.4 - 0.3(2.616667) + 0.2(-2.794524)}{10} = 7.005610$$

Note how $x_1^1$ was used immediately when computing $x_2^1$, and both $x_1^1$ and $x_2^1$ were used when computing $x_3^1$.

**Step 3 -- Second iteration** ($j = 2$):

$$x_1^2 = \frac{7.85 + 0.1(-2.794524) + 0.2(7.005610)}{3} = 2.990557$$

$$x_2^2 = \frac{-19.3 - 0.1(2.990557) + 0.3(7.005610)}{7} = -2.499625$$

$$x_3^2 = \frac{71.4 - 0.3(2.990557) + 0.2(-2.499625)}{10} = 7.000291$$

**Step 4 -- Compute relative errors** (for iteration 2):

$$\varepsilon_{a,1}^2 = \left| \frac{2.990557 - 2.616667}{2.990557} \right| \times 100\% = 12.5\%$$

The relative errors of $x_2$ and $x_3$ are $11.8\%$ and $0.0766\%$, respectively.

We repeat this procedure until **all** relative errors are smaller than $\varepsilon_s$, the user-specified tolerance.

---

<br>

## 3. Convergence and Diagonal Dominance

### 3.1 Diagonal Dominance Condition

The Gauss-Seidel method **may diverge** for some systems. However, if the system is **diagonally dominant**, it will **definitely converge**.

A matrix is diagonally dominant if, for every row $i$:

$$|a_{ii}| > \sum_{\substack{j=1 \\ j \neq i}}^{n} |a_{ij}|$$

In words: the absolute value of each diagonal element is **strictly greater than** the sum of the absolute values of all other elements in that row.

### 3.2 Example: 3x3 Diagonal Dominance Check

For a $3 \times 3$ system:

$$\begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \\ x_3 \end{Bmatrix} = \begin{Bmatrix} b_1 \\ b_2 \\ b_3 \end{Bmatrix}$$

The three conditions to verify are:

$$|a_{11}| > |a_{12}| + |a_{13}|$$

$$|a_{22}| > |a_{21}| + |a_{23}|$$

$$|a_{33}| > |a_{31}| + |a_{32}|$$

For Example 12.1: $|3| > |-0.1| + |-0.2|$, $|7| > |0.1| + |-0.3|$, $|10| > |0.3| + |-0.2|$. All conditions hold, so convergence is guaranteed.

### 3.3 Sufficient vs. Necessary Condition

Diagonal dominance is a **sufficient condition** but **not a necessary condition** for convergence.

- If the matrix is diagonally dominant, convergence is guaranteed
- If the matrix is **not** diagonally dominant, the method **may still converge** (it just isn't guaranteed)
- The Gauss-Seidel method also works for **symmetric positive definite** matrices

> **[Linear Algebra]** If the system is not diagonally dominant in its original form, it may be possible to **swap rows** (pivot) to achieve diagonal dominance. This is the strategy used in Example 12.2.

---

<br>

## 4. Gauss-Seidel Method in Matrix Form

### 4.1 Matrix Splitting: L, D, U Decomposition

The coefficient matrix $[A]$ can be decomposed as:

$$[A] = \tilde{L} + \tilde{D} + \tilde{U}$$

where:

- $\tilde{L}$ is the **strictly lower triangular** part of $[A]$ (below diagonal, zeros on and above diagonal)
- $\tilde{D}$ is the **diagonal** matrix of $[A]$ (diagonal elements only)
- $\tilde{U}$ is the **strictly upper triangular** part of $[A]$ (above diagonal, zeros on and below diagonal)

For the $3 \times 3$ case:

$$\tilde{L} = \begin{bmatrix} 0 & 0 & 0 \\ a_{21} & 0 & 0 \\ a_{31} & a_{32} & 0 \end{bmatrix}, \quad \tilde{D} = \begin{bmatrix} a_{11} & 0 & 0 \\ 0 & a_{22} & 0 \\ 0 & 0 & a_{33} \end{bmatrix}, \quad \tilde{U} = \begin{bmatrix} 0 & a_{12} & a_{13} \\ 0 & 0 & a_{23} \\ 0 & 0 & 0 \end{bmatrix}$$

> **[Linear Algebra]** This $\tilde{L} + \tilde{D} + \tilde{U}$ splitting is different from the LU factorization in Chapter 10. Here, $\tilde{L}$, $\tilde{D}$, $\tilde{U}$ are simply the parts of the original matrix $A$ -- no elimination is performed. The notation uses tildes to distinguish from the LU factors.

### 4.2 Matrix Form of the Iteration

Rearranging the Gauss-Seidel equations into matrix form:

$$[\tilde{L} + \tilde{D}]\{x\}^{j} = \{b\} - [\tilde{U}]\{x\}^{j-1}$$

Let $\tilde{L}_* = \tilde{L} + \tilde{D}$ (the lower triangular part of $A$ including the diagonal). The matrix-vector form is:

$$\{x\}^{j} = (\tilde{L}_*)^{-1}\left(\{b\} - [\tilde{U}]\{x\}^{j-1}\right)$$

The left-hand side involves a **lower triangular** system, which can be solved efficiently by **forward substitution** -- no explicit matrix inversion is needed.

---

<br>

## 5. Jacobi Iteration

### 5.1 Comparison with Gauss-Seidel

The **Jacobi iteration** is similar to the Gauss-Seidel method, except all variables in the $j$th iteration are updated using **only** the $(j-1)$th iteration values. No newly computed values are used within the same iteration.

| Feature | Gauss-Seidel | Jacobi |
|---|---|---|
| Update strategy | Uses latest values immediately | Uses only previous iteration values |
| Parallelism | Sequential (depends on update order) | **Naturally parallel** (all updates independent) |
| Convergence speed | Generally faster | Generally slower |
| Storage | One copy of solution vector | Two copies needed (old and new) |

The Jacobi method converges to the true solution if $[A]$ is diagonally dominant.

### 5.2 Component Form

For a $3 \times 3$ system, the Jacobi iteration computes:

$$x_1^{j} = \frac{b_1 - a_{12} x_2^{j-1} - a_{13} x_3^{j-1}}{a_{11}}$$

$$x_2^{j} = \frac{b_2 - a_{21} x_1^{j-1} - a_{23} x_3^{j-1}}{a_{22}}$$

$$x_3^{j} = \frac{b_3 - a_{31} x_1^{j-1} - a_{32} x_2^{j-1}}{a_{33}}$$

Note that **all** right-hand side values use the superscript $j-1$. This means each $x_i^j$ can be computed **independently**, making Jacobi ideal for parallel computation.

### 5.3 Matrix Form

$$[\tilde{D}]\{x\}^{j} = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}$$

$$\{x\}^{j} = [\tilde{D}]^{-1}\left(\{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}\right)$$

Since $\tilde{D}$ is diagonal, its inverse is trivial: $[\tilde{D}]^{-1}_{ii} = 1/a_{ii}$.

---

<br>

## 6. Relaxation Method

### 6.1 Relaxation Formula

Relaxation is a slight modification of the Gauss-Seidel method designed to **enhance convergence**. The idea is to combine the old value and the newly calculated Gauss-Seidel value using a weighting factor:

$$x_i^{\text{new}} = \lambda \, x_i^{\text{new(GS)}} + (1 - \lambda) \, x_i^{\text{old}}$$

where $\lambda$ is a weighting factor assigned a value between $0$ and $2$, and $x_i^{\text{new(GS)}}$ is the value that the standard Gauss-Seidel method would produce.

### 6.2 Types of Relaxation

| $\lambda$ Range | Name | Effect |
|---|---|---|
| $\lambda = 1$ | No modification | Standard Gauss-Seidel method |
| $0 < \lambda < 1$ | **Under-relaxation** | Dampens oscillations; may help non-convergent systems |
| $1 < \lambda \leq 2$ | **Over-relaxation** (SOR) | Accelerates convergence when Gauss-Seidel is slow |

> **[Linear Algebra]** The optimal value of $\lambda$ depends on the specific problem and is generally not known in advance. For the SOR (Successive Over-Relaxation) method, theoretical results exist for special matrix structures (e.g., tridiagonal systems), but in practice $\lambda$ is often determined experimentally.

### 6.3 Relaxation Applied to 3x3 System

For a $3 \times 3$ system, the relaxation update is:

$$x_1^{\text{new}} = \frac{b_1 - a_{12} x_2^{\text{old}} - a_{13} x_3^{\text{old}}}{a_{11}} \cdot \lambda + (1 - \lambda) \, x_1^{\text{old}}$$

$$x_2^{\text{new}} = \frac{b_2 - a_{21} x_1^{\text{new}} - a_{23} x_3^{\text{old}}}{a_{22}} \cdot \lambda + (1 - \lambda) \, x_2^{\text{old}}$$

$$x_3^{\text{new}} = \frac{b_3 - a_{31} x_1^{\text{new}} - a_{32} x_2^{\text{new}}}{a_{33}} \cdot \lambda + (1 - \lambda) \, x_3^{\text{old}}$$

Notice the Gauss-Seidel pattern is preserved: $x_1^{\text{new}}$ is used in the computation of $x_2^{\text{new}}$, and both are used for $x_3^{\text{new}}$.

### 6.4 Matrix Form of Relaxation (SOR)

Multiplying through and collecting terms, the relaxation method in matrix form is:

$$[\tilde{D} + \lambda \tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - \left(\lambda[\tilde{U}] + (\lambda - 1)[\tilde{D}]\right)\{x\}^{\text{old}}$$

**Derivation outline:**

Starting from the component form and multiplying by $a_{ii}$:

$$a_{ii} x_i^{\text{new}} = \left(b_i - \sum_{k<i} a_{ik} x_k^{\text{new}} - \sum_{k>i} a_{ik} x_k^{\text{old}}\right) \lambda + (1-\lambda) \, a_{ii} \, x_i^{\text{old}}$$

Expanding and rearranging into matrix notation yields:

$$a_{ii} x_i^{\text{new}} = \lambda b_i - \lambda \sum_{k>i} a_{ik} x_k^{\text{old}} + (1-\lambda) a_{ii} x_i^{\text{old}} - \lambda \sum_{k<i} a_{ik} x_k^{\text{new}}$$

Moving terms with "new" to the left and "old" to the right gives the final matrix form.

### 6.5 Example 12.2: Relaxation with lambda = 1.2

**Given:** $\lambda = 1.2$, $\varepsilon_s = 10\%$

**Original system:**

$$-3x_1 + 12x_2 = 9$$

$$10x_1 - 2x_2 = 8$$

**Step 1 -- Swap rows** for diagonal dominance:

$$10x_1 - 2x_2 = 8$$

$$-3x_1 + 12x_2 = 9$$

After swapping: $|10| > |-2|$ and $|12| > |-3|$. The system is now diagonally dominant.

Rearranged:

$$x_1^{j} = \frac{8 + 2x_2^{j-1}}{10}, \quad x_2^{j} = \frac{9 + 3x_1^{j}}{12}$$

**Iteration 1** (initial guess $x_1^0 = 0$, $x_2^0 = 0$):

Gauss-Seidel step:

$$x_1^1 = \frac{8 + 2(0)}{10} = 0.8$$

Apply relaxation:

$$x_1^{\text{new}} = \lambda \cdot x_1^1 + (1-\lambda) \cdot x_1^0 = 1.2(0.8) + (-0.2)(0) = 0.96$$

Gauss-Seidel step (using relaxed $x_1^{\text{new}} = 0.96$):

$$x_2^1 = \frac{9 + 3(0.96)}{12} = 0.99$$

Apply relaxation:

$$x_2^{\text{new}} = 1.2(0.99) + (-0.2)(0) = 1.188$$

**Iteration 2** (using $x_1^1 = 0.96$, $x_2^1 = 1.188$):

Gauss-Seidel step:

$$x_1^2 = \frac{8 + 2(1.188)}{10} = 1.0376$$

Apply relaxation:

$$x_1^{\text{new}} = 1.2(1.0376) + (-0.2)(0.96) = 1.05312$$

Gauss-Seidel step:

$$x_2^2 = \frac{9 + 3(1.05312)}{12} = 1.01328$$

Apply relaxation:

$$x_2^{\text{new}} = 1.2(1.01328) + (-0.2)(1.188) = 0.978336$$

Relative errors:

$$\varepsilon_{a,1} = \left|\frac{1.05312 - 0.96}{1.05312}\right| \times 100\% = 8.84\%$$

$$\varepsilon_{a,2} = \left|\frac{0.978336 - 1.188}{0.978336}\right| \times 100\% = 21.43\%$$

Since $\varepsilon_{a,2} > \varepsilon_s = 10\%$, we continue iterating.

**Iteration 3:**

$$x_1^3 = 0.984177, \quad \varepsilon_{a,1} = 7.01\%$$

$$x_2^3 = 0.999586, \quad \varepsilon_{a,2} = 2.13\%$$

Both errors are below $10\%$, so the iteration converges. The true solution is $x_1 = 1$, $x_2 = 1$.

---

<br>

## 7. Preconditioning

### 7.1 General Framework

The three iterative methods (Gauss-Seidel, Jacobi, Relaxation) are all **traditional iteration methods** to solve a linear system $A\mathbf{x} = \mathbf{b}$, labeled as equation (1).

We introduce a **preconditioner** $P$, which is an approximation to $A$. Adding $P\mathbf{x}$ to both sides of $A\mathbf{x} = \mathbf{b}$:

$$P\mathbf{x} + A\mathbf{x} = P\mathbf{x} + \mathbf{b}$$

$$P\mathbf{x} = (P - A)\mathbf{x} + \mathbf{b} \quad \cdots (2)$$

Iteratively updating the solution $\mathbf{x}^{k+1}$:

$$P\mathbf{x}^{k+1} = (P - A)\mathbf{x}^k + \mathbf{b} \quad \cdots (3)$$

The hope is that equation (3) is **easier to solve** than equation (1), because $P$ is chosen to be a matrix whose system is cheap to solve (e.g., triangular or diagonal).

### 7.2 Preconditioners for Each Method

Each iterative method corresponds to a specific choice of preconditioner $P$:

| Method | Preconditioner $P$ | $P - A$ |
|---|---|---|
| **Gauss-Seidel** | $P = \tilde{L} + \tilde{D}$ | $P - A = -\tilde{U}$ |
| **Jacobi** | $P = \tilde{D}$ | $P - A = -\tilde{L} - \tilde{U}$ |
| **Relaxation (SOR)** | $P = \tilde{D} + \lambda\tilde{L}$ | $P - A = (\lambda - 1)\tilde{L} - \tilde{U}$ |

The resulting iteration equations:

- **Gauss-Seidel:** $[\tilde{L} + \tilde{D}]\{x\}^j = \{b\} - [\tilde{U}]\{x\}^{j-1}$
- **Jacobi:** $[\tilde{D}]\{x\}^j = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}$
- **Relaxation:** $[\tilde{D} + \lambda\tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - (\lambda[\tilde{U}] + (\lambda-1)[\tilde{D}])\{x\}^{\text{old}}$

### 7.3 Error Analysis and Convergence

Define the error at iteration $k$ as:

$$\mathbf{e}^k = \mathbf{x} - \mathbf{x}^k$$

where $\mathbf{x}$ is the true solution. Subtracting equation (2) from equation (3):

$$P(\mathbf{x} - \mathbf{x}^{k+1}) = (P - A)(\mathbf{x} - \mathbf{x}^k)$$

$$P\mathbf{e}^{k+1} = (P - A)\mathbf{e}^k$$

$$\mathbf{e}^{k+1} = (I - P^{-1}A)\mathbf{e}^k$$

Define the **iteration matrix**:

$$M = I - P^{-1}A$$

Then each step multiplies the error vector by $M$:

$$\mathbf{e}^{k+1} = M\mathbf{e}^k$$

The convergence is governed by the **eigenvalues of $M$**. For convergence, every eigenvalue of $M$ must satisfy $|\lambda(M)| < 1$.

### 7.4 Spectral Radius

The largest eigenvalue (in absolute value) of $M$ is called the **spectral radius**:

$$\rho(M) = \max(|\lambda(M)|)$$

**Convergence requires that:**

$$\rho(M) < 1$$

When the error vector happens to be an eigenvector of $M$, the error at the next step becomes:

$$\mathbf{e}^{k+1} = M\mathbf{e}^k = \lambda \mathbf{e}^k$$

$$\Rightarrow \mathbf{e}^k = \lambda^k \mathbf{e}^0$$

So the error decays geometrically as $\lambda^k$. The smaller $\rho(M)$, the faster the convergence.

> **[Linear Algebra]** This is the fundamental reason why diagonal dominance ensures convergence. For a diagonally dominant matrix, the iteration matrix $M$ of both Gauss-Seidel and Jacobi methods has spectral radius less than 1. The choice of $\lambda$ in SOR aims to minimize $\rho(M)$.

### 7.5 Eigenvalue Relationship for I - A

If we don't use a preconditioner ($P = I$), then $M = I - A$.

**Theorem:** If $[A]$ has eigenvalues $\lambda_1, \lambda_2, \ldots, \lambda_n$, then $[I - A]$ has eigenvalues $1 - \lambda_1, 1 - \lambda_2, \ldots, 1 - \lambda_n$.

**Proof:**

Given $A\mathbf{v} = \lambda\mathbf{v}$ where $\mathbf{v}$ is an eigenvector:

$$(I - A)\mathbf{v} = \mathbf{v} - A\mathbf{v} = \mathbf{v} - \lambda\mathbf{v} = (1 - \lambda)\mathbf{v}$$

Therefore $(1 - \lambda)$ is an eigenvalue of $[I - A]$ with the same eigenvector $\mathbf{v}$.

**Consequence:** Without preconditioning, all eigenvalues of $A$ must be inside a unit circle centered at $1$ in the complex plane for the iteration to converge.

---

<br>

## 8. Python Implementations

### Gauss-Seidel Method

```python
import numpy as np

def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=100):
    """
    Gauss-Seidel iterative method for solving Ax = b.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Coefficient matrix.
    b : array_like, shape (n,)
        Right-hand side vector.
    x0 : array_like, shape (n,), optional
        Initial guess (default: zeros).
    tol : float
        Stopping tolerance for approximate relative error (%).
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    x : ndarray
        Solution vector.
    iterations : int
        Number of iterations performed.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)

    for iteration in range(1, max_iter + 1):
        x_old = x.copy()

        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (b[i] - sigma) / A[i][i]

        # Check convergence (max approximate relative error)
        errors = [abs((x[i] - x_old[i]) / x[i]) * 100
                  for i in range(n) if x[i] != 0]
        if max(errors) < tol:
            return x, iteration

    return x, max_iter


# Example 12.1
A = np.array([[3, -0.1, -0.2],
              [0.1, 7, -0.3],
              [0.3, -0.2, 10]], dtype=float)
b = np.array([7.85, -19.3, 71.4], dtype=float)

x, iters = gauss_seidel(A, b, tol=0.1)
print(f"Solution: {x}")
print(f"Iterations: {iters}")
```

### Jacobi Method

```python
def jacobi(A, b, x0=None, tol=1e-6, max_iter=100):
    """
    Jacobi iterative method for solving Ax = b.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Coefficient matrix.
    b : array_like, shape (n,)
        Right-hand side vector.
    x0 : array_like, shape (n,), optional
        Initial guess (default: zeros).
    tol : float
        Stopping tolerance for approximate relative error (%).
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    x : ndarray
        Solution vector.
    iterations : int
        Number of iterations performed.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)

    for iteration in range(1, max_iter + 1):
        x_new = np.zeros(n)

        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - sigma) / A[i][i]

        # Check convergence
        errors = [abs((x_new[i] - x[i]) / x_new[i]) * 100
                  for i in range(n) if x_new[i] != 0]
        if max(errors) < tol:
            return x_new, iteration

        x = x_new

    return x, max_iter
```

### Relaxation Method (SOR)

```python
def sor(A, b, lam=1.2, x0=None, tol=1e-6, max_iter=100):
    """
    Successive Over-Relaxation (SOR) method for solving Ax = b.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Coefficient matrix.
    b : array_like, shape (n,)
        Right-hand side vector.
    lam : float
        Relaxation factor (0 < lam < 2).
        lam = 1 => Gauss-Seidel,
        lam < 1 => under-relaxation,
        lam > 1 => over-relaxation.
    x0 : array_like, shape (n,), optional
        Initial guess (default: zeros).
    tol : float
        Stopping tolerance for approximate relative error (%).
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    x : ndarray
        Solution vector.
    iterations : int
        Number of iterations performed.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)

    for iteration in range(1, max_iter + 1):
        x_old = x.copy()

        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_gs = (b[i] - sigma) / A[i][i]          # Gauss-Seidel value
            x[i] = lam * x_gs + (1 - lam) * x_old[i]  # Relaxation

        # Check convergence
        errors = [abs((x[i] - x_old[i]) / x[i]) * 100
                  for i in range(n) if x[i] != 0]
        if max(errors) < tol:
            return x, iteration

    return x, max_iter


# Example 12.2
A2 = np.array([[10, -2],
               [-3, 12]], dtype=float)
b2 = np.array([8, 9], dtype=float)

x2, iters2 = sor(A2, b2, lam=1.2, tol=10)
print(f"Solution: {x2}")
print(f"Iterations: {iters2}")
```

### Diagonal Dominance Check

```python
def is_diagonally_dominant(A):
    """
    Check if a matrix is strictly diagonally dominant.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Matrix to check.

    Returns
    -------
    bool
        True if the matrix is strictly diagonally dominant.
    """
    n = len(A)
    for i in range(n):
        diag = abs(A[i][i])
        off_diag_sum = sum(abs(A[i][j]) for j in range(n) if j != i)
        if diag <= off_diag_sum:
            return False
    return True


# Check Example 12.1
print(is_diagonally_dominant(A))   # True
```

---

<br>

## Summary

| Topic | Key Formula / Concept |
|---|---|
| **General iteration** | $\mathbf{x}^{k+1} := \Psi(\mathbf{x}^k)$, $k \geq 0$ |
| **Gauss-Seidel** | $x_i^j = \frac{1}{a_{ii}}\left(b_i - \sum_{k<i} a_{ik} x_k^j - \sum_{k>i} a_{ik} x_k^{j-1}\right)$; uses latest values |
| **Gauss-Seidel (matrix)** | $[\tilde{L} + \tilde{D}]\{x\}^j = \{b\} - [\tilde{U}]\{x\}^{j-1}$ |
| **Jacobi** | $x_i^j = \frac{1}{a_{ii}}\left(b_i - \sum_{k \neq i} a_{ik} x_k^{j-1}\right)$; all from previous iteration |
| **Jacobi (matrix)** | $[\tilde{D}]\{x\}^j = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}$ |
| **Relaxation (SOR)** | $x_i^{\text{new}} = \lambda \, x_i^{\text{GS}} + (1-\lambda) \, x_i^{\text{old}}$ |
| **Relaxation (matrix)** | $[\tilde{D} + \lambda\tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - (\lambda[\tilde{U}] + (\lambda-1)[\tilde{D}])\{x\}^{\text{old}}$ |
| **Diagonal dominance** | $\|a_{ii}\| > \sum_{j \neq i} \|a_{ij}\|$ (sufficient, not necessary) |
| **Convergence criterion** | $\varepsilon_{a,i} = \left\|\frac{x_i^j - x_i^{j-1}}{x_i^j}\right\| \times 100\% \leq \varepsilon_s$ for all $i$ |
| **Preconditioning** | $P\mathbf{x}^{k+1} = (P-A)\mathbf{x}^k + \mathbf{b}$; iteration matrix $M = I - P^{-1}A$ |
| **Spectral radius** | $\rho(M) = \max(\|\lambda(M)\|) < 1$ required for convergence |
| **Eigenvalue shift** | If $A$ has $\lambda_i$, then $I-A$ has $1-\lambda_i$ |
