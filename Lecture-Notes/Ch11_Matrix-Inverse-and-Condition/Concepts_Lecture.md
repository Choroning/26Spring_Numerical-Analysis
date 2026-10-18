# Chapter 11 — Matrix Inverse and Condition Number

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 11

> **Learning Objectives:**
> 1. Compute matrix inverses numerically
> 2. Calculate and interpret condition numbers
> 3. Assess solution sensitivity to input perturbations

---

<br>

## Table of Contents

- [1. Matrix Inverse Definition](#1-matrix-inverse-definition)
- [2. Computing the Inverse via LU Decomposition](#2-computing-the-inverse-via-lu-decomposition)
  - [2.1 The Column-by-Column Method](#21-the-column-by-column-method)
  - [2.2 Example 11.1: 3x3 Inverse](#22-example-111-3x3-inverse)
- [3. Physical Interpretation: Mass-Spring System](#3-physical-interpretation-mass-spring-system)
- [4. Error Analysis: Three Checks for Ill-Conditioning](#4-error-analysis-three-checks-for-ill-conditioning)
- [5. Vector Norms](#5-vector-norms)
- [6. Matrix Norms](#6-matrix-norms)
- [7. Matrix Condition Number](#7-matrix-condition-number)
- [8. Error Bound Theorem](#8-error-bound-theorem)
- [9. Precision Rule](#9-precision-rule)
- [10. Example 11.3: Hilbert Matrix](#10-example-113-hilbert-matrix)
- [11. Python Implementation](#11-python-implementation)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Matrix Inverse Definition

The inverse of a square matrix $`[A]`$ is defined as the unique matrix $`[A]^{-1}`$ satisfying:

```math
[A][A]^{-1} = [A]^{-1}[A] = [I]
```

where $`[I]`$ is the identity matrix. The inverse exists if and only if:

- $`[A]`$ is **square** ($`n \times n`$).
- $`[A]`$ is **non-singular**, that is, $`\det(A) \neq 0`$.

If $`[A]`$ is singular, that is, its determinant is zero, the matrix has no inverse. Geometrically, a singular matrix collapses space. It maps at least one non-zero vector to the zero vector, so the transformation cannot be reversed.

Given $`[A]\{x\} = \{b\}`$, the formal solution is:

```math
\{x\} = [A]^{-1}\{b\}
```

However, in practice, computing the inverse explicitly is rarely the best way to solve a linear system. The inverse is primarily useful for theoretical analysis and physical interpretation.

---

<br>

## 2. Computing the Inverse via LU Decomposition

### 2.1 The Column-by-Column Method

The inverse can be computed by solving $`n`$ linear systems, one for each column of the identity matrix:

```math
[A]\{x_k\} = \{e_k\}, \quad k = 1, 2, \ldots, n
```

where $`\{e_k\}`$ is the $`k`$th column of the identity matrix (all zeros except a 1 in position $`k`$). Each solution $`\{x_k\}`$ becomes the $`k`$th column of $`[A]^{-1}`$.

The key advantage of using LU decomposition for this task is that the decomposition $`[A] = [L][U]`$ is performed **only once**. After that, each of the $`n`$ right-hand sides requires only a forward substitution ($`[L]\{d\} = \{e_k\}`$) and a back substitution ($`[U]\{x_k\} = \{d\}`$), both of which are $`O(n^2)`$ operations. The total cost is therefore:

- **One** LU decomposition: $`O(n^3)`$
- **$`n`$** forward and back substitutions: $`n \times O(n^2) = O(n^3)`$
- **Total**: $`O(n^3)`$, the same order as a single Gauss elimination.

### 2.2 Example 11.1: 3x3 Inverse

Consider the matrix from the textbook:

```math
[A] = \begin{bmatrix} 3 & -0.1 & -0.2 \\ 0.1 & 7 & -0.3 \\ 0.3 & -0.2 & 10 \end{bmatrix}
```

Using the LU decomposition obtained in Chapter 10, we solve for each column of the inverse.

**Column 1**: Solve $`[L]\{d\} = \{1, 0, 0\}^T`$, then $`[U]\{x_1\} = \{d\}`$.

Forward substitution ($`[L]\{d\} = \{1, 0, 0\}^T`$):

```math
d_1 = 1, \quad d_2 = 0 - (0.033333)(1) = -0.033333, \quad d_3 = 0 - (0.1)(1) - (-0.02713)(-0.033333) = -0.1009
```

Back substitution ($`[U]\{x_1\} = \{d\}`$):

```math
x_3 = \frac{-0.1009}{10.012} = -0.010078, \quad x_2 = \frac{-0.033333 - (-0.29333)(-0.010078)}{7.00333} = -0.005182
```

```math
x_1 = \frac{1 - (-0.1)(-0.005182) - (-0.2)(-0.010078)}{3} = 0.3325
```

**Column 2**: Repeating the same procedure with $`\{e_2\} = \{0, 1, 0\}^T`$ yields the second column.

**Column 3**: Repeating the same procedure with $`\{e_3\} = \{0, 0, 1\}^T`$ yields the third column.

**Result:**

```math
[A]^{-1} = \begin{bmatrix} 0.3325 & 0.0049 & 0.0068 \\ -0.0052 & 0.1429 & 0.0042 \\ -0.0101 & 0.0027 & 0.0999 \end{bmatrix}
```

The result can be verified by checking that $`[A][A]^{-1} \approx [I]`$ within rounding error.

---

<br>

## 3. Physical Interpretation: Mass-Spring System

The inverse provides direct physical insight when the system $`[A]\{x\} = \{b\}`$ has a physical meaning. The formal solution is:

```math
\{x\} = [A]^{-1}\{b\}
```

Each element $`a_{ij}^{-1}`$ of the inverse matrix has the interpretation:

```math
a_{ij}^{-1} = \text{contribution of force } b_j \text{ to displacement } x_i
```

In other words, $`a_{ij}^{-1}`$ tells us how much the $`j`$th input (force) affects the $`i`$th output (displacement). This is analogous to an **influence coefficient** or **transfer function** in engineering.

For a mass-spring system, reading the columns of $`[A]^{-1}`$ reveals how each individual force component propagates through the system. The $`j`$th column of the inverse shows the displacement pattern that results from a unit force applied at position $`j`$ and zero forces elsewhere.

> **Note:** **Flexibility matrix.** In structural engineering, the inverse of the stiffness matrix is called the **flexibility matrix** (or compliance matrix). Its entries directly measure how flexible the structure is: large values indicate that a small force produces a large displacement. This physical interpretation is one of the few situations where computing the full inverse is genuinely useful.

---

<br>

## 4. Error Analysis: Three Checks for Ill-Conditioning

When solving a linear system numerically, it is important to determine whether the results can be trusted. A system is **ill-conditioned** when small changes in the input data (coefficients or right-hand side) produce disproportionately large changes in the solution. Three practical checks can detect ill-conditioning:

1. **Scale and inspect the inverse:** Normalize the matrix so that the element of largest absolute value in each row is 1, then compute the inverse. If any elements of the inverse are **much larger than 1**, the system is likely ill-conditioned.

2. **Check if $`[A]^{-1}[A] = [I]`$:** Compute the product of the matrix and its computed inverse. If the result deviates significantly from the identity matrix, the inverse is inaccurate, which signals ill-conditioning.

3. **Check if $`([A]^{-1})^{-1} = [A]`$:** Compute the inverse of the inverse. If this does not recover the original matrix to adequate precision, the system is ill-conditioned.

These checks are useful but qualitative. A quantitative measure of conditioning requires the **condition number**, and defining the condition number requires the concept of matrix norms.

---

<br>

## 5. Vector Norms

A **norm** is a function that assigns a non-negative size, or length, to a vector. The general $`L_p`$ norm ($`p \geq 1`$) is defined as:

```math
\|x\|_p = \left(\sum_{i=1}^{n} |x_i|^p\right)^{1/p}
```

The three most commonly used norms are:

| Norm | Formula | Name |
|:-----|:--------|:-----|
| $`L_1`$ | $`\Vert x \Vert_1 = \sum_{i=1}^{n} \vert x_i \vert`$ | Manhattan (taxicab) norm |
| $`L_2`$ | $`\Vert x \Vert_2 = \sqrt{\sum_{i=1}^{n} x_i^2}`$ | Euclidean norm |
| $`L_\infty`$ | $`\Vert x \Vert_\infty = \max_{i} \vert x_i \vert`$ | Max (Chebyshev) norm |

*Table 1. Three commonly used vector norms*

**Example:** For $`x = [1, 2, 3]`$:

- $`\|x\|_1 = |1| + |2| + |3| = 6`$
- $`\|x\|_2 = \sqrt{1^2 + 2^2 + 3^2} = \sqrt{14} \approx 3.742`$
- $`\|x\|_\infty = \max(|1|, |2|, |3|) = 3`$

All norms satisfy four fundamental properties:

1. **Non-negativity:** $`\|x\| \geq 0`$
2. **Definiteness:** $`\|x\| = 0`$ if and only if $`x = 0`$.
3. **Homogeneity:** $`\|\alpha x\| = |\alpha| \cdot \|x\|`$
4. **Triangle inequality:** $`\|x + y\| \leq \|x\| + \|y\|`$

> **[Linear Algebra]** Norms generalize the concept of length to arbitrary vector spaces. All norms in finite dimensions are equivalent (same convergence), but they give different numerical values. The choice of norm therefore affects the value of the condition number.

---

<br>

## 6. Matrix Norms

Matrix norms extend the concept of vector norms to matrices. There are two main types:

**Element-based norm:**

- **Frobenius norm:** $`\|A\|_F = \sqrt{\sum_i \sum_j a_{ij}^2}`$. It treats the matrix as a long vector and takes its Euclidean norm.

**Induced (operator) norms:** They are defined as the maximum factor by which the matrix stretches a vector:

```math
\|A\|_a = \max_{\|x\|_a = 1} \|Ax\|_a
```

This measures the worst-case amplification: the largest factor by which $`A`$ can multiply the norm of any vector.

The induced norms corresponding to the common vector norms are:

| Matrix Norm | Formula | Computation |
|:------------|:--------|:------------|
| Column-sum ($`L_1`$) | $`\Vert A \Vert_1 = \max_j \sum_i \vert a_{ij} \vert`$ | Maximum absolute column sum |
| Spectral ($`L_2`$) | $`\Vert A \Vert_2 = \sqrt{\mu_{max}}`$ | $`\mu_{max}`$ = largest eigenvalue of $`A^T A`$ |
| Row-sum ($`L_\infty`$) | $`\Vert A \Vert_\infty = \max_i \sum_j \vert a_{ij} \vert`$ | Maximum absolute row sum |

*Table 2. Induced matrix norms corresponding to vector norms*

Matrix norms satisfy the same four properties as vector norms (non-negativity, definiteness, homogeneity, triangle inequality), and the induced norms and the Frobenius norm have an additional crucial property:

**Submultiplicativity:** $`\|AB\| \leq \|A\| \cdot \|B\|`$, and in particular $`\|Ax\| \leq \|A\| \cdot \|x\|`$ holds for induced norms.

This property is what makes matrix norms useful for bounding errors: it allows us to bound the norm of a product by the product of the norms.

> **[Linear Algebra]** The spectral norm (L2 matrix norm) equals the largest singular value of A. It measures the maximum stretching the matrix can apply to a unit vector.

---

<br>

## 7. Matrix Condition Number

The **condition number** of a matrix $`[A]`$ is defined as:

```math
\text{Cond}[A] = \|A\| \cdot \|A^{-1}\|
```

The condition number is also written $`\kappa(A)`$. Its value depends on which matrix norm is used (hence we may write $`\kappa_1`$, $`\kappa_2`$, $`\kappa_\infty`$ for different norms), but the qualitative interpretation is the same regardless of norm choice.

Key properties for induced norms:

- $`\text{Cond}[I] = 1`$. The identity matrix is perfectly conditioned (best possible).
- $`\text{Cond}[A] \geq 1`$ for any non-singular matrix.
- If $`\text{Cond}[A] \approx 1`$, the system is **well-conditioned**. Small perturbations in input produce proportionally small changes in output.
- If $`\text{Cond}[A] \gg 1`$, the system is **ill-conditioned**. Small perturbations in input may produce disproportionately large changes in output.
- If $`\text{Cond}[A] = \infty`$, the matrix is **singular** and has no inverse.

The condition number quantifies how close to singular a matrix is. A large condition number means the matrix is nearly singular: its rows (or columns) are nearly linearly dependent.

---

<br>

## 8. Error Bound Theorem

The condition number provides a rigorous bound on how errors in the coefficient matrix propagate to errors in the solution. If $`[A]`$ is perturbed by a small amount $`[\delta A]`$ and we solve $`(A + \delta A)(x + \delta x) = b`$, then for small perturbations:

```math
\frac{\|\delta x\|}{\|x\|} \leq \text{Cond}[A] \cdot \frac{\|\delta A\|}{\|A\|}
```

Strictly, the inequality holds exactly with $`\|x + \delta x\|`$ in the denominator of the left side, and for small perturbations the difference from $`\|x\|`$ is negligible. This states that the **relative error in the solution** is bounded by the condition number times the **relative error in the coefficients**.

The condition number acts as an **amplification factor**: it tells us the worst-case ratio by which relative errors in the input are magnified in the output. For a well-conditioned system ($`\text{Cond} \approx 1`$), errors are not amplified. For an ill-conditioned system ($`\text{Cond} \gg 1`$), even tiny perturbations in the coefficients can produce enormous changes in the solution.

Similarly, for perturbations in the right-hand side, $`A(x + \delta x) = b + \delta b`$:

```math
\frac{\|\delta x\|}{\|x\|} \leq \text{Cond}[A] \cdot \frac{\|\delta b\|}{\|b\|}
```

The bounds thus apply to perturbations in $`[A]`$ and in $`\{b\}`$ separately.

---

<br>

## 9. Precision Rule

The error bound theorem leads to a practical rule for estimating how many significant digits of accuracy are lost due to ill-conditioning.

If the coefficients of $`[A]`$ are known to $`t`$ significant digits, then the solution $`\{x\}`$ is valid to approximately:

```math
t - \log_{10}(\text{Cond}[A]) \quad \text{significant digits}
```

**Example:** If $`\text{Cond}[A] = 10^3`$ and the coefficients are known to 7 significant digits (the level of single-precision rounding), then the solution is valid to approximately $`7 - 3 = 4`$ significant digits.

This rule provides an immediate practical warning: if $`\log_{10}(\text{Cond}[A]) \geq t`$, then the computed solution may have **no** reliable significant digits. In double precision ($`t \approx 16`$), a condition number of $`10^{16}`$ or larger means the solution is essentially meaningless.

> **Note:** This is a worst-case bound. The actual error may be much smaller, but you cannot guarantee it without further analysis. When in doubt, check the residual $`\|Ax - b\|`$ as a practical accuracy indicator. When the condition number is large, however, even a small residual does not guarantee a small error.

---

<br>

## 10. Example 11.3: Hilbert Matrix

The **Hilbert matrix** is the canonical example of an ill-conditioned matrix. It is defined by:

```math
H_{ij} = \frac{1}{i + j - 1}
```

For $`n = 3`$:

```math
[H] = \begin{bmatrix} 1 & 1/2 & 1/3 \\ 1/2 & 1/3 & 1/4 \\ 1/3 & 1/4 & 1/5 \end{bmatrix}
```

**Computing the condition number using the row-sum norm ($`L_\infty`$):**

1. **Normalize** the matrix so the maximum element in each row is 1. Divide row 1 by 1, row 2 by 1/2, and row 3 by 1/3 to obtain the scaled matrix $`[A'] = \begin{bmatrix} 1 & 1/2 & 1/3 \\ 1 & 2/3 & 1/2 \\ 1 & 3/4 & 3/5 \end{bmatrix}`$.

2. **Compute $`\|A'\|_\infty`$:** Sum the absolute values of each row and take the maximum. The row sums are $`1.8333`$, $`2.1667`$, and $`2.35`$, so $`\|A'\|_\infty = 2.35`$.

3. **Compute the inverse** of the scaled matrix and find $`\|[A']^{-1}\|_\infty`$. The inverse is $`[A']^{-1} = \begin{bmatrix} 9 & -18 & 10 \\ -36 & 96 & -60 \\ 30 & -90 & 60 \end{bmatrix}`$, and its absolute row sums are $`37`$, $`192`$, and $`180`$. The maximum is row 2 with $`36 + 96 + 60 = 192`$, giving $`\|[A']^{-1}\|_\infty = 192`$.

4. **Condition number:** $`\text{Cond}[A'] = 2.35 \times 192 = 451.2`$

Since $`\text{Cond}[A'] = 451.2 \gg 1`$, the $`3 \times 3`$ Hilbert matrix is already poorly conditioned. This means solutions to systems involving this matrix will lose approximately $`\log_{10}(451) \approx 2.7`$ digits of accuracy.

The situation deteriorates rapidly with increasing size. For larger Hilbert matrices, the condition number grows approximately as $`\kappa(H_n) \sim e^{3.5n}`$, that is, in proportion to $`e^{3.5n}`$, making them practically unsolvable for $`n > 12`$ or so in double precision.

> **Note:** **Where Hilbert matrices arise.** Hilbert matrices are the canonical example of ill-conditioning. Their condition number grows exponentially with size: $`\kappa(H_n) \sim e^{3.5n}`$. They arise naturally in polynomial least-squares fitting.

---

<br>

## 11. Python Implementation

Python's NumPy provides built-in functions for computing matrix inverses, norms, and condition numbers:

```python
import numpy as np

# Matrix inverse
A_inv = np.linalg.inv(A)

# Condition number (2-norm by default)
kappa = np.linalg.cond(A)

# Condition number with specific norms
kappa_1   = np.linalg.cond(A, 1)       # 1-norm condition number
kappa_2   = np.linalg.cond(A, 2)       # 2-norm (default)
kappa_inf = np.linalg.cond(A, np.inf)  # infinity-norm condition number

# Norms
norm_1   = np.linalg.norm(A, ord=1)       # max abs column sum
norm_2   = np.linalg.norm(A, ord=2)       # spectral norm
norm_inf = np.linalg.norm(A, ord=np.inf)  # max abs row sum
norm_fro = np.linalg.norm(A, ord='fro')   # Frobenius norm
```

> **Note:** In practice, always use `np.linalg.solve(A, b)` instead of `np.linalg.inv(A) @ b` for solving linear systems. Direct solving is faster and more numerically accurate. The inverse should only be computed when the matrix entries themselves are needed (e.g., for physical interpretation or covariance analysis).

---

<br>

## Summary

| Concept | Key Summary |
|:------|:---------------------|
| Inverse definition | The inverse satisfies $`[A][A]^{-1} = [A]^{-1}[A] = [I]`$ and exists only for square non-singular matrices. |
| Inverse via LU | After factoring once, solving $`[A]\{x_k\} = \{e_k\}`$ for each column $`k`$ gives the inverse in $`O(n^3)`$. |
| Physical meaning of the inverse | $`a_{ij}^{-1}`$ measures how much the $`j`$th input affects the $`i`$th output. |
| $`L_1`$ norm | For vectors $`\Vert x \Vert_1 = \sum \vert x_i \vert`$, and the matrix norm $`\Vert A \Vert_1`$ is the maximum absolute column sum. |
| $`L_2`$ norm | For vectors $`\Vert x \Vert_2 = \sqrt{\sum x_i^2}`$, and the matrix norm $`\Vert A \Vert_2`$ is the largest singular value $`\sigma_{max}`$. |
| $`L_\infty`$ norm | For vectors $`\Vert x \Vert_\infty = \max \vert x_i \vert`$, and the matrix norm $`\Vert A \Vert_\infty`$ is the maximum absolute row sum. |
| Condition number | $`\kappa(A) = \Vert A \Vert \cdot \Vert A^{-1} \Vert`$ measures how close a matrix is to being singular. |
| Error bound | Since $`\Vert \delta x \Vert / \Vert x \Vert \leq \kappa(A) \cdot \Vert \delta A \Vert / \Vert A \Vert`$, the condition number is the amplification factor of relative errors. |
| Precision loss | If the coefficients are accurate to $`t`$ digits, the solution is valid to approximately $`t - \log_{10}\kappa`$ digits. |
| Hilbert matrix | It is defined by $`H_{ij} = 1/(i+j-1)`$, and its condition number grows exponentially with $`n`$. |

---

<br>

## Self-Check Questions

1. **Vector norms:** Find the $`L_1`$, $`L_2`$, and $`L_\infty`$ norms of $`x = (3, -4, 12)`$.

   <details>
   <summary>Answer</summary>

   $`\|x\|_1 = 3 + 4 + 12 = 19`$, $`\|x\|_2 = \sqrt{9 + 16 + 144} = \sqrt{169} = 13`$, and $`\|x\|_\infty = \max(3, 4, 12) = 12`$.

   </details>

2. **Matrix norms and condition number:** For $`A = \begin{bmatrix} 4 & 2 \\ 1 & 3 \end{bmatrix}`$, find $`\|A\|_1`$, $`\|A\|_\infty`$, $`A^{-1}`$, and the row-sum condition number $`\kappa_\infty(A)`$.

   <details>
   <summary>Answer</summary>

   The column sums are $`5, 5`$, so $`\|A\|_1 = 5`$, and the row sums are $`6, 4`$, so $`\|A\|_\infty = 6`$. Since $`\det(A) = 12 - 2 = 10`$, $`A^{-1} = \frac{1}{10}\begin{bmatrix} 3 & -2 \\ -1 & 4 \end{bmatrix} = \begin{bmatrix} 0.3 & -0.2 \\ -0.1 & 0.4 \end{bmatrix}`$. The absolute row sums of $`A^{-1}`$ are $`0.5, 0.5`$, so $`\|A^{-1}\|_\infty = 0.5`$ and $`\kappa_\infty(A) = 6 \times 0.5 = 3`$.

   </details>

3. **Precision rule:** If coefficients stored in double precision ($`t \approx 16`$) have a condition number of $`10^{10}`$, to about how many digits is the solution valid? What about a $`3 \times 3`$ Hilbert system ($`\text{Cond} = 451.2`$) whose coefficients have 7 significant digits?

   <details>
   <summary>Answer</summary>

   The first case gives $`16 - \log_{10}(10^{10}) = 16 - 10 = 6`$ digits. In the second case, $`\log_{10}(451.2) \approx 2.65`$, so $`7 - 2.65 \approx 4.3`$, that is, about 4 digits.

   </details>

4. **Error bound:** In a system with $`\text{Cond}[A] = 1000`$, the relative error in the right-hand side is $`\|\delta b\| / \|b\| = 10^{-6}`$. What is the largest possible relative error in the solution?

   <details>
   <summary>Answer</summary>

   $`\|\delta x\| / \|x\| \leq \text{Cond}[A] \cdot \|\delta b\| / \|b\| = 1000 \times 10^{-6} = 10^{-3}`$. The relative error in the solution is therefore at most 0.1%.

   </details>

5. **Inverse via LU and its meaning:** State the operation orders of the factorization and the substitutions when computing the inverse of an $`n \times n`$ matrix with LU decomposition, and explain what the element $`a_{23}^{-1}`$ of the inverse means for a mass-spring system.

   <details>
   <summary>Answer</summary>

   The LU factorization is done only once at $`O(n^3)`$. Each of the $`n`$ columns of the identity matrix needs a forward and a back substitution at $`O(n^2)`$, so all substitutions together also cost $`n \times O(n^2) = O(n^3)`$, and the total cost is $`O(n^3)`$. $`a_{23}^{-1}`$ is the displacement of mass 2 caused by a unit force at position 3, that is, the contribution of force $`b_3`$ to displacement $`x_2`$.

   </details>
