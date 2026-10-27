# Chapter 8 — Solving Linear Algebraic Equations

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-27
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 8, 9

> **Learning Objectives:**
> 1. Formulate engineering problems as linear algebraic equations
> 2. Distinguish the special types of matrices and apply the matrix operating rules, including multiplication, transposition, and inversion
> 3. Identify solution existence and uniqueness conditions
> 4. Interpret 2x2 linear systems graphically from the row view and the column view
> 5. Solve linear systems with naive Gauss elimination and analyze its computational cost

---

<br>

## Table of Contents

- [1. Motivation](#1-motivation)
- [2. Matrices and Matrix Operations](#2-matrices-and-matrix-operations)
  - [2.1 Matrix Notation](#21-matrix-notation)
  - [2.2 Special Types of Square Matrices](#22-special-types-of-square-matrices)
  - [2.3 Row and Column Vectors](#23-row-and-column-vectors)
  - [2.4 Equality, Addition, and Scalar Multiplication](#24-equality-addition-and-scalar-multiplication)
  - [2.5 Matrix Multiplication](#25-matrix-multiplication)
  - [2.6 Transpose](#26-transpose)
  - [2.7 Matrix Inverse](#27-matrix-inverse)
  - [2.8 Permutation Matrix](#28-permutation-matrix)
  - [2.9 Augmented Matrix](#29-augmented-matrix)
  - [2.10 Matrix Form of Linear Algebraic Equations](#210-matrix-form-of-linear-algebraic-equations)
  - [2.11 Stiffness Matrix of a Spring-Mass System](#211-stiffness-matrix-of-a-spring-mass-system)
  - [2.12 Solving Linear Equations with Python](#212-solving-linear-equations-with-python)
- [3. Small Systems (2x2): Graphical Method](#3-small-systems-2x2-graphical-method)
  - [3.1 Two Equations as Two Lines](#31-two-equations-as-two-lines)
  - [3.2 Three Possibilities](#32-three-possibilities)
  - [3.3 Column View](#33-column-view)
- [4. Determinant and Solvability](#4-determinant-and-solvability)
  - [4.1 Determinant of a 2x2 System](#41-determinant-of-a-2x2-system)
  - [4.2 Singular and Ill-Conditioned Systems](#42-singular-and-ill-conditioned-systems)
- [5. Cramer's Rule](#5-cramers-rule)
  - [5.1 Derivation for 2x2](#51-derivation-for-2x2)
  - [5.2 General Form](#52-general-form)
  - [5.3 Limitation](#53-limitation)
- [6. Naive Gauss Elimination](#6-naive-gauss-elimination)
  - [6.1 Forward Elimination](#61-forward-elimination)
  - [6.2 Back Substitution](#62-back-substitution)
- [7. Computational Cost](#7-computational-cost)
- [8. Example: 3x3 System (Complete Step-by-Step)](#8-example-3x3-system-complete-step-by-step)
  - [8.1 Problem Setup](#81-problem-setup)
  - [8.2 Forward Elimination: Step 1](#82-forward-elimination-step-1)
  - [8.3 Forward Elimination: Step 2](#83-forward-elimination-step-2)
  - [8.4 Back Substitution](#84-back-substitution)
- [9. Pitfalls of Naive Gauss Elimination](#9-pitfalls-of-naive-gauss-elimination)
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

## 2. Matrices and Matrix Operations

The methods of this chapter operate on whole arrays of coefficients, so they are stated in matrix notation. This section defines the notation, the special types of matrices that recur in numerical linear algebra, and the operating rules for matrices, and then writes a system of linear algebraic equations in matrix form.

### 2.1 Matrix Notation

A **matrix** consists of a rectangular array of elements and is denoted by a single symbol in brackets, such as $`[A]`$. The element in row $`i`$ and column $`j`$ is written $`a_{ij}`$, where $`i`$ is the **row index** and $`j`$ is the **column index**. A horizontal set of elements is a **row**, and a vertical set is a **column**.

A matrix with $`m`$ rows and $`n`$ columns has the **dimension** $`m \times n`$ and is called an $`m`$ by $`n`$ matrix:

```math
[A] = \begin{bmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{bmatrix}
```

For example, the following matrix is $`2 \times 3`$, and its elements include $`a_{12} = -7`$, $`a_{21} = 8`$, and $`a_{23} = -6`$:

```math
[A] = \begin{bmatrix} 2 & -7 & 4 \\ 8 & 3 & -6 \end{bmatrix}
```

A matrix with $`m = n`$ is a **square matrix**. The elements $`a_{11}, a_{22}, \ldots, a_{nn}`$ of a square matrix form its **main diagonal**, and all other elements are **off-diagonal** elements. A system with as many equations as unknowns has a square coefficient matrix.

### 2.2 Special Types of Square Matrices

Square matrices are classified by relations among their elements or by the positions in which their nonzero elements can appear:

- **Symmetric matrix:** $`a_{ij} = a_{ji}`$ for all $`i`$ and $`j`$, so the elements are mirrored across the main diagonal.
- **Diagonal matrix:** a square matrix whose off-diagonal elements are all zero.
- **Identity matrix:** a diagonal matrix whose diagonal elements are all 1. The $`n \times n`$ identity matrix is written $`[I_n]`$, or $`[I]`$ when the size is clear from the context.
- **Upper triangular matrix:** a square matrix whose elements below the main diagonal are all zero ($`a_{ij} = 0`$ for $`i > j`$).
- **Lower triangular matrix:** a square matrix whose elements above the main diagonal are all zero ($`a_{ij} = 0`$ for $`i < j`$).
- **Banded matrix:** a square matrix whose nonzero elements are confined to a band of diagonals centered on the main diagonal, with all elements outside the band equal to zero.
- **Tridiagonal matrix:** a banded matrix whose band consists of the main diagonal and the diagonals immediately above and below it ($`a_{ij} = 0`$ for $`|i - j| > 1`$).

Figure 1 shows, for $`n = 5`$, the positions of the elements that can be nonzero in diagonal, triangular, tridiagonal, and banded matrices. An asterisk marks an element that can be nonzero, and 0 marks an element that is zero by definition. The banded example has a band of five diagonals, two on each side of the main diagonal.

```text
Diagonal           Upper triangular   Lower triangular   Tridiagonal        Banded
* 0 0 0 0          * * * * *          * 0 0 0 0          * * 0 0 0          * * * 0 0
0 * 0 0 0          0 * * * *          * * 0 0 0          * * * 0 0          * * * * 0
0 0 * 0 0          0 0 * * *          * * * 0 0          0 * * * 0          * * * * *
0 0 0 * 0          0 0 0 * *          * * * * 0          0 0 * * *          0 * * * *
0 0 0 0 *          0 0 0 0 *          * * * * *          0 0 0 * *          0 0 * * *
```

*Figure 1. Positions of the elements that can be nonzero in 5x5 diagonal, triangular, tridiagonal, and banded matrices*

For example, the first of the following matrices is symmetric ($`a_{12} = a_{21} = 1`$, $`a_{13} = a_{31} = 2`$, $`a_{23} = a_{32} = 7`$), the second is diagonal, and the third is the $`3 \times 3`$ identity matrix:

```math
\begin{bmatrix} 5 & 1 & 2 \\ 1 & 3 & 7 \\ 2 & 7 & 8 \end{bmatrix}, \quad \begin{bmatrix} 4 & 0 & 0 \\ 0 & -2 & 0 \\ 0 & 0 & 7 \end{bmatrix}, \quad [I_3] = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}
```

The identity matrix plays the role of the number 1 in matrix multiplication (Section 2.5): for a square matrix $`[A]`$ of the same size, $`[A][I] = [I][A] = [A]`$. Triangular matrices are the target forms of elimination. Gauss elimination (Section 6) reduces $`[A]`$ to upper triangular form, and LU factorization (Chapter 10) writes $`[A]`$ as the product of a lower and an upper triangular matrix. Tridiagonal and other banded matrices arise when each unknown is coupled only to its neighbors, and Chapter 9 solves tridiagonal systems with the Thomas algorithm.

### 2.3 Row and Column Vectors

A matrix with only one row is a **row vector**, and a matrix with only one column is a **column vector**. A row vector with $`n`$ elements has the dimension $`1 \times n`$, and a column vector with $`m`$ elements has the dimension $`m \times 1`$. Row vectors are written with floor brackets $`\lfloor \ \rfloor`$ and column vectors with braces $`\{\ \}`$:

```math
[b] = \begin{bmatrix} b_1 & b_2 & \cdots & b_n \end{bmatrix} = \lfloor b \rfloor
```
```math
[c] = \begin{bmatrix} c_1 \\ c_2 \\ \vdots \\ c_m \end{bmatrix} = \{c\}
```

The unknown vector $`\{x\}`$ and the right-hand side vector $`\{b\}`$ of Section 1 are column vectors.

### 2.4 Equality, Addition, and Scalar Multiplication

Two matrices of the same dimension are **equal** if and only if all their corresponding elements are equal:

```math
[A] = [B] \iff a_{ij} = b_{ij} \quad \forall\, i, j
```

**Addition** and **subtraction** are defined for matrices of the same dimension and are carried out element by element:

```math
[C] = [A] + [B], \quad c_{ij} = a_{ij} + b_{ij} \quad \forall\, i, j
```
```math
[D] = [E] - [F], \quad d_{ij} = e_{ij} - f_{ij} \quad \forall\, i, j
```

**Multiplication by a scalar** $`g`$ multiplies every element of the matrix by $`g`$:

```math
g[A] = \begin{bmatrix} ga_{11} & ga_{12} & \cdots & ga_{1n} \\ ga_{21} & ga_{22} & \cdots & ga_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ ga_{m1} & ga_{m2} & \cdots & ga_{mn} \end{bmatrix}
```

For example, take the following two $`2 \times 2`$ matrices:

```math
[A] = \begin{bmatrix} 1 & 2 \\ 3 & 4 \end{bmatrix}, \quad [B] = \begin{bmatrix} 5 & -1 \\ 0 & 2 \end{bmatrix}
```

Their sum, their difference, and the scalar multiple $`3[A]`$ are:

```math
[A] + [B] = \begin{bmatrix} 6 & 1 \\ 3 & 6 \end{bmatrix}, \quad [A] - [B] = \begin{bmatrix} -4 & 3 \\ 3 & 2 \end{bmatrix}, \quad 3[A] = \begin{bmatrix} 3 & 6 \\ 9 & 12 \end{bmatrix}
```

Because each element is added as an ordinary number, matrix addition is commutative and associative (Table 1).

### 2.5 Matrix Multiplication

The product $`[C] = [A][B]`$ is defined only when the number of columns of $`[A]`$ equals the number of rows of $`[B]`$. If $`[A]`$ is $`m \times n`$ and $`[B]`$ is $`n \times p`$, the product is $`m \times p`$:

```math
\underset{m \times p}{[C]} = \underset{m \times n}{[A]}\,\underset{n \times p}{[B]}
```

Each element $`c_{ij}`$ is the sum of the products of the elements of row $`i`$ of $`[A]`$ and the corresponding elements of column $`j`$ of $`[B]`$:

```math
c_{ij} = \sum_{k=1}^{n} a_{ik}b_{kj} = a_{i1}b_{1j} + a_{i2}b_{2j} + \cdots + a_{in}b_{nj}
```

For instance, $`c_{11} = a_{11}b_{11} + a_{12}b_{21} + \cdots + a_{1n}b_{n1} = \sum_{k=1}^{n} a_{1k}b_{k1}`$ combines the first row of $`[A]`$ with the first column of $`[B]`$.

**Example:** Multiply the $`3 \times 2`$ matrix $`[X]`$ by the $`2 \times 2`$ matrix $`[Y]`$:

```math
[X] = \begin{bmatrix} 3 & 1 \\ 8 & 6 \\ 0 & 4 \end{bmatrix}, \quad [Y] = \begin{bmatrix} 5 & 9 \\ 7 & 2 \end{bmatrix}
```

The inner dimensions agree ($`[X]`$ has 2 columns and $`[Y]`$ has 2 rows), so $`[Z] = [X][Y]`$ exists and is $`3 \times 2`$. Each element combines a row of $`[X]`$ with a column of $`[Y]`$:

- $`z_{11} = (3)(5) + (1)(7) = 22`$, $`z_{12} = (3)(9) + (1)(2) = 29`$
- $`z_{21} = (8)(5) + (6)(7) = 82`$, $`z_{22} = (8)(9) + (6)(2) = 84`$
- $`z_{31} = (0)(5) + (4)(7) = 28`$, $`z_{32} = (0)(9) + (4)(2) = 8`$

```math
[Z] = [X][Y] = \begin{bmatrix} 22 & 29 \\ 82 & 84 \\ 28 & 8 \end{bmatrix}
```

The reverse product $`[Y][X]`$ is not defined, because $`[Y]`$ has 2 columns but $`[X]`$ has 3 rows. A matrix-vector product is the special case $`p = 1`$, for example:

```math
\begin{bmatrix} 2 & 1 \\ 1 & -1 \end{bmatrix} \begin{bmatrix} 2 \\ 1 \end{bmatrix} = \begin{bmatrix} (2)(2) + (1)(1) \\ (1)(2) + (-1)(1) \end{bmatrix} = \begin{bmatrix} 5 \\ 1 \end{bmatrix}
```

Matrix multiplication is associative and distributive whenever the dimensions allow the products, but it is **not commutative** in general: even when both $`[A][B]`$ and $`[B][A]`$ exist, they usually differ. For the matrix $`[A]`$ of Section 2.4 and the matrix $`[S]`$ below:

```math
[S] = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}, \quad [A][S] = \begin{bmatrix} 2 & 1 \\ 4 & 3 \end{bmatrix}, \quad [S][A] = \begin{bmatrix} 3 & 4 \\ 1 & 2 \end{bmatrix}
```

Multiplying by $`[S]`$ on the right swaps the two columns of $`[A]`$, while multiplying on the left swaps the two rows, so $`[A][S] \neq [S][A]`$. The order of the factors must therefore be kept, and multiplying $`[A]`$ on the left (**premultiplication**) is distinguished from multiplying it on the right (**postmultiplication**). Table 1 collects the operating rules.

| Rule | Statement | Condition |
|:-----|:-----|:-----|
| Addition is commutative | $`[A] + [B] = [B] + [A]`$ | Same dimension |
| Addition is associative | $`([A] + [B]) + [C] = [A] + ([B] + [C])`$ | Same dimension |
| Multiplication is associative | $`([A][B])[C] = [A]([B][C])`$ | Inner dimensions agree |
| Multiplication is distributive | $`[A]([B] + [C]) = [A][B] + [A][C]`$, $`([A] + [B])[C] = [A][C] + [B][C]`$ | Every sum and product is defined |
| Multiplication is not commutative | $`[A][B] \neq [B][A]`$ in general | Both products are defined |
| Identity matrix | $`[A][I] = [I][A] = [A]`$ | $`[A]`$ is square and of the same size as $`[I]`$ |

*Table 1. Operating rules for matrix addition and multiplication*

### 2.6 Transpose

The **transpose** $`[A]^T`$ of $`[A]`$ turns the rows of $`[A]`$ into columns, that is, $`a^T_{ji} = a_{ij}`$ for all $`i`$ and $`j`$. If $`[A]`$ is $`m \times n`$, then $`[A]^T`$ is $`n \times m`$. For a $`3 \times 3`$ matrix:

```math
[A] = \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix}, \quad [A]^T = \begin{bmatrix} a_{11} & a_{21} & a_{31} \\ a_{12} & a_{22} & a_{32} \\ a_{13} & a_{23} & a_{33} \end{bmatrix}
```

The $`2 \times 3`$ example of Section 2.1 has a $`3 \times 2`$ transpose, and transposing a column vector gives a row vector:

```math
\begin{bmatrix} 2 & -7 & 4 \\ 8 & 3 & -6 \end{bmatrix}^T = \begin{bmatrix} 2 & 8 \\ -7 & 3 \\ 4 & -6 \end{bmatrix}, \quad \{c\} = \begin{bmatrix} c_1 \\ c_2 \\ c_3 \end{bmatrix} \implies \{c\}^T = \lfloor c_1 \quad c_2 \quad c_3 \rfloor
```

Since a symmetric matrix is defined by $`a_{ij} = a_{ji}`$, a square matrix is symmetric exactly when $`[A]^T = [A]`$.

### 2.7 Matrix Inverse

The **inverse** $`[A]^{-1}`$ of a square matrix $`[A]`$ is the matrix that satisfies:

```math
[A]^{-1}[A] = [A][A]^{-1} = [I]
```

The inverse of $`[A]`$ exists when $`[A]`$ is square and **nonsingular**, that is, when its determinant is not zero (Section 4). A matrix that has an inverse is called **invertible**. For a $`2 \times 2`$ matrix, the inverse is:

```math
[A] = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix}, \quad [A]^{-1} = \frac{1}{\det(A)} \begin{bmatrix} a_{22} & -a_{12} \\ -a_{21} & a_{11} \end{bmatrix}, \quad \det(A) = \begin{vmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{vmatrix} = a_{11}a_{22} - a_{12}a_{21}
```

That is, the two diagonal elements are exchanged, the signs of the two off-diagonal elements are reversed, and the result is divided by the determinant. **Example:** For $`[A] = \begin{bmatrix} 4 & 7 \\ 2 & 6 \end{bmatrix}`$, $`\det(A) = (4)(6) - (7)(2) = 10`$, so:

```math
[A]^{-1} = \frac{1}{10}\begin{bmatrix} 6 & -7 \\ -2 & 4 \end{bmatrix} = \begin{bmatrix} 0.6 & -0.7 \\ -0.2 & 0.4 \end{bmatrix}
```

Multiplying the inverse by the original matrix satisfies the definition: the elements of $`[A][A]^{-1}`$ are $`(4)(0.6) + (7)(-0.2) = 1`$, $`(4)(-0.7) + (7)(0.4) = 0`$, $`(2)(0.6) + (6)(-0.2) = 0`$, and $`(2)(-0.7) + (6)(0.4) = 1`$, so $`[A][A]^{-1} = [I]`$. For larger matrices, the inverse is computed column by column with the Gauss-Jordan method (Chapter 9) or with LU factorization (Chapter 11).

### 2.8 Permutation Matrix

A **permutation matrix** $`[P]`$ is an identity matrix whose rows (equivalently, columns) have been interchanged. Interchanging rows 1 and 3 of $`[I_3]`$ gives:

```math
[P] = \begin{bmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{bmatrix}
```

Premultiplying by $`[P]`$ selects rows: row $`i`$ of $`[P][A]`$ is the row of $`[A]`$ whose number is the column position of the 1 in row $`i`$ of $`[P]`$. For example:

```math
[P][A] = \begin{bmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{bmatrix} \begin{bmatrix} 2 & -7 & 4 \\ 8 & 3 & -6 \\ 5 & 1 & 9 \end{bmatrix} = \begin{bmatrix} 5 & 1 & 9 \\ 8 & 3 & -6 \\ 2 & -7 & 4 \end{bmatrix}
```

The first element, for instance, is $`(0)(2) + (0)(8) + (1)(5) = 5`$, because the 1 in column 3 of the first row of $`[P]`$ picks the third row of $`[A]`$. Postmultiplying by $`[P]`$ selects columns in the same way: column $`j`$ of $`[A][P]`$ is the column of $`[A]`$ whose number is the row position of the 1 in column $`j`$ of $`[P]`$.

```math
[A][P] = \begin{bmatrix} 2 & -7 & 4 \\ 8 & 3 & -6 \\ 5 & 1 & 9 \end{bmatrix} \begin{bmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{bmatrix} = \begin{bmatrix} 4 & -7 & 2 \\ -6 & 3 & 8 \\ 9 & 1 & 5 \end{bmatrix}
```

Thus $`[P][A]`$ interchanges rows 1 and 3 of $`[A]`$, and $`[A][P]`$ interchanges columns 1 and 3. The row interchanges of partial pivoting (Chapter 9) are recorded in a permutation matrix of this kind, which leads to the factorization $`[P][A] = [L][U]`$ (Chapter 10).

### 2.9 Augmented Matrix

A matrix is **augmented** by appending one or more columns to it. For example, augmenting a $`3 \times 3`$ matrix $`[A]`$ with the $`3 \times 3`$ identity matrix gives a $`3 \times 6`$ matrix:

```math
[A|I] = \left[\begin{array}{ccc|ccc} a_{11} & a_{12} & a_{13} & 1 & 0 & 0 \\ a_{21} & a_{22} & a_{23} & 0 & 1 & 0 \\ a_{31} & a_{32} & a_{33} & 0 & 0 & 1 \end{array}\right]
```

The vertical line separates the original columns from the appended ones. A row operation on $`[A|I]`$ acts on both parts at once, and the Gauss-Jordan method turns $`[A|I]`$ into $`[I|A^{-1}]`$ (Chapter 9). Appending the right-hand side vector instead gives the augmented matrix $`[A|b]`$, which the rank conditions of Section 3.2 and Gauss elimination (Section 6) work with.

### 2.10 Matrix Form of Linear Algebraic Equations

Finding $`x`$ such that $`f(x) = 0`$ for a single equation is the root problem of Chapters 5 and 6. For example, $`x^2 + 2x + 1 = 0`$ is equivalent to $`(x + 1)^2 = 0`$, so $`x = -1`$. With two unknowns, the problem is to find $`x`$ and $`y`$ such that $`f(x, y) = 0`$ and $`g(x, y) = 0`$. For $`f(x, y) = y - x`$ and $`g(x, y) = y + x`$, the equations $`y - x = 0`$ and $`y + x = 0`$ describe two lines through the origin, and their intersection gives the solution $`x = 0`$, $`y = 0`$. In general, $`n`$ unknowns $`x_1, x_2, \ldots, x_n`$ must satisfy $`n`$ equations simultaneously:

```math
f_1(x_1, x_2, \ldots, x_n) = 0, \quad f_2(x_1, x_2, \ldots, x_n) = 0, \quad \ldots, \quad f_n(x_1, x_2, \ldots, x_n) = 0
```

When every $`f_i`$ is linear in $`x_1, x_2, \ldots, x_n`$, the system is a system of linear algebraic equations:

```math
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1
```
```math
a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2
```
```math
\vdots
```
```math
a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nn}x_n = b_n
```

By the definition of matrix multiplication (Section 2.5), the left-hand side of equation $`i`$ is element $`i`$ of the product $`[A]\{x\}`$. The whole system is therefore a single matrix equation whose dimensions agree, $`(n \times n)(n \times 1) = n \times 1`$:

```math
\underset{n \times n}{[A]}\,\underset{n \times 1}{\{x\}} = \underset{n \times 1}{\{b\}}
```

For three equations in three unknowns, it reads:

```math
\begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} b_1 \\ b_2 \\ b_3 \end{bmatrix}
```

If $`[A]`$ is invertible, premultiplying both sides by $`[A]^{-1}`$ and using $`[A]^{-1}[A]\{x\} = [I]\{x\} = \{x\}`$ gives the solution:

```math
[A]^{-1}[A]\{x\} = [A]^{-1}\{b\} \quad \Longrightarrow \quad \{x\} = [A]^{-1}\{b\}
```

For the matrix of Section 2.7 and $`\{b\} = \lfloor 1 \quad 2 \rfloor^T`$, the elements of $`\{x\} = [A]^{-1}\{b\}`$ are $`(0.6)(1) + (-0.7)(2) = -0.8`$ and $`(-0.2)(1) + (0.4)(2) = 0.6`$. Substituting back gives $`4(-0.8) + 7(0.6) = 1`$ and $`2(-0.8) + 6(0.6) = 2`$. This formula expresses the solution compactly, but computing the inverse takes more work than solving the system directly, so Gauss elimination (Section 6) and the methods of Chapters 9 and 10 solve the system without forming $`[A]^{-1}`$.

When the coefficient matrix is not square, let $`[A]`$ be $`m \times n`$, so that the system has $`m`$ equations and $`n`$ unknowns. The system is **overdetermined** if $`m > n`$ (more equations than unknowns) and **underdetermined** if $`m < n`$ (fewer equations than unknowns). For example, the first system below is overdetermined ($`m = 3 > n = 2`$) and the second is underdetermined ($`m = 1 < n = 2`$):

```math
\begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \\ a_{31} & a_{32} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} b_1 \\ b_2 \\ b_3 \end{bmatrix}, \qquad \lfloor a_{11} \quad a_{12} \rfloor \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \{b_1\}
```

An overdetermined system in general has no vector $`\{x\}`$ that satisfies all of its equations exactly, and the regression methods of Chapters 14 and 15 instead choose the $`\{x\}`$ that minimizes the sum of the squared residuals. An underdetermined system, if it has a solution at all, has infinitely many: in the example, the single equation $`a_{11}x_1 + a_{12}x_2 = b_1`$ (with $`a_{11}`$ and $`a_{12}`$ not both zero) holds at every point of a line.

### 2.11 Stiffness Matrix of a Spring-Mass System

Figure 2 shows two masses connected in series. Mass $`m_1`$ is connected to a stationary wall by a spring with stiffness $`k_1`$, and mass $`m_2`$ is connected to $`m_1`$ by a spring with stiffness $`k_2`$. Initially, the springs are unstretched and no forces act. Steady forces $`f_1`$ and $`f_2`$ are then applied to $`m_1`$ and $`m_2`$, and the problem is to find the equilibrium displacements $`x_1`$ and $`x_2`$ of the two masses. The displacements are measured from the initial positions and are positive in the direction of the forces.

```text
          f1 →        f2 →
||--k1--[ m1 ]--k2--[ m2 ]
          x1 →        x2 →
```

*Figure 2. Two masses connected in series to a wall by springs with stiffness k1 and k2*

Spring 1 is stretched by $`x_1`$ and pulls $`m_1`$ back toward the wall with the force $`k_1x_1`$ (Hooke's law). Spring 2 is stretched by $`x_2 - x_1`$, so it pulls $`m_1`$ forward and $`m_2`$ back, each with the force $`k_2(x_2 - x_1)`$. The force balances on $`m_1`$ and $`m_2`$ are, respectively:

```math
f_1 + k_2(x_2 - x_1) - k_1x_1 = 0
```
```math
f_2 - k_2(x_2 - x_1) = 0
```

Collecting the terms in each unknown gives two equations in two unknowns:

```math
(k_1 + k_2)x_1 - k_2x_2 = f_1
```
```math
-k_2x_1 + k_2x_2 = f_2
```

With the unknown displacement vector $`\{x\}`$ and the known force vector $`\{f\}`$, the system takes the matrix form $`[K]\{x\} = \{f\}`$:

```math
\underbrace{\begin{bmatrix} k_1 + k_2 & -k_2 \\ -k_2 & k_2 \end{bmatrix}}_{[K]} \underbrace{\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}}_{\{x\}} = \underbrace{\begin{bmatrix} f_1 \\ f_2 \end{bmatrix}}_{\{f\}}
```

The stiffness matrix $`[K]`$ is symmetric, since $`k_{12} = k_{21} = -k_2`$, and it has one more property. A square matrix $`[A]`$ is **positive definite** if the scalar $`\{v\}^T[A]\{v\}`$, a product with the dimensions $`(1 \times n)(n \times n)(n \times 1) = 1 \times 1`$, is positive for every nonzero vector $`\{v\}`$ in $`\mathbb{R}^n`$:

```math
\{v\}^T[A]\{v\} > 0, \quad \{v\} \in \mathbb{R}^n, \ \{v\} \neq \{0\}
```

To check this property for $`[K]`$, the product is expanded with $`\{v\}^T = \lfloor v_1 \quad v_2 \rfloor`$:

```math
\{v\}^T[K]\{v\} = \lfloor v_1 \quad v_2 \rfloor \begin{bmatrix} k_1 + k_2 & -k_2 \\ -k_2 & k_2 \end{bmatrix} \begin{bmatrix} v_1 \\ v_2 \end{bmatrix} = \lfloor v_1 \quad v_2 \rfloor \begin{bmatrix} (k_1 + k_2)v_1 - k_2v_2 \\ -k_2v_1 + k_2v_2 \end{bmatrix}
```
```math
= (k_1 + k_2)v_1^2 - k_2v_1v_2 - k_2v_1v_2 + k_2v_2^2 = k_1v_1^2 + k_2(v_1^2 - 2v_1v_2 + v_2^2) = k_1v_1^2 + k_2(v_1 - v_2)^2
```

If $`k_1 > 0`$ and $`k_2 > 0`$, both terms are nonnegative, and their sum is zero only when $`v_1 = 0`$ and $`v_1 = v_2`$, that is, only when $`\{v\} = \{0\}`$. Hence $`\{v\}^T[K]\{v\} > 0`$ for every nonzero $`\{v\}`$, and $`[K]`$ is **symmetric positive definite** (SPD). The Cholesky factorization of Chapter 10 relies on this property.

The product $`\{v\}^T[K]\{v\}`$ has a physical meaning. With $`v_1 = x_1`$ and $`v_2 = x_2`$, the term $`\frac{1}{2}k_1x_1^2`$ is the energy stored in spring 1 and $`\frac{1}{2}k_2(x_2 - x_1)^2`$ is the energy stored in spring 2. Their sum is the potential energy (PE) of the spring system:

```math
\text{PE} = \frac{1}{2}k_1x_1^2 + \frac{1}{2}k_2(x_2 - x_1)^2 = \frac{1}{2}\lfloor x_1 \quad x_2 \rfloor [K] \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \frac{1}{2}\{x\}^T[K]\{x\}
```

The positive definiteness of $`[K]`$ thus states that every nonzero displacement stores a positive amount of energy in the springs.

**Example:** For $`k_1 = 100`$ N/m, $`k_2 = 200`$ N/m, $`f_1 = 10`$ N, and $`f_2 = 20`$ N, the system and its determinant are:

```math
\begin{bmatrix} 300 & -200 \\ -200 & 200 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 10 \\ 20 \end{bmatrix}, \quad \det(K) = (300)(200) - (-200)(-200) = 20000
```

The $`2 \times 2`$ inverse formula of Section 2.7 then gives:

```math
\{x\} = [K]^{-1}\{f\} = \frac{1}{20000}\begin{bmatrix} 200 & 200 \\ 200 & 300 \end{bmatrix} \begin{bmatrix} 10 \\ 20 \end{bmatrix} = \frac{1}{20000}\begin{bmatrix} 6000 \\ 8000 \end{bmatrix} = \begin{bmatrix} 0.3 \\ 0.4 \end{bmatrix}
```

The displacements are $`x_1 = 0.3`$ m and $`x_2 = 0.4`$ m. They agree with a direct force argument: spring 2 carries $`f_2 = 20`$ N and stretches by $`20/200 = 0.1`$ m $`= x_2 - x_1`$, and spring 1 carries $`f_1 + f_2 = 30`$ N and stretches by $`30/100 = 0.3`$ m $`= x_1`$. The stored energy is $`\frac{1}{2}(100)(0.3)^2 + \frac{1}{2}(200)(0.1)^2 = 4.5 + 1.0 = 5.5`$ J, a positive value, as positive definiteness requires.

### 2.12 Solving Linear Equations with Python

NumPy offers two ways to solve $`[A]\{x\} = \{b\}`$. The first implements $`\{x\} = [A]^{-1}\{b\}`$ literally, and the second calls the `solve` function:

```python
import numpy as np

A = np.array([[4.0, 7.0],
              [2.0, 6.0]])
b = np.array([1.0, 2.0])

x = np.linalg.inv(A) @ b   # form the inverse, then multiply by b
print(x)
x = np.linalg.solve(A, b)  # solve without forming the inverse
print(x)
```

Both `print` calls output `[-0.8  0.6]`, the solution found in Section 2.10. `np.linalg.inv(A)` computes the entire inverse, which amounts to solving one system for each column of $`[I]`$, and then the result is multiplied by `b`. `np.linalg.solve(A, b)` factors $`[A]`$ once by LU factorization with partial pivoting (Chapter 10) and obtains $`\{x\}`$ by forward and back substitution, so it needs fewer operations and is the more efficient choice. Avoiding the explicit inverse also generally reduces the roundoff error.

> **Caution:** For NumPy arrays, `*` multiplies element by element. `np.linalg.inv(A) * b` therefore multiplies each row of the inverse element by element by `b` and returns a $`2 \times 2`$ array instead of the matrix-vector product. The matrix product of arrays is written with `@` (or `np.dot`).

---

<br>

## 3. Small Systems (2x2): Graphical Method

### 3.1 Two Equations as Two Lines

Consider a system of two equations in two unknowns:

```math
a_{11}x_1 + a_{12}x_2 = b_1
```
```math
a_{21}x_1 + a_{22}x_2 = b_2
```

Each equation represents a straight line in the $`x_1`$-$`x_2`$ plane. Solving the system means finding the point(s) that satisfy both line equations simultaneously.

This is the **row view** of the system, in which each row of $`[A]\{x\} = \{b\}`$ is read as one line. When $`a_{12} \neq 0`$ and $`a_{22} \neq 0`$, solving each equation for $`x_2`$ gives the slope-intercept forms:

```math
x_2 = -\frac{a_{11}}{a_{12}}x_1 + \frac{b_1}{a_{12}}
```
```math
x_2 = -\frac{a_{21}}{a_{22}}x_1 + \frac{b_2}{a_{22}}
```

The two lines have the slopes $`-a_{11}/a_{12}`$ and $`-a_{21}/a_{22}`$, and the solution $`(x_1, x_2)`$ is their point of intersection.

### 3.2 Three Possibilities

There are exactly three possible geometric configurations for two lines in 2D:

1. **Unique solution**: The two lines intersect at exactly one point. The system has a single solution $`(x_1, x_2)`$. This is the typical case when the equations are independent and the determinant is nonzero.

2. **No solution**: The two lines are parallel (same slope, different intercepts). The system is **inconsistent**: no point satisfies both equations simultaneously. This occurs when the rows of $`[A]`$ are proportional but $`\{b\}`$ is not similarly proportional.

3. **Infinite solutions**: The two lines are identical (one equation is a scalar multiple of the other). Every point on the line is a solution. The system is **dependent**: the equations carry redundant information.

These three cases correspond precisely to the rank conditions of the augmented matrix $`[A|b]`$. In the unique-solution case, $`\text{rank}([A]) = \text{rank}([A|b]) = n`$. When the system is inconsistent, $`\text{rank}([A]) < \text{rank}([A|b])`$. When the system is dependent, $`\text{rank}([A]) = \text{rank}([A|b]) < n`$.

### 3.3 Column View

The same system can also be read by columns. With the column vectors $`\vec{a}_1`$ and $`\vec{a}_2`$ of $`[A]`$ and the right-hand side $`\vec{b}`$, the product $`[A]\{x\}`$ is a linear combination of the columns:

```math
\begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix} \iff x_1\underbrace{\begin{bmatrix} a_{11} \\ a_{21} \end{bmatrix}}_{\vec{a}_1} + x_2\underbrace{\begin{bmatrix} a_{12} \\ a_{22} \end{bmatrix}}_{\vec{a}_2} = \underbrace{\begin{bmatrix} b_1 \\ b_2 \end{bmatrix}}_{\vec{b}}
```

In this **column view**, solving the system means finding the multiples $`x_1`$ and $`x_2`$ for which the linear combination $`x_1\vec{a}_1 + x_2\vec{a}_2`$ of the column vectors equals $`\vec{b}`$. Geometrically, $`\vec{b}`$ is the diagonal of the parallelogram whose sides are $`x_1\vec{a}_1`$ and $`x_2\vec{a}_2`$. For the system $`2x_1 + x_2 = 5`$, $`x_1 - x_2 = 1`$, the solution $`x_1 = 2`$, $`x_2 = 1`$ combines the columns as follows:

```math
2\begin{bmatrix} 2 \\ 1 \end{bmatrix} + 1\begin{bmatrix} 1 \\ -1 \end{bmatrix} = \begin{bmatrix} 4 + 1 \\ 2 - 1 \end{bmatrix} = \begin{bmatrix} 5 \\ 1 \end{bmatrix}
```

Every $`\vec{b}`$ can be written as exactly one such combination if $`\vec{a}_1 \times \vec{a}_2 \neq 0`$, that is, if the two column vectors are not parallel. Treating the plane vectors as vectors in space with a zero third component, the cross product has only a third component, $`a_{11}a_{22} - a_{21}a_{12}`$, which is the determinant of $`[A]`$ (Section 4.1). When $`\vec{a}_1`$ and $`\vec{a}_2`$ are parallel, their combinations reach only the points of the line along them: if $`\vec{b}`$ lies on that line, there are infinitely many solutions, and otherwise there is no solution. For example, with $`\vec{a}_1 = (1, 2)`$ and $`\vec{a}_2 = (2, 4)`$, the right-hand side $`\vec{b} = (3, 6)`$ lies on the common line and $`x_1 + 2x_2 = 3`$ has infinitely many solutions, whereas $`\vec{b} = (3, 7)`$ does not and the system has no solution.

---

<br>

## 4. Determinant and Solvability

### 4.1 Determinant of a 2x2 System

For the $`2 \times 2`$ coefficient matrix:

```math
[A] = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix}
```

the **determinant** is:

```math
D = a_{11}a_{22} - a_{12}a_{21}
```

The determinant provides a single scalar that characterizes the solvability of the system.

The determinant arises directly from the row view of Section 3.1. When $`a_{12} \neq 0`$ and $`a_{22} \neq 0`$, the two lines intersect in exactly one point if and only if their slopes differ:

```math
\frac{a_{11}}{a_{12}} \neq \frac{a_{21}}{a_{22}} \iff \frac{a_{11}}{a_{12}} - \frac{a_{21}}{a_{22}} \neq 0 \iff a_{11}a_{22} - a_{12}a_{21} \neq 0
```

The last step multiplies by $`a_{12}a_{22} \neq 0`$. The final condition contains no division, and it characterizes a unique solution for every $`2 \times 2`$ system, including those in which one of the lines is vertical ($`a_{12} = 0`$ or $`a_{22} = 0`$). For a $`2 \times 2`$ matrix, the statements $`D \neq 0`$, the two rows of $`[A]`$ are linearly independent, $`[A]`$ is nonsingular, and $`[A]`$ is invertible are therefore equivalent, and each of them means that the system has exactly one solution.

### 4.2 Singular and Ill-Conditioned Systems

- **$`D = 0`$**: The system is **singular**, and no unique solution exists. The rows of $`[A]`$ are linearly dependent (one row is a scalar multiple of the other), and the system is either inconsistent or dependent.

- **$`D`$ very small relative to the coefficients**: The system is **ill-conditioned**. A unique solution technically exists, but it is extremely sensitive to small perturbations in the coefficients or the right-hand side. Tiny roundoff errors in the input data or during computation can produce wildly inaccurate results. Because the magnitude of the determinant changes when an equation is merely multiplied by a constant, the comparison is meaningful only after each equation is normalized.

> **[Linear Algebra]** The absolute value of the determinant is the area of the parallelogram (2D) or the volume of the parallelepiped (3D) formed by the row (or column) vectors. A zero determinant means the vectors are linearly dependent and the area or volume collapses to zero; the redundant rows of the coefficient matrix in the singular systems above are exactly this case.

Geometrically, ill-conditioning for a $`2 \times 2`$ system means the two lines intersect at a very shallow angle. In such cases, small changes in the coefficients (shifting a line slightly) can move the intersection point dramatically. This sensitivity is quantified by the **condition number**, which in the 2-norm is the ratio of the largest to the smallest singular value of $`[A]`$. The condition number is discussed further in Chapter 11.

---

<br>

## 5. Cramer's Rule

### 5.1 Derivation for 2x2

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

### 5.2 General Form

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

### 5.3 Limitation

Cramer's rule requires evaluating $`n + 1`$ determinants, and computing each one by cofactor expansion takes about $`n!`$ operations. The total cost is therefore about $`(n+1) \cdot n! = (n+1)!`$ operations, which grows factorially, that is, astronomically fast. The table below lists the operation count for one determinant:

| $`n`$ | $`n!`$ operations |
|:----|:----------------|
| 3 | 6 |
| 5 | 120 |
| 10 | 3,628,800 |
| 20 | $`\approx 2.4 \times 10^{18}`$ |

*Table 2. Operation count for one determinant by cofactor expansion*

For $`n > 3`$, Cramer's rule is computationally impractical. Even a computer performing $`10^9`$ operations per second would need approximately **77 years** to evaluate a single $`20 \times 20`$ determinant ($`20! \approx 2.4 \times 10^{18}`$ operations), and Cramer's rule needs 21 such determinants. This motivates the development of more efficient algorithms.

---

<br>

## 6. Naive Gauss Elimination

Gauss elimination is the workhorse algorithm for solving linear systems. The method consists of two phases. **Forward elimination** reduces the system to upper triangular form, and **back substitution** solves the triangular system from bottom to top.

### 6.1 Forward Elimination

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

**Continue through column $`n-1`$:** At each step $`i`$, use $`a_{ii}^{(i-1)}`$ as the pivot element and eliminate all entries below it. For a $`3 \times 3`$ augmented matrix, the process proceeds as in Figure 3.

```text
[ * * * | * ]       [ * * * | * ]       [ * * * | * ]
[ * * * | * ]  -->  [ 0 * * | * ]  -->  [ 0 * * | * ]
[ * * * | * ]       [ 0 * * | * ]       [ 0 0 * | * ]
```

*Figure 3. A 3x3 augmented matrix reaching upper triangular form through the two steps of forward elimination*

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

### 6.2 Back Substitution

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

## 7. Computational Cost

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

*Table 3. Dominant terms of the operation count in each phase of Gauss elimination*

This is a dramatic improvement over Cramer's rule. For $`n = 20`$:

- Cramer's rule: about $`2.4 \times 10^{18}`$ operations for just one determinant.
- Gauss elimination: about $`(2/3)(20)^3 = 5333`$ operations.

---

<br>

## 8. Example: 3x3 System (Complete Step-by-Step)

### 8.1 Problem Setup

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

### 8.2 Forward Elimination: Step 1

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

### 8.3 Forward Elimination: Step 2

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

### 8.4 Back Substitution

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

> **Note:** $`x_3`$ comes out as 7.00003 instead of the exact value 7 because of roundoff error accumulated by rounding to six significant digits at each step. The exact solution of this system is $`x_1 = 3`$, $`x_2 = -2.5`$, $`x_3 = 7`$, and the computed solution is very close to it.

---

<br>

## 9. Pitfalls of Naive Gauss Elimination

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

**2. Near-zero pivots amplify roundoff error:**

When a pivot element is very small but nonzero, the computed factor $`a_{ji}/a_{ii}`$ becomes very large. Multiplying row $`i`$ by this large factor and subtracting magnifies any roundoff errors present in row $`i`$, potentially corrupting the result severely.

**3. Large factor values magnify errors:**

Related to the previous point, when the factor is much larger than 1, the operation $`\text{Row}_j - \text{factor} \times \text{Row}_i`$ subtracts a large multiple of one row from another. This can lead to catastrophic cancellation and significant loss of significant digits.

These three pitfalls all stem from the same root cause: **the pivot element is too small relative to other entries in the column**. The remedy is **partial pivoting**. At each step, swap rows so that the available entry with the largest absolute value in the current column becomes the pivot. This strategy, covered in Chapter 9, effectively mitigates all three problems and is used in virtually all production-quality implementations of Gauss elimination.

---

<br>

## Summary

| Concept | Key Summary |
|:------|:----------|
| System form | It has the form $`[A]\{x\} = \{b\}`$ and arises in spring, truss, and circuit analysis. |
| Matrix Types | Square matrices are classified as symmetric ($`a_{ij} = a_{ji}`$), diagonal, identity, upper or lower triangular, banded, and tridiagonal by the pattern of their elements. |
| Matrix Multiplication | An $`m \times n`$ matrix times an $`n \times p`$ matrix gives an $`m \times p`$ matrix with $`c_{ij} = \sum_{k} a_{ik}b_{kj}`$, and the product is associative and distributive but not commutative in general. |
| Transpose and Inverse | The transpose satisfies $`a^T_{ji} = a_{ij}`$, and the inverse satisfies $`[A]^{-1}[A] = [A][A]^{-1} = [I]`$ and exists when $`[A]`$ is square and nonsingular. |
| Permutation and Augmented Matrices | A permutation matrix interchanges rows in $`[P][A]`$ and columns in $`[A][P]`$, and an augmented matrix appends columns such as $`\{b\}`$ or $`[I]`$ to a matrix. |
| Matrix Form | A linear system is $`[A]\{x\} = \{b\}`$ with $`\{x\} = [A]^{-1}\{b\}`$ for invertible $`[A]`$, and an $`m \times n`$ system is overdetermined if $`m > n`$ and underdetermined if $`m < n`$. |
| Stiffness Matrix | The stiffness matrix of a spring-mass system is symmetric positive definite, and $`\frac{1}{2}\{x\}^T[K]\{x\}`$ is the potential energy of the springs. |
| Python Solution | `np.linalg.solve(A, b)` solves the system with fewer operations than `np.linalg.inv(A) @ b`. |
| Graphical (2x2) | Two lines meeting at one point give a unique solution, parallel lines give no solution, and identical lines give infinitely many solutions. |
| Column View | It finds the combination $`x_1\vec{a}_1 + x_2\vec{a}_2`$ of the columns that equals $`\vec{b}`$, which is unique when the columns are not parallel, that is, when $`\det(A) \neq 0`$. |
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

1. **Matrix multiplication:** For $`[A] = \begin{bmatrix} 1 & 2 \\ 3 & 4 \end{bmatrix}`$ and $`[B] = \begin{bmatrix} 2 & 0 \\ 1 & 3 \end{bmatrix}`$, compute $`[A][B]`$ and $`[B][A]`$. Is the multiplication commutative for these two matrices?

   <details>
   <summary>Answer</summary>

   $`[A][B] = \begin{bmatrix} (1)(2) + (2)(1) & (1)(0) + (2)(3) \\ (3)(2) + (4)(1) & (3)(0) + (4)(3) \end{bmatrix} = \begin{bmatrix} 4 & 6 \\ 10 & 12 \end{bmatrix}`$ and $`[B][A] = \begin{bmatrix} (2)(1) + (0)(3) & (2)(2) + (0)(4) \\ (1)(1) + (3)(3) & (1)(2) + (3)(4) \end{bmatrix} = \begin{bmatrix} 2 & 4 \\ 10 & 14 \end{bmatrix}`$. The two products differ, so the multiplication is not commutative for these matrices, as is usual for matrices.

   </details>

2. **Stiffness matrix:** For the spring-mass system of Section 2.11 with $`k_1 = k_2 = 100`$ N/m, write $`[K]`$, find $`[K]^{-1}`$, and compute the displacements for $`f_1 = 0`$ and $`f_2 = 50`$ N. Also evaluate $`\{v\}^T[K]\{v\}`$ for $`\{v\} = \lfloor 1 \quad 1 \rfloor^T`$.

   <details>
   <summary>Answer</summary>

   $`[K] = \begin{bmatrix} 200 & -100 \\ -100 & 100 \end{bmatrix}`$ and $`\det(K) = (200)(100) - (-100)(-100) = 10000`$, so $`[K]^{-1} = \frac{1}{10000}\begin{bmatrix} 100 & 100 \\ 100 & 200 \end{bmatrix} = \begin{bmatrix} 0.01 & 0.01 \\ 0.01 & 0.02 \end{bmatrix}`$. Then $`x_1 = (0.01)(0) + (0.01)(50) = 0.5`$ m and $`x_2 = (0.01)(0) + (0.02)(50) = 1.0`$ m, which agrees with the fact that both springs carry 50 N and each stretches by 0.5 m. Finally, $`\{v\}^T[K]\{v\} = k_1v_1^2 + k_2(v_1 - v_2)^2 = (100)(1)^2 + (100)(0)^2 = 100 > 0`$, as positive definiteness requires.

   </details>

3. **Cramer's rule:** Solve $`2x_1 + x_2 = 5`$, $`x_1 - x_2 = 1`$ by Cramer's rule.

   <details>
   <summary>Answer</summary>

   $`D = (2)(-1) - (1)(1) = -3`$. Then $`x_1 = \frac{a_{22}b_1 - a_{12}b_2}{D} = \frac{(-1)(5) - (1)(1)}{-3} = 2`$ and $`x_2 = \frac{a_{11}b_2 - a_{21}b_1}{D} = \frac{(2)(1) - (1)(5)}{-3} = 1`$. The two lines meet at the point $`(2, 1)`$.

   </details>

4. **Solvability:** Find the determinant of $`x_1 + 2x_2 = 3`$, $`2x_1 + 4x_2 = 7`$ and decide how many solutions the system has. What happens if the right-hand side of the second equation is 6?

   <details>
   <summary>Answer</summary>

   $`D = (1)(4) - (2)(2) = 0`$, so the system is singular. The coefficient rows are proportional (by a factor of 2), but $`7 \neq 2 \times 3`$, so the lines are parallel and there is no solution (inconsistent). If the right-hand side is 6, the second equation is twice the first, the lines coincide, and there are infinitely many solutions (dependent).

   </details>

5. **Naive Gauss elimination:** Solve the following system by naive Gauss elimination: $`2x_1 + x_2 - x_3 = 8`$, $`-3x_1 - x_2 + 2x_3 = -11`$, $`-2x_1 + x_2 + 2x_3 = -3`$.

   <details>
   <summary>Answer</summary>

   Eliminating column 1 with factors $`-3/2 = -1.5`$ and $`-2/2 = -1`$ turns row 2 into $`0.5x_2 + 0.5x_3 = 1`$ and row 3 into $`2x_2 + x_3 = 5`$. Eliminating column 2 with the factor $`2/0.5 = 4`$ turns row 3 into $`-x_3 = 1`$. Back substitution gives $`x_3 = -1`$, $`x_2 = (1 - 0.5(-1))/0.5 = 3`$, and $`x_1 = (8 - x_2 + x_3)/2 = (8 - 3 - 1)/2 = 2`$.

   </details>

6. **Growth of the cost:** If $`n`$ grows tenfold from 100 to 1000, by roughly what factor does the operation count of Gauss elimination grow? Also give the operation count in each case.

   <details>
   <summary>Answer</summary>

   The operation count is about $`(2/3)n^3`$, so it is about $`6.7 \times 10^5`$ for $`n = 100`$ and about $`6.7 \times 10^8`$ for $`n = 1000`$. A tenfold increase in $`n`$ multiplies the cost by $`10^3 = 1000`$. The $`n^2`$ term of back substitution, $`10^4`$ and $`10^6`$ respectively, is small compared with forward elimination.

   </details>
