# HW3 Assignment Handout — LU Factorization and Matrix Inverse

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-16

This is the Week 4 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab04-LU.ipynb` distributed by the instructor and run it. The notebook has five parts: LU factorization, multiple right-hand sides, Cholesky factorization, the matrix inverse, and the condition number with ill-conditioning. There are six blanks in total, each between two comment lines marked `COMPLETE`. The first code cell of the notebook imports `numpy` as `np`, `scipy` as `sc`, `matplotlib.pyplot` as `plt`, and `time`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The distributed blank notebook is [Lab04-LU.ipynb](Lab04-LU.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blanks and running it.

---

<br>

## Problem 1. Doolittle's Algorithm (Notebook Section 1.2)

Gaussian elimination turns $`A`$ into an upper triangular matrix $`U`$ through row operations that subtract a multiple of an upper row from a lower row. The multipliers $`\ell_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}`$ used there form a lower triangular matrix $`L`$ with 1s on the diagonal, and $`A = LU`$ holds. Once $`L`$ and $`U`$ are stored, $`A\mathbf{x} = \mathbf{b}`$ for a new right-hand side is solved with two $`O(n^2)`$ triangular solves instead of repeating the $`O(n^3)`$ factorization.

Doolittle's method computes $`L`$ and $`U`$ directly, without row operations. With the convention $`\ell_{ii} = 1`$ for all $`i`$, matching $`A = LU`$ entry by entry gives the following for $`k = 0, 1, \ldots, n-1`$.

```math
u_{kj} = a_{kj} - \sum_{m=0}^{k-1} \ell_{km}\, u_{mj} \quad (j = k, \ldots, n-1), \qquad
\ell_{ik} = \frac{a_{ik} - \displaystyle\sum_{m=0}^{k-1} \ell_{im}\, u_{mk}}{u_{kk}} \quad (i = k+1, \ldots, n-1)
```

Both sums are empty when $`k = 0`$, so $`u_{0j} = a_{0j}`$ and $`\ell_{i0} = a_{i0}/u_{00}`$. The notebook shows the process step by step on the following matrix.

```math
A = \begin{bmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{bmatrix}
```

- **$`k = 0`$:** Row 0 of $`U`$ is $`u_{00} = 2`$, $`u_{01} = 1`$, $`u_{02} = 1`$, and column 0 of $`L`$ is $`\ell_{10} = 4/2 = 2`$, $`\ell_{20} = 8/2 = 4`$.
- **$`k = 1`$:** Row 1 of $`U`$ is $`u_{11} = 3 - 2 \cdot 1 = 1`$, $`u_{12} = 3 - 2 \cdot 1 = 1`$, and column 1 of $`L`$ is $`\ell_{21} = (7 - 4 \cdot 1)/1 = 3`$.
- **$`k = 2`$:** Row 2 of $`U`$ is $`u_{22} = 9 - 4 \cdot 1 - 3 \cdot 1 = 2`$.

```math
L = \begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 4 & 3 & 1 \end{bmatrix},\quad
U = \begin{bmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{bmatrix}
```

Fill in the blanks of the function `lu_doolittle` to compute the rows of $`U`$ and the columns of $`L`$.

```python
def lu_doolittle(A):
    """
    Compute the LU factorization A = L @ U using Doolittle's algorithm.
    L has 1s on the diagonal (unit lower triangular).
    No pivoting is performed.

    Parameters
    ----------
    A : (n, n) array-like

    Returns
    -------
    L : (n, n) ndarray  — unit lower triangular
    U : (n, n) ndarray  — upper triangular
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    L = np.eye(n)          # start with identity (diagonal = 1)
    U = np.zeros((n, n))

    for k in range(n):
        # ── Row k of U: U[k, j] = A[k, j] - sum_{m<k} L[k,m]*U[m,j]
        for j in range(k, n):
            # --- COMPLETE ---
            U[k, j] =
            # --- COMPLETE ---

        # ── Column k of L: L[i, k] = (A[i,k] - sum_{m<k} L[i,m]*U[m,k]) / U[k,k]
        for i in range(k + 1, n):
            # --- COMPLETE ---
            L[i, k] =
            # --- COMPLETE ---

    return L, U
```

The cell that follows compares the result of the function with SciPy's `sc.linalg.lu`.

```python
# Verify on the 3x3 example traced above
A_test = np.array([[8, 7, 9],
                   [2, 1, 1],
                   [4, 3, 3]], dtype=float)

L, U = lu_doolittle(A_test)

print('L =\n', L)
print('U =\n', U)
print('LU =\n', L@U)
# print('\nSanity check  L @ U - A =\n', L @ U - A_test)

# Compare with scipy 
P, L_sc, U_sc = sc.linalg.lu(A_test)
print('L_sc =\n', L_sc)
print('U_sc =\n', U_sc)
print('L_scU_sc =\n', P@L_sc@U_sc)
```

## Problem 2. Solving with L and U (Notebook Section 1.4)

With $`A = LU`$, the system $`A\mathbf{x} = \mathbf{b}`$ splits into two triangular systems. Without pivoting the right-hand side stays $`\mathbf{b}`$, and if pivoting was used and the matrix was factored as $`A = PLU`$ following SciPy's convention, the right-hand side is $`\mathbf{b}' = P^T \mathbf{b}`$.

**Step 1 (forward substitution):** Solve $`L\mathbf{d} = \mathbf{b}`$.

```math
d_0 = b_0, \qquad d_i = b_i - \sum_{j=0}^{i-1} \ell_{ij}\, d_j, \quad i = 1, \ldots, n-1
```

**Step 2 (backward substitution):** Solve $`U\mathbf{x} = \mathbf{d}`$.

```math
x_{n-1} = \frac{d_{n-1}}{u_{n-1,n-1}}, \qquad
x_i = \frac{d_i - \displaystyle\sum_{j=i+1}^{n-1} u_{ij}\, x_j}{u_{ii}}, \quad i = n-2, \ldots, 0
```

Fill in the blanks of the functions `forward_sub` and `backward_sub`. The function `lu_solve` chains SciPy's LU factorization with the two functions.

```python
def forward_sub(L, b):
    """
    Solve L d = b by forward substitution.
    L must be unit lower triangular (diagonal = 1).

    Returns
    -------
    d : (n,) ndarray
    """
    n = len(b)
    d = np.zeros(n)
    for i in range(n):
        # --- COMPLETE ---
        d[i] =
        # --- COMPLETE ---
    return d


def backward_sub(U, d):
    """
    Solve U x = d by backward substitution.
    U must be upper triangular.

    Returns
    -------
    x : (n,) ndarray
    """
    n = len(d)
    x = np.zeros(n)
    x[-1] = d[-1] / U[-1, -1]
    for i in range(n - 2, -1, -1):
        # --- COMPLETE ---
        x[i] =
        # --- COMPLETE ---
    return x


def lu_solve(A, b):
    """Solve Ax = b using scipy's LU (with pivoting) + triangular substitutions."""
    P, L, U = sc.linalg.lu(A)
    b_perm = P.T @ b          # apply permutation: P^T b
    d = forward_sub(L, b_perm)
    x = backward_sub(U, d)
    return x
```

The cell that follows tests the functions on a $`3 \times 3`$ system with a known solution.

```python
# Test on a 3x3 system with known solution
A_sys = np.array([[ 3, -0.1, -0.2],
                  [ 0.1,  7, -0.3],
                  [ 0.3, -0.2, 10]], dtype=float)
b_sys = np.array([7.85, -19.3, 71.4])

x_lu    = lu_solve(A_sys, b_sys)
x_ref   = np.linalg.solve(A_sys, b_sys)

print('Solution via LU solve :', np.round(x_lu, 6))
print('Solution via np.linalg:', np.round(x_ref, 6))
print('Max difference        :', np.max(np.abs(x_lu - x_ref)))
print('Residual ||Ax - b||   :', np.linalg.norm(A_sys @ x_lu - b_sys))
```

## Problem 3. Spring-Mass System under Multiple Load Cases (Notebook Section 2.2)

Suppose $`m`$ systems with the same matrix $`A`$ and different right-hand sides, $`A\mathbf{x}_1 = \mathbf{b}_1, \ldots, A\mathbf{x}_m = \mathbf{b}_m`$, must be solved. Gaussian elimination for each system costs $`m \times O(n^3/3)`$, while factoring once and solving only the triangular systems costs $`O(n^3/3) + m \times O(n^2)`$.

The notebook reuses the three-mass spring system of Lab 03 (HW2).

```math
K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix}
```

The goal is to find the displacements under three load cases, gravity, wind, and the combined loading of both, without refactoring $`K`$ each time. Fill in the blanks so that the cell factors $`K`$ once with `scipy.linalg.lu_factor` and solves all three load cases with `scipy.linalg.lu_solve`.

```python
# Spring constants (N/m)
k1, k2, k3, k4 = 40.0, 30.0, 20.0, 10.0

K_spring = np.array([
    [k1+k2,   -k2,      0],
    [  -k2, k2+k3,    -k3],
    [    0,   -k3,  k3+k4]
])

# Three load cases (N): gravity, wind, combined
F_gravity  = np.array([10.0, 10.0, 10.0])
F_wind     = np.array([50.0, 20.0,  0.0])
F_combined = F_gravity + F_wind

# --- COMPLETE: factor K once and solve for all three load cases ---
lu_K = sc.linalg.lu_factor(K_spring)                # factor K

x_gravity  =
x_wind     =
x_combined =
# --- COMPLETE ---

labels = ['Mass 1', 'Mass 2', 'Mass 3']
cases  = [('Gravity', x_gravity), ('Wind', x_wind), ('Combined', x_combined)]

print(f'{"Load Case":<12}  {"x1 (m)":>10}  {"x2 (m)":>10}  {"x3 (m)":>10}')
print('-' * 48)
for name, x in cases:
    print(f'{name:<12}  {x[0]:>10.5f}  {x[1]:>10.5f}  {x[2]:>10.5f}')

# Verify superposition: x_combined should equal x_gravity + x_wind
print(f'\nSuperposition check ||x_comb - (x_grav + x_wind)||: '
      f'{np.linalg.norm(x_combined - (x_gravity + x_wind)):.2e}')
```

The distributed notebook shows the following output for this cell.

```text
Load Case         x1 (m)      x2 (m)      x3 (m)
------------------------------------------------
Gravity          0.52000     0.88000     0.92000
Wind             1.46000     1.74000     1.16000
Combined         1.98000     2.62000     2.08000

Superposition check ||x_comb - (x_grav + x_wind)||: 4.97e-16
```

## Problem 4. Solving with the Cholesky Factorization (Notebook Section 3.2)

A matrix $`A`$ is called **symmetric positive definite (SPD)** if $`A = A^T`$ and $`\mathbf{v}^T A \mathbf{v} > 0`$ for every nonzero vector $`\mathbf{v}`$. This is equivalent to all eigenvalues being positive and to all leading principal minors being positive, and SPD matrices arise as stiffness matrices in structural mechanics, as the normal equations $`A^T A`$ of least squares, and as covariance matrices in statistics.

An SPD matrix factors as $`A = LL^T`$ with a lower triangular $`L`$ whose diagonal entries are positive (the Cholesky factorization). Compared with general LU, it applies only to SPD matrices, but the factorization costs about $`n^3/6`$, half as much, stores a single matrix, and is stable without pivoting. The formulas, column by column for $`k = 0, \ldots, n-1`$, are as follows.

```math
\ell_{kk} = \sqrt{a_{kk} - \sum_{m=0}^{k-1} \ell_{km}^2}, \qquad
\ell_{ik} = \frac{a_{ik} - \displaystyle\sum_{m=0}^{k-1} \ell_{im}\,\ell_{km}}{\ell_{kk}}, \quad i > k
```

Fill in the blanks of the function `cholesky_solve`. As the comments instruct, use `forward_sub` defined in Problem 2, and note that $`L^T`$ is upper triangular.

```python
def cholesky_solve(A, b):
    """
    Solve A x = b for an SPD matrix A using Cholesky factorization.

    Algorithm:
      1. Factor: A = L L^T
      2. Forward sub: L d = b
      3. Backward sub: L^T x = d

    Returns
    -------
    x : (n,) ndarray
    """
    # Step 1: Cholesky factorization
    L = sc.linalg.cholesky(A, lower=True)

    # Step 2: Forward substitution  L d = b
    # --- COMPLETE: use forward_sub defined in Section 1 ---
    d =
    # --- COMPLETE ---

    # Step 3: Backward substitution  L^T x = d
    # --- COMPLETE: note that L^T is upper triangular ---
    x =
    # --- COMPLETE ---

    return x
```

The cell that follows applies the function to an FEM spring network whose stiffness matrix is SPD and compares the time of SciPy's SPD solver with that of a general LU solve. The function `is_spd`, defined in an earlier cell, decides whether a matrix is SPD from its symmetry and the signs of its eigenvalues.

```python
# Application: FEM spring network with SPD stiffness matrix
# 4-node bar element: K x = f, where K is SPD
K_fem = np.array([
    [ 2, -1,  0,  0],
    [-1,  2, -1,  0],
    [ 0, -1,  2, -1],
    [ 0,  0, -1,  1]
], dtype=float) * 1e4   # stiffness in N/m

f_ext = np.array([0.0, 0.0, 500.0, 1000.0])   # external forces in N

print('Is K_fem SPD?', is_spd(K_fem))

x_chol = cholesky_solve(K_fem, f_ext)
x_ref  = np.linalg.solve(K_fem, f_ext)

print('\nDisplacements (m):')
for i, xi in enumerate(x_chol):
    print(f'  x{i+1} = {xi*1000:.4f} mm')
print('\nMax error vs np.linalg.solve:', np.max(np.abs(x_chol - x_ref)))

# Efficiency: Cholesky vs LU for SPD systems
n = 500
# Build a large SPD matrix: A = R^T R for random R
R = np.random.randn(n, n)
A_large = R.T @ R + n * np.eye(n)
b_large = np.random.randn(n)

t0 = time.perf_counter(); sc.linalg.solve(A_large, b_large, assume_a='pos'); t_chol = time.perf_counter() - t0
t0 = time.perf_counter(); np.linalg.solve(A_large, b_large);                 t_lu   = time.perf_counter() - t0

print(f'\nSolve time (n={n}):  Cholesky = {t_chol:.4f}s,  LU = {t_lu:.4f}s,  speedup = {t_lu/t_chol:.1f}x')
```

## Problem 5. Matrix Inverse via LU (Notebook Section 4.2)

The inverse of a square matrix $`A`$ is the unique matrix $`A^{-1}`$ with $`AA^{-1} = A^{-1}A = I`$, and it exists only when $`A`$ is nonsingular ($`\det(A) \neq 0`$). Since $`AX = I`$ for $`X = A^{-1}`$, the inverse can be found by solving $`n`$ systems with the $`j`$-th column $`\mathbf{e}_j`$ of the identity matrix.

```math
A \mathbf{x}_j = \mathbf{e}_j, \quad j = 0, 1, \ldots, n-1
```

The algorithm is as follows.

1. Factor $`A = LU`$ once ($`O(n^3/3)`$).
2. For each $`j = 0, \ldots, n-1`$, solve $`LU\mathbf{x}_j = \mathbf{e}_j`$ ($`O(n^2)`$ per column).
3. Assemble $`A^{-1} = [\mathbf{x}_0 \mid \mathbf{x}_1 \mid \cdots \mid \mathbf{x}_{n-1}]`$.

Fill in the blank of the function `lu_inverse` to compute each column of $`A^{-1}`$.

```python
def lu_inverse(A):
    """
    Compute A^{-1} by solving A X = I column by column using LU factorization.

    Returns
    -------
    A_inv : (n, n) ndarray
    """
    n = A.shape[0]
    I = np.eye(n)

    # Factor once
    lu_obj = sc.linalg.lu_factor(A)

    # --- COMPLETE: solve A x_j = e_j for each column j of I ---
    A_inv = np.zeros((n, n))
    for j in range(n):
        A_inv[:, j] =
    # --- COMPLETE ---

    return A_inv
```

The cell that follows tests the function on the matrix `A_sys` of Problem 2 and shows what happens when the inverse of a singular matrix is requested.

```python
# Test lu_inverse
A_inv_lu  = lu_inverse(A_sys)
A_inv_ref = np.linalg.inv(A_sys)

print('A^{-1} via lu_inverse =\n', np.round(A_inv_lu, 6))
print('\nA^{-1} via np.linalg.inv =\n', np.round(A_inv_ref, 6))
print('\nMax error:', np.max(np.abs(A_inv_lu - A_inv_ref)))
print('Sanity check  A @ A_inv - I  max =', np.max(np.abs(A_sys @ A_inv_lu - np.eye(3))))

# Singular matrix — inverse does not exist
A_sing = np.array([[1.0, 2.0], [2.0, 4.0]])   # rank 1: row 2 = 2 * row 1
print('\nDeterminant of singular matrix:', np.linalg.det(A_sing))
try:
    np.linalg.inv(A_sing)
except np.linalg.LinAlgError as e:
    print('Error:', e)
```

## Problem 6. Diagnosing an Ill-Conditioned System (Notebook Section 5.4)

The condition number of a matrix $`A`$ is $`\kappa(A) = \lVert A \rVert \lVert A^{-1} \rVert`$, and it measures how much a small perturbation $`\delta\mathbf{b}`$ in the right-hand side is amplified into the error $`\delta\mathbf{x}`$ of the solution.

```math
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa(A)\,\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}
```

| $`\kappa(A)`$ | Interpretation |
|:-----|:-----|
| $`\approx 1`$ | Well-conditioned; full floating-point accuracy |
| $`10^6`$ | The solution may lose about 6 significant digits |
| $`\sim 1/\varepsilon_{\text{mach}} \approx 10^{16}`$ | Effectively singular; the solution is unreliable |

*Table 1. Interpretation of the size of the condition number*

As a rule of thumb, if $`\kappa(A) \approx 10^k`$, about $`k`$ significant digits are lost. In NumPy, `np.linalg.cond(A)` returns the 2-norm condition number by default, which is the ratio of the largest to the smallest singular value, $`\kappa_2(A) = \sigma_{\max} / \sigma_{\min}`$.

The exercise reads as follows. A structural engineer reports that two different programs give very different displacements for the same frame stiffness problem. Diagnose the matrix given below using the condition number and suggest a remedy. Fill in the blanks of the next cell to compute the condition number, the number of significant digits lost, and the solution.

```python
# Suspicious stiffness matrix (values in N/m)
K_susp = np.array([
    [1e6,  1e6 - 1,       0],
    [1e6 - 1, 1e6,        0],
    [0,        0,      1e-2]
], dtype=float)

f_susp = np.array([1000.0, 1000.0, 0.001])

# --- COMPLETE ---
# 1. Compute the condition number
kappa =
# 2. Estimate how many significant digits are lost
digits_lost =
# 3. Solve the system
x_susp =
# --- COMPLETE ---

print(f'Condition number kappa = {kappa:.2e}')
print(f'Significant digits lost ≈ {digits_lost}')
print(f'Solution x = {x_susp}')
print(f'Residual ||Kx - f|| = {np.linalg.norm(K_susp @ x_susp - f_susp):.2e}')
print()
print('Diagnosis: The upper-left 2x2 block has nearly identical rows')
print('(rows differ by only 1 N/m out of 1e6), making the system near-singular.')
print('The solution in the first two DOFs is unreliable.')
```

The distributed notebook shows the following output for this cell.

```text
Condition number kappa = 2.00e+08
Significant digits lost ≈ 8
Solution x = [0.0005 0.0005 0.1   ]
Residual ||Kx - f|| = 0.00e+00

Diagnosis: The upper-left 2x2 block has nearly identical rows
(rows differ by only 1 N/m out of 1e6), making the system near-singular.
The solution in the first two DOFs is unreliable.
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab04-LU_Completed.ipynb) for the full code.
