# HW2 과제 안내문 — 연립방정식과 가우스 소거법

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-16

수치해석 3주차 과제이다. 교수자가 배포한 실습 노트북 `Lab03-LinearSystem.ipynb`의 빈칸을 채워 실행한다. 노트북은 $`A\mathbf{x} = \mathbf{b}`$의 설정과 풀이, 단순 가우스 소거법, 부분 피봇팅, 삼중대각 시스템과 토마스 알고리즘의 네 부분으로 이루어지며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이에 있고 각 부분에 하나씩 모두 4곳이다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고 `time`을 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 빈칸 옆에 남은 주석은 배포된 노트북에 적혀 있던 설명이다. 배포된 빈 노트북은 [Lab03-LinearSystem.ipynb](Lab03-LinearSystem.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

---

<br>

## 문제 1. 질량-스프링 시스템의 우변 벡터 (노트북 1.3절)

고정된 두 벽 사이에 질량 세 개가 스프링 네 개로 연결되어 있다.

```text
Wall ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ Wall
              ↓F₁         ↓F₂         ↓F₃
```

*그림 1. 두 벽 사이에 스프링 네 개로 연결된 질량 세 개와 각 질량에 작용하는 힘*

오른쪽을 $`x`$의 양의 방향으로 하여 각 질량에서 힘의 평형을 적용하면 다음 식을 얻는다.

```math
\begin{aligned}
m_1: &\quad (k_1 + k_2)\,x_1 - k_2\,x_2 = F_1 \\
m_2: &\quad -k_2\,x_1 + (k_2 + k_3)\,x_2 - k_3\,x_3 = F_2 \\
m_3: &\quad -k_3\,x_2 + (k_3 + k_4)\,x_3 = F_3
\end{aligned}
```

이를 행렬 형태 $`A\mathbf{x} = \mathbf{b}`$로 쓰면 다음과 같다.

```math
\begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} =
\begin{bmatrix} F_1 \\ F_2 \\ F_3 \end{bmatrix}
```

행렬 $`A`$는 대칭이고 삼중대각이며, 이 구조는 노트북 4절에서 중요해진다. 다음 셀의 빈칸을 채워 우변 벡터 $`\mathbf{b} = [F_1, F_2, F_3]`$를 만들어라.

```python
# Spring constants (N/m)
k1, k2, k3, k4 = 40.0, 30.0, 20.0, 10.0

# Applied forces (N)
F1, F2, F3 = 20.0, 40.0, 30.0

# Coefficient matrix (derived from equilibrium equations above)
A = np.array([
    [k1+k2,   -k2,      0],
    [  -k2, k2+k3,    -k3],
    [    0,   -k3,  k3+k4]
])

# --- COMPLETE: assemble the right-hand side vector b = [F1, F2, F3] ---
b = np.array([ ,  ,  ])
# --- COMPLETE ---

print('Coefficient matrix A:')
print(A)
print('\nRight-hand side b:', b)
```

이어지는 셀은 NumPy로 시스템을 풀고 잔차로 해를 검증한다.

```python
# Solve with NumPy
x = np.linalg.solve(A, b)

print('Displacements:')
for i, xi in enumerate(x):
    print(f'  x{i+1} = {xi:.6f} m')

# Verify: compute the residual r = A*x - b (should be ~0)
residual = A @ x - b
print(f'\nResidual A*x - b = {residual}')
print(f'Max |residual|   = {np.max(np.abs(residual)):.2e}  (machine precision: ~2e-16)')
```

## 문제 2. 단순 가우스 소거법 (노트북 2.4절)

가우스 소거법은 기본 행 연산으로 $`A\mathbf{x} = \mathbf{b}`$를 동치인 상삼각 시스템으로 바꾼 뒤 후진 대입으로 푼다. 먼저 $`A`$와 $`\mathbf{b}`$를 확대 행렬 $`[A \mid \mathbf{b}]`$로 합친다.

```math
\left[\begin{array}{ccc|c}
a_{11} & a_{12} & a_{13} & b_1 \\
a_{21} & a_{22} & a_{23} & b_2 \\
a_{31} & a_{32} & a_{33} & b_3
\end{array}\right]
```

**전진 소거:** 각 피봇 열 $`k = 0, 1, \ldots, n-2`$와 그 아래의 각 행 $`i = k+1, \ldots, n-1`$에 대해 다음을 수행한다.

1. 소거 계수 $`f_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}`$를 계산한다.
2. 행 $`i`$를 $`\text{row}_i \leftarrow \text{row}_i - f_{ik} \cdot \text{row}_k`$로 갱신한다.

이렇게 하면 $`(i, k)`$ 성분이 0이 되어 상삼각 행렬이 남는다.

**후진 대입:** 마지막 식에서 시작해 위로 올라간다.

```math
x_{n-1} = \frac{b_{n-1}^{(n)}}{a_{n-1,n-1}^{(n)}}, \qquad
x_i = \frac{b_i^{(i+1)} - \displaystyle\sum_{j=i+1}^{n-1} a_{ij}^{(i+1)}\, x_j}{a_{ii}^{(i+1)}}, \quad i = n-2, \ldots, 0
```

노트북은 다음 시스템으로 과정을 단계별로 추적한다. 참해는 $`(x_1, x_2, x_3) = (2, 3, -1)`$이다.

```math
\begin{bmatrix} 2 & 1 & -1 \\ -3 & -1 & 2 \\ -2 & 1 & 2 \end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} =
\begin{bmatrix} 8 \\ -11 \\ -3 \end{bmatrix}
```

- **1단계 (피봇 열 $`k = 0`$):** $`f_{10} = -3/2`$, $`f_{20} = -2/2 = -1`$로 0행 아래를 소거하면 확대 행렬의 1, 2행은 $`[0, 1/2, 1/2 \mid 1]`$, $`[0, 2, 1 \mid 5]`$가 된다.
- **2단계 (피봇 열 $`k = 1`$):** $`f_{21} = 2/(1/2) = 4`$로 1행 아래를 소거하면 2행은 $`[0, 0, -1 \mid 1]`$이 된다.
- **후진 대입:** $`x_3 = 1/(-1) = -1`$, $`x_2 = (1 - (1/2)(-1))/(1/2) = 3`$, $`x_1 = (8 - (1)(3) - (-1)(-1))/2 = 2`$이다.

함수 `gauss_naive`의 빈칸을 채워 전진 소거와 후진 대입을 완성하라.

```python
def gauss_naive(A, b):
    """
    Solve Ax = b using Naive Gaussian Elimination (no pivoting).

    Parameters
    ----------
    A : (n, n) array-like  — coefficient matrix
    b : (n,) array-like    — right-hand side vector

    Returns
    -------
    x   : (n,) ndarray — solution vector
    Aug : (n, n+1) ndarray — augmented matrix after elimination
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).reshape(-1, 1)
    n  = A.shape[0]
    nb = n + 1                         # number of columns in augmented matrix

    # Build augmented matrix [A | b]
    Aug = np.hstack((A, b))

    # ── 1. Forward Elimination ────────────────────────────────────────────────
    for k in range(n - 1):             # pivot column k = 0, 1, ..., n-2
        for i in range(k + 1, n):     # rows to eliminate below pivot
            # --- COMPLETE: compute elimination factor and update row i ---
            factor =
            Aug[i, k:nb] =                                  # row i -= factor * row k 
            # --- COMPLETE ---

    # ── 2. Back Substitution ──────────────────────────────────────────────────
    x = np.zeros(n)
    # --- COMPLETE: solve for x[n-1] first, then work backwards ---
    x[n-1] =
    for i in range(n - 2, -1, -1):
        x[i] =
    # --- COMPLETE ---

    return x, Aug
```

이어지는 두 셀은 위의 $`3 \times 3`$ 예제와 문제 1의 질량-스프링 시스템으로 함수를 시험한다.

```python
# Test gauss_naive on the 3x3 example traced above
A_test = np.array([[ 2,  1, -1],
                   [-3, -1,  2],
                   [-2,  1,  2]], dtype=float)
b_test = np.array([8.0, -11.0, -3.0])

x_sol, Aug_result = gauss_naive(A_test, b_test)

print('Upper triangular augmented matrix after elimination:')
print(np.round(Aug_result, 6))
print()
print(f'Solution x  = {x_sol}')
print(f'Expected    = [2, 3, -1]')
print(f'Residual    = {A_test @ x_sol - b_test}')
```

```python
# Apply gauss_naive to the spring-mass system
x_gauss, _ = gauss_naive(A, b)
x_numpy     = np.linalg.solve(A, b)

print('Spring-mass displacements:')
print(f'  Gauss Elimination : {x_gauss}')
print(f'  np.linalg.solve   : {x_numpy}')
print(f'  Max difference    : {np.max(np.abs(x_gauss - x_numpy)):.2e}')
```

## 문제 3. 부분 피봇팅 (노트북 3.2절)

단순 가우스 소거법은 각 단계에서 피봇 원소 $`a_{kk}^{(k)}`$로 나눈다. 피봇이 0이면 알고리즘이 실패하고, 매우 작으면 재앙적 소거와 같은 방식으로 반올림 오차가 증폭된다. 예를 들어 다음 시스템은 첫 피봇이 $`a_{11} = 0`$이어서 단순 소거법이 첫 단계에서 0으로 나누게 된다.

```math
\begin{bmatrix} 0 & -3 & 7 \\ 1 & 2 & -1 \\ 5 & -2 & 0 \end{bmatrix}\mathbf{x} =
\begin{bmatrix} 4 \\ 0 \\ 3 \end{bmatrix}
```

```python
A_bad = np.array([[0.0, -3.0,  7.0],
                  [1.0,  2.0, -1.0],
                  [5.0, -2.0,  0.0]])
b_bad = np.array([4.0, 0.0, 3.0])

print('True solution (via numpy):')
print('  x =', np.linalg.solve(A_bad, b_bad))

print('\nNaive Gaussian Elimination result:')
x_fail, Aug_fail = gauss_naive(A_bad, b_bad)
print('  x =', x_fail)
print('Augmented matrix after elimination:')
print(Aug_fail)
```

**부분 피봇팅**은 각 소거 단계 전에 행의 순서를 바꾼다. 피봇 단계 $`k`$에서 다음을 수행한다.

1. 열 $`k`$에서 절댓값이 가장 큰 행 $`r \geq k`$를 찾는다.
2. 확대 행렬에서 행 $`k`$와 행 $`r`$을 맞바꾼다.
3. 전진 소거를 진행한다.

이렇게 하면 시스템의 해가 유일한 한 피봇이 0이 되지 않고, 소거 계수가 $`|f_{ik}| = |a_{ik}| / |a_{kk}| \leq 1`$이 되어 오차 증폭을 막는다. 이 알고리즘을 부분 피봇팅을 적용한 가우스 소거법(GEPP)이라 한다. 함수 `gauss_pivot`의 빈칸을 채워 행 `k`와 행 `imax`를 맞바꾸게 하라.

```python
def gauss_pivot(A, b):
    """
    Solve Ax = b using Gaussian Elimination with Partial Pivoting.

    Parameters
    ----------
    A : (n, n) array-like  — coefficient matrix
    b : (n,) array-like    — right-hand side vector

    Returns
    -------
    x   : (n,) ndarray — solution vector
    Aug : (n, n+1) ndarray — augmented matrix after elimination
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).reshape(-1, 1)
    n  = A.shape[0]
    nb = n + 1
    Aug = np.hstack((A, b))

    for k in range(n - 1):
        # ── Partial pivoting ────────────────────────────────────────────────
        # Find the row index (>= k) with the largest |value| in column k
        imax = np.argmax(np.abs(Aug[k:, k])) + k
        if imax != k:
            # --- COMPLETE: swap row k and row imax in the augmented matrix ---
            Aug[[k, imax]] =                        # swap rows
            # --- COMPLETE ---

        # ── Forward elimination (same as naive) ─────────────────────────────
        for i in range(k + 1, n):
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]

    # ── Back substitution (same as naive) ───────────────────────────────────
    x = np.zeros(n)
    x[n-1] = Aug[n-1, nb-1] / Aug[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, nb-1] - Aug[i, i+1:n] @ x[i+1:n]) / Aug[i, i]

    return x, Aug
```

이어지는 셀은 앞에서 실패한 예제로 함수를 시험한다.

```python
# Test gauss_pivot on the previously-failing zero-pivot example
x_pivot, Aug_pivot = gauss_pivot(A_bad, b_bad)
x_ref              = np.linalg.solve(A_bad, b_bad)

print('Gauss with partial pivoting:')
print(f'  x       = {x_pivot}')
print(f'  x (ref) = {x_ref}')
print(f'  Max error: {np.max(np.abs(x_pivot - x_ref)):.2e}')
print()
print('Upper triangular augmented matrix after pivoting:')
print(np.round(Aug_pivot, 6))
```

## 문제 4. 토마스 알고리즘 (노트북 4.2절)

1차원 열전도, 보의 처짐, 유한 차분 편미분방정식 같은 많은 공학 문제는 주대각선과 그 바로 아래, 바로 위의 대각선만 0이 아닌 삼중대각 시스템을 만든다.

```math
\begin{bmatrix}
f_0 & g_0 &        &        &     \\
e_1 & f_1 & g_1    &        &     \\
    & e_2 & f_2    & g_2    &     \\
    &     & \ddots & \ddots & \ddots \\
    &     &        & e_{n-1} & f_{n-1}
\end{bmatrix}\mathbf{x} = \mathbf{r}
```

이 구조에 가우스 소거법을 적용하면 $`O(n^3)`$ 대신 $`O(n)`$ 번의 연산으로 끝나는 알고리즘을 얻는다. 토마스 알고리즘의 전진 소거는 $`k = 1, 2, \ldots, n-1`$에 대해 아래 대각선($`e`$ 성분)을 소거한다.

```math
\text{factor} = \frac{e_k}{f_{k-1}}, \qquad
f_k \leftarrow f_k - \text{factor} \cdot g_{k-1}, \qquad
r_k \leftarrow r_k - \text{factor} \cdot r_{k-1}
```

후진 대입은 $`k = n-2, n-3, \ldots, 0`$에 대해 다음과 같다.

```math
x_{n-1} = \frac{r_{n-1}}{f_{n-1}}, \qquad x_k = \frac{r_k - g_k\, x_{k+1}}{f_k}
```

함수 `thomas`의 빈칸을 채워 전진 소거를 완성하라.

```python
def thomas(e, f, g, r):
    """
    Solve a tridiagonal system using the Thomas Algorithm.

    Parameters (all 1-D arrays of length n)
    ----------
    e : subdiagonal    e[0] unused, e[1]..e[n-1] active
    f : main diagonal  f[0]..f[n-1]  (modified internally)
    g : superdiagonal  g[0]..g[n-2] active, g[n-1] unused
    r : right-hand side r[0]..r[n-1]  (modified internally)

    Returns
    -------
    x : (n,) ndarray — solution vector
    """
    f = np.array(f, dtype=float)   # work on copies
    r = np.array(r, dtype=float)
    n = len(f)

    # ── Forward sweep: eliminate subdiagonal ─────────────────────────────────
    for k in range(1, n):
        # --- COMPLETE: compute factor and update f[k] and r[k] ---
        factor =
        f[k] =
        r[k] =
        # --- COMPLETE ---

    # ── Back substitution ────────────────────────────────────────────────────
    x = np.zeros(n)
    x[-1] = r[-1] / f[-1]
    for k in range(n - 2, -1, -1):
        x[k] = (r[k] - g[k] * x[k+1]) / f[k]

    return x
```

노트북 4.3절은 이 함수를 1차원 정상상태 열전도에 적용한다. 길이 $`L = 1\,\text{m}`$인 막대의 양 끝 온도가 $`T(0) = T_L`$, $`T(L) = T_R`$이다. 막대를 같은 간격 $`\Delta x = L/(n+1)`$의 $`n+1`$개 구간으로 나누고 $`d^2T/dx^2 \approx (T_{i-1} - 2T_i + T_{i+1})/\Delta x^2 = 0`$을 적용하면 $`T_{i-1} - 2T_i + T_{i+1} = 0`$ $`(i = 1, \ldots, n)`$을 얻는다. 이는 $`e_i = g_i = 1`$, $`f_i = -2`$인 삼중대각 시스템이며, 경계 조건은 $`r_1 = -T_L`$, $`r_n = -T_R`$로 들어간다(나머지 $`r_i = 0`$).

```python
def heat_conduction_1d(n_interior, T_left=50.0, T_right=150.0):
    """Solve steady-state 1D heat conduction on n_interior interior nodes."""
    n = n_interior
    f =  -2.0 * np.ones(n)   # main diagonal
    e =   1.0 * np.ones(n)   # subdiagonal  (e[0] unused)
    g =   1.0 * np.ones(n)   # superdiagonal (g[n-1] unused)
    r =   np.zeros(n)

    # Boundary conditions modify the first and last RHS entries
    r[0]  -= T_left
    r[-1] -= T_right

    T_interior = thomas(e, f, g, r)

    x_all = np.linspace(0, 1, n + 2)
    T_all = np.concatenate([[T_left], T_interior, [T_right]])
    return x_all, T_all


# Exact solution: linear profile (for zero internal heat generation)
T_left, T_right = 50.0, 150.0
x_exact = np.linspace(0, 1, 500)
T_exact = T_left + (T_right - T_left) * x_exact

plt.figure(figsize=(8, 4))
for n in [3, 6, 15]:
    x, T = heat_conduction_1d(n)
    plt.plot(x, T, 'o-', ms=6, label=f'Thomas, {n} interior nodes')
plt.plot(x_exact, T_exact, 'k--', lw=2, label='Exact (linear)')
plt.xlabel('Position $x$ (m)')
plt.ylabel('Temperature $T$ (°C)')
plt.title('Steady-State 1D Heat Conduction — Thomas Algorithm')
plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab03-LinearSystem_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
