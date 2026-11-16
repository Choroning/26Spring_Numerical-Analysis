# HW3 과제 안내문 — LU 분해와 역행렬

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-16

수치해석 4주차 과제이다. 교수자가 배포한 실습 노트북 `Lab04-LU.ipynb`의 빈칸을 채워 실행한다. 노트북은 LU 분해, 다중 우변, 촐레스키 분해, 역행렬, 조건수와 악조건의 다섯 부분으로 이루어지며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이에 모두 6곳 있다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `scipy`를 `sc`로, `matplotlib.pyplot`을 `plt`로 불러오고 `time`을 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab04-LU.ipynb](Lab04-LU.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

---

<br>

## 문제 1. Doolittle 알고리즘 (노트북 1.2절)

가우스 소거법은 아래 행에서 위 행의 배수를 빼는 행 연산으로 $`A`$를 상삼각 행렬 $`U`$로 바꾼다. 이때 쓰는 소거 계수 $`\ell_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}`$를 모으면 대각 성분이 1인 하삼각 행렬 $`L`$이 되고, $`A = LU`$가 성립한다. $`L`$과 $`U`$를 저장해 두면 새로운 우변에 대해 $`O(n^3)`$의 분해를 반복하지 않고 $`O(n^2)`$의 삼각 시스템 풀이 두 번으로 $`A\mathbf{x} = \mathbf{b}`$를 풀 수 있다.

Doolittle 방법은 행 연산 없이 $`L`$과 $`U`$를 직접 계산한다. 모든 $`i`$에 대해 $`\ell_{ii} = 1`$로 정하고 $`A = LU`$를 성분별로 맞추면 $`k = 0, 1, \ldots, n-1`$에 대해 다음을 얻는다.

```math
u_{kj} = a_{kj} - \sum_{m=0}^{k-1} \ell_{km}\, u_{mj} \quad (j = k, \ldots, n-1), \qquad
\ell_{ik} = \frac{a_{ik} - \displaystyle\sum_{m=0}^{k-1} \ell_{im}\, u_{mk}}{u_{kk}} \quad (i = k+1, \ldots, n-1)
```

$`k = 0`$이면 두 합이 모두 비어 있으므로 $`u_{0j} = a_{0j}`$, $`\ell_{i0} = a_{i0}/u_{00}`$이다. 노트북은 다음 행렬로 과정을 단계별로 보인다.

```math
A = \begin{bmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{bmatrix}
```

- **$`k = 0`$:** $`U`$의 0행은 $`u_{00} = 2`$, $`u_{01} = 1`$, $`u_{02} = 1`$이고, $`L`$의 0열은 $`\ell_{10} = 4/2 = 2`$, $`\ell_{20} = 8/2 = 4`$이다.
- **$`k = 1`$:** $`U`$의 1행은 $`u_{11} = 3 - 2 \cdot 1 = 1`$, $`u_{12} = 3 - 2 \cdot 1 = 1`$이고, $`L`$의 1열은 $`\ell_{21} = (7 - 4 \cdot 1)/1 = 3`$이다.
- **$`k = 2`$:** $`U`$의 2행은 $`u_{22} = 9 - 4 \cdot 1 - 3 \cdot 1 = 2`$이다.

```math
L = \begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 4 & 3 & 1 \end{bmatrix},\quad
U = \begin{bmatrix} 2 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{bmatrix}
```

함수 `lu_doolittle`의 빈칸을 채워 $`U`$의 행과 $`L`$의 열을 계산하라.

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

이어지는 셀은 함수의 결과를 SciPy의 `sc.linalg.lu`와 비교한다.

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

## 문제 2. L과 U를 이용한 풀이 (노트북 1.4절)

$`A = LU`$이면 $`A\mathbf{x} = \mathbf{b}`$는 두 삼각 시스템으로 나뉜다. 피봇팅이 없으면 우변은 $`\mathbf{b}`$ 그대로이고, 피봇팅을 써서 SciPy의 관례대로 $`A = PLU`$로 분해했으면 우변은 $`\mathbf{b}' = P^T \mathbf{b}`$이다.

**1단계 (전진 대입):** $`L\mathbf{d} = \mathbf{b}`$를 푼다.

```math
d_0 = b_0, \qquad d_i = b_i - \sum_{j=0}^{i-1} \ell_{ij}\, d_j, \quad i = 1, \ldots, n-1
```

**2단계 (후진 대입):** $`U\mathbf{x} = \mathbf{d}`$를 푼다.

```math
x_{n-1} = \frac{d_{n-1}}{u_{n-1,n-1}}, \qquad
x_i = \frac{d_i - \displaystyle\sum_{j=i+1}^{n-1} u_{ij}\, x_j}{u_{ii}}, \quad i = n-2, \ldots, 0
```

함수 `forward_sub`와 `backward_sub`의 빈칸을 채워라. 함수 `lu_solve`는 SciPy의 LU 분해와 두 함수를 이어 붙인 것이다.

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

이어지는 셀은 해가 알려진 $`3 \times 3`$ 시스템으로 함수를 시험한다.

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

## 문제 3. 다중 하중 케이스의 질량-스프링 시스템 (노트북 2.2절)

같은 행렬 $`A`$에 우변만 다른 $`m`$개의 시스템 $`A\mathbf{x}_1 = \mathbf{b}_1, \ldots, A\mathbf{x}_m = \mathbf{b}_m`$을 풀어야 하는 경우를 생각한다. 시스템마다 가우스 소거를 하면 $`m \times O(n^3/3)`$이 들지만, LU 분해를 한 번 하고 삼각 시스템만 풀면 $`O(n^3/3) + m \times O(n^2)`$으로 끝난다.

노트북은 Lab 03(HW2)의 질량 세 개짜리 스프링 시스템을 다시 쓴다.

```math
K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix}
```

중력, 바람, 두 하중을 합한 결합 하중의 세 가지 하중 케이스에서 $`K`$를 매번 다시 분해하지 않고 변위를 구하려 한다. `scipy.linalg.lu_factor`로 한 번 분해하고 `scipy.linalg.lu_solve`로 세 하중 케이스를 모두 풀도록 빈칸을 채워라.

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

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Load Case         x1 (m)      x2 (m)      x3 (m)
------------------------------------------------
Gravity          0.52000     0.88000     0.92000
Wind             1.46000     1.74000     1.16000
Combined         1.98000     2.62000     2.08000

Superposition check ||x_comb - (x_grav + x_wind)||: 4.97e-16
```

## 문제 4. 촐레스키 분해를 이용한 풀이 (노트북 3.2절)

행렬 $`A`$가 $`A = A^T`$이고 모든 영이 아닌 벡터 $`\mathbf{v}`$에 대해 $`\mathbf{v}^T A \mathbf{v} > 0`$이면 **대칭 양정치 행렬(symmetric positive definite, SPD)** 이라 한다. 모든 고유값이 양수인 것, 모든 선행 주소행렬식이 양수인 것과 동치이며, 구조 역학의 강성 행렬, 최소제곱법의 정규 방정식 $`A^T A`$, 통계의 공분산 행렬에서 나타난다.

SPD 행렬은 대각 성분이 양수인 하삼각 행렬 $`L`$로 $`A = LL^T`$와 같이 분해된다(촐레스키 분해). 일반 LU와 비교하면 SPD 행렬에만 쓸 수 있지만 분해 비용이 약 $`n^3/6`$로 절반이고, 행렬 하나만 저장하며, 피봇팅 없이 안정하다. 공식은 열 단위로 $`k = 0, \ldots, n-1`$에 대해 다음과 같다.

```math
\ell_{kk} = \sqrt{a_{kk} - \sum_{m=0}^{k-1} \ell_{km}^2}, \qquad
\ell_{ik} = \frac{a_{ik} - \displaystyle\sum_{m=0}^{k-1} \ell_{im}\,\ell_{km}}{\ell_{kk}}, \quad i > k
```

함수 `cholesky_solve`의 빈칸을 채워라. 주석의 지시대로 문제 2에서 정의한 `forward_sub`를 쓰고, $`L^T`$가 상삼각 행렬임에 유의한다.

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

이어지는 셀은 강성 행렬이 SPD인 FEM 스프링 네트워크에 함수를 적용하고, SciPy의 SPD 전용 풀이와 일반 LU 풀이의 시간을 비교한다. 함수 `is_spd`는 앞 셀에서 대칭성과 고유값의 부호로 SPD 여부를 판정하도록 정의되어 있다.

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

## 문제 5. LU 분해를 이용한 역행렬 (노트북 4.2절)

정방 행렬 $`A`$의 역행렬은 $`AA^{-1} = A^{-1}A = I`$를 만족하는 유일한 행렬 $`A^{-1}`$이며, $`A`$가 정칙($`\det(A) \neq 0`$)일 때만 존재한다. $`X = A^{-1}`$이면 $`AX = I`$이므로, 단위 행렬의 $`j`$번째 열 $`\mathbf{e}_j`$에 대해 $`n`$개의 시스템을 풀어 역행렬을 구할 수 있다.

```math
A \mathbf{x}_j = \mathbf{e}_j, \quad j = 0, 1, \ldots, n-1
```

알고리즘은 다음과 같다.

1. $`A = LU`$로 한 번 분해한다($`O(n^3/3)`$).
2. 각 $`j = 0, \ldots, n-1`$에 대해 $`LU\mathbf{x}_j = \mathbf{e}_j`$를 푼다(열마다 $`O(n^2)`$).
3. $`A^{-1} = [\mathbf{x}_0 \mid \mathbf{x}_1 \mid \cdots \mid \mathbf{x}_{n-1}]`$로 모은다.

함수 `lu_inverse`의 빈칸을 채워 $`A^{-1}`$의 각 열을 구하라.

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

이어지는 셀은 문제 2의 행렬 `A_sys`로 함수를 시험하고, 특이 행렬의 역행렬을 구하려 할 때의 동작을 보인다.

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

## 문제 6. 악조건 시스템의 진단 (노트북 5.4절)

행렬 $`A`$의 조건수는 $`\kappa(A) = \lVert A \rVert \lVert A^{-1} \rVert`$이며, 우변의 작은 섭동 $`\delta\mathbf{b}`$가 해의 오차 $`\delta\mathbf{x}`$로 얼마나 증폭되는지를 나타낸다.

```math
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa(A)\,\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}
```

| $`\kappa(A)`$ | 해석 |
|:-----|:-----|
| $`\approx 1`$ | 조건이 좋다. 부동 소수점의 정확도를 모두 얻는다 |
| $`10^6`$ | 해의 유효 숫자를 약 6자리 잃을 수 있다 |
| $`\sim 1/\varepsilon_{\text{mach}} \approx 10^{16}`$ | 사실상 특이하다. 해를 믿을 수 없다 |

*표 1. 조건수의 크기에 따른 해석*

경험 규칙으로, $`\kappa(A) \approx 10^k`$이면 유효 숫자를 약 $`k`$자리 잃는다고 본다. NumPy에서 `np.linalg.cond(A)`는 기본값으로 2-노름 조건수, 곧 최대 특이값과 최소 특이값의 비 $`\kappa_2(A) = \sigma_{\max} / \sigma_{\min}`$을 준다.

연습 문제는 다음과 같다. 한 구조 엔지니어가 같은 골조 강성 문제에 대해 두 프로그램이 매우 다른 변위를 낸다고 보고한다. 아래에 주어진 행렬을 조건수로 진단하고 해결책을 제시하라. 다음 셀의 빈칸을 채워 조건수, 잃는 유효 숫자의 개수, 해를 구한다.

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

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

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

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab04-LU_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
