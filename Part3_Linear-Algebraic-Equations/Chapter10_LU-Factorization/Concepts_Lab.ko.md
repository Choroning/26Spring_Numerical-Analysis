# 제10장 실습 — LU 분해와 촐레스키(LU Factorization and Cholesky)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 10

> **선수 지식**: [프로그래밍언어] MATLAB/Python. [선형대수학] 가우스 소거법 (제8-9장).
>
> **학습 목표**:
> 1. 행렬을 L과 U 인수로 분해할 수 있다
> 2. 효율적인 다중 풀이를 위해 LU 분해를 적용할 수 있다
> 3. 대칭 행렬에 대해 촐레스키 분해를 구현할 수 있다

---

<br>

## 목차

- [1. Doolittle 알고리즘](#1-doolittle-알고리즘)
  - [1.1 구현](#11-구현)
  - [1.2 SciPy와 비교](#12-scipy와-비교)
- [2. L과 U를 이용한 풀이](#2-l과-u를-이용한-풀이)
  - [2.1 전진 대입(Forward Substitution)](#21-전진-대입forward-substitution)
  - [2.2 후진 대입(Backward Substitution)](#22-후진-대입backward-substitution)
  - [2.3 피봇팅을 포함한 완전한 LU 풀이](#23-피봇팅을-포함한-완전한-lu-풀이)
- [3. 열저항 네트워크(Thermal Resistor Network)](#3-열저항-네트워크thermal-resistor-network)
- [4. 다중 우변 — LU가 빛나는 이유](#4-다중-우변--lu가-빛나는-이유)
  - [4.1 성능 비교](#41-성능-비교)
  - [4.2 다중 하중 케이스의 질량-스프링 시스템](#42-다중-하중-케이스의-질량-스프링-시스템)
- [5. 촐레스키 분해(Cholesky Factorization)](#5-촐레스키-분해cholesky-factorization)
  - [5.1 SPD 확인](#51-spd-확인)
  - [5.2 촐레스키 분해 및 풀이](#52-촐레스키-분해-및-풀이)
  - [5.3 FEM 스프링 네트워크 응용](#53-fem-스프링-네트워크-응용)
- [요약](#요약)

---

<br>

## 1. Doolittle 알고리즘

### 1.1 구현

Doolittle 알고리즘은 $U$의 행과 $L$의 열을 체계적으로 구축하여 LU 분해를 계산합니다. 각 피봇 인덱스 $k$에 대해, 알고리즘은 먼저 이전에 계산된 $L$과 $U$의 성분을 사용하여 $U$의 $k$번째 행을 계산한 다음, 피봇 $U_{kk}$로 나누어 $L$의 $k$번째 열을 계산합니다.

```python
import numpy as np
import scipy as sc
from scipy.linalg import lu

def lu_doolittle(A):
    """LU factorization using Doolittle's algorithm (no pivoting)."""
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros((n, n))
    A = A.astype(float).copy()

    for k in range(n):
        # Compute U row k
        for j in range(k, n):
            U[k, j] = A[k, j] - L[k, :k] @ U[:k, j]
        # Compute L column k
        for i in range(k + 1, n):
            L[i, k] = (A[i, k] - L[i, :k] @ U[:k, k]) / U[k, k]

    return L, U

# Test
A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
L, U = lu_doolittle(A)
print(f'||A - LU|| = {np.linalg.norm(A - L @ U):.2e}')  # ≈ 0
```

핵심 아이디어는 $U$의 $k$행에 대해:

$$U_{kj} = A_{kj} - \sum_{s=1}^{k-1} L_{ks} U_{sj}, \quad j = k, k+1, \ldots, n$$

그리고 $L$의 $k$열에 대해:

$$L_{ik} = \frac{A_{ik} - \sum_{s=1}^{k-1} L_{is} U_{sk}}{U_{kk}}, \quad i = k+1, k+2, \ldots, n$$

이 공식들은 $A$의 성분과 곱 $LU$의 대응하는 성분을 같다고 놓는 것에서 직접 도출됩니다.

### 1.2 SciPy와 비교

SciPy의 `lu` 함수는 부분 피봇팅을 적용한 LU 분해를 수행하며, 치환 행렬 $P$, 하삼각 행렬 $L$, 상삼각 행렬 $U$를 반환합니다:

```python
# Compare with scipy
P, L_sp, U_sp = sc.linalg.lu(A)
print(f'scipy: P = \n{P}\nL = \n{L_sp}\nU = \n{U_sp}')
```

> **참고:** SciPy는 $PA = LU$를 만족하는 $P$를 반환합니다. 행 교환이 필요하지 않은 경우(이 양조건 행렬처럼), $P$는 단위 행렬이며 결과는 우리의 Doolittle 구현과 일치합니다.

---

<br>

## 2. L과 U를 이용한 풀이

### 2.1 전진 대입(Forward Substitution)

전진 대입은 위에서 아래로 각 $d_i$를 순차적으로 계산하여 하삼각 시스템 $Ld = b$를 풉니다. $L$이 단위 하삼각(대각 성분이 1)이므로 대각으로 나눌 필요가 없습니다:

```python
def forward_sub(L, b):
    """Solve Ld = b by forward substitution."""
    n = len(b)
    d = np.zeros(n)
    for i in range(n):
        d[i] = b[i] - L[i, :i] @ d[:i]
    return d
```

표현식 `L[i, :i] @ d[:i]`는 NumPy의 슬라이스 표기법을 사용하여 $\sum_{j=0}^{i-1} L_{ij} d_j$를 계산합니다. $i = 0$일 때 슬라이스가 비어 있어 내적은 0이므로, 기대대로 $d_0 = b_0$이 됩니다.

### 2.2 후진 대입(Backward Substitution)

후진 대입은 아래에서 위로 각 $x_i$를 계산하여 상삼각 시스템 $Ux = d$를 풉니다. 각 단계에서 대각 성분 $U_{ii}$로 나누어야 합니다:

```python
def backward_sub(U, d):
    """Solve Ux = d by backward substitution."""
    n = len(d)
    x = np.zeros(n)
    x[n-1] = d[n-1] / U[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - U[i, i+1:] @ x[i+1:]) / U[i, i]
    return x
```

### 2.3 피봇팅을 포함한 완전한 LU 풀이

분해와 두 대입 단계를 피봇팅을 처리하는 단일 함수로 결합합니다:

```python
def lu_solve(A, b):
    """Solve Ax = b using LU factorization with pivoting."""
    P, L, U = sc.linalg.lu(A)
    b_perm = P.T @ b        # Apply permutation
    d = forward_sub(L, b_perm)
    x = backward_sub(U, d)
    return x

# Test
b = np.array([7.85, -19.3, 71.4])
x = lu_solve(A, b)
print(f'Solution: {x}')  # [3.0, -2.5, 7.0003]
print(f'Residual: {np.linalg.norm(A @ x - b):.2e}')
```

치환은 `P.T @ b`로 적용됩니다. 이는 SciPy가 $PA = LU$를 만족하도록 $P$를 정의하므로, $Ax = b$가 $LUx = PAx = Pb$가 되기 때문입니다. 그러나 $P$는 좌측 곱 행렬로 반환되므로, 행 교환에 따라 $b$를 재배열하려면 $P^T b$가 필요합니다. 잔차 $\|Ax - b\|$는 기계 엡실론 근처여야 하며, 해의 정확도를 확인합니다.

---

<br>

## 3. 열저항 네트워크(Thermal Resistor Network)

연립방정식의 실용적 응용은 열저항 네트워크의 해석입니다. 경계 노드의 온도가 고정되고 내부 노드에 열원이 가해지는 4-노드 네트워크를 고려합니다:

- $T_1 = 100°C$ (고정), $T_4 = 20°C$ (고정)
- 열원: $Q_2 = 50\text{W}$, $Q_3 = 30\text{W}$
- 열저항: $R_{12} = 10$, $R_{23} = 20$, $R_{34} = 15$ (W/K)

각 내부 노드에서의 에너지 평형:

$$\frac{T_2 - T_1}{R_{12}} + \frac{T_2 - T_3}{R_{23}} = Q_2$$

$$\frac{T_3 - T_2}{R_{23}} + \frac{T_3 - T_4}{R_{34}} = Q_3$$

행렬 형태 $A \mathbf{T} = \mathbf{b}$로 재배열:

```python
# 4-node network: T1=100 (fixed), T4=20 (fixed)
# Heat sources: Q2=50W, Q3=30W
# Resistances: R12=10, R23=20, R34=15 (W/K)

# Energy balance: (T2-T1)/R12 + (T2-T3)/R23 = Q2
#                  (T3-T2)/R23 + (T3-T4)/R34 = Q3

A_net = np.array([
    [1/10 + 1/20, -1/20],
    [-1/20, 1/20 + 1/15]
])
b_net = np.array([50 + 100/10, 30 + 20/15])
T = np.linalg.solve(A_net, b_net)
print(f'T2 = {T[0]:.1f}°C, T3 = {T[1]:.1f}°C')
```

$A_\text{net}$의 대각 성분은 각 노드에 연결된 전도도(저항의 역수)의 합이고, 비대각 성분은 노드 간의 음의 전도도입니다. 우변에는 열원 항과 고정 경계 온도의 기여가 모두 포함됩니다.

---

<br>

## 4. 다중 우변 — LU가 빛나는 이유

### 4.1 성능 비교

LU 분해의 주요 이점은 동일한 계수 행렬을 여러 다른 우변 벡터와 사용해야 할 때 분명해집니다. $A$를 한 번 분해하고 인자를 재사용하는 것이 매번 처음부터 푸는 것보다 훨씬 저렴합니다:

```python
import time

n = 500
A = np.random.rand(n, n) + n * np.eye(n)
m = 100  # number of different RHS

# Method 1: Repeated np.linalg.solve
start = time.time()
for _ in range(m):
    b = np.random.rand(n)
    x = np.linalg.solve(A, b)
t1 = time.time() - start

# Method 2: Factor once, solve many
start = time.time()
lu_obj = sc.linalg.lu_factor(A)
for _ in range(m):
    b = np.random.rand(n)
    x = sc.linalg.lu_solve(lu_obj, b)
t2 = time.time() - start

print(f'Repeated solve: {t1:.3f}s')
print(f'LU factor + solve: {t2:.3f}s')
print(f'Speedup: {t1/t2:.1f}x')
```

$n = 500$이고 $m = 100$개의 우변에 대해:

- **반복 풀이** 는 `np.linalg.solve`를 100번 호출하며, 각각 내부적으로 전체 $O(n^3)$ 분해를 수행합니다. 총 비용: $100 \times O(n^3) = O(100 \, n^3)$.
- **LU 분해 + 풀이** 는 한 번의 $O(n^3)$ 분해와 100번의 전진/후진 대입을 각 $O(n^2)$로 수행합니다. 총 비용: $O(n^3) + O(100 \, n^2) \approx O(n^3)$ (큰 $n$에서).

속도 향상은 $m$에 따라 선형으로 증가합니다: 우변이 많을수록 사전 분해의 이점이 커집니다.

### 4.2 다중 하중 케이스의 질량-스프링 시스템

질량-스프링 시스템은 여러 하중 케이스가 동일한 강성 행렬을 공유하는 구체적인 공학 예시를 제공합니다:

```python
K = np.array([
    [300, -200, 0],
    [-200, 350, -150],
    [0, -150, 400]
], dtype=float)

F_gravity = np.array([0, 0, -9.81 * 10])  # 10 kg mass at node 3
F_wind = np.array([50, 30, 0])
F_combined = F_gravity + F_wind

lu_K = sc.linalg.lu_factor(K)
x_grav = sc.linalg.lu_solve(lu_K, F_gravity)
x_wind = sc.linalg.lu_solve(lu_K, F_wind)
x_comb = sc.linalg.lu_solve(lu_K, F_combined)

# Superposition check
print(f'||x_comb - (x_grav + x_wind)|| = {np.linalg.norm(x_comb - (x_grav + x_wind)):.2e}')
# Should be ~1e-16 (machine epsilon)
```

> **[물리]** 중첩 원리(principle of superposition)는 선형 시스템에서 결합된 하중에 대한 응답이 개별 응답의 합과 같음을 명시합니다. 이는 선형성의 직접적 결과입니다: $A(x_1 + x_2) = Ax_1 + Ax_2 = b_1 + b_2$.

중첩 검사는 이 성질을 수치적으로 확인합니다: 결합 해와 개별 해의 합 사이의 차이가 기계 엡실론 수준이며, 이는 연립방정식이 (부동소수점 정밀도까지) 중첩을 정확히 보존함을 보여줍니다.

---

<br>

## 5. 촐레스키 분해(Cholesky Factorization)

### 5.1 SPD 확인

촐레스키 분해를 적용하기 전에 행렬이 대칭 양정치(SPD)인지 확인해야 합니다. 행렬이 대칭($A = A^T$)이고 모든 고유값이 엄격히 양수이면 SPD입니다:

```python
def is_spd(A):
    """Check if A is symmetric positive definite."""
    if not np.allclose(A, A.T):
        return False
    eigenvalues = np.linalg.eigvalsh(A)
    return np.all(eigenvalues > 0)

# Example SPD matrix
A_spd = np.array([[6, 15, 55], [15, 55, 225], [55, 225, 979]], dtype=float)
print(f'Is SPD: {is_spd(A_spd)}')
```

`np.linalg.eigvalsh`는 대칭 행렬에 최적화되어 있으므로(`h`는 에르미트(Hermitian)를 의미) `np.linalg.eigvals` 대신 사용됩니다. 실수 고유값만 반환하며 이 목적에서 더 빠르고 수치적으로 더 안정적입니다.

### 5.2 촐레스키 분해 및 풀이

SciPy는 `sc.linalg.cholesky`를 통해 촐레스키 분해를 제공합니다. 기본적으로 상삼각 인자를 반환하며; `lower=True`를 전달하면 하삼각 인자를 반환합니다:

```python
# Cholesky
L_chol = sc.linalg.cholesky(A_spd, lower=True)
print(f'L = \n{L_chol}')
print(f'||A - L L^T|| = {np.linalg.norm(A_spd - L_chol @ L_chol.T):.2e}')

def cholesky_solve(A, b):
    """Solve Ax = b using Cholesky factorization."""
    L = sc.linalg.cholesky(A, lower=True)
    d = forward_sub(L, b)
    x = backward_sub(L.T, d)
    return x
```

풀이 절차는 일반 LU를 반영합니다: 먼저 전진 대입으로 $Ld = b$를 풀고, 후진 대입으로 $L^T x = d$를 풉니다. 핵심 차이는 하나의 삼각 인자만 계산하면 된다는 것($U = L^T$이므로)으로, 계산과 저장 모두 절반으로 줄입니다.

> **참고:** 2절에서 정의한 `forward_sub` 함수는 단위 하삼각 행렬(대각 성분이 1)을 가정합니다. 촐레스키에서는 $L$의 대각 성분이 1이 아니므로, 전진 대입 공식은 $d_i = (b_i - \sum_{j=1}^{i-1} L_{ij} d_j) / L_{ii}$가 됩니다. 실무에서는 SciPy의 `cho_solve`가 이를 자동으로 처리합니다.

### 5.3 FEM 스프링 네트워크 응용

유한요소법(FEM)으로 조립된 스프링 네트워크의 강성 행렬은 항상 SPD(유효하고 완전히 구속된 스프링 시스템의 경우)이므로, 촐레스키가 이상적인 분해 방법입니다:

```python
# Spring constants
k = [100, 200, 150, 250]

# Assembly (stiffness matrix is always SPD for valid spring systems)
K_fem = np.array([
    [k[0]+k[1], -k[1], 0],
    [-k[1], k[1]+k[2], -k[2]],
    [0, -k[2], k[2]+k[3]]
], dtype=float)

F_fem = np.array([50.0, 0.0, -30.0])
print(f'K is SPD: {is_spd(K_fem)}')
x_fem = cholesky_solve(K_fem, F_fem)
print(f'Displacements: {x_fem}')
```

강성 행렬 $K$는 각 스프링 요소의 기여를 합산하여 조립됩니다. 노드 $i$와 $i+1$을 연결하는 요소 $k_i$는 대각 성분 $(i,i)$와 $(i+1,i+1)$에 $+k_i$를, 비대각 성분 $(i,i+1)$과 $(i+1,i)$에 $-k_i$를 기여합니다. 결과 행렬은 항상 대칭(구성에 의해)이고 양정치(에너지 $\frac{1}{2} x^T K x > 0$이 영이 아닌 모든 변위 $x$에 대해 성립)입니다.

---

<br>

## 요약

| 방법 | 함수 | 사용 시기 |
|:-----|:-----|:----------|
| Doolittle LU | `lu_doolittle(A)` | 알고리즘 이해 |
| SciPy LU | `sc.linalg.lu(A)` | 일반 분해 |
| LU 분해 + 풀이 | `lu_factor(A)` + `lu_solve(lu, b)` | 다중 우변 (한 번 분해!) |
| 촐레스키(Cholesky) | `sc.linalg.cholesky(A)` | SPD 행렬 (2배 빠름) |
| `np.linalg.solve` | 직접 풀이 | 단일 우변, 실무 사용 |

**실용적 지침:**

1. 단일 우변에는 항상 `np.linalg.solve` 사용 — 수동으로 역행렬 계산하지 않기
2. 같은 $A$로 다중 우변: `lu_factor` + `lu_solve` 사용
3. $A$가 SPD이면: 촐레스키 사용 (더 빠르고 더 안정적)
4. 결과를 신뢰하기 전에 조건수 확인
5. 잔차로 해 검증: $\|Ax - b\|$

---
