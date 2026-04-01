# 제11장 실습 — 역행렬과 조건수

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 11

> **선수 지식**: [프로그래밍언어] MATLAB/Python. [선형대수학] LU 분해 (제8-10장).
>
> **학습 목표**:
> 1. 행렬의 역행렬을 수치적으로 계산할 수 있다
> 2. 조건수를 계산하고 해석할 수 있다
> 3. 입력 변동에 대한 해의 민감도를 평가할 수 있다

---

<br>

## 목차

- [1. 역행렬](#1-역행렬)
  - [1.1 LU 분해를 이용한 역행렬 계산](#11-lu-분해를-이용한-역행렬-계산)
  - [1.2 황금 규칙: 역행렬 대신 직접 풀기!](#12-황금-규칙-역행렬-대신-직접-풀기)
- [2. 조건수와 비정상 조건](#2-조건수와-비정상-조건)
  - [2.1 조건수 계산](#21-조건수-계산)
  - [2.2 기하학적 해석](#22-기하학적-해석)
  - [2.3 힐베르트 행렬 — 지수적 비정상 조건](#23-힐베르트-행렬--지수적-비정상-조건)
  - [2.4 오차 증폭 시연](#24-오차-증폭-시연)
  - [2.5 진단 실습: 의심스러운 강성행렬](#25-진단-실습-의심스러운-강성행렬)
- [요약](#요약)

---

<br>

## 1. 역행렬

### 1.1 LU 분해를 이용한 역행렬 계산

행렬 $[A]$의 역행렬은 LU 분해를 사용하여 열별로 계산할 수 있습니다: $[A]$를 한 번 분해한 후, 단위행렬의 각 열 $\{e_k\}$에 대해 $[A]\{x_k\} = \{e_k\}$를 풉니다. 역행렬에 대한 유용한 항등식들:

- $(AB)^{-1} = B^{-1}A^{-1}$
- $(A^T)^{-1} = (A^{-1})^T$
- $2 \times 2$ 행렬의 경우: $\begin{bmatrix} a & b \\ c & d \end{bmatrix}^{-1} = \frac{1}{ad - bc} \begin{bmatrix} d & -b \\ -c & a \end{bmatrix}$

```python
import numpy as np
import scipy as sc

def lu_inverse(A):
    """Compute A^{-1} using LU factorization."""
    n = A.shape[0]
    A_inv = np.zeros((n, n))
    lu_obj = sc.linalg.lu_factor(A)
    I = np.eye(n)

    for j in range(n):
        A_inv[:, j] = sc.linalg.lu_solve(lu_obj, I[:, j])

    return A_inv

A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
A_inv = lu_inverse(A)
print(f'A^{{-1}} = \n{A_inv}')
print(f'||A A^{{-1}} - I|| = {np.linalg.norm(A @ A_inv - np.eye(3)):.2e}')
```

LU 기반 접근법은 $[A]$를 하삼각행렬과 상삼각행렬 인수로 한 번 분해한 후, 모든 $n$개의 우변에 대해 재사용합니다. 이는 $n$개의 독립적인 가우스 소거를 계산하는 것보다 효율적입니다.

---

<br>

### 1.2 황금 규칙: 역행렬 대신 직접 풀기!

$[A]\{x\} = \{b\}$를 풀기 위해 $[A]^{-1}\{b\}$를 명시적으로 계산하는 것은 거의 항상 **잘못된** 접근법입니다. LU 분해를 통한 직접 풀기가 더 빠르고 수치적으로 더 정확합니다:

```python
import time

n = 500
A = np.random.rand(n, n) + n * np.eye(n)
b = np.random.rand(n)

# Method 1: Invert then multiply (SLOW, LESS ACCURATE)
start = time.time()
for _ in range(100):
    x_inv = np.linalg.inv(A) @ b
t_inv = (time.time() - start) / 100

# Method 2: Direct solve (FAST, MORE ACCURATE)
start = time.time()
for _ in range(100):
    x_solve = np.linalg.solve(A, b)
t_solve = (time.time() - start) / 100

print(f'Inverse + multiply: {t_inv:.6f}s')
print(f'Direct solve:       {t_solve:.6f}s')
print(f'Speedup: {t_inv/t_solve:.1f}x')
print(f'Error (inv):   {np.linalg.norm(A @ x_inv - b):.2e}')
print(f'Error (solve): {np.linalg.norm(A @ x_solve - b):.2e}')
```

> **참고:** $Ax = b$를 직접 풀 수 있을 때 $A^{-1}b$를 계산하지 마십시오. 직접 풀기가 약 3배 빠르고 더 정확한 결과를 생성합니다. 전체 역행렬을 계산하는 유일한 이유는 행렬 원소 자체가 필요한 경우입니다 (예: 물리적 해석이나 공분산 분석).

---

<br>

## 2. 조건수와 비정상 조건

### 2.1 조건수 계산

조건수 $\kappa(A) = \|A\| \cdot \|A^{-1}\|$는 다양한 노름을 사용하여 계산할 수 있습니다. NumPy는 모든 변형에 대해 단일 함수를 제공합니다:

```python
# Computing condition number
A = np.array([[1, 2], [1.0001, 2]])
print(f'κ_2(A) = {np.linalg.cond(A, 2):.2e}')
print(f'κ_1(A) = {np.linalg.cond(A, 1):.2e}')
print(f'κ_∞(A) = {np.linalg.cond(A, np.inf):.2e}')

# 2-norm condition = ratio of singular values
U, S, Vt = np.linalg.svd(A)
print(f'σ_max/σ_min = {S[0]/S[-1]:.2e}')  # equals κ_2
```

2-노름 조건수는 특이값 분해(SVD)를 통해 깔끔한 해석을 갖습니다: $\kappa_2(A) = \sigma_{max} / \sigma_{min}$, 즉 최대 특이값과 최소 특이값의 비율입니다. 매우 작은 특이값을 가진 행렬은 거의 특이하므로 큰 조건수를 갖게 됩니다.

---

<br>

### 2.2 기하학적 해석

비정상 조건은 $2 \times 2$ 시스템에서 명확한 기하학적 의미를 갖습니다. 각 방정식은 평면에서 직선을 정의합니다. 정상 조건 시스템은 직선들이 급한 각도로 교차하므로 교점이 잘 정의됩니다. 비정상 조건 시스템은 거의 평행한 직선을 가지므로 작은 섭동만으로도 교점이 극적으로 이동합니다:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Well-conditioned
x_range = np.linspace(-2, 4, 100)
axes[0].plot(x_range, (5 - 2*x_range), label='2x + y = 5')
axes[0].plot(x_range, (x_range - 1), label='x - y = 1')
axes[0].set_title(f'Well-conditioned (κ ≈ {np.linalg.cond(np.array([[2,1],[1,-1]])):.1f})')
axes[0].legend()
axes[0].grid(True)

# Ill-conditioned
axes[1].plot(x_range, (5 - 2*x_range), label='2x + y = 5')
axes[1].plot(x_range, (5.001 - 2*x_range), label='2x + y = 5.001')
axes[1].set_title(f'Ill-conditioned (nearly parallel)')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()
```

정상 조건의 경우, 두 직선은 서로 다른 기울기를 가지며 명확한 각도로 교차합니다. 교점은 작은 섭동에도 안정적입니다. 비정상 조건의 경우, 두 직선은 거의 동일하므로(거의 평행) 계수의 작은 변화가 교점을 그래프의 한쪽 끝에서 다른 쪽 끝으로 이동시킵니다.

---

<br>

### 2.3 힐베르트 행렬 — 지수적 비정상 조건

힐베르트 행렬 $H_{ij} = 1/(i + j - 1)$은 비정상 조건의 표준 교과서 예제입니다. 조건수가 크기에 따라 **지수적으로** 증가하여 적당한 차원을 넘어서면 사실상 풀 수 없게 됩니다:

```python
from scipy.linalg import hilbert

print(f'{"n":>4} {"κ(H_n)":>15} {"digits lost":>12}')
print('-' * 35)
for n in range(2, 16):
    H = hilbert(n)
    kappa = np.linalg.cond(H)
    digits = int(np.log10(kappa)) if kappa > 1 else 0
    print(f'{n:4d} {kappa:15.2e} {digits:12d}')
```

배정밀도($\approx 16$개 유효 자릿수)에서 표를 보면, $n = 12$일 때 힐베르트 행렬의 조건수가 약 $10^{16}$이 되어 계산된 해가 본질적으로 신뢰할 수 있는 자릿수가 **전혀 없음** 을 알 수 있습니다. 힐베르트 행렬은 다항식 최소제곱 적합에서 자연스럽게 나타나며, 이것이 적당한 차수를 넘는 다항식 적합이 수치적으로 위험한 이유 중 하나입니다.

---

<br>

### 2.4 오차 증폭 시연

다음 실험은 조건수가 오차 증폭을 어떻게 제어하는지 보여줍니다. $\{x_{true}\} = \{1, 1, \ldots, 1\}$인 $[H]\{x\} = \{b\}$를 풀고, $\{b\}$를 아주 작은 양($10^{-10}$)만큼 섭동한 후 다시 풉니다:

```python
for n in [4, 8, 12]:
    H = hilbert(n)
    x_true = np.ones(n)
    b = H @ x_true

    # Solve
    x_computed = np.linalg.solve(H, b)

    # Perturb b
    b_pert = b + 1e-10 * np.random.randn(n)
    x_pert = np.linalg.solve(H, b_pert)

    rel_err = np.linalg.norm(x_computed - x_true) / np.linalg.norm(x_true)
    pert_err = np.linalg.norm(x_pert - x_computed) / np.linalg.norm(x_computed)

    print(f'n={n:2d}: κ={np.linalg.cond(H):.2e}, '
          f'solve error={rel_err:.2e}, perturbation effect={pert_err:.2e}')
```

결과는 명확한 패턴을 보여줍니다: $n$이 증가함에 따라 조건수가 지수적으로 커지고, 풀기 오차(정확한 산술을 의도했음에도)와 섭동 민감도 모두 비례하여 증가합니다. $n = 12$인 경우, $\{b\}$의 $10^{-10}$ 섭동이 $\{x\}$에서 $10^6$ 이상의 변화를 일으킬 수 있으며 — 이는 조건수에 맞먹는 $10^{16}$의 증폭 계수입니다.

---

<br>

### 2.5 진단 실습: 의심스러운 강성행렬

실제로 비정상 조건은 경계조건이 부적절하게 적용되거나 메시에 거의 퇴화된 요소가 있을 때 유한요소법(FEM) 강성행렬에서 자주 나타납니다. 이 예제는 의심스러운 행렬을 진단하는 방법을 보여줍니다:

```python
# A suspicious stiffness matrix from a FEM simulation
K_susp = np.array([
    [1e6, -1e6, 0, 0],
    [-1e6, 2e6, -1e6, 0],
    [0, -1e6, 2e6, -1e6],
    [0, 0, -1e6, 1e6+1e-2]
])

f_susp = np.array([0, 100, 0, 0])

kappa = np.linalg.cond(K_susp)
digits_lost = int(np.log10(kappa))
x_susp = np.linalg.solve(K_susp, f_susp)

print(f'Condition number: {kappa:.2e}')
print(f'Digits of accuracy lost: ~{digits_lost}')
print(f'Solution: {x_susp}')
print(f'Residual: {np.linalg.norm(K_susp @ x_susp - f_susp):.2e}')
```

> **참고:** 큰 조건수가 반드시 계산된 답이 틀렸음을 의미하지는 않습니다 — 주의를 기울이고 검증해야 한다는 것을 의미합니다. 잔차 $\|Ax - b\|$를 확인하십시오: 높은 $\kappa$에서도 작은 잔차는 주어진 데이터에 대해 상당히 좋은 해임을 나타냅니다.

---

<br>

## 요약

| 개념 | 함수 / 공식 |
|:--------|:------------------|
| 역행렬 계산 | `np.linalg.inv(A)` 또는 `lu_inverse(A)` |
| 풀기 (권장!) | `np.linalg.solve(A, b)` |
| 조건수 | `np.linalg.cond(A)`, `np.linalg.cond(A, 2)` |
| 노름 | `np.linalg.norm(A, ord)` (ord=1, 2, np.inf, 'fro') |
| 특이값 | `np.linalg.svd(A)` — $\kappa_2 = \sigma_{max} / \sigma_{min}$ |

핵심 규칙:

1. **역행렬 대신 직접 풀기** — `solve(A, b)`가 `inv(A) @ b`보다 빠르고 정확
2. 결과를 신뢰하기 전에 항상 $\kappa(A)$를 확인
3. 손실 자릿수 $\approx \log_{10}\kappa(A)$
4. 힐베르트 행렬: 조건수가 지수적으로 증가 — 비정상 조건의 전형적 예
5. 작은 잔차 $\|Ax - b\|$가 $\kappa$가 클 때 작은 오차 $\|x - x_{true}\|$를 보장하지 않음

---
