# 제9장 실습 — 피봇팅, 조건수, 삼중대각 시스템(Pivoting, Condition Number, and Tridiagonal Systems)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 9

> **선수 지식**: [프로그래밍언어] MATLAB/Python. [선형대수학] 행렬 연산 (제8장).
>
> **학습 목표**:
> 1. 부분 피봇팅을 포함한 가우스 소거법을 수행할 수 있다
> 2. 소거법의 계산 복잡도를 분석할 수 있다
> 3. 후진 대입 알고리즘을 구현할 수 있다

---

<br>

## 목차

- [1. 부분 피봇팅 구현](#1-부분-피봇팅-구현)
  - [1.1 단순 소거법의 문제점](#11-단순-소거법의-문제점)
  - [1.2 부분 피봇팅을 적용한 가우스 소거법](#12-부분-피봇팅을-적용한-가우스-소거법)
- [2. 조건수(Condition Number)](#2-조건수condition-number)
  - [2.1 양조건 시스템 vs. 악조건 시스템](#21-양조건-시스템-vs-악조건-시스템)
  - [2.2 악조건 시스템의 섭동 민감도](#22-악조건-시스템의-섭동-민감도)
- [3. 토마스 알고리즘과 열전도](#3-토마스-알고리즘과-열전도)
  - [3.1 토마스 알고리즘 구현](#31-토마스-알고리즘-구현)
  - [3.2 1차원 정상상태 열전도 응용](#32-1차원-정상상태-열전도-응용)
  - [3.3 성능 비교: 토마스 vs. 일반 솔버](#33-성능-비교-토마스-vs-일반-솔버)
- [요약](#요약)

---

<br>

## 1. 부분 피봇팅 구현

### 1.1 단순 소거법의 문제점

단순 가우스 소거법(Naive Gauss elimination)은 피봇 요소가 0일 때 실패합니다. $(1,1)$ 성분이 0인 다음 시스템을 고려합니다:

```python
import numpy as np

# Problem: naive Gauss fails when pivot is zero
A_bad = np.array([[0, 2, 1], [1, 1, 2], [2, 1, 1]], dtype=float)
b_bad = np.array([1, 1, 1], dtype=float)

# gauss_naive(A_bad, b_bad) would fail — division by zero!
```

첫 번째 피봇 $a_{11} = 0$이 즉시 실패를 유발합니다. 피봇이 정확히 0이 아니더라도 매우 작은 경우, 결과적으로 큰 승수가 반올림 오차를 재앙적으로 증폭합니다.

### 1.2 부분 피봇팅을 적용한 가우스 소거법

해결책은 현재 열에서 (피봇 행 이하에서) 가장 큰 요소를 찾아 각 소거 단계 전에 행을 교환하는 것입니다:

```python
def gauss_pivot(A, b):
    """Solve Ax = b using Gaussian elimination with partial pivoting."""
    A = A.astype(float).copy()
    b = b.astype(float).copy()
    n = len(b)
    Aug = np.hstack([A, b.reshape(-1, 1)])
    nb = n + 1

    for k in range(n - 1):
        # Partial pivoting: find row with max |value| in column k
        imax = np.argmax(np.abs(Aug[k:, k])) + k
        if imax != k:
            Aug[[k, imax]] = Aug[[imax, k]]  # swap rows

        # Forward elimination (same as naive)
        for i in range(k + 1, n):
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]

    # Back substitution
    x = np.zeros(n)
    x[n-1] = Aug[n-1, nb-1] / Aug[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, nb-1] - Aug[i, i+1:n] @ x[i+1:n]) / Aug[i, i]

    return x

# Test: now works!
x = gauss_pivot(A_bad, b_bad)
print(f'Solution: {x}')
print(f'Residual: {np.linalg.norm(A_bad @ x - b_bad):.2e}')
```

> **참고:** 핵심 줄은 `imax = np.argmax(np.abs(Aug[k:, k])) + k`입니다 — 이는 현재 피봇 행 이하의 행만 검색합니다. `+ k` 오프셋은 지역 인덱스를 전역 행 인덱스로 변환합니다. 이 오프셋이 없으면 잘못된 행이 교환됩니다.

행 교환 `Aug[[k, imax]] = Aug[[imax, k]]`는 NumPy의 고급 인덱싱을 사용하여 한 번의 연산으로 두 행을 교환합니다. 전진 소거와 후진 대입은 단순 버전과 동일합니다 — 유일한 추가 사항은 피봇 탐색과 교환입니다.

---

<br>

## 2. 조건수(Condition Number)

### 2.1 양조건 시스템 vs. 악조건 시스템

행렬의 **조건수(condition number)** $\kappa(A)$는 $Ax = b$의 해가 입력 데이터의 섭동에 얼마나 민감한지를 정량화합니다. 큰 조건수는 $A$ 또는 $b$의 작은 변화가 $x$의 큰 변화를 유발할 수 있음을 의미합니다:

```python
# Well-conditioned system
A_good = np.array([[1, 0], [0, 1]], dtype=float)
print(f'Identity κ = {np.linalg.cond(A_good):.1f}')  # 1.0

# Ill-conditioned: Hilbert matrix
from scipy.linalg import hilbert
H = hilbert(5)
print(f'Hilbert(5) κ = {np.linalg.cond(H):.1f}')  # ~476607
```

단위 행렬(identity matrix)의 조건수는 $\kappa = 1$ (가능한 최선)이며, 섭동이 전혀 증폭되지 않음을 의미합니다. 힐베르트 행렬(Hilbert matrix)은 악조건 행렬의 고전적 예시입니다 — $5 \times 5$ 크기에서도 조건수가 거의 $500{,}000$입니다.

> **[선형대수]** 조건수 $\kappa(A)$는 최대 특이값과 최소 특이값의 비율입니다. 이는 최악의 오차 증폭을 제한합니다: 입력의 상대 섭동 $\varepsilon$이 출력에서 최대 $\kappa(A) \cdot \varepsilon$의 상대 섭동을 유발할 수 있습니다.

### 2.2 악조건 시스템의 섭동 민감도

다음 실험은 악조건의 실질적 영향을 보여줍니다. 우변 $b$의 미세한 섭동이 해에 불균형적으로 큰 변화를 유발합니다:

```python
# Perturbation sensitivity
n = 10
H10 = hilbert(n)
b_true = H10 @ np.ones(n)
x_exact = np.linalg.solve(H10, b_true)

# Perturb b slightly
b_perturbed = b_true + 1e-10 * np.random.randn(n)
x_perturbed = np.linalg.solve(H10, b_perturbed)

print(f'Hilbert(10) κ = {np.linalg.cond(H10):.2e}')
print(f'Input perturbation: {np.linalg.norm(b_perturbed - b_true):.2e}')
print(f'Output perturbation: {np.linalg.norm(x_perturbed - x_exact):.2e}')
# Output perturbation >> input perturbation!
```

출력 섭동은 입력 섭동보다 수 차수 크며, 증폭 계수는 조건수에 의해 제한됩니다. `hilbert(10)`의 경우 $\kappa \approx 10^{13}$이므로, $b$의 $10^{-10}$ 섭동이 $x$에서 최대 $10^{3}$의 변화를 유발할 수 있습니다. 이는 시스템을 배정밀도(double precision)에서 본질적으로 풀 수 없게 만듭니다.

> **참고:** $\kappa(A) \cdot \varepsilon_{\text{mach}} \gtrsim 1$일 때, 시스템은 너무 악조건이어서 표준 배정밀도에서 어떤 수치 방법도 의미 있는 해를 생성할 수 없습니다. 크기 $n \geq 13$인 힐베르트 행렬이 이 임계값에 도달합니다.

---

<br>

## 3. 토마스 알고리즘과 열전도

### 3.1 토마스 알고리즘 구현

**토마스 알고리즘(Thomas Algorithm)**은 밴드 구조를 활용하여 삼중대각 시스템을 $O(n)$ 시간에 풉니다. 세 개의 대각선만 저장하고 처리하면 됩니다:

```python
def thomas(e, f, g, r):
    """Solve tridiagonal system using Thomas algorithm.

    Args:
        e: sub-diagonal (n elements, e[0] unused)
        f: main diagonal (n elements)
        g: super-diagonal (n elements, g[n-1] unused)
        r: right-hand side (n elements)
    Returns:
        x: solution vector
    """
    e = e.astype(float).copy()
    f = f.astype(float).copy()
    g = g.astype(float).copy()
    r = r.astype(float).copy()
    n = len(f)

    # Forward sweep
    for k in range(1, n):
        factor = e[k] / f[k-1]
        f[k] = f[k] - factor * g[k-1]
        r[k] = r[k] - factor * r[k-1]

    # Back substitution
    x = np.zeros(n)
    x[n-1] = r[n-1] / f[n-1]
    for k in range(n - 2, -1, -1):
        x[k] = (r[k] - g[k] * x[k+1]) / f[k]

    return x
```

`.copy()` 호출은 필수입니다 — 이것이 없으면 NumPy 배열은 참조로 전달되므로 함수가 호출자의 배열을 수정하게 됩니다. 전진 소거는 하대각 성분을 하나씩 제거하고, 후진 대입은 아래에서 위로 해를 복원합니다.

> **[Python]** `.astype(float).copy()` 패턴은 두 가지를 보장합니다: (1) 정수 배열이 전달되어도 부동소수점 산술을 사용하고, (2) 원본 배열이 수정되지 않습니다. 이는 수치 Python에서 흔한 방어적 프로그래밍 패턴입니다.

### 3.2 1차원 정상상태 열전도 응용

삼중대각 시스템의 고전적 응용은 1차원 정상상태 열 방정식(heat equation)입니다:

$$\frac{d^2T}{dx^2} = 0 \quad \text{with} \quad T(0) = T_L, \quad T(L) = T_R$$

균일 격자 간격 $\Delta x$에서 중심 유한 차분(central finite differences)으로 이산화하면:

$$\frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2} = 0$$

이는 삼중대각 시스템 $T_{i-1} - 2T_i + T_{i+1} = 0$으로 단순화되며, 경계 조건은 우변을 통해 들어갑니다:

```python
import matplotlib.pyplot as plt

def heat_conduction_1d(T_left, T_right, n_interior):
    """Solve 1D steady heat conduction with Dirichlet BCs."""
    n = n_interior

    e = np.ones(n)        # sub-diagonal
    f = -2 * np.ones(n)   # main diagonal
    g = np.ones(n)        # super-diagonal
    r = np.zeros(n)       # RHS

    # Boundary conditions
    r[0] = -T_left
    r[-1] = -T_right

    T_interior = thomas(e, f, g, r)

    # Full temperature profile
    x = np.linspace(0, 1, n + 2)
    T = np.concatenate([[T_left], T_interior, [T_right]])
    return x, T

# Solve for different mesh sizes
for n in [3, 6, 15]:
    x, T = heat_conduction_1d(100, 20, n)
    plt.plot(x, T, 'o-', label=f'n={n}')

plt.plot([0, 1], [100, 20], 'k--', label='Exact (linear)')
plt.xlabel('Position x')
plt.ylabel('Temperature T')
plt.title('1D Steady Heat Conduction')
plt.legend()
plt.grid(True)
plt.show()
```

라플라스 방정식($d^2T/dx^2 = 0$)의 경우 정확한 해는 선형입니다: $T(x) = T_L + (T_R - T_L) \cdot x / L$. 수치 해는 중심 차분 근사가 선형 함수를 정확히 표현하므로 어떤 격자 크기에서도 (반올림까지) 정확히 일치합니다.

> **[미분방정식]** 중심 차분 근사 $T''(x_i) \approx (T_{i-1} - 2T_i + T_{i+1}) / (\Delta x)^2$는 2차 정확도를 가집니다. 라플라스 방정식에 대해서는 절단 오차가 $T^{(4)}(x)$를 포함하는데, 이는 선형(및 2차) 함수에 대해 0이므로 정확한 결과를 줍니다.

### 3.3 성능 비교: 토마스 vs. 일반 솔버

토마스 알고리즘과 일반 가우스 소거법 사이의 $O(n)$ 대 $O(n^3)$ 차이는 큰 문제 크기에서 극적입니다:

```python
import time

sizes = [100, 1000, 10000, 100000]
t_thomas = []
t_numpy = []

for n in sizes:
    e = np.ones(n)
    f = -2 * np.ones(n)
    g = np.ones(n)
    r = np.random.rand(n)

    # Thomas
    start = time.time()
    for _ in range(100):
        thomas(e, f, g, r)
    t_thomas.append((time.time() - start) / 100)

    # numpy (build full matrix)
    A = np.diag(f) + np.diag(g[:-1], 1) + np.diag(e[1:], -1)
    start = time.time()
    for _ in range(10):
        np.linalg.solve(A, r)
    t_numpy.append((time.time() - start) / 10)

print(f'{"n":>8} {"Thomas (s)":>12} {"NumPy (s)":>12} {"Speedup":>10}')
for i, n in enumerate(sizes):
    speedup = t_numpy[i] / t_thomas[i] if t_thomas[i] > 0 else float('inf')
    print(f'{n:8d} {t_thomas[i]:12.6f} {t_numpy[i]:12.6f} {speedup:10.1f}x')
```

$n = 100{,}000$에서 토마스 알고리즘은 밀리초 단위로 완료되지만, `np.linalg.solve` (전체 $100{,}000 \times 100{,}000$ 행렬을 구축하고 분해)는 테라바이트의 메모리와 수 시간의 계산이 필요합니다. 이는 대규모 과학 계산에서 행렬 구조를 활용하는 것이 왜 필수적인지를 보여줍니다.

> **참고:** 실무에서는 삼중대각 시스템에 `scipy.linalg.solve_banded()`를 사용하십시오 — 순수 Python 토마스 알고리즘보다 빠르면서 $O(n)$ 복잡도를 유지하는 최적화된 LAPACK 구현(DGTSV)을 제공합니다.

---

<br>

## 요약

| 방법 | 구현 | 사용 시기 |
|:-----|:-----|:----------|
| `gauss_pivot()` | 부분 피봇팅 가우스 | 일반 밀집 시스템 |
| `thomas()` | 토마스 알고리즘 $O(n)$ | 삼중대각 시스템 |
| `np.linalg.solve()` | 최적화된 LAPACK | 실무 사용 (항상) |
| `np.linalg.cond()` | 조건수 | 악조건 진단 |

핵심 구현 요점:

- **부분 피봇팅**은 최소한의 코드(피봇 탐색 + 행 교환)를 추가하지만 0으로 나눗셈을 방지하고 오차 성장을 제어
- **조건수**는 해를 신뢰하기 전에 항상 확인해야 함 — $\kappa(A) \cdot \varepsilon_{\text{mach}} \approx 1$이면 결과를 신뢰할 수 없음
- **토마스 알고리즘**은 삼중대각 구조를 활용하여 $O(n)$을 달성 — 행렬이 특수 구조를 가질 때는 항상 특화된 솔버 사용

---
