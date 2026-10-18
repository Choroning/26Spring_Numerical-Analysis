# 9장 실습 — 피봇팅, 조건수, 삼중대각 시스템(Pivoting, Condition Number, and Tridiagonal Systems)

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **최종 수정일:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 9

> **학습 목표:**
> 1. 부분 피봇팅을 포함한 가우스 소거법을 수행할 수 있다
> 2. 소거법의 계산 복잡도를 분석할 수 있다
> 3. 후진 대입 알고리즘을 구현할 수 있다
> 4. 조건수로 악조건 시스템을 진단할 수 있다
> 5. 토마스 알고리즘으로 삼중대각 시스템을 풀 수 있다

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
- [점검 문제](#점검-문제)

---

<br>

## 1. 부분 피봇팅 구현

### 1.1 단순 소거법의 문제점

단순 가우스 소거법(Naive Gauss elimination)은 피봇 요소가 0이면 실패한다. $`(1,1)`$ 성분이 0인 다음 시스템을 생각한다:

```python
import numpy as np

# Problem: naive Gauss fails when pivot is zero
A_bad = np.array([[0, 2, 1], [1, 1, 2], [2, 1, 1]], dtype=float)
b_bad = np.array([1, 1, 1], dtype=float)

# gauss_naive(A_bad, b_bad) would fail — division by zero!
```

첫 번째 피봇 $`a_{11} = 0`$이 즉시 실패를 일으킨다. 피봇이 정확히 0이 아니더라도 매우 작으면, 그 결과로 생기는 큰 승수가 반올림 오차를 재앙적으로 증폭한다.

### 1.2 부분 피봇팅을 적용한 가우스 소거법

해결책은 각 소거 단계 전에 현재 열에서 (피봇 행과 그 아래에서) 절댓값이 가장 큰 요소를 찾아 행을 교환하는 것이다:

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

> **참고:** 핵심 줄은 `imax = np.argmax(np.abs(Aug[k:, k])) + k`이다. 이 줄은 현재 피봇 행과 그 아래의 행만 검색한다. `+ k` 오프셋은 부분 배열의 지역 인덱스를 전체 행 인덱스로 바꾼다. 이 오프셋이 없으면 잘못된 행이 교환된다.

행 교환 `Aug[[k, imax]] = Aug[[imax, k]]`는 NumPy의 고급 인덱싱으로 두 행을 한 번에 교환한다. 전진 소거와 후진 대입은 단순 버전과 같으며, 추가된 것은 피봇 탐색과 교환뿐이다. 이 시스템의 해는 $`x = (0.2, 0.4, 0.2)`$이고 잔차는 0이다.

---

<br>

## 2. 조건수(Condition Number)

### 2.1 양조건 시스템 vs. 악조건 시스템

행렬의 **조건수(condition number)** $`\kappa(A)`$는 $`Ax = b`$의 해가 입력 데이터의 섭동에 얼마나 민감한지를 정량화한다. 조건수가 크면 $`A`$나 $`b`$의 작은 변화가 $`x`$의 큰 변화를 일으킬 수 있다:

```python
# Well-conditioned system
A_good = np.array([[1, 0], [0, 1]], dtype=float)
print(f'Identity κ = {np.linalg.cond(A_good):.1f}')  # 1.0

# Ill-conditioned: Hilbert matrix
from scipy.linalg import hilbert
H = hilbert(5)
print(f'Hilbert(5) κ = {np.linalg.cond(H):.1f}')  # ~476607
```

단위 행렬(identity matrix)의 조건수는 가능한 최솟값인 $`\kappa = 1`$이며, 섭동이 전혀 증폭되지 않는다는 뜻이다. 힐베르트 행렬(Hilbert matrix, $`H_{ij} = 1/(i + j - 1)`$)은 악조건 행렬의 고전적 예시로, $`5 \times 5`$ 크기에서도 조건수가 약 $`4.8 \times 10^5`$으로 거의 $`500{,}000`$에 이른다.

> **참고:** **조건수와 특이값.** 2-노름 기준의 조건수 $`\kappa(A)`$는 최대 특이값과 최소 특이값의 비율이다. 조건수는 최악의 오차 증폭을 제한한다. 입력의 상대 섭동 $`\varepsilon`$은 출력에서 최대 $`\kappa(A) \cdot \varepsilon`$의 상대 섭동을 일으킬 수 있다.

### 2.2 악조건 시스템의 섭동 민감도

다음 실험은 악조건의 실제 영향을 보여 준다. 우변 $`b`$의 미세한 섭동이 해에 불균형적으로 큰 변화를 일으킨다:

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

출력 섭동은 입력 섭동보다 몇 자릿수 크며, 증폭 정도의 상한은 조건수가 정한다. `hilbert(10)`은 $`\kappa \approx 1.6 \times 10^{13}`$이므로 $`b`$의 $`10^{-10}`$ 섭동이 $`x`$에서 최대 $`10^{3}`$ 정도의 변화를 일으킬 수 있다. 따라서 데이터에 $`10^{-10}`$ 수준의 오차만 있어도 이 시스템의 해는 사실상 의미를 잃는다.

> **참고:** $`\kappa(A) \cdot \varepsilon_{\text{mach}} \gtrsim 1`$이면 시스템이 너무 악조건이어서, 표준 배정밀도에서는 어떤 수치 방법도 의미 있는 해를 내지 못한다. 힐베르트 행렬은 $`n \geq 12`$에서 이 임계값에 도달한다($`\kappa(H_{12}) \approx 1.7 \times 10^{16}`$).

---

<br>

## 3. 토마스 알고리즘과 열전도

### 3.1 토마스 알고리즘 구현

**토마스 알고리즘(Thomas Algorithm)** 은 밴드 구조를 활용하여 삼중대각 시스템을 $`O(n)`$ 시간에 푼다. 세 대각선만 저장하고 처리하면 된다. 아래 구현은 세 대각선을 길이 $`n`$인 배열 `e`, `f`, `g`에 그림 1과 같이 저장한다. 행 $`i`$(0부터 셈)의 방정식은 $`e_i x_{i-1} + f_i x_i + g_i x_{i+1} = r_i`$이며, `e[0]`과 `g[n-1]`은 쓰지 않는다.

```text
index     0       1       2      ...     n-1
e         -      e[1]    e[2]    ...    e[n-1]
f        f[0]    f[1]    f[2]    ...    f[n-1]
g        g[0]    g[1]    g[2]    ...      -
```

*그림 1. 토마스 알고리즘에서 세 대각선을 배열에 저장하는 방식*

구현은 다음과 같다:

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

입력을 복사하는 것은 필수이다. 복사하지 않으면 NumPy 배열이 함수에 참조로 전달되므로 함수가 호출자의 배열을 수정하게 된다. 다만 `astype(float)`는 기본적으로 새 배열을 반환하므로 여기서 `.copy()`는 중복이며, 복사 의도를 드러내는 역할을 한다. 전진 소거는 하대각 성분을 하나씩 제거하고, 후진 대입은 아래에서 위로 해를 복원한다.

> **[컴퓨터프로그래밍]** `.astype(float).copy()` 패턴은 두 가지를 보장한다. (1) 정수 배열이 전달되어도 부동소수점 산술을 사용하고, (2) 원본 배열이 수정되지 않는다. 이는 수치 계산용 Python 코드에서 흔한 방어적 프로그래밍 패턴이다.

### 3.2 1차원 정상상태 열전도 응용

삼중대각 시스템의 고전적인 응용은 1차원 정상상태 열 방정식(heat equation)이다:

```math
\frac{d^2T}{dx^2} = 0 \quad \text{with} \quad T(0) = T_L, \quad T(L) = T_R
```

균일한 격자 간격 $`\Delta x`$에서 중심 유한 차분(central finite differences)으로 이산화하면 다음과 같다:

```math
\frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2} = 0
```

이 식은 삼중대각 시스템 $`T_{i-1} - 2T_i + T_{i+1} = 0`$으로 단순해지며, 경계 조건은 우변을 통해 들어간다. 예를 들어 첫 내부점의 식은 $`-2T_1 + T_2 = -T_L`$이 된다:

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

라플라스 방정식($`d^2T/dx^2 = 0`$)의 정확한 해는 선형 함수 $`T(x) = T_L + (T_R - T_L) \cdot x / L`$이다. 중심 차분 근사는 선형 함수를 정확히 표현하므로, 수치 해는 격자 크기와 관계없이 (반올림 오차를 제외하면) 정확히 일치한다.

> **참고:** **중심 차분의 정확도.** 중심 차분 근사 $`T''(x_i) \approx (T_{i-1} - 2T_i + T_{i+1}) / (\Delta x)^2`$는 2차 정확도를 가진다. 절단 오차가 $`T^{(4)}(x)`$를 포함하는데, 이 값은 3차 이하의 다항식(선형 함수 포함)에서 0이므로 라플라스 방정식에 대해서는 정확한 결과를 준다.

### 3.3 성능 비교: 토마스 vs. 일반 솔버

토마스 알고리즘과 일반 가우스 소거법의 $`O(n)`$ 대 $`O(n^3)`$ 차이는 문제 크기가 클수록 극적이다:

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

$`n = 100{,}000`$에서 순수 Python 토마스 알고리즘은 수십 밀리초 안에 끝나지만, `np.linalg.solve`는 $`100{,}000 \times 100{,}000`$ 밀집 행렬을 만들고 분해해야 하므로 행렬 하나에만 약 80 GB의 메모리와 수 시간의 계산이 필요하다. 따라서 위 코드를 그대로 실행하면 $`n = 100{,}000`$에서 밀집 행렬을 만들다가 메모리가 부족해진다. 이는 대규모 과학 계산에서 행렬 구조를 활용하는 것이 왜 필수적인지 보여 준다.

> **참고:** 실무에서는 삼중대각 시스템에 `scipy.linalg.solve_banded()`를 사용하는 것이 좋다. 이 함수는 $`O(n)`$ 복잡도를 유지하면서 순수 Python 토마스 알고리즘보다 빠른 최적화된 LAPACK 구현(DGTSV)을 사용한다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:-----|:-----|
| `gauss_pivot()` | 부분 피봇팅 가우스 소거법을 구현한 함수로, 일반 밀집 시스템에 쓴다. |
| `thomas()` | $`O(n)`$ 토마스 알고리즘을 구현한 함수로, 삼중대각 시스템에 쓴다. |
| `np.linalg.solve()` | 최적화된 LAPACK 루틴으로, 실무에서 일반 시스템을 풀 때 쓴다. |
| `np.linalg.cond()` | 조건수를 계산하여 악조건을 진단한다. |
| 부분 피봇팅 | 피봇 탐색과 행 교환이라는 최소한의 코드만 더하지만 0으로 나눗셈을 막고 오차 성장을 제어한다. |
| 조건수 | 해를 신뢰하기 전에 항상 확인해야 하며, $`\kappa(A) \cdot \varepsilon_{\text{mach}} \approx 1`$이면 결과를 신뢰할 수 없다. |
| 토마스 알고리즘 | 삼중대각 구조를 활용하여 $`O(n)`$을 달성하며, 행렬에 특수 구조가 있으면 항상 특화된 솔버를 쓴다. |

---

<br>

## 점검 문제

1. **피봇 행 인덱스:** `k = 1`이고 `Aug[1:, 1]`이 `[0.5, -3.0, 2.0]`일 때 `imax`는 얼마인가? `+ k`를 빼면 어떻게 되는가?

   <details>
   <summary>정답</summary>

   `np.abs(Aug[1:, 1])`은 `[0.5, 3.0, 2.0]`이므로 `np.argmax`는 지역 인덱스 1을 반환하고 `imax = 1 + 1 = 2`이다. `+ k`가 없으면 `imax = 1`이 되어 현재 피봇 행(값 0.5)을 자기 자신과 교환하므로 피봇팅이 일어나지 않는다.

   </details>

2. **조건수와 신뢰할 수 있는 자릿수:** $`\kappa(A) = 10^6`$인 시스템에서 $`b`$의 상대 오차가 $`10^{-8}`$이면 $`x`$의 상대 오차의 상한은 얼마인가? 배정밀도($`\varepsilon_{\text{mach}} \approx 2.2 \times 10^{-16}`$)에서 반올림 오차만 있다면 대략 몇 자리를 신뢰할 수 있는가?

   <details>
   <summary>정답</summary>

   상한은 $`\kappa(A) \times 10^{-8} = 10^{-2}`$, 즉 1%이다. 반올림 오차만 있으면 $`\kappa(A) \cdot \varepsilon_{\text{mach}} \approx 2.2 \times 10^{-10}`$이므로, 약 16자리 중 6자리를 잃고 약 9~10자리를 신뢰할 수 있다.

   </details>

3. **열전도 이산화:** $`T_L = 100`$, $`T_R = 20`$, 내부점 3개($`L = 1`$, $`\Delta x = 0.25`$)일 때 삼중대각 시스템을 쓰고 해를 구하라.

   <details>
   <summary>정답</summary>

   시스템은 $`-2T_1 + T_2 = -100`$, $`T_1 - 2T_2 + T_3 = 0`$, $`T_2 - 2T_3 = -20`$이다. 정확한 해 $`T(x) = 100 - 80x`$를 $`x = 0.25, 0.5, 0.75`$에서 구하면 $`T_1 = 80`$, $`T_2 = 60`$, $`T_3 = 40`$이다. 대입하면 $`-160 + 60 = -100`$, $`80 - 120 + 40 = 0`$, $`60 - 80 = -20`$으로 세 식을 모두 만족하며, `thomas()`도 같은 값을 낸다.

   </details>

4. **메모리 비교:** $`n = 10{,}000`$인 삼중대각 시스템을 배정밀도 밀집 행렬로 저장할 때와 세 대각선 배열로 저장할 때의 메모리를 비교하라.

   <details>
   <summary>정답</summary>

   밀집 행렬은 원소 $`10^8`$개에 원소당 8바이트이므로 $`8 \times 10^8`$바이트, 약 800 MB이다. 세 대각선은 원소 $`3 \times 10^4`$개이므로 240,000바이트, 약 240 KB로 밀집 행렬의 약 3,300분의 1이다.

   </details>
