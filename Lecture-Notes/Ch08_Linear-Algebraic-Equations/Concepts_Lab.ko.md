# 8장 실습 — 연립방정식과 가우스 소거법

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **최종 수정일:** 2026-10-27
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 8, 9

> **학습 목표:**
> 1. 공학 문제를 연립일차방정식으로 공식화할 수 있다
> 2. 해의 존재성과 유일성 조건을 식별할 수 있다
> 3. 연립방정식의 그래프적 해석을 적용할 수 있다
> 4. 단순 가우스 소거법을 구현하고 계산 비용을 실험으로 확인할 수 있다

---

<br>

## 목차

- [1. $`Ax = b`$ 설정하기](#1-ax--b-설정하기)
  - [1.1 행렬 형태](#11-행렬-형태)
  - [1.2 기하학적 해석 (2x2)](#12-기하학적-해석-2x2)
- [2. 질량-스프링 시스템 응용](#2-질량-스프링-시스템-응용)
  - [2.1 문제 설정](#21-문제-설정)
  - [2.2 풀이 및 검증](#22-풀이-및-검증)
- [3. 단순 가우스 소거법 구현](#3-단순-가우스-소거법-구현)
  - [3.1 알고리즘 구현](#31-알고리즘-구현)
  - [3.2 교과서 예제로 테스트](#32-교과서-예제로-테스트)
- [4. 계산 비용 분석](#4-계산-비용-분석)
  - [4.1 경험적 시간 측정](#41-경험적-시간-측정)
  - [4.2 $`O(n^3)`$ 스케일링 검증](#42-on3-스케일링-검증)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. $`Ax = b`$ 설정하기

### 1.1 행렬 형태

$`n`$개의 미지수에 대한 $`n`$개의 연립 일차방정식은 행렬-벡터 형태로 간결하게 나타낼 수 있다:

```math
[A]\{x\} = \{b\}
```

여기서 각 기호는 다음과 같다:
- $`[A]`$는 알려진 계수를 담은 $`n \times n`$ **계수 행렬(coefficient matrix)** 이다.
- $`\{x\}`$는 구하려는 변수를 담은 $`n \times 1`$ **미지수 벡터(unknown vector)** 이다.
- $`\{b\}`$는 알려진 상수를 담은 $`n \times 1`$ **우변 벡터(right-hand side, RHS vector)** 이다.

예를 들어 다음 시스템

```math
2x_1 + x_2 = 5
```
```math
x_1 - x_2 = 1
```

은 다음과 같이 나타낸다:

```math
\begin{bmatrix} 2 & 1 \\ 1 & -1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 5 \\ 1 \end{bmatrix}
```

### 1.2 기하학적 해석 (2x2)

$`2 \times 2`$ 시스템에서는 각 방정식이 평면에서 직선 하나를 정의한다. 따라서 유일한 해, 해 없음, 무한히 많은 해라는 세 가지 결과를 직접 시각화할 수 있다:

```python
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Case 1: Unique solution (2x + y = 5, x - y = 1)
x = np.linspace(-1, 4, 100)
axes[0].plot(x, 5 - 2*x, label='2x + y = 5')
axes[0].plot(x, x - 1, label='x - y = 1')
axes[0].plot(2, 1, 'ko', markersize=8)
axes[0].set_title('Unique Solution')
axes[0].legend()
axes[0].grid(True)

# Case 2: No solution (parallel: y = 2x + 1, y = 2x + 3)
axes[1].plot(x, 2*x + 1, label='y = 2x + 1')
axes[1].plot(x, 2*x + 3, label='y = 2x + 3')
axes[1].set_title('No Solution (Parallel)')
axes[1].legend()
axes[1].grid(True)

# Case 3: Infinite solutions (y = 2x + 1, 2y = 4x + 2)
axes[2].plot(x, 2*x + 1, label='y = 2x + 1', linewidth=3)
axes[2].plot(x, 2*x + 1, '--', label='2y = 4x + 2', linewidth=1)
axes[2].set_title('Infinite Solutions')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.show()
```

유일한 해의 경우 교차점 $`(2, 1)`$이 검은 점으로 표시된다. 평행선의 경우 기울기는 같고 절편이 다른 두 직선이 그려지며, 이 직선들은 만나지 않는다. 무한히 많은 해의 경우 두 방정식이 정확히 같은 직선을 나타낸다(한 방정식이 다른 방정식의 두 배).

---

<br>

## 2. 질량-스프링 시스템 응용

### 2.1 문제 설정

두 고정벽 사이에 스프링 4개로 연결된 질량 3개로 이루어진 시스템을 생각한다(그림 1). 스프링 상수는 $`k_1 = 100`$, $`k_2 = 200`$, $`k_3 = 150`$, $`k_4 = 250`$ N/m이고, 가해진 힘은 $`F_1 = 50`$, $`F_2 = 0`$, $`F_3 = -30`$ N이다. 변위 $`x_i`$와 힘 $`F_i`$는 같은 방향을 양으로 잡는다.

```text
||--k1--[m1]--k2--[m2]--k3--[m3]--k4--||
```

*그림 1. 두 고정벽 사이에 스프링 4개로 연결된 질량 3개*

각 질량에 훅의 법칙($`F = kx`$)과 힘 평형을 적용하면 강성 행렬 방정식 $`[K]\{x\} = \{F\}`$를 얻는다:

```math
\begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} F_1 \\ F_2 \\ F_3 \end{bmatrix}
```

각 대각 성분 $`(k_i + k_{i+1})`$은 질량 $`i`$에 작용하는 총 강성을 나타낸다. 비대각 성분 $`(-k_j)`$는 공유 스프링을 통한 인접 질량 사이의 결합을 나타낸다.

### 2.2 풀이 및 검증

```python
import numpy as np

# 3 masses, 4 springs between walls
# k1=100, k2=200, k3=150, k4=250 (N/m)
# F1=50, F2=0, F3=-30 (N)

# Stiffness matrix
K = np.array([
    [100+200, -200, 0],
    [-200, 200+150, -150],
    [0, -150, 150+250]
], dtype=float)

F = np.array([50.0, 0.0, -30.0])

x = np.linalg.solve(K, F)
print(f'Displacements: x = {x}')

# Verify: residual should be near zero
residual = K @ x - F
print(f'Residual norm: {np.linalg.norm(residual):.2e}')
```

잔차(residual) $`\|[K]\{x\} - \{F\}\|`$는 머신 엡실론 수준이어야 하며, 이는 계산된 해가 원래 방정식을 부동 소수점 정밀도 한계까지 만족함을 확인해 준다. 이 시스템의 변위는 $`x_1 \approx 0.2584`$, $`x_2 \approx 0.1377`$, $`x_3 \approx -0.0234`$ m이고, 잔차 노름은 약 $`4 \times 10^{-16}`$이다.

> **참고:** **강성 행렬의 성질.** 모든 스프링 상수가 양수이고 시스템이 벽에 고정되어 있으면 강성 행렬 $`K`$는 대칭 양정치(SPD, Symmetric Positive Definite)이다. 각 대각 성분은 해당 질량에 연결된 스프링 상수의 합이다. 비대각 성분은 공유 스프링을 통한 질량 사이의 결합을 나타낸다.

SPD 성질은 시스템이 항상 유일한 해를 가지고, 유한한 힘에 대해 변위가 유한하며, 특화된 솔버(예: 촐레스키 분해(Cholesky factorization))가 대칭성을 활용해 효율을 높일 수 있음을 보장한다.

---

<br>

## 3. 단순 가우스 소거법 구현

### 3.1 알고리즘 구현

다음 함수는 전진 소거와 후진 대입으로 이루어진 단순 가우스 소거법을 구현한다. $`A`$와 $`b`$를 따로 관리하지 않도록 **확대 행렬(augmented matrix)** $`[A|b]`$에 대해 연산한다:

```python
import numpy as np

def gauss_naive(A, b):
    """Solve Ax = b using naive Gaussian elimination.

    Args:
        A: n×n coefficient matrix (not modified)
        b: n×1 right-hand side vector (not modified)
    Returns:
        x: solution vector
    """
    A = A.astype(float).copy()
    b = b.astype(float).copy()
    n = len(b)

    # Augmented matrix [A|b]
    Aug = np.hstack([A, b.reshape(-1, 1)])
    nb = n + 1

    # Forward Elimination
    for k in range(n - 1):
        for i in range(k + 1, n):
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]

    # Back Substitution
    x = np.zeros(n)
    x[n-1] = Aug[n-1, nb-1] / Aug[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, nb-1] - Aug[i, i+1:n] @ x[i+1:n]) / Aug[i, i]

    return x
```

주요 구현 세부사항은 다음과 같다:

- **입력 복사**: `A.astype(float).copy()`는 원본 배열이 수정되지 않게 하고, 정수 배열을 실수형으로 바꾸어 나눗셈이 올바르게 수행되게 한다.
- **확대 행렬**: $`[A|b]`$를 배열 하나로 합쳐 행 연산이 계수 행렬과 우변에 일관되게 적용되도록 한다.
- **벡터화된 내부 루프**: `Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]` 줄은 NumPy 슬라이싱으로 행 전체를 한 번에 갱신하여 열에 대한 명시적인 내부 루프를 피한다.

### 3.2 교과서 예제로 테스트

8장 강의 정리본 8절에서 손으로 푼 $`3 \times 3`$ 시스템으로 구현을 시험한다:

```python
# Test with the 3x3 system from the lecture
A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
b = np.array([7.85, -19.3, 71.4])
x = gauss_naive(A, b)
print(f'Solution: {x}')  # [ 3.  -2.5  7. ]

# Compare with NumPy's built-in solver
x_np = np.linalg.solve(A, b)
print(f'NumPy solution: {x_np}')
print(f'Difference: {np.abs(x - x_np)}')
```

배정밀도 계산에서는 손 계산의 7.00003과 달리 머신 엡실론 수준의 오차 안에서 정확한 해 $`(3, -2.5, 7)`$이 나온다. 두 방법은 사실상 같은 결과를 내며, 머신 엡실론 수준의 작은 차이는 내부 구현의 차이에서 비롯된다. `np.linalg.solve`는 부분 피봇팅을 포함하는 LAPACK의 최적화된 루틴을 사용하지만, 이 단순 구현은 피봇팅을 하지 않는다.

---

<br>

## 4. 계산 비용 분석

### 4.1 경험적 시간 측정

크기가 커지는 행렬에 대해 `np.linalg.solve`의 실행 시간을 측정하여 $`O(n^3)`$ 스케일링을 실험으로 확인할 수 있다:

```python
import numpy as np
import time

sizes = [50, 100, 200, 400, 800]
times = []

for n in sizes:
    A = np.random.rand(n, n)
    b = np.random.rand(n)

    start = time.time()
    for _ in range(10):
        np.linalg.solve(A, b)
    elapsed = (time.time() - start) / 10
    times.append(elapsed)
    print(f'n={n:4d}: {elapsed:.6f}s')
```

### 4.2 $`O(n^3)`$ 스케일링 검증

로그-로그 그래프에서 $`O(n^3)`$ 알고리즘은 기울기가 약 3인 직선으로 나타난다:

```python
import matplotlib.pyplot as plt

plt.loglog(sizes, times, 'bo-')
plt.xlabel('Matrix size n')
plt.ylabel('Time (s)')
plt.title('np.linalg.solve: O(n³) scaling')
plt.grid(True)
plt.show()
```

로그-로그 그래프의 기울기는 임의의 두 데이터 점으로 추정할 수 있다. 크기 $`n_1`$에서의 시간이 $`t_1`$이고 크기 $`n_2`$에서의 시간이 $`t_2`$이면 다음과 같다:

```math
\text{slope} = \frac{\log(t_2/t_1)}{\log(n_2/n_1)} \approx 3
```

작은 행렬에서는 오버헤드(함수 호출, 메모리 할당)가 지배적일 수 있어 기울기가 3보다 작게 나타날 수 있다. 매우 큰 행렬에서는 3차 항이 모든 저차 항을 지배하므로 기울기가 3에 가까워진다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:-----|:-----|
| 행렬 형태 | $`n`$개의 연립 일차방정식은 $`[A]\{x\} = \{b\}`$로 간결하게 쓴다. |
| 2x2 시스템의 기하 | 두 직선이 한 점에서 만나면 유일한 해, 평행하면 해 없음, 일치하면 무한히 많은 해이다. |
| 강성 행렬 | 질량-스프링 시스템의 강성 행렬은 대칭 양정치이므로 시스템은 유일한 해를 가진다. |
| 잔차 검증 | 잔차가 머신 엡실론 수준이면 계산된 해가 원래 방정식을 만족함을 확인할 수 있다. |
| 크래머 법칙(Cramer's Rule) | 비용이 $`O((n+1)!)`$이고 피봇팅과 무관하며, 매우 작은 시스템에만 쓴다. |
| 단순 가우스(Naive Gauss) | 비용은 $`O(n^3)`$이고 피봇팅을 하지 않으므로 피봇이 0이면 실패할 수 있다. |
| `np.linalg.solve` | 비용은 $`O(n^3)`$이고 부분 피봇팅을 하는 최적화된 LAPACK 루틴이며, 실무에서 권장된다. |
| 스케일링 검증 | 로그-로그 그래프에서 실행 시간의 기울기가 약 3이면 $`O(n^3)`$ 스케일링이 확인된다. |

---

<br>

## 점검 문제

1. **강성 행렬 구성:** 두 고정벽 사이에 스프링 3개($`k_1 = k_2 = k_3 = 100`$ N/m)로 연결된 질량 2개에 $`F_1 = 10`$ N, $`F_2 = 0`$ N이 작용한다. 강성 행렬을 쓰고 변위를 구하라.

   <details>
   <summary>정답</summary>

   $`[K] = \begin{bmatrix} 200 & -100 \\ -100 & 200 \end{bmatrix}`$이다. $`\det [K] = 200^2 - 100^2 = 30000`$이므로 크래머 법칙으로 $`x_1 = (200 \times 10 - (-100) \times 0)/30000 \approx 0.0667`$ m, $`x_2 = (200 \times 0 - (-100) \times 10)/30000 \approx 0.0333`$ m이다. 힘이 가해진 질량 1이 더 많이 움직이고, 질량 2는 스프링 2를 통해 끌려간다.

   </details>

2. **로그-로그 기울기:** `np.linalg.solve`의 실행 시간이 $`n = 400`$에서 0.010초, $`n = 800`$에서 0.075초였다. 스케일링 기울기를 추정하라.

   <details>
   <summary>정답</summary>

   $`\text{slope} = \log(0.075/0.010)/\log(800/400) = \log 7.5/\log 2 \approx 2.015/0.693 \approx 2.91`$이다. 3에 가까우므로 $`O(n^3)`$ 스케일링과 일치한다.

   </details>

3. **피봇이 0인 경우:** $`A = \begin{bmatrix} 0 & 1 \\ 1 & 1 \end{bmatrix}`$, $`b = (1, 2)`$를 `gauss_naive`에 넣으면 어떻게 되는가? 이 시스템의 해는 무엇인가?

   <details>
   <summary>정답</summary>

   첫 피봇 `Aug[0, 0]`이 0이므로 승수를 계산할 때 0으로 나누게 된다. NumPy는 경고와 함께 `inf`와 `nan`을 만들고, 결과는 `[nan nan]`이 된다. 그러나 시스템 자체는 정칙이며 해는 $`x_2 = 1`$, $`x_1 = 2 - 1 = 1`$이다. 두 행을 바꾸면 풀 수 있으며, 이를 체계화한 것이 9장의 부분 피봇팅이다.

   </details>

4. **잔차의 한계:** 잔차가 머신 엡실론 수준으로 작으면 해의 오차도 항상 작다고 할 수 있는가?

   <details>
   <summary>정답</summary>

   그렇지 않다. 잔차는 계산된 해가 방정식을 얼마나 잘 만족하는지를 보여 줄 뿐이다. 악조건 시스템에서는 해의 상대 오차가 상대 잔차의 최대 조건수배까지 커질 수 있으므로, 잔차가 작아도 해의 오차는 클 수 있다. 조건수는 9장 실습과 11장에서 다룬다.

   </details>
