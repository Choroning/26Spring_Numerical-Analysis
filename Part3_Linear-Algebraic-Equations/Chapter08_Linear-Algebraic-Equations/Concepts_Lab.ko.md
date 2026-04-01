# 제8장 실습 — 연립방정식과 가우스 소거법(Linear Systems and Gaussian Elimination)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 8

> **선수 지식**: [프로그래밍언어] MATLAB/Python. [선형대수학] 행렬 연산. [미적분학] 선형 시스템 (제1-7장).
>
> **학습 목표**:
> 1. 공학 문제를 연립일차방정식으로 공식화할 수 있다
> 2. 해의 존재성과 유일성 조건을 식별할 수 있다
> 3. 연립방정식의 그래프적 해석을 적용할 수 있다

---

<br>

## 목차

- [1. $Ax = b$ 설정하기](#1-ax--b-설정하기)
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
  - [4.2 O(n^3) 스케일링 검증](#42-on3-스케일링-검증)
- [요약](#요약)

---

<br>

## 1. $Ax = b$ 설정하기

### 1.1 행렬 형태

$n$개의 미지수에 대한 $n$개의 연립 일차방정식은 행렬-벡터 형태로 간결하게 표현할 수 있습니다:

$$[A]\{x\} = \{b\}$$

여기서:
- $[A]$는 알려진 계수를 포함하는 $n \times n$ **계수 행렬(coefficient matrix)**
- $\{x\}$는 구하고자 하는 변수를 포함하는 $n \times 1$ **미지수 벡터(unknown vector)**
- $\{b\}$는 알려진 상수를 포함하는 $n \times 1$ **우변 벡터(right-hand side, RHS vector)**

예를 들어, 다음 시스템:

$$2x_1 + x_2 = 5$$
$$x_1 - x_2 = 1$$

은 다음과 같이 표현됩니다:

$$\begin{bmatrix} 2 & 1 \\ 1 & -1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 5 \\ 1 \end{bmatrix}$$

### 1.2 기하학적 해석 (2x2)

$2 \times 2$ 시스템의 경우, 각 방정식은 평면에서 하나의 직선을 정의합니다. 유일한 해, 해 없음, 무한한 해 — 세 가지 가능한 결과를 직접 시각화할 수 있습니다:

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

유일한 해의 경우 교차점 $(2, 1)$이 검은 점으로 표시됩니다. 평행선의 경우 기울기는 같지만 절편이 다른 두 직선이 표시됩니다 — 이들은 결코 만나지 않습니다. 무한한 해의 경우 정확히 같은 직선을 나타내는 두 방정식이 표시됩니다(하나는 단순히 다른 것의 두 배).

---

<br>

## 2. 질량-스프링 시스템 응용

### 2.1 문제 설정

두 고정벽 사이에 4개의 스프링으로 연결된 3개의 질량으로 구성된 시스템을 고려합니다. 스프링 상수는 $k_1 = 100$, $k_2 = 200$, $k_3 = 150$, $k_4 = 250$ N/m이고, 가해진 힘은 $F_1 = 50$, $F_2 = 0$, $F_3 = -30$ N입니다.

각 질량에 훅의 법칙($F = kx$)과 힘 평형을 적용하면 강성 행렬 방정식 $[K]\{x\} = \{F\}$가 됩니다:

$$\begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} F_1 \\ F_2 \\ F_3 \end{bmatrix}$$

각 대각 성분 $(k_i + k_{i+1})$은 질량 $i$에 작용하는 총 강성을 나타냅니다. 비대각 성분 $(-k_j)$은 공유 스프링을 통한 인접 질량 간의 결합을 나타냅니다.

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

잔차(residual) $\|[K]\{x\} - \{F\}\|$는 기계 엡실론(machine epsilon) 근처여야 하며, 이는 계산된 해가 원래 방정식을 완전한 부동소수점 정밀도로 만족함을 확인합니다.

> **[물리]** 강성 행렬 K는 물리적으로 유효한 스프링 시스템에 대해 대칭 양정치(SPD, Symmetric Positive Definite)입니다. 각 대각 성분은 해당 질량에 연결된 스프링 상수의 합입니다. 비대각 성분은 공유 스프링을 통한 질량 간의 결합을 나타냅니다.

SPD 속성은 시스템이 항상 유일한 해를 가지며, 유한한 힘에 대해 변위가 유한하고, 특화된 솔버(예: 촐레스키 분해(Cholesky decomposition))가 대칭성을 활용하여 효율성을 높일 수 있음을 보장합니다.

---

<br>

## 3. 단순 가우스 소거법 구현

### 3.1 알고리즘 구현

다음 함수는 전진 소거와 후진 대입을 포함한 단순 가우스 소거법을 구현합니다. $A$와 $b$를 별도로 관리하지 않도록 **확대 행렬(augmented matrix)** $[A|b]$에 대해 연산합니다:

```python
import numpy as np

def gauss_naive(A, b):
    """Solve Ax = b using naive Gaussian elimination.

    Args:
        A: n×n coefficient matrix (will be modified)
        b: n×1 right-hand side vector (will be modified)
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

주요 구현 세부사항:

- **입력 복사**: `A.astype(float).copy()`는 원본 배열이 수정되지 않도록 하고, 정수 배열이 올바른 나눗셈을 위해 실수형으로 승격되도록 합니다.
- **확대 행렬**: $[A|b]$를 단일 배열로 결합하여 행 연산이 계수 행렬과 우변에 일관되게 적용되도록 합니다.
- **벡터화된 내부 루프**: `Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]` 줄은 NumPy 슬라이싱을 사용하여 전체 행을 한 번에 갱신하며, 열에 대한 명시적 내부 루프를 피합니다.

### 3.2 교과서 예제로 테스트

```python
# Test with the 3x3 system from the lecture
A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
b = np.array([7.85, -19.3, 71.4])
x = gauss_naive(A, b)
print(f'Solution: {x}')  # [3.0, -2.5, 7.0003...]

# Compare with NumPy's built-in solver
x_np = np.linalg.solve(A, b)
print(f'NumPy solution: {x_np}')
print(f'Difference: {np.abs(x - x_np)}')
```

두 방법 모두 본질적으로 같은 결과를 산출해야 합니다. 기계 엡실론 수준의 작은 수치적 차이는 내부 구현의 차이에서 비롯됩니다 — `np.linalg.solve`는 부분 피봇팅을 포함하는 LAPACK의 최적화된 루틴을 사용하지만, 우리의 단순 구현은 그렇지 않습니다.

---

<br>

## 4. 계산 비용 분석

### 4.1 경험적 시간 측정

크기가 증가하는 행렬에 대해 `np.linalg.solve`의 시간을 측정하여 $O(n^3)$ 스케일링을 경험적으로 검증할 수 있습니다:

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

### 4.2 O(n^3) 스케일링 검증

로그-로그 그래프에서 $O(n^3)$ 알고리즘은 기울기가 약 3인 직선으로 나타납니다:

```python
import matplotlib.pyplot as plt

plt.loglog(sizes, times, 'bo-')
plt.xlabel('Matrix size n')
plt.ylabel('Time (s)')
plt.title('np.linalg.solve: O(n³) scaling')
plt.grid(True)
plt.show()
```

로그-로그 그래프의 기울기는 임의의 두 데이터 점에서 추정할 수 있습니다. 크기 $n_1$에서의 시간이 $t_1$이고 크기 $n_2$에서의 시간이 $t_2$이면:

$$\text{slope} = \frac{\log(t_2/t_1)}{\log(n_2/n_1)} \approx 3$$

작은 행렬의 경우 오버헤드 비용(함수 호출, 메모리 할당)이 지배적일 수 있어 기울기가 3보다 작게 나타날 수 있습니다. 매우 큰 행렬의 경우 3차 항이 모든 저차 항을 지배하므로 기울기가 3에 접근합니다.

---

<br>

## 요약

| 방법 | 복잡도 | 피봇팅 | 비고 |
|:-----|:-------|:-------|:-----|
| 크래머 법칙(Cramer's Rule) | $O(n!)$ | 해당 없음 | 매우 작은 시스템에만 사용 |
| 단순 가우스(Naive Gauss) | $O(n^3)$ | 없음 | 피봇이 0이면 실패 가능 |
| `np.linalg.solve` | $O(n^3)$ | 있음 | 최적화된 LAPACK, 권장 |

---
