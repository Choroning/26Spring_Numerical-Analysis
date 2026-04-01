# 제24장 강의 — 경계값 문제(Boundary-Value Problems)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 24

> **선수 지식**: [미적분학] ODE 방법 (제22-23장). [선형대수학] 연립방정식.
>
> **학습 목표**:
> 1. ODE의 경계값 문제를 공식화할 수 있다
> 2. 사격법과 유한 차분법을 적용할 수 있다
> 3. BVP 풀이의 수렴을 분석할 수 있다

---

<br>

## 목차

- [1. 경계값 문제(BVP) 소개](#1-경계값-문제bvp-소개)
  - [1.1 동기 부여 예제: 가열된 막대](#11-동기-부여-예제-가열된-막대)
  - [1.2 지배 방정식](#12-지배-방정식)
  - [1.3 BVP vs. IVP](#13-bvp-vs-ivp)
- [2. 해석해](#2-해석해)
  - [2.1 지배 방정식 재정리](#21-지배-방정식-재정리)
  - [2.2 동차해(Homogeneous Solution)](#22-동차해homogeneous-solution)
  - [2.3 특수해(Particular Solution)](#23-특수해particular-solution)
  - [2.4 전체(일반)해](#24-전체일반해)
  - [2.5 경계 조건으로부터 상수 결정](#25-경계-조건으로부터-상수-결정)
- [3. 사격법(Shooting Method, 절 24.2)](#3-사격법shooting-method-절-242)
  - [3.1 핵심 아이디어: BVP를 IVP로 변환](#31-핵심-아이디어-bvp를-ivp로-변환)
  - [3.2 2차 ODE를 두 개의 1차 ODE로 변환](#32-2차-ode를-두-개의-1차-ode로-변환)
  - [3.3 추정 절차](#33-추정-절차)
  - [3.4 알고리즘 요약](#34-알고리즘-요약)
  - [3.5 Python 구현](#35-python-구현)
- [4. 유한 차분법(Finite Difference Methods, 절 24.3)](#4-유한-차분법finite-difference-methods-절-243)
  - [4.1 핵심 아이디어: 영역의 이산화](#41-핵심-아이디어-영역의-이산화)
  - [4.2 2차 도함수 근사](#42-2차-도함수-근사)
  - [4.3 유한 차분 방정식 생성](#43-유한-차분-방정식-생성)
  - [4.4 미지수와 방정식 개수 세기](#44-미지수와-방정식-개수-세기)
  - [4.5 삼중대각 시스템 조립](#45-삼중대각-시스템-조립)
  - [4.6 행렬 공식화](#46-행렬-공식화)
  - [4.7 Python 구현](#47-python-구현)
- [요약](#요약)

---

<br>

## 1. 경계값 문제(BVP) 소개

### 1.1 동기 부여 예제: 가열된 막대

두 개의 일정 온도 벽 사이에 놓인 길고 가는 막대의 **정상 상태 온도 분포** 를 고려하자:

- 왼쪽 벽은 온도 $T_a$로 유지
- 오른쪽 벽은 온도 $T_b$로 유지
- 주변 공기 온도는 $T_\infty$
- 막대의 길이는 $L$

막대에 두 가지 열전달 메커니즘이 작용한다:

| 메커니즘 | 방향 | 설명 |
|:----------|:----------|:------------|
| **전도(Conduction)** | 막대 방향 (수평) | 막대 재료를 통해 한쪽 끝에서 다른 쪽으로 열이 흐름 |
| **대류(Convection)** | 막대에 수직 (수직) | 막대 표면과 $T_\infty$ 온도의 주변 공기 사이에서 열이 교환됨 |

기본적인 질문은: **막대를 따른 온도 분포 $T(x)$는 무엇인가?**

온도 $T(x)$는 $x = 0$에서 $T_a$부터 $x = L$에서 $T_b$까지 변하며, 이 끝점들 사이의 프로파일은 전도와 대류의 균형에 의존한다.

> **[물리]** 전도는 푸리에 법칙에 의해 지배되며 막대 재료의 열전도도에 의존한다. 대류는 막대 표면과 주변 유체 사이의 열전달 계수에 의존한다. 지배 방정식의 매개변수 $h'$는 대류적 효과 대 전도적 효과의 비율을 캡슐화한다.

### 1.2 지배 방정식

막대를 따른 열전도 및 열확산의 지배 방정식은:

$$0 = \frac{d^2 T}{dx^2} + h'(T_\infty - T) \quad \cdots (**)$$

**경계 조건**:

$$T(x = 0) = T_a, \quad T(x = L) = T_b$$

여기서 $h'$는 대류 열전달 계수, 막대 기하학 및 열전도도와 관련된 매개변수이다.

영역의 끝점에서 두 개의 경계값이 주어지면, 온도 분포 $T(x)$를 풀 수 있다.

### 1.3 BVP vs. IVP

| 특성 | 초기값 문제(IVP) | 경계값 문제(BVP) |
|:--------|:---------------------------|:-----------------------------|
| **조건이 주어지는 곳** | 단일 점 ($x = x_0$) | 둘 이상의 서로 다른 점 |
| **전형적인 조건** | $y(x_0) = y_0$, $y'(x_0) = y'_0$ | $y(a) = \alpha$, $y(b) = \beta$ |
| **풀이 접근법** | $x_0$에서 앞으로 전진 | 단순히 전진 불가; 전역 전략 필요 |
| **차수** | $n$차 ODE에 한 점에서 $n$개의 조건 필요 | 조건이 경계에 걸쳐 분산됨 |

> **[미적분]** 2차 ODE는 유일한 해를 결정하기 위해 정확히 두 개의 조건이 필요하다. IVP에서는 두 조건이 같은 점에서 주어지고 (예: $t = 0$에서 위치와 속도), BVP에서는 다른 공간 위치에서 주어진다 (예: 막대 양쪽 끝의 온도).

---

<br>

## 2. 해석해

### 2.1 지배 방정식 재정리

방정식 $(**)$를 다음과 같이 표현한다:

$$\frac{d^2 T}{dx^2} - h' T = -h' T_\infty$$

이것은 **상수 계수를 가진 2차 선형 ODE** 이며, 우변에 비동차항 $-h' T_\infty$가 있다.

전체 해는 다음의 합이다:
1. **동차** (일반)해
2. **특수** 해

### 2.2 동차해(Homogeneous Solution)

동차 방정식을 풀자:

$$\frac{d^2 T}{dx^2} - h' T = 0 \quad \cdots (***)$$

$T = C e^{\lambda x}$ 형태의 시행해를 가정한다.

$(***)$에 대입하면:

$$\lambda^2 C e^{\lambda x} - h' C e^{\lambda x} = (\lambda^2 - h') C e^{\lambda x} = 0$$

$T = C e^{\lambda x}$는 임의의 (비자명) 값이므로, 특성 방정식이 만족되어야 한다:

$$\lambda^2 - h' = 0$$

$$\lambda = \pm \sqrt{h'}$$

따라서 **동차(일반)해** 는:

$$T_h = C_1 e^{\lambda x} + C_2 e^{-\lambda x}, \quad \lambda = \sqrt{h'}$$

> **[미적분]** 특성 방정식 $\lambda^2 - h' = 0$은 $h' > 0$일 때 (물리적 열전달에서 항상 성립) 두 개의 서로 다른 실근을 갖는다. 이것은 두 개의 선형 독립인 해 $e^{\lambda x}$와 $e^{-\lambda x}$를 제공하며, 이들의 선형결합이 완전한 동차해를 형성한다.

### 2.3 특수해(Particular Solution)

특수해는 완전한 비동차 방정식을 만족하는 상수 함수이다. 직관적으로, $T = T_\infty$이면:

$$\frac{d^2 T_\infty}{dx^2} - h' T_\infty = 0 - h' T_\infty = -h' T_\infty \quad \checkmark$$

따라서 특수해는:

$$T_p = T_\infty$$

### 2.4 전체(일반)해

전체 해는 동차해와 특수해의 합이다:

$$T(x) = T_\infty + C_1 e^{\lambda x} + C_2 e^{-\lambda x}, \quad \lambda = \sqrt{h'}$$

### 2.5 경계 조건으로부터 상수 결정

$C_1$과 $C_2$를 결정하기 위해 두 경계 조건을 적용한다:

**$x = 0$에서:**

$$T_a = T(0) = T_\infty + C_1 + C_2$$

**$x = L$에서:**

$$T_b = T(L) = T_\infty + C_1 e^{\lambda L} + C_2 e^{-\lambda L}$$

이것은 두 미지수 ($C_1$, $C_2$)에 대한 두 방정식 시스템이다.

**$C_2$ 풀기:** 첫 번째 방정식에 $e^{\lambda L}$을 곱하고 두 번째를 빼면:

$$e^{\lambda L} T_a - T_b = T_\infty (e^{\lambda L} - 1) + C_2 (e^{\lambda L} - e^{-\lambda L})$$

$$\therefore \; C_2 = \frac{e^{\lambda L}(T_a - T_\infty) + (T_\infty - T_b)}{e^{\lambda L} - e^{-\lambda L}}$$

**$C_1$ 풀기:** 첫 번째 방정식에 $e^{-\lambda L}$을 곱하고 두 번째를 빼면:

$$e^{-\lambda L} T_a - T_b = e^{-\lambda L} T_\infty - T_\infty + C_1 (e^{-\lambda L} - e^{\lambda L})$$

$$\therefore \; C_1 = \frac{e^{-\lambda L}(T_a - T_\infty) + (T_\infty - T_b)}{e^{-\lambda L} - e^{\lambda L}}$$

> **[미적분]** 분모 $e^{\lambda L} - e^{-\lambda L} = 2\sinh(\lambda L)$은 $\lambda > 0$이고 $L > 0$일 때 항상 0이 아니므로, BVP의 유일한 해를 보장한다. 이것은 온도 분포가 경계 온도에 의해 유일하게 결정된다는 물리적 기대를 확인한다.

---

<br>

## 3. 사격법(Shooting Method, 절 24.2)

### 3.1 핵심 아이디어: BVP를 IVP로 변환

사격법은 경계값 문제를 초기값 문제로 변환한다:

1. 한쪽 끝의 알려진 경계 조건을 초기 조건으로 **사용**
2. 누락된 초기 조건 (그 끝에서의 도함수)을 **추정**
3. 임의의 IVP 솔버 (예: 오일러, RK4)를 사용하여 한 경계에서 다른 경계까지 ODE를 **적분**
4. 먼 경계에서의 계산 결과를 알려진 경계값과 **비교**
5. 추정을 **조정** 하고 먼 경계 조건이 만족될 때까지 반복

"사격법"이라는 이름은 발사체를 쏘는 비유에서 왔다: 발사체가 목표 (먼 경계 조건)에 맞을 때까지 발사각 (초기 기울기 추정)을 조정한다.

### 3.2 2차 ODE를 두 개의 1차 ODE로 변환

지배 방정식에서 출발한다:

$$0 = \frac{d^2 T}{dx^2} + h'(T_\infty - T)$$

새 변수 $z$를 도입한다:

$$z := \frac{dT}{dx}$$

그러면 하나의 2차 ODE가 **두 개의 1차 ODE 시스템** 이 된다:

$$\begin{cases} \dfrac{dT}{dx} = z \\[10pt] \dfrac{dz}{dx} = h'(T - T_\infty) \end{cases}$$

> **[미적분]** 이것은 일반적인 기법이다: 임의의 $n$차 ODE를 $n-1$차까지의 각 도함수에 대한 보조 변수를 도입하여 $n$개의 1차 ODE 시스템으로 변환할 수 있다. 대부분의 수치 ODE 솔버 (오일러, 룽게-쿠타 등)는 1차 시스템용으로 설계되었기 때문에 이것은 필수적이다.

### 3.3 추정 절차

$z(x = 0)$와 $T(x = 0)$ 모두 알면, ODE 시스템을 $x = 0$에서 $x = L$까지 적분할 수 있다.

그러나 BVP에서:
- $T(x = 0) = T_a$는 **알려져 있음** (경계 조건)
- $z(x = 0) = \frac{dT}{dx}\big|_{x=0}$은 **미지**

**절차:**

1. $z(x = 0) = z_{a,1}$을 **추정** (초기 기울기의 첫 번째 추정)
2. 수치 방법 (예: RK4)을 사용하여 1차 ODE 쌍을 $x = 0$에서 $x = L$까지 **적분**
3. $T(x = L) = T_{b,1}$ (오른쪽 경계에서의 계산된 온도) 획득
4. $T_{b,1}$을 알려진 $T_b$와 **비교**
5. 일반적으로, $T_{b,1} \neq T_b$
6. 추정을 **조정**: $z(x = 0) = z_{a,2}$ (새 추정) 설정
7. $T(L) \approx T_b$가 원하는 허용 오차 내에 들 때까지 적분과 비교를 **반복**

추정의 조정은 근 찾기 방법을 사용하여 체계적으로 수행할 수 있다:

- **이분법(Bisection method)**: 두 추정 사이에서 올바른 기울기를 구간으로 잡음
- **할선법(Secant method) / 선형 보간**: 두 추정과 그 결과를 사용하여 다음 추정을 보간

$$z_{a,\text{new}} = z_{a,1} + \frac{(T_b - T_{b,1})(z_{a,2} - z_{a,1})}{T_{b,2} - T_{b,1}}$$

### 3.4 알고리즘 요약

```
1. Set T(0) = Ta  (known boundary condition)
2. Guess z(0) = z_a1
3. Integrate the system of ODEs from x = 0 to x = L
4. Compute T(L) = T_b1
5. If |T_b1 - Tb| < tolerance: DONE
6. Else: adjust z(0) using root-finding, go to step 3
```

### 3.5 Python 구현

```python
import numpy as np

def shooting_method(h_prime, T_inf, Ta, Tb, L, n_steps=100, tol=1e-6, max_iter=50):
    """
    Solve the BVP: d^2T/dx^2 + h'(T_inf - T) = 0
    with T(0) = Ta, T(L) = Tb using the shooting method.

    Parameters
    ----------
    h_prime : float  - heat transfer parameter
    T_inf   : float  - ambient temperature
    Ta      : float  - temperature at x = 0
    Tb      : float  - temperature at x = L
    L       : float  - length of the rod
    n_steps : int    - number of integration steps
    tol     : float  - convergence tolerance
    max_iter: int    - maximum number of shooting iterations
    """
    dx = L / n_steps

    def integrate(z0):
        """Integrate the system using RK4 from x=0 to x=L."""
        T = Ta
        z = z0
        T_vals = [T]
        x_vals = [0.0]

        for i in range(n_steps):
            x = i * dx

            # RK4 for the system dT/dx = z, dz/dx = h'(T - T_inf)
            k1_T = z
            k1_z = h_prime * (T - T_inf)

            k2_T = z + 0.5 * dx * k1_z
            k2_z = h_prime * ((T + 0.5 * dx * k1_T) - T_inf)

            k3_T = z + 0.5 * dx * k2_z
            k3_z = h_prime * ((T + 0.5 * dx * k2_T) - T_inf)

            k4_T = z + dx * k3_z
            k4_z = h_prime * ((T + dx * k3_T) - T_inf)

            T = T + (dx / 6) * (k1_T + 2*k2_T + 2*k3_T + k4_T)
            z = z + (dx / 6) * (k1_z + 2*k2_z + 2*k3_z + k4_z)

            T_vals.append(T)
            x_vals.append((i + 1) * dx)

        return np.array(x_vals), np.array(T_vals), T  # T at x = L

    # Two initial guesses for z(0)
    z1 = 0.0
    _, _, Tb1 = integrate(z1)

    z2 = 1.0
    _, _, Tb2 = integrate(z2)

    # Secant method to find the correct z(0)
    for iteration in range(max_iter):
        if abs(Tb2 - Tb1) < 1e-15:
            break

        # Linear interpolation for next guess
        z_new = z1 + (Tb - Tb1) * (z2 - z1) / (Tb2 - Tb1)
        x_vals, T_vals, Tb_new = integrate(z_new)

        if abs(Tb_new - Tb) < tol:
            print(f"Converged in {iteration + 1} iterations, z(0) = {z_new:.6f}")
            return x_vals, T_vals

        # Update for next iteration
        z1, Tb1 = z2, Tb2
        z2, Tb2 = z_new, Tb_new

    print(f"Warning: did not converge after {max_iter} iterations")
    return x_vals, T_vals
```

---

<br>

## 4. 유한 차분법(Finite Difference Methods, 절 24.3)

### 4.1 핵심 아이디어: 영역의 이산화

BVP를 IVP로 변환하는 대신, 유한 차분(FD) 방법은 격자점 집합에서 지배 방정식을 **직접 이산화** 하여 동시에 풀 수 있는 대수 방정식 시스템을 생성한다.

출발점:

$$\frac{d^2 T}{dx^2} - h' T = -h' T_\infty$$

우변 $-h' T_\infty$는 **알려진 값** 임에 유의하라.

### 4.2 2차 도함수 근사

영역 $[0, L]$을 간격 $\Delta x = L / n$인 $n$개의 동일한 부분 구간으로 나누어, 노드 $x_0, x_1, \ldots, x_n$을 생성한다.

노드 $x_i$에서의 2차 도함수에 대한 **중심 차분 근사** 는:

$$\left. \frac{d^2 T}{dx^2} \right|_{x_i} \approx \frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2}$$

> **[미적분]** 이 근사는 $x_i$ 주변의 전진 및 후진 테일러 전개를 더하여 유도된다. 선행 오차 항은 $O((\Delta x)^2)$이며, 이는 근사가 2차 정확도를 가짐을 의미한다.

### 4.3 유한 차분 방정식 생성

중심 차분 근사를 노드 $x_i$에서의 지배 방정식에 대입하면:

$$\frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2} - h' T_i = -h' T_\infty$$

이 방정식은 각 **내부** 노드 $i = 1, 2, \ldots, n-1$에서 성립한다.

경계 조건과 결합하면:
- $T_0 = T_a$ (왼쪽 경계)
- $T_n = T_b$ (오른쪽 경계)

이것은 $(n+1)$개의 미지수 $T_0, T_1, \ldots, T_n$에 대한 총 $(n+1)$개의 방정식을 준다.

### 4.4 미지수와 방정식 개수 세기

| 항목 | 개수 |
|:-----|:------|
| 총 미지수: $T_0, T_1, \ldots, T_n$ | $n + 1$ |
| 경계 조건: $T_0 = T_a$, $T_n = T_b$ | $2$ |
| 내부 FD 방정식 ($x_1, x_2, \ldots, x_{n-1}$에서) | $n - 1$ |
| **총 방정식** | $n + 1$ |

$T_0$과 $T_n$은 경계 조건에서 알려져 있으므로, 시스템은 $(n - 1)$개의 미지수 $T_1, T_2, \ldots, T_{n-1}$에 대한 $(n - 1)$개의 방정식으로 축소된다.

### 4.5 삼중대각 시스템 조립

편의를 위해 다음 상수를 정의한다:

$$G = 2 + (\Delta x)^2 h', \quad F = h' T_\infty (\Delta x)^2$$

각 노드에서의 FD 방정식에 $(\Delta x)^2$를 곱하고 재정리하면:

**i) 첫 번째 내부 노드 ($i = 1$):**

$$T_0 - 2T_1 + T_2 - (\Delta x)^2 h' T_1 = -h' T_\infty (\Delta x)^2$$

$T_0$은 알려져 있으므로, 우변으로 이동:

$$(-2 - (\Delta x)^2 h') T_1 + T_2 = -h' T_\infty (\Delta x)^2 - T_0$$

$$-G \cdot T_1 + T_2 = -F - T_0$$

**ii) 일반 내부 노드 ($i = 2, 3, \ldots, n-2$):**

$$T_{i-1} - 2T_i + T_{i+1} - (\Delta x)^2 h' T_i = -h' T_\infty (\Delta x)^2$$

$$T_{i-1} - (2 + (\Delta x)^2 h') T_i + T_{i+1} = -h' T_\infty (\Delta x)^2$$

$$T_{i-1} - G \cdot T_i + T_{i+1} = -F$$

**iii) 마지막 내부 노드 ($i = n-1$):**

$$T_{n-2} - 2T_{n-1} + T_n - (\Delta x)^2 h' T_{n-1} = -h' T_\infty (\Delta x)^2$$

$T_n$은 알려져 있으므로, 우변으로 이동:

$$T_{n-2} - (2 + (\Delta x)^2 h') T_{n-1} = -h' T_\infty (\Delta x)^2 - T_n$$

$$T_{n-2} - G \cdot T_{n-1} = -F - T_n$$

### 4.6 행렬 공식화

모든 내부 방정식을 행렬 형태로 결합하면:

$$\begin{bmatrix} -G & 1 & & & \\ 1 & -G & 1 & & \\ & 1 & -G & 1 & \\ & & \ddots & \ddots & \ddots \\ & & & 1 & -G \end{bmatrix} \begin{bmatrix} T_1 \\ T_2 \\ T_3 \\ \vdots \\ T_{n-1} \end{bmatrix} = \begin{bmatrix} -F \\ -F \\ -F \\ \vdots \\ -F \end{bmatrix} - \begin{bmatrix} T_0 \\ 0 \\ 0 \\ \vdots \\ T_n \end{bmatrix}$$

이것은 다음을 사용하여 효율적으로 풀 수 있는 $(n-1)$개 방정식의 **삼중대각 시스템(tridiagonal system)** 이다:
- **토마스 알고리즘(Thomas algorithm)** (삼중대각 행렬 알고리즘): $O(n)$ 시간
- **NumPy의 `linalg.solve`** 또는 표준 선형대수 솔버

> **[선형대수]** 계수 행렬은 삼중대각이고, 대칭이며, 대각 우세 ($G = 2 + (\Delta x)^2 h' > 2$이고 비대각 항은 $\pm 1$)이다. 이것은 시스템이 유일한 해를 가지며 가우스-자이델과 같은 반복법이 수렴함을 보장한다.

### 4.7 Python 구현

```python
import numpy as np

def finite_difference_bvp(h_prime, T_inf, Ta, Tb, L, n=10):
    """
    Solve the BVP: d^2T/dx^2 + h'(T_inf - T) = 0
    with T(0) = Ta, T(L) = Tb using finite differences.

    Parameters
    ----------
    h_prime : float  - heat transfer parameter
    T_inf   : float  - ambient temperature
    Ta      : float  - temperature at x = 0 (T_0)
    Tb      : float  - temperature at x = L (T_n)
    L       : float  - length of the rod
    n       : int    - number of subintervals
    """
    dx = L / n
    G = 2 + dx**2 * h_prime       # diagonal coefficient
    F = h_prime * T_inf * dx**2   # forcing term

    # Build the (n-1) x (n-1) tridiagonal coefficient matrix
    size = n - 1
    A = np.zeros((size, size))
    b = np.full(size, -F)

    for i in range(size):
        A[i, i] = -G                          # main diagonal
        if i > 0:
            A[i, i - 1] = 1.0                 # lower diagonal
        if i < size - 1:
            A[i, i + 1] = 1.0                 # upper diagonal

    # Adjust RHS for boundary conditions
    b[0] -= Ta        # T_0 = Ta is known
    b[-1] -= Tb       # T_n = Tb is known

    # Solve the tridiagonal system
    T_interior = np.linalg.solve(A, b)

    # Assemble full solution including boundary nodes
    T_full = np.concatenate(([Ta], T_interior, [Tb]))
    x_full = np.linspace(0, L, n + 1)

    return x_full, T_full
```

**예제 사용: 두 방법을 해석해와 비교:**

```python
import matplotlib.pyplot as plt

# Problem parameters
h_prime = 0.01    # heat transfer parameter
T_inf = 20.0      # ambient temperature
Ta = 40.0         # left boundary temperature
Tb = 200.0        # right boundary temperature
L = 10.0          # rod length

# Analytical solution
lam = np.sqrt(h_prime)
C1 = (np.exp(-lam*L) * (Ta - T_inf) + (T_inf - Tb)) / (np.exp(-lam*L) - np.exp(lam*L))
C2 = (np.exp(lam*L) * (Ta - T_inf) + (T_inf - Tb)) / (np.exp(lam*L) - np.exp(-lam*L))

x_exact = np.linspace(0, L, 200)
T_exact = T_inf + C1 * np.exp(lam * x_exact) + C2 * np.exp(-lam * x_exact)

# Shooting method
x_shoot, T_shoot = shooting_method(h_prime, T_inf, Ta, Tb, L, n_steps=100)

# Finite difference method
x_fd, T_fd = finite_difference_bvp(h_prime, T_inf, Ta, Tb, L, n=10)

# Plot comparison
plt.figure(figsize=(10, 6))
plt.plot(x_exact, T_exact, 'k-', label='Analytical', linewidth=2)
plt.plot(x_shoot, T_shoot, 'b--o', label='Shooting method', markersize=4)
plt.plot(x_fd, T_fd, 'r--s', label='Finite difference', markersize=6)
plt.xlabel('x')
plt.ylabel('T(x)')
plt.title('BVP Solution: Heated Rod')
plt.legend()
plt.grid(True)
plt.show()
```

---

<br>

## 요약

| 주제 | 핵심 사항 |
|:------|:-----------|
| **경계값 문제** | 둘 이상의 서로 다른 점에서 조건이 주어진 ODE. 모든 조건이 한 점에 있는 IVP와 다름 |
| **가열된 막대 BVP** | $\dfrac{d^2 T}{dx^2} + h'(T_\infty - T) = 0$, $T(0) = T_a$, $T(L) = T_b$ |
| **해석해** | $T(x) = T_\infty + C_1 e^{\lambda x} + C_2 e^{-\lambda x}$, $\lambda = \sqrt{h'}$이며 상수는 경계 조건에서 결정 |
| **사격법** | BVP를 IVP로 변환: 미지의 초기 기울기 $z(0) = dT/dx\|_{x=0}$를 추정; 전진 적분; $T(L) = T_b$가 될 때까지 근 찾기로 추정 조정 |
| **ODE 변환** | 단일 $n$차 ODE를 보조 변수 도입으로 $n$개의 1차 ODE로 변환 가능 |
| **유한 차분법** | 영역을 격자로 이산화; $d^2T/dx^2 \approx \frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2}$로 근사; 결과 삼중대각 선형 시스템을 풀음 |
| **FD 행렬 구조** | 대각 $-G = -(2 + (\Delta x)^2 h')$, 비대각 $1$인 삼중대각; 토마스 알고리즘으로 $O(n)$에 풀이 가능 |
| **핵심 차이** | 사격법은 반복적 (추정-확인); FD 방법은 직접 대수 시스템 생성 (선형 BVP에 대해 한 번에 풀이) |
