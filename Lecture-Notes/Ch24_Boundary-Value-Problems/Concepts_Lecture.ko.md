# 24장 — 경계값 문제

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 24

> **학습 목표:**
> 1. ODE의 경계값 문제를 공식화할 수 있다
> 2. 선형 BVP의 해석해를 동차해와 특수해의 합으로 구할 수 있다
> 3. 사격법과 유한 차분법을 적용할 수 있다
> 4. BVP 풀이의 정확도와 수렴을 분석할 수 있다

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
- [3. 사격법(Shooting Method, 24.2절)](#3-사격법shooting-method-242절)
  - [3.1 핵심 아이디어: BVP를 IVP로 변환](#31-핵심-아이디어-bvp를-ivp로-변환)
  - [3.2 2차 ODE를 두 개의 1차 ODE로 변환](#32-2차-ode를-두-개의-1차-ode로-변환)
  - [3.3 추정 절차](#33-추정-절차)
  - [3.4 알고리즘 요약](#34-알고리즘-요약)
  - [3.5 Python 구현](#35-python-구현)
- [4. 유한 차분법(Finite Difference Methods, 24.3절)](#4-유한-차분법finite-difference-methods-243절)
  - [4.1 핵심 아이디어: 영역의 이산화](#41-핵심-아이디어-영역의-이산화)
  - [4.2 2계도함수 근사](#42-2계도함수-근사)
  - [4.3 유한 차분 방정식 생성](#43-유한-차분-방정식-생성)
  - [4.4 미지수와 방정식 개수 세기](#44-미지수와-방정식-개수-세기)
  - [4.5 삼중대각 시스템 조립](#45-삼중대각-시스템-조립)
  - [4.6 행렬 공식화](#46-행렬-공식화)
  - [4.7 Python 구현](#47-python-구현)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 경계값 문제(BVP) 소개

### 1.1 동기 부여 예제: 가열된 막대

두 개의 일정 온도 벽 사이에 놓인 길고 가는 막대의 **정상 상태 온도 분포**를 고려한다:

- 왼쪽 벽은 온도 $`T_a`$로 유지된다.
- 오른쪽 벽은 온도 $`T_b`$로 유지된다.
- 주변 공기 온도는 $`T_\infty`$이다.
- 막대의 길이는 $`L`$이다.

```text
                 주변 공기 T∞
 ┃   ↕    ↕    ↕    ↕    ↕    ↕    ↕    ↕   ┃   ↕ 대류
 ┃██████████████████████████████████████████┃   ◀─▶ 전도(막대 방향)
T_a                                        T_b
 ├──────────────────────────────────────────┤──▶ x
 0                                          L
```

*그림 1. 두 벽 사이에 놓인 가열된 막대의 전도와 대류*

막대에 두 가지 열전달 메커니즘이 작용한다:

| 메커니즘 | 방향 | 설명 |
|:----------|:----------|:------------|
| **전도(Conduction)** | 막대 방향 (수평) | 막대 재료를 통해 한쪽 끝에서 다른 쪽으로 열이 흐름 |
| **대류(Convection)** | 막대에 수직 (수직) | 막대 표면과 $`T_\infty`$ 온도의 주변 공기 사이에서 열이 교환됨 |

*표 1. 가열된 막대에 작용하는 두 열전달 메커니즘*

기본적인 질문은 다음과 같다: **막대를 따른 온도 분포 $`T(x)`$는 무엇인가?**

온도 $`T(x)`$는 $`x = 0`$에서 $`T_a`$부터 $`x = L`$에서 $`T_b`$까지 변하며, 이 끝점들 사이의 프로파일은 전도와 대류의 균형에 의존한다.

> **참고:** **물리 배경.** 전도는 푸리에 법칙을 따르며 막대 재료의 열전도도에 의존한다. 대류는 막대 표면과 주변 유체 사이의 열전달 계수에 의존한다. 지배 방정식의 매개변수 $`h'`$는 대류 효과와 전도 효과의 상대적 크기를 나타낸다.

### 1.2 지배 방정식

막대를 따른 열전도와 주변 공기로의 대류를 함께 나타내는 지배 방정식은 다음과 같다:

```math
0 = \frac{d^2 T}{dx^2} + h'(T_\infty - T) \quad \cdots (**)
```

**경계 조건**은 다음과 같다:

```math
T(x = 0) = T_a, \quad T(x = L) = T_b
```

여기서 $`h'`$는 대류 열전달 계수, 막대의 기하 형상, 열전도도와 관련된 매개변수이다.

영역의 끝점에서 두 개의 경계값이 주어지면, 온도 분포 $`T(x)`$를 풀 수 있다.

### 1.3 BVP vs. IVP

| 특성 | 초기값 문제(IVP) | 경계값 문제(BVP) |
|:--------|:---------------------------|:-----------------------------|
| **조건이 주어지는 곳** | 단일 점 ($`x = x_0`$) | 둘 이상의 서로 다른 점 |
| **전형적인 조건** | $`y(x_0) = y_0`$, $`y'(x_0) = y'_0`$ | $`y(a) = \alpha`$, $`y(b) = \beta`$ |
| **풀이 접근법** | $`x_0`$에서 앞으로 전진 | 한쪽 끝에서 단순히 전진할 수 없음; 반복 추정(사격법)이나 전역 풀이(유한 차분법) 필요 |
| **차수** | $`n`$차 ODE에 한 점에서 $`n`$개의 조건 필요 | 조건이 경계에 걸쳐 분산됨 |

*표 2. 초기값 문제와 경계값 문제의 비교*

> **[미적분학]** 2차 ODE의 해를 하나로 정하려면 일반적으로 두 개의 조건이 필요하다. IVP에서는 두 조건이 같은 점에서 주어지고(예: $`t = 0`$에서 위치와 속도), BVP에서는 다른 공간 위치에서 주어진다(예: 막대 양쪽 끝의 온도). 계수가 연속인 선형 IVP는 해가 하나로 정해지지만, BVP는 조건이 두 개여도 해가 없거나 무수히 많을 수 있다. 예를 들어 $`y'' + y = 0`$, $`y(0) = 0`$, $`y(\pi) = 0`$은 모든 상수 $`C`$에 대해 $`y = C \sin x`$를 해로 가진다.

---

<br>

## 2. 해석해

### 2.1 지배 방정식 재정리

방정식 $`(**)`$를 다음과 같이 표현한다:

```math
\frac{d^2 T}{dx^2} - h' T = -h' T_\infty
```

이것은 **상수 계수를 가진 2차 선형 ODE**이며, 우변에 비동차항 $`-h' T_\infty`$가 있다.

전체 해는 다음의 합이다:
1. **동차해**
2. **특수해**

### 2.2 동차해(Homogeneous Solution)

먼저 동차 방정식을 푼다:

```math
\frac{d^2 T}{dx^2} - h' T = 0 \quad \cdots (***)
```

$`T = C e^{\lambda x}`$ 형태의 시행해를 가정한다.

$`(***)`$에 대입하면:

```math
\lambda^2 C e^{\lambda x} - h' C e^{\lambda x} = (\lambda^2 - h') C e^{\lambda x} = 0
```

비자명해($`C \neq 0`$)를 구하고 있고 $`e^{\lambda x}`$는 0이 되지 않으므로, 특성 방정식이 만족되어야 한다:

```math
\lambda^2 - h' = 0
```

```math
\lambda = \pm \sqrt{h'}
```

따라서 **동차해**는:

```math
T_h = C_1 e^{\lambda x} + C_2 e^{-\lambda x}, \quad \lambda = \sqrt{h'}
```

> **[미적분학]** 특성 방정식 $`\lambda^2 - h' = 0`$은 $`h' > 0`$일 때(물리적 열전달에서 항상 성립) 두 개의 서로 다른 실근을 갖는다. 이것은 두 개의 선형독립인 해 $`e^{\lambda x}`$와 $`e^{-\lambda x}`$를 제공하며, 이들의 선형결합이 완전한 동차해를 형성한다.

### 2.3 특수해(Particular Solution)

특수해는 비동차 방정식 전체를 만족하는 해 하나이다. 우변이 상수이므로 상수 함수를 시도할 수 있으며, $`T = T_\infty`$이면 다음이 성립한다:

```math
\frac{d^2 T_\infty}{dx^2} - h' T_\infty = 0 - h' T_\infty = -h' T_\infty \quad \checkmark
```

따라서 특수해는:

```math
T_p = T_\infty
```

### 2.4 전체(일반)해

전체 해는 동차해와 특수해의 합이다:

```math
T(x) = T_\infty + C_1 e^{\lambda x} + C_2 e^{-\lambda x}, \quad \lambda = \sqrt{h'}
```

### 2.5 경계 조건으로부터 상수 결정

$`C_1`$과 $`C_2`$를 결정하기 위해 두 경계 조건을 적용한다:

**$`x = 0`$에서:**

```math
T_a = T(0) = T_\infty + C_1 + C_2
```

**$`x = L`$에서:**

```math
T_b = T(L) = T_\infty + C_1 e^{\lambda L} + C_2 e^{-\lambda L}
```

이것은 두 미지수($`C_1`$, $`C_2`$)에 대한 두 방정식 시스템이다.

**$`C_2`$ 풀기:** 첫 번째 방정식에 $`e^{\lambda L}`$을 곱하고 두 번째를 빼면:

```math
e^{\lambda L} T_a - T_b = T_\infty (e^{\lambda L} - 1) + C_2 (e^{\lambda L} - e^{-\lambda L})
```

```math
\therefore \; C_2 = \frac{e^{\lambda L}(T_a - T_\infty) + (T_\infty - T_b)}{e^{\lambda L} - e^{-\lambda L}}
```

**$`C_1`$ 풀기:** 첫 번째 방정식에 $`e^{-\lambda L}`$을 곱하고 두 번째를 빼면:

```math
e^{-\lambda L} T_a - T_b = e^{-\lambda L} T_\infty - T_\infty + C_1 (e^{-\lambda L} - e^{\lambda L})
```

```math
\therefore \; C_1 = \frac{e^{-\lambda L}(T_a - T_\infty) + (T_\infty - T_b)}{e^{-\lambda L} - e^{\lambda L}}
```

> **[미적분학]** 분모 $`e^{\lambda L} - e^{-\lambda L} = 2\sinh(\lambda L)`$은 $`\lambda > 0`$이고 $`L > 0`$일 때 항상 0이 아니므로, BVP의 유일한 해를 보장한다. 이것은 온도 분포가 경계 온도에 의해 유일하게 결정된다는 물리적 기대를 확인한다.

---

<br>

## 3. 사격법(Shooting Method, 24.2절)

### 3.1 핵심 아이디어: BVP를 IVP로 변환

사격법은 다음 과정으로 경계값 문제를 초기값 문제로 변환한다:

1. 한쪽 끝의 알려진 경계 조건을 초기 조건으로 **사용한다**.
2. 누락된 초기 조건(그 끝에서의 도함수)을 **추정한다**.
3. 임의의 IVP 솔버(예: 오일러, RK4)를 사용하여 한 경계에서 다른 경계까지 ODE를 **적분한다**.
4. 먼 경계에서의 계산 결과를 알려진 경계값과 **비교한다**.
5. 추정을 **조정하고**, 먼 경계 조건이 만족될 때까지 반복한다.

"사격법"이라는 이름은 발사체를 쏘는 비유에서 왔다: 발사체가 목표(먼 경계 조건)에 맞을 때까지 발사각(초기 기울기 추정)을 조정한다.

### 3.2 2차 ODE를 두 개의 1차 ODE로 변환

지배 방정식에서 출발한다:

```math
0 = \frac{d^2 T}{dx^2} + h'(T_\infty - T)
```

새 변수 $`z`$를 도입한다:

```math
z := \frac{dT}{dx}
```

그러면 하나의 2차 ODE가 **두 개의 1차 ODE 시스템**이 된다:

```math
\begin{cases} \dfrac{dT}{dx} = z \\[10pt] \dfrac{dz}{dx} = h'(T - T_\infty) \end{cases}
```

> **[미적분학]** 이것은 일반적인 기법이다: 임의의 $`n`$차 ODE를 $`n-1`$차까지의 각 도함수에 대한 보조 변수를 도입하여 $`n`$개의 1차 ODE 시스템으로 변환할 수 있다. 대부분의 수치 ODE 솔버(오일러, 룽게-쿠타 등)는 1차 시스템용으로 설계되었기 때문에 이것은 필수적이다.

### 3.3 추정 절차

$`z(x = 0)`$와 $`T(x = 0)`$를 모두 알면, ODE 시스템을 $`x = 0`$에서 $`x = L`$까지 적분할 수 있다.

그러나 BVP에서는 다음과 같다:
- $`T(x = 0) = T_a`$는 경계 조건으로 **알려져 있다**.
- $`z(x = 0) = \frac{dT}{dx}\big|_{x=0}`$은 **알려져 있지 않다**.

**절차:**

1. 초기 기울기를 $`z(x = 0) = z_{a,1}`$로 **추정한다**(첫 번째 추정).
2. 수치 방법(예: RK4)을 사용하여 1차 ODE 쌍을 $`x = 0`$에서 $`x = L`$까지 **적분한다**.
3. 오른쪽 경계에서 계산된 온도 $`T(x = L) = T_{b,1}`$을 얻는다.
4. $`T_{b,1}`$을 알려진 $`T_b`$와 **비교한다**.
5. 일반적으로 $`T_{b,1} \neq T_b`$이다.
6. 추정을 $`z(x = 0) = z_{a,2}`$(새 추정)로 **조정한다**.
7. $`T(L) \approx T_b`$가 원하는 허용 오차 안에 들 때까지 적분과 비교를 **반복한다**.

추정의 조정은 근 찾기 방법을 사용하여 체계적으로 수행할 수 있다:

- **이분법(Bisection method)**: 올바른 기울기를 두 추정값 사이의 구간에 가두고 그 구간을 줄여 나간다.
- **할선법(Secant method) / 선형 보간**: 두 추정값과 그 결과를 사용하여 다음 추정값을 보간한다.

```math
z_{a,\text{new}} = z_{a,1} + \frac{(T_b - T_{b,1})(z_{a,2} - z_{a,1})}{T_{b,2} - T_{b,1}}
```

> **참고:** 이 예제처럼 ODE가 선형이면 $`T(L)`$은 $`z(0)`$의 일차 함수이다. 따라서 두 번의 시행 결과를 선형 보간한 $`z_{a,\text{new}}`$는 적분의 수치 오차를 제외하면 정확한 초기 기울기이며, 더 반복하지 않아도 경계 조건을 만족한다. 4.7절의 예제 값으로 3.5절의 구현을 실행하면 첫 보간에서 수렴한다. 비선형 ODE에서는 이 관계가 성립하지 않으므로 보간을 반복해야 한다.

### 3.4 알고리즘 요약

```text
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

## 4. 유한 차분법(Finite Difference Methods, 24.3절)

### 4.1 핵심 아이디어: 영역의 이산화

BVP를 IVP로 변환하는 대신, 유한 차분(FD) 방법은 격자점 집합에서 지배 방정식을 **직접 이산화**하여 동시에 풀 수 있는 대수 방정식 시스템을 생성한다.

출발점:

```math
\frac{d^2 T}{dx^2} - h' T = -h' T_\infty
```

우변 $`-h' T_\infty`$는 **알려진 값**임에 유의하라.

### 4.2 2계도함수 근사

영역 $`[0, L]`$을 간격 $`\Delta x = L / n`$인 $`n`$개의 동일한 부분 구간으로 나누어, 노드 $`x_0, x_1, \ldots, x_n`$을 생성한다.

```text
T_0 = T_a                                         T_n = T_b
  ●───────○───────○── ... ──○───────○───────○── ... ──●
 x_0     x_1     x_2     x_{i-1}   x_i   x_{i+1}     x_n
  ├── Δx ─┤
 ● 경계 노드(값을 앎)   ○ 내부 노드(미지수)
```

*그림 2. 유한 차분 격자의 경계 노드와 내부 노드*

노드 $`x_i`$에서의 2계도함수에 대한 **중심 차분 근사**는:

```math
\left. \frac{d^2 T}{dx^2} \right|_{x_i} \approx \frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2}
```

> **참고:** 이 근사는 4장 3.4절과 21장 3.1절에서 유도한 2계도함수의 중심 차분 공식으로, $`x_i`$ 주변의 전진 및 후진 테일러 전개를 더하여 얻는다. 선행 오차 항은 $`O((\Delta x)^2)`$이며, 이는 근사가 2차 정확도를 가짐을 의미한다.

### 4.3 유한 차분 방정식 생성

중심 차분 근사를 노드 $`x_i`$에서의 지배 방정식에 대입하면:

```math
\frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2} - h' T_i = -h' T_\infty
```

이 방정식은 각 **내부** 노드 $`i = 1, 2, \ldots, n-1`$에서 성립한다.

경계 조건과 결합하면:
- $`T_0 = T_a`$ (왼쪽 경계)
- $`T_n = T_b`$ (오른쪽 경계)

이것은 $`(n+1)`$개의 미지수 $`T_0, T_1, \ldots, T_n`$에 대한 총 $`(n+1)`$개의 방정식을 준다.

### 4.4 미지수와 방정식 개수 세기

| 항목 | 개수 |
|:-----|:------|
| 총 미지수: $`T_0, T_1, \ldots, T_n`$ | $`n + 1`$ |
| 경계 조건: $`T_0 = T_a`$, $`T_n = T_b`$ | $`2`$ |
| 내부 FD 방정식 ($`x_1, x_2, \ldots, x_{n-1}`$에서) | $`n - 1`$ |
| **총 방정식** | $`n + 1`$ |

*표 3. 미지수와 방정식의 개수*

$`T_0`$과 $`T_n`$은 경계 조건에서 알려져 있으므로, 시스템은 $`(n - 1)`$개의 미지수 $`T_1, T_2, \ldots, T_{n-1}`$에 대한 $`(n - 1)`$개의 방정식으로 축소된다.

### 4.5 삼중대각 시스템 조립

편의를 위해 다음 상수를 정의한다:

```math
G = 2 + (\Delta x)^2 h', \quad F = h' T_\infty (\Delta x)^2
```

각 노드에서의 FD 방정식에 $`(\Delta x)^2`$를 곱하고 재정리하면:

**i) 첫 번째 내부 노드 ($`i = 1`$):**

```math
T_0 - 2T_1 + T_2 - (\Delta x)^2 h' T_1 = -h' T_\infty (\Delta x)^2
```

$`T_0`$은 알려져 있으므로, 우변으로 이동:

```math
(-2 - (\Delta x)^2 h') T_1 + T_2 = -h' T_\infty (\Delta x)^2 - T_0
```

```math
-G \cdot T_1 + T_2 = -F - T_0
```

**ii) 일반 내부 노드 ($`i = 2, 3, \ldots, n-2`$):**

```math
T_{i-1} - 2T_i + T_{i+1} - (\Delta x)^2 h' T_i = -h' T_\infty (\Delta x)^2
```

```math
T_{i-1} - (2 + (\Delta x)^2 h') T_i + T_{i+1} = -h' T_\infty (\Delta x)^2
```

```math
T_{i-1} - G \cdot T_i + T_{i+1} = -F
```

**iii) 마지막 내부 노드 ($`i = n-1`$):**

```math
T_{n-2} - 2T_{n-1} + T_n - (\Delta x)^2 h' T_{n-1} = -h' T_\infty (\Delta x)^2
```

$`T_n`$은 알려져 있으므로, 우변으로 이동:

```math
T_{n-2} - (2 + (\Delta x)^2 h') T_{n-1} = -h' T_\infty (\Delta x)^2 - T_n
```

```math
T_{n-2} - G \cdot T_{n-1} = -F - T_n
```

$`n = 2`$이면 내부 노드가 $`T_1`$ 하나뿐이므로 i)과 iii)이 같은 식이 되어 $`-G \cdot T_1 = -F - T_0 - T_2`$가 된다. 일반 내부 노드의 식 ii)는 $`n \geq 4`$일 때만 나타난다.

### 4.6 행렬 공식화

모든 내부 방정식을 행렬 형태로 결합하면:

```math
\begin{bmatrix} -G & 1 & & & \\ 1 & -G & 1 & & \\ & 1 & -G & 1 & \\ & & \ddots & \ddots & \ddots \\ & & & 1 & -G \end{bmatrix} \begin{bmatrix} T_1 \\ T_2 \\ T_3 \\ \vdots \\ T_{n-1} \end{bmatrix} = \begin{bmatrix} -F \\ -F \\ -F \\ \vdots \\ -F \end{bmatrix} - \begin{bmatrix} T_0 \\ 0 \\ 0 \\ \vdots \\ T_n \end{bmatrix}
```

이것은 $`(n-1)`$개 방정식의 **삼중대각 시스템(tridiagonal system)** 이며, 다음 방법으로 풀 수 있다:
- **토마스 알고리즘(Thomas algorithm)** (삼중대각 행렬 알고리즘): $`O(n)`$ 시간에 효율적으로 푼다.
- **NumPy의 `linalg.solve`** 또는 표준 선형대수 솔버: 사용은 간단하지만 삼중대각 구조를 활용하지 않고 밀집 행렬로 풀므로 $`O(n^3)`$ 시간이 걸린다.

> **참고:** 계수 행렬은 삼중대각이고, 대칭이며, 엄밀한 대각 우위이다. $`h' > 0`$이면 각 대각 원소의 절댓값 $`G = 2 + (\Delta x)^2 h'`$가 2보다 크고, 각 행의 비대각 원소는 값이 1인 원소가 많아야 두 개이기 때문이다. 엄밀한 대각 우위 행렬은 정칙이므로 시스템은 유일한 해를 가지며, 12장 3절의 수렴 조건에 따라 가우스-자이델과 같은 반복법이 수렴한다.

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

**사용 예:** 두 방법을 해석해와 비교한다.

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

이 값으로 실행하면 해석해의 상수는 $`C_1 \approx 73.4523`$, $`C_2 \approx -53.4523`$이다. 사격법은 첫 보간에서 $`z(0) \approx 12.6905`$로 수렴하고, $`n = 10`$인 유한 차분 해는 $`x = 5`$에서 $`108.6904`$(해석해 $`108.6819`$)이며 모든 노드에서 오차가 $`0.009`$ 미만이다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:------|:-----------|
| 경계값 문제 | 둘 이상의 서로 다른 점에서 조건이 주어지는 ODE 문제이며, 모든 조건이 한 점에서 주어지는 IVP와 다르다. |
| 가열된 막대 BVP | 지배 방정식은 $`\dfrac{d^2 T}{dx^2} + h'(T_\infty - T) = 0`$이고, 경계 조건은 $`T(0) = T_a`$, $`T(L) = T_b`$이다. |
| 해석해 | $`T(x) = T_\infty + C_1 e^{\lambda x} + C_2 e^{-\lambda x}`$($`\lambda = \sqrt{h'}`$)이며, 상수 $`C_1`$, $`C_2`$는 경계 조건에서 결정한다. |
| 사격법 | 미지의 초기 기울기 $`z(0) = dT/dx\vert_{x=0}`$를 추정해 BVP를 IVP로 바꾸고, 전진 적분한 뒤 $`T(L) = T_b`$가 될 때까지 근 찾기로 추정을 조정한다. |
| ODE 변환 | 단일 $`n`$차 ODE는 보조 변수를 도입하여 $`n`$개의 1차 ODE 시스템으로 바꿀 수 있다. |
| 유한 차분법 | 영역을 격자로 이산화하고 $`d^2T/dx^2 \approx \frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2}`$로 근사한 뒤, 그 결과로 생기는 삼중대각 선형 시스템을 푼다. |
| FD 행렬 구조 | 대각 원소가 $`-G = -(2 + (\Delta x)^2 h')`$이고 비대각 원소가 $`1`$인 삼중대각 행렬이며, 토마스 알고리즘으로 $`O(n)`$에 풀 수 있다. |
| 핵심 차이 | 사격법은 추정과 확인을 반복하는 방법이고, FD 방법은 대수 시스템을 직접 만들어 선형 BVP를 한 번에 푼다. |

---

<br>

## 점검 문제

1. **BVP와 IVP:** BVP를 IVP처럼 한쪽 끝에서 바로 전진 적분할 수 없는 이유를 설명하고, 사격법이 이 문제를 어떻게 해결하는지 쓰라.

   <details>
   <summary>정답</summary>

   2차 ODE를 전진 적분하려면 시작점에서 $`T`$와 $`dT/dx`$를 모두 알아야 한다. BVP에서는 $`x = 0`$의 $`T`$와 $`x = L`$의 $`T`$가 주어지므로 시작점의 기울기 $`z(0)`$를 모른다. 사격법은 $`z(0)`$를 추정하여 IVP로 적분하고, 계산된 $`T(L)`$과 주어진 $`T_b`$의 차이가 허용 오차 안에 들 때까지 근 찾기 방법으로 추정을 조정한다.

   </details>

2. **해석해의 상수:** $`h' = 0.01`$, $`L = 10`$, $`T_\infty = 20`$, $`T_a = 40`$, $`T_b = 200`$일 때 $`\lambda`$, $`C_1`$, $`C_2`$를 구하고 $`T(5)`$를 계산하라.

   <details>
   <summary>정답</summary>

   $`\lambda = \sqrt{0.01} = 0.1`$이고 $`\lambda L = 1`$이다. $`e \approx 2.718282`$, $`e^{-1} \approx 0.367879`$를 대입하면 $`C_1 = (0.367879 \times 20 - 180)/(0.367879 - 2.718282) = (-172.6424)/(-2.350402) \approx 73.4523`$, $`C_2 = (2.718282 \times 20 - 180)/(2.718282 - 0.367879) = (-125.6344)/2.350402 \approx -53.4523`$이다. 검산하면 $`C_1 + C_2 = 20 = T_a - T_\infty`$이다. 따라서 $`T(5) = 20 + 73.4523 e^{0.5} - 53.4523 e^{-0.5} \approx 20 + 121.1023 - 32.4204 = 108.6819`$이다.

   </details>

3. **사격법의 선형 보간:** 2번 문제의 조건에서 $`z(0) = 10`$으로 적분하면 $`T(10) \approx 168.3817`$, $`z(0) = 20`$으로 적분하면 $`T(10) \approx 285.9019`$를 얻는다. 다음 추정값 $`z_{a,\text{new}}`$를 구하고, 이 값이 2번 문제의 해석해에서 얻은 $`dT/dx`$의 $`x = 0`$ 값과 일치하는 이유를 설명하라.

   <details>
   <summary>정답</summary>

   $`z_{a,\text{new}} = 10 + (200 - 168.3817)(20 - 10)/(285.9019 - 168.3817) = 10 + 316.183/117.5202 \approx 12.6905`$이다. 해석해에서 $`dT/dx\big|_{x=0} = \lambda(C_1 - C_2) = 0.1 \times (73.4523 + 53.4523) \approx 12.6905`$이므로 두 값이 같다. 이 ODE는 선형이어서 $`T(L) = T_\infty + (T_a - T_\infty)\cosh(\lambda L) + (z(0)/\lambda)\sinh(\lambda L)`$처럼 $`T(L)`$이 $`z(0)`$의 일차 함수이다. 따라서 두 점을 지나는 직선으로 보간한 값이 정확한 초기 기울기가 된다.

   </details>

4. **유한 차분 계산:** 2번 문제의 조건에서 $`n = 2`$($`\Delta x = 5`$)로 유한 차분법을 적용하여 $`T_1 \approx T(5)`$를 구하고, 해석해 및 $`n = 10`$일 때의 결과($`108.6904`$)와 비교하라.

   <details>
   <summary>정답</summary>

   $`G = 2 + 25 \times 0.01 = 2.25`$, $`F = 0.01 \times 20 \times 25 = 5`$이다. 내부 노드가 하나뿐이므로 $`-2.25 T_1 = -5 - 40 - 200 = -245`$이고, $`T_1 = 245/2.25 \approx 108.8889`$이다. 해석해 $`108.6819`$와의 오차는 약 $`0.2070`$이고, $`n = 10`$일 때의 오차는 약 $`0.0085`$이다. $`\Delta x`$가 5에서 1로 1/5이 되자 오차가 약 1/24로 줄었으며, 이는 오차가 $`(\Delta x)^2`$에 비례하는 2차 정확도(1/25)와 거의 맞는다.

   </details>

5. **계수 행렬의 성질:** $`n = 10`$일 때 유한 차분법이 푸는 선형 시스템의 크기를 쓰고, 4.6절의 시스템이 유일한 해를 가지며 가우스-자이델 반복법이 수렴하는 이유를 설명하라.

   <details>
   <summary>정답</summary>

   $`T_0`$과 $`T_{10}`$은 경계 조건으로 알려져 있으므로 미지수는 $`T_1, \ldots, T_9`$이고, $`9 \times 9`$ 시스템을 푼다. 각 행에서 대각 원소의 절댓값은 $`G = 2 + (\Delta x)^2 h'`$로 $`h' > 0`$이면 2보다 크고, 비대각 원소의 절댓값 합은 많아야 2이다. 따라서 행렬은 엄밀한 대각 우위이며, 엄밀한 대각 우위 행렬은 정칙이므로 해가 유일하고 가우스-자이델 반복법이 수렴한다.

   </details>
