# 제6장 강의 — 근: 개방법(Open Methods)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 6

> **선수 지식**: [미적분학] 미분 (제1-5장).
>
> **학습 목표**:
> 1. 근 구하기를 위해 뉴턴-랩슨 방법을 적용할 수 있다
> 2. 시컨트 방법과 수정 뉴턴 방법을 구현할 수 있다
> 3. 개방형 vs 브래킷 방법의 수렴을 비교할 수 있다

---

<br>

## 목차

1. [목표](#1-목표)
2. [개방법 vs 구간법](#2-개방법-vs-구간법)
3. [고정점 반복(Fixed-Point Iteration)](#3-고정점-반복fixed-point-iteration)
   - 3.1 [개념과 정의](#31-개념과-정의)
   - 3.2 [알고리즘](#32-알고리즘)
   - 3.3 [고정점과 근의 관계](#33-고정점과-근의-관계)
   - 3.4 [수렴 분석](#34-수렴-분석)
   - 3.5 [평균값 정리를 통한 수렴 증명](#35-평균값-정리를-통한-수렴-증명)
   - 3.6 [예제 6.1 — 수렴하는 경우](#36-예제-61--수렴하는-경우)
   - 3.7 [예제 6.1 — 발산하는 경우](#37-예제-61--발산하는-경우)
   - 3.8 [예제 6.1 — 순환하는 경우](#38-예제-61--순환하는-경우)
4. [Wegstein 방법](#4-wegstein-방법)
   - 4.1 [가위치법에서의 동기](#41-가위치법에서의-동기)
   - 4.2 [기하학적 아이디어](#42-기하학적-아이디어)
   - 4.3 [유도](#43-유도)
   - 4.4 [일반 공식](#44-일반-공식)
5. [뉴턴-랩슨법(Newton-Raphson Method)](#5-뉴턴-랩슨법newton-raphson-method)
   - 5.1 [기하학적 유도](#51-기하학적-유도)
   - 5.2 [2단계 접근법](#52-2단계-접근법)
   - 5.3 [테일러 전개로부터의 유도](#53-테일러-전개로부터의-유도)
   - 5.4 [이차 수렴(Quadratic Convergence)](#54-이차-수렴quadratic-convergence)
   - 5.5 [이차 수렴 증명](#55-이차-수렴-증명)
   - 5.6 [예제 6.5 — 느리게 수렴하는 함수](#56-예제-65--느리게-수렴하는-함수)
6. [할선법(Secant Method)](#6-할선법secant-method)
   - 6.1 [표준 할선법](#61-표준-할선법)
   - 6.2 [수정 할선법(Modified Secant Method)](#62-수정-할선법modified-secant-method)
7. [역 이차 보간(Inverse Quadratic Interpolation)](#7-역-이차-보간inverse-quadratic-interpolation)
   - 7.1 [개념](#71-개념)
   - 7.2 [라그랑주 다항식 공식화](#72-라그랑주-다항식-공식화)
   - 7.3 [할선법과의 비교](#73-할선법과의-비교)
8. [브렌트 방법(Brent's Method)](#8-브렌트-방법brents-method)
   - 8.1 [개요](#81-개요)
   - 8.2 [SciPy를 이용한 Python 구현](#82-scipy를-이용한-python-구현)
9. [요약 표](#9-요약-표)

<br>

---

<br>

## 1. 목표

- **고정점 반복(Fixed-point iteration)** 을 이해한다
- **Wegstein 방법** 을 구현한다
- 이차 수렴하는 **뉴턴-랩슨법(Newton-Raphson method)** 을 이해한다
- **할선법(Secant method)** 을 이해한다
- **브렌트 방법(Brent's method)** 을 이해한다

> **선수 과목:** 제5장 — 근: 구간 방법(이분법, 가위치법). 제5장에서 다룬 구간, 오차 허용 범위, 반복 수렴의 개념이 여기서 직접적으로 확장된다.

<br>

---

<br>

## 2. 개방법 vs 구간법

**구간법**(이분법, 가위치법)과 달리, 개방법은 근을 둘러싸는 두 개의 초기 추정값을 **필요로 하지 않는다**. 주요 특징:

| 속성 | 구간법(Bracketing Methods) | 개방법(Open Methods) |
|---|---|---|
| 초기 추정값 | $f(x_l) \cdot f(x_u) < 0$인 두 점 $x_l, x_u$ | 하나 또는 두 개의 초기 추정값 (부호 변화 요구 없음) |
| 수렴 보장 | 항상 수렴 (보장됨) | **발산** 하거나 근에서 멀어질 수 있음 |
| 속도 | 느림 (선형 수렴) | 수렴할 때 **더 빠른** 수렴 |
| 구간 동작 | 매 반복마다 구간이 축소 | 단일 점(들)이 근을 향해 반복 |

> **핵심 통찰:** 개방법은 수렴 보장을 속도와 교환한다. 수렴할 경우, 구간법보다 훨씬 빠르게 근에 접근한다.

<br>

---

<br>

## 3. 고정점 반복(Fixed-Point Iteration)

### 3.1 개념과 정의

다른 이름:
- **단순 1점 반복법(Simple one-point iteration)**
- **연속 대입법(Method of successive substitutions)**

$f: [a, b] \to \mathbb{R}$인 방정식 $f(x) = 0$이 주어지면, 이를 다음 형태로 재배열한다:

$$f(x) = x - g(x) = 0 \quad \Longrightarrow \quad x = g(x)$$

반복 $x_{i+1} = g(x_i)$를 **고정점 반복** 이라 한다.

**정의:** 함수 $g$의 **고정점(fixed point)** 은 다음을 만족하는 실수 $x^*$이다:

$$g(x^*) = x^*$$

> **선수 과목 연결:** 구간법(제5장)은 $f(x) = 0$을 직접 다루는 반면, 고정점 반복은 문제를 $x = g(x)$로 재공식화하며, 이는 $y = g(x)$와 $y = x$의 교점을 찾는 것으로 볼 수 있다.

### 3.2 알고리즘

근 $x^* \in (a, b)$를 찾기 위해:

1. 시작점 $x_0 \in (a, b)$를 선택
2. 다음을 사용하여 $x_1, x_2, \ldots, x_{i+1}$을 순차적으로 계산:

$$x_{i+1} = g(x_i)$$

3. $x_{i+1} \approx x_i$가 될 때까지 반복 (즉, 수렴이 달성될 때까지)

이것은 **순환 계산(circular calculation)** 이다: 한 단계의 출력 $x_{i+1}$이 다음 단계의 입력 $x_i$가 된다.

```
input x_i  -->  [ g ]  -->  output x_{i+1}
     ^                            |
     |____________________________|
```

### 3.3 고정점과 근의 관계

$g$를 연속 함수라 하자. 수열 $\{x_i : i \geq 0\}$이 고정점 반복 $x_{i+1} = g(x_i)$로 생성된다고 하자.

$\displaystyle\lim_{i \to \infty} x_i = x^*$이면, $x^*$는 $g$의 고정점이다.

따라서 $x^*$는 $f$의 **근** 이다:

$$f(x^*) = x^* - g(x^*) = 0$$

### 3.4 수렴 분석

수렴은 고정점에서의 **$g$의 기울기** 에 의존한다. 평균값 정리를 사용하면:

$$\text{err}_{i+1} = \text{err}_i \cdot |g'(\xi)|$$

여기서 $\xi$는 $x_i$와 $x^*$ 사이의 어떤 점이다.

**수렴 기준:**

| 조건 | 결과 |
|---|---|
| $\|g'(\xi)\| < 1$ | $x_{i+1}$이 **수렴** |
| $\|g'(\xi)\| > 1$ | $x_{i+1}$이 **발산** |
| $\|g'(\xi)\| = 1$ | 판정 불가 (순환할 수 있음) |

### 3.5 평균값 정리를 통한 수렴 증명

$x^* \in (a, b)$를 $g$의 고정점, $x_i \in (a, b)$라 하자.

**평균값 정리** 에 의해, $\xi_i \in (x_i, x^*)$ 또는 $(x^*, x_i)$인 $\xi_i$가 존재하여:

$$g(x_i) - g(x^*) = g'(\xi_i)(x_i - x^*)$$

$x_i = g(x_{i-1})$이고 $x^* = g(x^*)$이므로:

$$|x_i - x^*| = |g(x_{i-1}) - g(x^*)| = |g'(\xi_{i-1})| \cdot |x_{i-1} - x^*|$$

$(a, b)$의 모든 $\xi$에 대해 $|g'(\xi)| \leq q < 1$이면:

$$|x_i - x^*| \leq q \cdot |x_{i-1} - x^*| \leq q^2 \cdot |x_{i-2} - x^*| \leq \cdots \leq q^i \cdot |x_0 - x^*|$$

$q < 1$이므로, $i \to \infty$일 때 $q^i \to 0$이 되어:

$$|x_i - x^*| \to 0$$

> **핵심 결과:** 근과 반복값을 포함하는 구간에서 **축소 상수(contraction constant)** $q = \max|g'(\xi)| < 1$이면 고정점 반복은 수렴한다.

### 3.6 예제 6.1 — 수렴하는 경우

**문제:** $f(x) = x - e^{-x}$의 근을 추정하라.

**설정:** $f(x) = 0$에서:

$$x = e^{-x} = g(x)$$

$x_0 = 0$으로 시작:

| 반복 | $x_i$ | $g(x_i) = e^{-x_i}$ |
|---|---|---|
| 0 | 0 | $e^{0} = 1.0000$ |
| 1 | 1.0000 | $e^{-1} = 0.3679$ |
| 2 | 0.3679 | $e^{-0.3679} = 0.6922$ |
| 3 | 0.6922 | $\cdots$ |
| $\vdots$ | $\vdots$ | $x^* \approx 0.5671$로 수렴 |

**수렴하는 이유:** $g(x) = e^{-x}$이므로 $g'(x) = -e^{-x}$. 근 $x^* \approx 0.5671$에서 $|g'(x^*)| = e^{-0.5671} \approx 0.5671 < 1$.

그래프적으로, 이것은 **나선형 수렴(spiral convergence)** 이다 -- 반복값이 $y = g(x) = e^{-x}$와 $y = x$의 교점을 향해 안쪽으로 나선을 그린다.

```python
import numpy as np

def fixed_point_ex1(x0, tol=1e-6, max_iter=100):
    """Fixed-point iteration for f(x) = x - e^(-x)"""
    g = lambda x: np.exp(-x)
    x = x0
    for i in range(max_iter):
        x_new = g(x)
        if abs(x_new - x) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.6f}")
            return x_new
        x = x_new
    print("Did not converge")
    return x
```

### 3.7 예제 6.1 — 발산하는 경우

**문제:** $f(x) = x + \ln x$의 근을 추정하라.

**설정:** $f(x) = 0$에서:

$$x = -\ln(x) = g(x)$$

**결과:** 반복이 **발산** 한다 (나선이 바깥으로 향함).

**발산하는 이유:** $g(x) = -\ln(x)$이므로 $g'(x) = -1/x$. 근 $x^* \approx 0.5671$에서 $|g'(x^*)| = 1/0.5671 \approx 1.763 > 1$.

> **교훈:** 같은 근이라도 $f(x) = 0$의 한 재배열로는 찾을 수 있지만 다른 재배열로는 찾을 수 없다. $g(x)$의 선택이 수렴을 결정적으로 좌우한다.

### 3.8 예제 6.1 — 순환하는 경우

$g(x) = \dfrac{1}{x}$에서 $x_0 = 4$인 경우를 고려:

| 반복 | $x_i$ | $g(x_i) = 1/x_i$ |
|---|---|---|
| 0 | 4 | $1/4 = 0.25$ |
| 1 | 0.25 | $1/0.25 = 4$ |
| 2 | 4 | $1/4 = 0.25$ |

반복이 두 점 $x_0$과 $x_1$ 사이를 **무한히 순환** 하며, 고정점 $x^* = 1$에 수렴하지 않는다.

**이유:** $g'(x) = -1/x^2$이므로 $|g'(1)| = 1$ (경계 경우).

<br>

---

<br>

## 4. Wegstein 방법

### 4.1 가위치법에서의 동기

**가위치법(False Position method)**(제5장)에서는 $f(x)$ 위의 두 추정값 사이의 직선이 x축과 만나는 점이 다음 추정값이 되었음을 상기하자.

Wegstein 방법은 유사한 선형 보간 아이디어를 적용하되, $f(x)$ 대신 함수 $g(x)$에 적용한다.

### 4.2 기하학적 아이디어

Wegstein 방법에서는 곡선 $y = g(x)$ 위의 두 추정값 사이의 직선이 **직선 $y = x$와 교차** 하는 점으로 다음 해 추정값을 결정한다.

곡선 $g(x)$ 위의 두 점 $(x_0, g(x_0))$과 $(x_1, g(x_1))$이 주어지면:
- 이 두 점을 지나는 직선을 그린다
- 이 직선이 $y = x$와 만나는 점을 찾는다
- 그 교점이 다음 추정값 $x_2$가 된다

### 4.3 유도

$(x_0, g(x_0))$과 $(x_1, g(x_1))$을 지나는 직선의 기울기는 $(x_1, g(x_1))$과 $(x_2, g(x_2))$를 지나는 직선의 기울기와 같다:

$$\text{slope} = \frac{g(x_1) - g(x_0)}{x_1 - x_0} = \frac{g(x_2) - g(x_1)}{x_2 - x_1}$$

교점이 직선 $y = x$ 위에 있으므로 $g(x_2) = x_2$이다. 따라서:

$$\frac{g(x_1) - g(x_0)}{x_1 - x_0} = \frac{x_2 - g(x_1)}{x_2 - x_1}$$

교차 곱하고 $x_2$에 대해 풀면:

$$x_2(g(x_1) - g(x_0)) - x_1(g(x_1) - g(x_0)) = (x_1 - x_0)x_2 - g(x_1)(x_1 - x_0)$$

$$x_2 = \frac{x_1 g(x_0) - x_0 g(x_1)}{x_1 - x_0 - g(x_1) + g(x_0)}$$

### 4.4 일반 공식

일반 Wegstein 반복 공식은:

$$\boxed{x_{i+1} = \frac{x_i \, g(x_{i-1}) - x_{i-1} \, g(x_i)}{x_i - x_{i-1} - g(x_i) + g(x_{i-1})}}$$

> **핵심 특성:** 마지막 두 추정값을 사용하여 다음 추정값을 결정한다. 이는 **두 개의 초기 추정값** 이 필요하지만(할선법과 유사), 근을 둘러쌀 필요는 없다.

> **선수 과목 연결:** Wegstein 방법은 본질적으로 가위치법(제5장)의 아이디어를 $f(x)$가 $y = 0$과 만나는 것이 아닌 $g(x)$가 $y = x$와 만나는 것에 적용한 것이다.

```python
def wegstein(g, x0, x1, tol=1e-6, max_iter=100):
    """Wegstein method for solving x = g(x)"""
    for i in range(max_iter):
        g0, g1 = g(x0), g(x1)
        denom = x1 - x0 - g1 + g0
        if abs(denom) < 1e-15:
            print("Denominator near zero")
            return x1
        x_new = (x1 * g0 - x0 * g1) / denom
        if abs(x_new - x1) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.6f}")
            return x_new
        x0, x1 = x1, x_new
    print("Did not converge")
    return x1
```

<br>

---

<br>

## 5. 뉴턴-랩슨법(Newton-Raphson Method)

### 5.1 기하학적 유도

$f(x) = 0$이 주어지면, 현재 추정값에서의 **1차 도함수**(접선)를 사용한다.

점 $(x_i, f(x_i))$에서 접선의 기울기는:

$$\text{slope} = f'(x_i) = \frac{0 - f(x_i)}{x_{i+1} - x_i}$$

$x_{i+1}$에 대해 풀면:

$$\boxed{x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}}$$

### 5.2 2단계 접근법

뉴턴-랩슨법은 2단계 과정으로 볼 수 있다:

**단계 1:** 보정값(뉴턴 스텝)을 계산:

$$\delta x_i = -\frac{f(x_i)}{f'(x_i)}$$

**단계 2:** 추정값을 갱신:

$$x_{i+1} = x_i + \delta x_i$$

> **선수 과목 연결:** 뉴턴-랩슨법은 **하나의** 초기 추정값과 도함수 $f'(x)$만 필요한 개방법이다. 구간법과 달리 수렴이 보장되지 않지만, 수렴할 때 이차적으로 수렴한다.

### 5.3 테일러 전개로부터의 유도

뉴턴-랩슨 공식을 **테일러 전개(Taylor expansion)** 로부터 유도할 수 있다:

$$f(x_{i+1}) = f(x_i + \delta x_i) = f(x_i) + f'(x_i)\,\delta x_i + \frac{f''(x_i)}{2}(\delta x_i)^2 + \cdots$$

1차 항 이후를 절단하면:

$$f(x_{i+1}) \approx f(x_i) + f'(x_i)\,\delta x_i$$

$f(x_{i+1}) = 0$이 되길 **기대** 한다 (즉, $x_{i+1}$이 근 $x^*$). 위 식을 영으로 설정하면:

$$0 = f(x_i) + f'(x_i)\,\delta x_i$$

$$\delta x_i = -\frac{f(x_i)}{f'(x_i)}$$

이로써 뉴턴-랩슨 공식이 복원된다.

### 5.4 이차 수렴(Quadratic Convergence)

뉴턴-랩슨법은 **이차 수렴(quadratic convergence)** 을 가진다, 즉:

$$|x_{i+1} - x^*| \lesssim |x_i - x^*|^2$$

**직관적 의미:** 현재 오차가 $\sim 10^{-2}$이면, 다음 오차는 $\sim 10^{-4}$가 된다. 정확한 자릿수가 매 반복마다 대략 **두 배** 가 된다.

**정형화된 서술:** $f$가 $(a, b)$에서 두 번 연속 미분 가능하다고 가정하자. $f'(x^*) \neq 0$인 $x^* \in (a, b)$가 존재한다고 하자. 뉴턴-랩슨 반복을 정의:

$$x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}, \quad i = 1, 2, \ldots$$

$i \to \infty$일 때 $x_i \to x^*$라 가정하면, 충분히 큰 $i$에 대해:

$$|x_{i+1} - x^*| \leq M |x_i - x^*|^2 \quad \text{where} \quad M \geq \frac{|f''(x^*)|}{2|f'(x^*)|}$$

### 5.5 이차 수렴 증명

**증명:**

$e_i = x_i - x^*$ (단계 $i$에서의 오차)라 하자, 따라서 $x^* = x_i - e_i$.

$x_i - e_i$에서 $f$의 테일러 전개를 고려:

$$f(x_i - e_i) = f(x_i) - f'(x_i)\,e_i + f''(\xi_i)\frac{e_i^2}{2}$$

여기서 $\xi_i$는 $x_i$와 $x^*$ 사이에 있다 ($f''$에 대한 평균값 정리 적용).

$x^* = x_i - e_i$이고 $f(x^*) = 0$이므로:

$$0 = f(x_i) - f'(x_i)\,e_i + f''(\xi_i)\frac{e_i^2}{2}$$

$f'(x_i) \neq 0$ ($x_i$가 $x^*$에 충분히 가까운 한)이므로, $f'(x_i)$로 나누면:

$$0 = \frac{f(x_i)}{f'(x_i)} - e_i + \frac{f''(\xi_i)}{f'(x_i)}\frac{e_i^2}{2}$$

뉴턴법의 정의에 의해 $x_{i+1} = x_i - \dfrac{f(x_i)}{f'(x_i)}$이므로:

$$x_{i+1} - x^* = x_i - x^* - \frac{f(x_i)}{f'(x_i)} = e_i - \frac{f(x_i)}{f'(x_i)}$$

위 식으로부터:

$$0 = (x^* - x_{i+1}) + \frac{f''(\xi_i)}{f'(x_i)} \cdot \frac{e_i^2}{2}$$

따라서:

$$x_{i+1} - x^* = \frac{f''(\xi_i)}{f'(x_i)} \cdot \frac{1}{2}(x_i - x^*)^2$$

절댓값을 취하면:

$$|x_{i+1} - x^*| = \frac{|f''(\xi_i)|}{2|f'(x_i)|} \cdot |x_i - x^*|^2$$

$\dfrac{|f''(\xi_i)|}{2|f'(x_i)|} \leq M$인 상수 $M$이 존재하므로:

$$\boxed{|x_{i+1} - x^*| \leq M |x_i - x^*|^2}$$

이로써 **이차 수렴** 이 증명된다. $\blacksquare$

```python
def newton_raphson(f, df, x0, tol=1e-6, max_iter=100):
    """Newton-Raphson method"""
    x = x0
    for i in range(max_iter):
        fx = f(x)
        dfx = df(x)
        if abs(dfx) < 1e-15:
            print("Derivative near zero — method fails")
            return x
        delta_x = -fx / dfx            # Step 1: Newton step
        x = x + delta_x                # Step 2: Update
        if abs(delta_x) < tol:
            print(f"Converged in {i+1} iterations: x* = {x:.10f}")
            return x
    print("Did not converge")
    return x
```

### 5.6 예제 6.5 — 느리게 수렴하는 함수

**문제:** $f(x) = x^{60} - 1$

이것은 근에서 멀리 떨어진 $x$ 값에서 기울기가 거의 영이므로 **느리게 수렴하는 함수** 이다.

- $x_0 = 0.5$인 경우: $x = 0.5$에서의 기울기가 매우 작아 ($f'(0.5) = 60 \cdot 0.5^{59} \approx 0$) 뉴턴 스텝이 매우 커져 발산할 수 있다.
- $x_0 = 0.9$인 경우: $x_0$이 근 $x^* = 1$에 충분히 가까우므로 **3번 반복** 이내에 수렴한다.

> **교훈:** 뉴턴-랩슨법은 다음 경우에 실패하거나 느리게 수렴할 수 있다:
> - $f'(x_i) \approx 0$ (거의 수평인 접선이 과도한 이동을 유발)
> - 초기 추정값이 근에서 먼 경우
> - 근 근처에서 함수가 매우 평탄한 경우

<br>

---

<br>

## 6. 할선법(Secant Method)

### 6.1 표준 할선법

할선법은 뉴턴-랩슨법의 **유한 차분 근사(finite difference approximation)** 이다. **해석적 도함수를 구할 수 없을 때** 유용하다.

뉴턴-랩슨법에서 출발:

$$x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}$$

$f'(x_i)$를 이전 두 점을 사용한 유한 차분 근사로 대체:

$$f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}}$$

대입하면:

$$\boxed{x_{i+1} = x_i - \frac{f(x_i)(x_i - x_{i-1})}{f(x_i) - f(x_{i-1})}}$$

**특성:**
- Wegstein 방법처럼 **두 개의 초기 추정값** $x_0$과 $x_1$이 필요하지만, 근을 둘러쌀 필요는 없다
- 두 추정값을 사용한다는 점에서 Wegstein 방법과 유사
- 수렴 차수는 **초선형(superlinear)** ($\approx 1.618$, 황금비)으로, 선형보다 빠르지만 이차보다 느리다

> **선수 과목 연결:** 할선법 공식은 형태상 가위치법 공식(제5장)과 동일하지만, 할선법은 구간을 **유지하지 않는다**. 항상 가장 최근의 두 점을 사용하는 반면, 가위치법은 구간을 유지한다.

```python
def secant_method(f, x0, x1, tol=1e-6, max_iter=100):
    """Standard secant method"""
    for i in range(max_iter):
        f0, f1 = f(x0), f(x1)
        if abs(f1 - f0) < 1e-15:
            print("Division by zero in secant method")
            return x1
        x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x_new - x1) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.10f}")
            return x_new
        x0, x1 = x1, x_new
    print("Did not converge")
    return x1
```

### 6.2 수정 할선법(Modified Secant Method)

두 점 $x_{i-1}$과 $x_i$가 **충분히 가깝지 않으면**, **수정 할선법** 을 사용할 수 있다.

두 개의 별도 점을 사용하는 대신, $x_i$를 작은 분율 $\delta x_i$만큼 교란한다:

$$f'(x_i) \approx \frac{f(x_i + \delta x_i) - f(x_i)}{\delta x_i}$$

뉴턴-랩슨에 대입하면:

$$\boxed{x_{i+1} = x_i - \frac{f(x_i) \cdot \delta x_i}{f(x_i + \delta x_i) - f(x_i)}}$$

여기서 $\delta x_i$는 작은 교란값이다 (예: 작은 $\epsilon$에 대해 $\delta x_i = \epsilon \cdot x_i$).

> **장점:** 뉴턴-랩슨법처럼 **하나의** 초기 추정값만 필요하며 해석적 도함수 계산을 피한다. 교란 $\delta x_i$가 유한 차분의 스텝 크기를 제어한다.

```python
def modified_secant(f, x0, delta=1e-4, tol=1e-6, max_iter=100):
    """Modified secant method using perturbation"""
    x = x0
    for i in range(max_iter):
        fx = f(x)
        dx = delta * x if abs(x) > 1e-10 else delta
        f_perturbed = f(x + dx)
        if abs(f_perturbed - fx) < 1e-15:
            print("Denominator near zero")
            return x
        x_new = x - fx * dx / (f_perturbed - fx)
        if abs(x_new - x) < tol:
            print(f"Converged in {i+1} iterations: x* = {x_new:.10f}")
            return x_new
        x = x_new
    print("Did not converge")
    return x
```

<br>

---

<br>

## 7. 역 이차 보간(Inverse Quadratic Interpolation)

### 7.1 개념

**세 개의 점** $(x_{i-2}, y_{i-2})$, $(x_{i-1}, y_{i-1})$, $(x_i, y_i)$ ($y_k = f(x_k)$)가 있다고 하자.

핵심 아이디어:
- 세 점을 지나는 순방향 이차식 $y = f(x)$는 **x축과 만나지 않을 수 있다**
- 그러나 이 세 점을 지나는 **역함수** $x = g(y)$는 x축과 **만난다** (즉, $g(0)$을 평가하여 근 추정값을 구할 수 있다)

$y$의 **이차 함수** 를 적합한다:

$$x = ay^2 + by + c$$

미지수 3개 $(a, b, c)$와 방정식 3개(각 점에 대해 하나)로 시스템이 유일하게 결정된다.

### 7.2 라그랑주 다항식 공식화

**라그랑주 다항식(Lagrange polynomials)** 을 사용하면, 역 이차 보간은:

$$g(y) = \frac{(y - y_{i-1})(y - y_i)}{(y_{i-2} - y_{i-1})(y_{i-2} - y_i)} x_{i-2} + \frac{(y - y_{i-2})(y - y_i)}{(y_{i-1} - y_{i-2})(y_{i-1} - y_i)} x_{i-1} + \frac{(y - y_{i-2})(y - y_{i-1})}{(y_i - y_{i-2})(y_i - y_{i-1})} x_i$$

이를 간결하게 쓰면:

$$g(y) = \ell_{i-2}(y) \cdot x_{i-2} + \ell_{i-1}(y) \cdot x_{i-1} + \ell_i(y) \cdot x_i$$

여기서 $\ell_k(y)$는 라그랑주 기저 다항식이다.

다음 근 추정값은 $x_{i+1} = g(0)$이다, 즉 위 공식에서 $y = 0$으로 설정한다.

> **대체 방법:** $y_{i-2} = y_{i-1}$이면(두 $y$ 값이 일치) 라그랑주 다항식이 퇴화한다. 이 경우 두 점을 사용하는 **할선법** 으로 대체한다.

### 7.3 할선법과의 비교

| 특성 | 할선법(Secant Method) | 역 이차 보간(Inverse Quadratic Interpolation) |
|---|---|---|
| 사용하는 점 수 | 2 | 3 |
| 보간 방식 | 선형 ($f$ 위의 2점을 지나는 직선) | 이차 (역함수 위의 3점을 지나는 포물선) |
| 수렴 차수 | $\approx 1.618$ (황금비) | $\approx 1.839$ |
| 견고성 | 더 견고 | $y$ 값이 일치하면 실패할 수 있음 |

<br>

---

<br>

## 8. 브렌트 방법(Brent's Method)

### 8.1 개요

**브렌트 방법**(Richard Brent, 1973)은 다음을 결합한 하이브리드 근 찾기 알고리즘이다:
- **이분법(Bisection method)** (구간법, 수렴 보장)
- **할선법(Secant method)** (빠른 개방법)
- **역 이차 보간(Inverse quadratic interpolation)** (적용 가능할 때 더 빠름)

$$\text{브렌트 방법} = \text{구간법} + \text{개방법}$$

**전략:** 더 빠른 개방법(할선법 또는 IQI)이 현재 구간 내의 결과를 생성할 때 이를 사용하고, 실패하거나 구간 밖의 결과를 생성할 때 이분법으로 대체한다.

> **핵심 장점:** 브렌트 방법은 이분법의 **수렴 보장** 을 상속하면서도 유리한 경우에 개방법의 **속도** 를 달성한다. 많은 수치 라이브러리에서 기본 근 찾기 알고리즘으로 사용된다.

### 8.2 SciPy를 이용한 Python 구현

SciPy는 `scipy.optimize.brentq`를 통해 브렌트 방법을 제공한다:

```python
from scipy import optimize

def f(x):
    return (x**2 - 1)

# Find root in [-2, 0]
root = optimize.brentq(f, -2, 0)
print(root)  # -1.0

# Find root in [0, 2]
root = optimize.brentq(f, 0, 2)
print(root)  # 1.0
```

**함수 시그니처:** `brentq(f, xl, xu)`
- `f` : 풀고자 하는 방정식의 함수 이름
- `xl` : 하한(왼쪽) 구간 값
- `xu` : 상한(오른쪽) 구간 값
- $f(x_l) \cdot f(x_u) < 0$ 필요 (이분법과 같은 부호 변화)

> **실용적 참고:** 실제 응용에서 브렌트 방법은 신뢰성과 속도의 조합 때문에 순수한 이분법, 뉴턴-랩슨법, 또는 할선법보다 거의 항상 선호된다.

<br>

---

<br>

## 9. 요약 표

| 방법 | 초기 추정값 | 도함수 필요 | 수렴 차수 | 수렴 보장 | 핵심 공식 |
|---|---|---|---|---|---|
| **고정점 반복(Fixed-Point Iteration)** | 1 | 아니오 | 선형 ($\sim q^i$) | $\|g'(x^*)\| < 1$일 때만 | $x_{i+1} = g(x_i)$ |
| **Wegstein** | 2 | 아니오 | 초선형(Superlinear) | 아니오 | $x_{i+1} = \frac{x_i g(x_{i-1}) - x_{i-1} g(x_i)}{x_i - x_{i-1} - g(x_i) + g(x_{i-1})}$ |
| **뉴턴-랩슨(Newton-Raphson)** | 1 | 예 ($f'$) | 이차 ($\sim e_i^2$) | 아니오 | $x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}$ |
| **할선법(Secant)** | 2 | 아니오 | 초선형 ($\approx 1.618$) | 아니오 | $x_{i+1} = x_i - \frac{f(x_i)(x_i - x_{i-1})}{f(x_i) - f(x_{i-1})}$ |
| **수정 할선법(Modified Secant)** | 1 | 아니오 | 초선형 | 아니오 | $x_{i+1} = x_i - \frac{f(x_i)\,\delta x_i}{f(x_i + \delta x_i) - f(x_i)}$ |
| **역 이차 보간(Inverse Quadratic Interp.)** | 3 | 아니오 | $\approx 1.839$ | 아니오 | $y = 0$에서의 라그랑주 역다항식 |
| **브렌트 방법(Brent's Method)** | 2 (구간) | 아니오 | 초선형 (하이브리드) | 예 (이분법 대체를 통해) | 적응형: 이분법 + 할선법 + IQI |
