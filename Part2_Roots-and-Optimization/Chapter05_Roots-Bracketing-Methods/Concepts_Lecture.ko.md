# 제5장 강의 — 근: 구간 방법(Bracketing Methods)

> **최종 수정일:** 2026-03-26

---

<br>

## 목차

1. [동기: 이차 공식](#1-동기-이차-공식)
2. [근 문제란?](#2-근-문제란)
3. [그래프 방법](#3-그래프-방법)
   - 3.1 [적금 예제](#31-적금-예제)
   - 3.2 [그래프를 이용한 근 추정](#32-그래프를-이용한-근-추정)
4. [구간법 vs. 개방법](#4-구간법-vs-개방법)
5. [증분 탐색(Incremental Search)](#5-증분-탐색incremental-search)
   - 5.1 [부호 변화 원리](#51-부호-변화-원리)
   - 5.2 [알고리즘과 단점](#52-알고리즘과-단점)
6. [이분법(Bisection Method, 반구간법)](#6-이분법bisection-method-반구간법)
   - 6.1 [알고리즘](#61-알고리즘)
   - 6.2 [단계별 절차](#62-단계별-절차)
7. [이분법의 수렴](#7-이분법의-수렴)
   - 7.1 [오차 한계 유도](#71-오차-한계-유도)
   - 7.2 [사전 오차 추정(A Priori Error Estimation)](#72-사전-오차-추정a-priori-error-estimation)
   - 7.3 [필요한 반복 횟수](#73-필요한-반복-횟수)
8. [가위치법(False Position Method, Regula Falsi)](#8-가위치법false-position-method-regula-falsi)
   - 8.1 [공식 유도](#81-공식-유도)
   - 8.2 [수렴 특성](#82-수렴-특성)
   - 8.3 [편측성 문제(One-Sidedness Problem)](#83-편측성-문제one-sidedness-problem)
9. [수정 가위치법(Modified False Position Method)](#9-수정-가위치법modified-false-position-method)
   - 9.1 [스케일링 매개변수 도입](#91-스케일링-매개변수-도입)
   - 9.2 [연속 두 점이 같은 쪽에 있는 경우](#92-연속-두-점이-같은-쪽에-있는-경우)
10. [요약 표](#10-요약-표)

---

<br>

## 1. 동기: 이차 공식

이차 방정식에 대해:

$$f(x) = ax^2 + bx + c = 0$$

완전제곱식을 통해 근을 해석적으로 유도할 수 있다:

> **[미적분]** 완전제곱식(Completing the square)은 이차식을 꼭짓점 형태로 변환하는 기본적인 대수 기법으로, 이차 공식으로 가는 직접적인 경로를 제공한다.

**단계 1.** $a$를 인수로 분리:

$$a\left(x^2 + \frac{b}{a}x\right) + c = 0$$

**단계 2.** 완전제곱식 완성:

$$a\left(x + \frac{b}{2a}\right)^2 - \frac{b^2}{4a} + c = 0$$

**단계 3.** $x$에 대해 풀기 ($a \neq 0$이라 가정):

$$\left(x + \frac{b}{2a}\right)^2 = \frac{b^2 - 4ac}{4a^2}$$

$$x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$$

근은 $f(x) = 0$인 $x$ 값이다. 기하학적으로 근은 곡선이 **x축과 만나는** 점이다.

그러나 고차 다항식과 초월 방정식의 경우, 일반적으로 **닫힌 형태의 공식이 존재하지 않는다**. 근을 찾기 위해 **수치적 방법**이 필요하다:

- **제5장**: 구간법(Bracketing methods) (수렴 보장, 느림)
- **제6장**: 개방법(Open methods) (발산 가능, 빠름)

---

<br>

## 2. 근 문제란?

"근" 문제는 어떤 함수 $f$가 하나 이상의 종속 변수 $x$로 표현되고, $f(x) = 0$의 해가 문제의 해를 제공하는 경우에 발생한다.

이러한 문제는 **설계 문제**에서 필요한 매개변수에 대한 **음함수 방정식**이 나타날 때 자주 발생한다.

> **[미적분]** 중간값 정리(Intermediate Value Theorem, IVT)는 모든 구간법의 이론적 기초이다: $f$가 $[a, b]$에서 연속이고 $f(a) \cdot f(b) < 0$이면, $f(x^*) = 0$인 $x^* \in (a, b)$가 적어도 하나 존재한다.

---

<br>

## 3. 그래프 방법

### 3.1 적금 예제

**문제 설정:** Shin은 매달 $p$를 정기적으로 입금하여 저축한다. 연이율은 $r$이다. $n$번 납입 후 총 금액 $A$는:

$$A = p\left(1 + \frac{r}{12}\right) + p\left(1 + \frac{r}{12}\right)^2 + \cdots + p\left(1 + \frac{r}{12}\right)^n$$

> **[미적분]** 이것은 첫째항이 $p(1 + r/12)$이고 공비가 $(1 + r/12)$인 등비급수이다.

등비급수 공식을 적용하면:

$$A = \sum_{i=1}^{n} p\left(1 + \frac{r}{12}\right)^i = p \cdot \frac{\left(1 + \frac{r}{12}\right)^n - 1}{\left(1 + \frac{r}{12}\right) - 1} = \frac{12p}{r}\left[\left(1 + \frac{r}{12}\right)^n - 1\right]$$

**근 찾기 공식화:** $A$가 목표 금액이라 하자. Shin은 $n$개월 이내에 $A$를 모으려 한다. 이 목표를 달성하기 위한 이율 $r$은?

$$\text{Find } r: \quad f(r) = \frac{12p}{r}\left[\left(1 + \frac{r}{12}\right)^n - 1\right] - A = 0$$

이 방정식은 $r$에 대해 **해석적으로 풀 수 없다** -- 이것은 근 찾기 문제이다.

### 3.2 그래프를 이용한 근 추정

그래프 접근법을 사용하여 $g(x) = f(r)$을 그리고, 곡선이 x축을 지나는 점($g(x) = 0$)을 시각적으로 확인함으로써 근을 추정할 수 있다.

```python
# ex
def func(interest, duration, deposit, target):
    ratio = 1. + interest / 12.
    return 12.0 * deposit / interest * (ratio**duration - 1.0) - target

from functools import partial

duration = 24
deposit = 60.
target = 1500.
g = partial(func, duration=duration, deposit=deposit, target=target)
```

점점 좁은 구간에서 그래프를 그려 보면:

| 확대 단계 | 이율 범위 | 관찰 결과 |
|:---:|:---:|:---|
| 1차 | 1% -- 10% | 근이 약 ~4% 부근 |
| 2차 | 4.0% -- 5.0% | 근이 약 ~4.2% 부근 |
| 3차 | 4.20% -- 4.30% | 근이 약 **4.23%** |

> **[수치해석]** 그래프 방법은 대략적인 초기 추정치를 제공하지만 정밀하지 않다. 일반적으로 더 체계적인 수치 방법의 출발점으로 사용된다.

---

<br>

## 4. 구간법 vs. 개방법

근 찾기 방법은 두 가지 계열로 분류된다:

### 구간법(Bracketing Methods)

- 근을 **둘러싸는**(bracket) **두 개의 초기 추정값**에 기반
- 근은 구간 $[x_l, x_u]$ 내에 위치
- 해로의 **수렴이 보장**됨
- 수렴 속도가 느림

### 개방법(Open Methods)

- **하나 이상의 초기 추정값**을 사용할 수 있지만, **구간이 없음**
- **발산 가능** (수렴이 보장되지 않음)
- 수렴할 경우 수렴 속도가 **빠름**

> **[수치해석]** 수렴 보장과 속도 사이의 상충 관계는 반복되는 주제이다. 실무에서는 구간법으로 시작하여 근에 가까워진 후 개방법으로 전환하는 하이브리드 전략이 자주 사용된다.

---

<br>

## 5. 증분 탐색(Incremental Search)

### 5.1 부호 변화 원리

$f$가 $x_l$에서 $x_u$까지의 구간에서 **실수이고 연속**이며, 근의 양쪽에서 $f$의 부호가 바뀌면:

$$f(x_l) \cdot f(x_u) < 0$$

이는 $x_l$(하한)과 $x_u$(상한) 사이 어딘가에 근이 존재함을 의미한다.

> **[미적분]** 이것은 중간값 정리의 직접적인 적용이다. $[x_l, x_u]$에서 $f$의 연속성은 필요 조건이다.

### 5.2 알고리즘과 단점

**알고리즘:**

1. 구간 $[x_l, x_u]$를 더 작은 부분 구간으로 분할
2. 각 부분 구간을 순회하며 근이 존재하는지 확인 (부호 변화 검사: $f(x_i) \cdot f(x_{i+1}) < 0$)
3. 부호 변화가 감지되면, 해당 부분 구간에 근이 존재

**단점:** **증분 길이**(스텝 크기)에 **민감하다**.

- 스텝 크기가 **너무 크면**, 두 연속 점 사이에 짝수 개의 근이 있을 수 있으므로 근을 포함하는 구간이 **누락**될 수 있다
- 스텝 크기가 **너무 작으면**, 방법이 **계산 비용이 많이** 든다

> **[수치해석]** $\sin(x^2)$과 같은 함수는 매우 촘촘하게 배치된 근을 가질 수 있다. 증분이 인접한 근 사이의 간격보다 크면, 함수가 같은 부호로 돌아오기 때문에 근의 쌍이 완전히 누락된다.

---

<br>

## 6. 이분법(Bisection Method, 반구간법)

이분법은 가장 간단한 구간법이다. 매 반복마다 탐색 구간을 체계적으로 반으로 나눈다.

### 6.1 알고리즘

**핵심 아이디어:**

1. 매 반복마다 구간의 **중점**을 선택하여 탐색 공간을 반으로 분할
2. 두 절반 중 어느 쪽에 근이 있는지 판별
3. 탐색 구간이 허용 오차 $\varepsilon$보다 작아질 때까지 단계 1과 2를 반복

### 6.2 단계별 절차

**단계 0 (초기화):** $f(x_l) \cdot f(x_u) < 0$인 구간 $[x_l, x_u]$로 시작.

$a_0 = x_l$, $b_0 = x_u$로 설정.

중점을 계산:

$$c_0 = \frac{x_l + x_u}{2} = \frac{a_0 + b_0}{2}$$

**단계 1 (어느 절반인지 판별):** 부호 변화 확인:

- $f(a_0) \cdot f(c_0) < 0$이면: 근은 **왼쪽** 구간 $[a_0, c_0]$에 있음
  - $a_1 = a_0$, $b_1 = c_0$으로 설정
- $f(c_0) \cdot f(b_0) < 0$이면: 근은 **오른쪽** 구간 $[c_0, b_0]$에 있음
  - $a_1 = c_0$, $b_1 = b_0$으로 설정

**단계 2 (수렴 확인):** $(b_1 - a_1) < \varepsilon$이면, 출력은:

$$x^* \approx \frac{a_1 + b_1}{2}$$

그렇지 않으면, 구간 $[a_1, b_1]$로 탐색을 계속한다.

**일반 반복:** 반복 $n$에서:

$$c_n = \frac{a_n + b_n}{2}$$

이후 부호 변화 검사에 따라 구간을 갱신하고 반복한다.

```python
def bisection(f, a, b, tol=1e-6, max_iter=100):
    """
    Bisection method for finding a root of f in [a, b].
    Assumes f(a) * f(b) < 0.
    """
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs")

    for n in range(max_iter):
        c = (a + b) / 2.0
        if (b - a) < tol:
            return c, n
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2.0, max_iter
```

---

<br>

## 7. 이분법의 수렴

### 7.1 오차 한계 유도

매 반복마다 오차(참 근 $x^*$와 중점 추정값 $c_n$ 사이의 거리)는 현재 구간 길이의 절반으로 한정된다.

**1번째 반복 후:**

$$\text{err}_1 = |x^* - c_0| \leq \frac{1}{2}(b_0 - a_0)$$

**2번째 반복 후:**

$$\text{err}_2 = |x^* - c_1| \leq \frac{1}{2}(b_1 - a_1)$$

**3번째 반복 후:**

$$\text{err}_3 = |x^* - c_2| \leq \frac{1}{2}(b_2 - a_2)$$

**4번째 반복 후:**

$$\text{err}_4 = |x^* - c_3| \leq \frac{1}{2}(b_3 - a_3)$$

각 구간은 이전 구간의 절반이므로:

$$\frac{1}{2}(b_3 - a_3) = \frac{1}{2^2}(b_2 - a_2) = \frac{1}{2^3}(b_1 - a_1) = \frac{1}{2^4}(b_0 - a_0)$$

### 7.2 사전 오차 추정(A Priori Error Estimation)

$n$번 반복 후, 오차의 상한은:

$$\boxed{\text{err}_n = |x^* - c_{n-1}| \leq \frac{\Delta x}{2^n}}$$

여기서 $\Delta x = b_0 - a_0$는 **초기 구간 폭**이다.

> **[수치해석]** 이것은 **사전 오차 추정(a priori error estimation)**이라 불리는데, 알고리즘을 실행하기 *전에* 초기 구간과 반복 횟수만으로 최대 오차를 예측할 수 있기 때문이다.

**실용적 예시:**
- **10번 반복 후**: 원래 구간 불확실성이 $2^{10} = 1024$배 감소 (약 3자리 정밀도)
- **20번 반복 후**: $2^{20} \approx 10^6$배 감소 (약 6자리 정밀도 -- "백만 배" 개선)

> **[미적분]** 수렴은 **선형(linear)**이다 -- 각 반복이 대략 1비트의 정확도를 추가한다. 이는 이차적으로 수렴하는 고차 방법(예: 뉴턴법)과 대조적이다.

### 7.3 필요한 반복 횟수

**질문:** 주어진 오차 허용 범위 $\varepsilon_s$에 대해 몇 번의 반복이 필요한가?

다음이 필요하다:

$$\text{err}_n \leq \frac{\Delta x}{2^n} < \varepsilon_s$$

정리하면:

$$2^n > \frac{\Delta x}{\varepsilon_s}$$

양변에 $\log_2$를 취하면:

$$\boxed{n > \log_2\left(\frac{\Delta x}{\varepsilon_s}\right)}$$

**예시:** $\Delta x = 1$이고 $\varepsilon_s = 10^{-5}$이면:

$$n > \log_2\left(\frac{1}{10^{-5}}\right) = \log_2(10^5) \approx 16.6$$

따라서 **17번 반복**으로 허용 오차를 보장할 수 있다.

```python
import math

def bisection_iterations_needed(delta_x, tolerance):
    """Calculate minimum number of bisection iterations needed."""
    return math.ceil(math.log2(delta_x / tolerance))

# Example
n = bisection_iterations_needed(1.0, 1e-5)
print(f"Iterations needed: {n}")  # Output: 17
```

---

<br>

## 8. 가위치법(False Position Method, Regula Falsi)

### 8.1 공식 유도

가위치법(FP)은 이분법과 유사하지만, **중점** 대신 **선형 보간(linear interpolation)**을 사용하여 새로운 근 추정값을 제안한다.

$(x_l, f(x_l))$과 $(x_u, f(x_u))$ 사이의 직선이 x축과 만나는 점 $c$를 취한다.

> **[선형대수]** 이 방법은 두 점을 지나는 직선(할선)의 방정식을 구성하고 그 x절편을 찾는다. 이것은 다항식 보간의 가장 간단한 형태이다.

**유도:** $(a, f(a))$와 $(b, f(b))$를 지나는 직선은:

$$y - f(b) = \frac{f(b) - f(a)}{b - a}(x - b)$$

x절편을 구하기 위해 $y = 0$으로 설정:

$$-f(b) \cdot \frac{b - a}{f(b) - f(a)} = x - b$$

$$x = b + \frac{-b \cdot f(b) + a \cdot f(b)}{f(b) - f(a)}$$

$$x = \frac{a \cdot f(b) - b \cdot f(a)}{f(b) - f(a)}$$

따라서 반복 $n$에서의 **가위치법 공식**은:

$$\boxed{c_n = \frac{a_n \cdot f(b_n) - b_n \cdot f(a_n)}{f(b_n) - f(a_n)}}$$

알고리즘은 이분법처럼 진행된다: 어느 부분 구간 $[a_n, c_n]$ 또는 $[c_n, b_n]$에 근이 있는지 확인하고, 구간을 갱신하며, $|f(c_n)| < \varepsilon$이 될 때까지 반복한다.

```python
def false_position(f, a, b, tol=1e-6, max_iter=100):
    """
    False Position (Regula Falsi) method.
    Assumes f(a) * f(b) < 0.
    """
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs")

    for n in range(max_iter):
        fa, fb = f(a), f(b)
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)

        if abs(fc) < tol:
            return c, n

        if fa * fc < 0:
            b = c
        else:
            a = c
    return c, max_iter
```

### 8.2 수렴 특성

- 가위치법은 이분법보다 빠르게 수렴할 **수도 있고 그렇지 않을 수도** 있다
- 주어진 반복 횟수 내에 기준을 충족한다는 보장을 **제공하지 않는다**
- 반복 횟수를 사용하여 **오차를 예측할 수 없다** (이분법의 사전 한계와 달리)

> **[수치해석]** 이분법과 달리, 가위치법에는 $\text{err}_n \leq \Delta x / 2^n$과 같은 간단한 공식이 없다. 수렴 속도는 근 근처에서 $f$의 형태(곡률)에 크게 의존한다.

### 8.3 편측성 문제(One-Sidedness Problem)

가위치법의 주요 약점은 **편측성 문제**이다: 구간 끝점 중 하나가 **고정**(이동하지 않음)된 채로 남는다. 이로 인해 가위치법은 **느린 수렴**을 보인다.

**예시:** $[-1, 1]$에서 $f(x) = 2x^3 - 4x^2 + 3x$.

함수의 곡률이 큰 경우, 한쪽 끝점은 "고정"되고 다른 쪽은 한 방향에서 천천히 근에 접근한다. 구간이 대칭적으로 줄어들지 않는다.

---

<br>

## 9. 수정 가위치법(Modified False Position Method)

### 9.1 스케일링 매개변수 도입

근을 포함하는 구간의 길이가 영으로 수렴하도록 보장하기 위해, Regula-Falsi 방법을 스케일링 매개변수 $\alpha$와 $\beta$ ($\alpha \cdot \beta = 1$)를 도입하여 **수정**할 수 있다:

**원래 공식:**

$$c_n = \frac{a_n \cdot f(b_n) - b_n \cdot f(a_n)}{f(b_n) - f(a_n)}$$

**수정 공식:**

$$c_n = \frac{a_n \cdot \beta \cdot f(b_n) - b_n \cdot \alpha \cdot f(a_n)}{\beta \cdot f(b_n) - \alpha \cdot f(a_n)}$$

### 9.2 연속 두 점이 같은 쪽에 있는 경우

수정은 **연속 두 중점 추정값** $c_0$과 $c_1$이 근의 **같은 쪽**에 나타날 때 ($f(c_0) \cdot f(c_1) > 0$) 트리거된다.

**경우 1: 연속 두 점이 근의 오른쪽에 있는 경우** (즉, $f(c_0) \cdot f(c_1) > 0$이고 고정 끝점이 $a_n$인 경우):

$\alpha = \frac{1}{2}$로 설정하여 **고정** 끝점의 함수값을 절반으로 줄인다:

$$c = \frac{a \cdot f(b) - b \cdot \frac{1}{2} f(a)}{f(b) - \frac{1}{2} f(a)}$$

$f(a)$를 절반으로 줄이면, 할선이 기울어져 다음 추정값이 고정된 쪽에 더 가까워지고, 따라서 실제 근에 더 가까워진다.

**경우 2: 연속 두 점이 근의 왼쪽에 있는 경우** (고정 끝점이 $b_n$인 경우):

$\beta = \frac{1}{2}$로 설정하여 고정 끝점의 함수값을 절반으로 줄인다:

$$c = \frac{a \cdot \frac{1}{2} f(b) - b \cdot f(a)}{\frac{1}{2} f(b) - f(a)}$$

> **[수치해석]** 수정 가위치법은 구간 폭이 영으로 줄어듦을 보장하여, 표준 Regula Falsi의 편측성 문제를 극복한다. Illinois 알고리즘은 이 절반 전략을 사용하는 잘 알려진 변형이다.

```python
def modified_false_position(f, a, b, tol=1e-6, max_iter=100):
    """
    Modified False Position method (Illinois algorithm variant).
    Halves the function value at the stationary endpoint when
    two consecutive estimates appear on the same side.
    """
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs")

    fa, fb = f(a), f(b)
    side = 0  # Track which side was updated last

    for n in range(max_iter):
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)

        if abs(fc) < tol:
            return c, n

        if fa * fc < 0:
            # Root is in [a, c]
            if side == -1:
                # Two consecutive on the right side: halve f(a)
                fa *= 0.5
            b = c
            fb = fc
            side = -1
        else:
            # Root is in [c, b]
            if side == 1:
                # Two consecutive on the left side: halve f(b)
                fb *= 0.5
            a = c
            fa = fc
            side = 1

    return c, max_iter
```

---

<br>

## 10. 요약 표

| 속성 | 증분 탐색(Incremental Search) | 이분법(Bisection) | 가위치법(False Position) | 수정 가위치법(Modified False Position) |
|:---|:---:|:---:|:---:|:---:|
| **유형** | 구간법 | 구간법 | 구간법 | 구간법 |
| **초기 요구 사항** | 구간 $[x_l, x_u]$ | $f(x_l) \cdot f(x_u) < 0$ | $f(x_l) \cdot f(x_u) < 0$ | $f(x_l) \cdot f(x_u) < 0$ |
| **중점 공식** | 균등 스텝 | $c = \frac{a + b}{2}$ | $c = \frac{a f(b) - b f(a)}{f(b) - f(a)}$ | FP의 스케일 버전 |
| **수렴 보장** | 스텝 크기에 따라 다름 | 항상 (연속이면) | 항상 (연속이면) | 항상 (연속이면) |
| **사전 오차 한계** | 없음 | $\frac{\Delta x}{2^n}$ | 없음 | 없음 |
| **허용 오차에 필요한 반복 횟수** | 알 수 없음 | $n > \log_2(\Delta x / \varepsilon_s)$ | 예측 불가 | 예측 불가 |
| **편측성 문제** | 해당 없음 | 없음 | 있음 | 없음 (수정됨) |
| **수렴 속도** | 해당 없음 | 선형 | 다양함 (정체 가능) | 일반적으로 FP보다 우수 |
| **주요 장점** | 구현이 간단 | 예측 가능하고 견고함 | 이분법보다 자주 빠름 | FP 정체 문제 해결 |
| **주요 단점** | 스텝 크기에 민감 | 느림 (선형 수렴) | 정체 가능 (편측성) | 약간 더 복잡 |
