# 제3장 실습 — 수치 오차의 원인

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 3

> **선수 지식**: [프로그래밍언어] MATLAB/Python. [미적분학] 테일러 급수. [프로그래밍언어] 기본 프로그래밍 (제1-2장).
>
> **학습 목표**:
> 1. 유효 숫자와 정확도를 구분할 수 있다
> 2. 부동소수점 연산에서 반올림 오차를 분석할 수 있다
> 3. IEEE 754 부동소수점 표현을 적용할 수 있다

---

<br>

## 목차

- [1. 오차 정의와 지표](#1-오차-정의와-지표)
  - [1.1 스카보로 기준 표](#11-스카보로-기준-표)
  - [1.2 sqrt(2) 근사 — 바빌로니아/헤론 방법](#12-sqrt2-근사--바빌로니아헤론-방법)
- [2. 절삭 오차와 테일러 급수](#2-절삭-오차와-테일러-급수)
  - [2.1 e^x의 매클로린 급수 — 항별 계산](#21-ex의-매클로린-급수--항별-계산)
  - [2.2 중지 기준을 사용한 반복 계산](#22-중지-기준을-사용한-반복-계산)
  - [2.3 절삭 오차 수렴의 시각화](#23-절삭-오차-수렴의-시각화)
- [3. 반올림 오차와 부동 소수점](#3-반올림-오차와-부동-소수점)
  - [3.1 머신 엡실론](#31-머신-엡실론)
  - [3.2 소수의 이진 표현](#32-소수의-이진-표현)
  - [3.3 반올림 오차의 누적](#33-반올림-오차의-누적)
  - [3.4 카한 합산 알고리즘](#34-카한-합산-알고리즘)
  - [3.5 재앙적 소거 — 이차 방정식 공식](#35-재앙적-소거--이차-방정식-공식)
- [요약](#요약)

---

<br>

## 1. 오차 정의와 지표

### 1.1 스카보로 기준 표

**스카보로 기준(Scarborough Criterion)**은 원하는 유효 숫자 수로부터 필요한 근사 백분율 상대 오차로의 직접적인 매핑을 제공한다. 이 기준은 $|\varepsilon_a| < (0.5 \times 10^{2-n})\%$이면 결과가 최소 $n$개의 유효 숫자까지 정확하다고 말한다.

간단한 람다 함수를 사용하여 각 $n$ 값에 대한 임계값을 계산할 수 있다:

```python
import numpy as np
import matplotlib.pyplot as plt
import math

# Scarborough criterion: minimum epsilon_a for n significant figures
scarborough = lambda n: 0.5 * 10**(2 - n)

# Table of criteria
for n in range(1, 8):
    print(f'{n} sig. figs → |ε_a| < {scarborough(n):.4e}%')
# 1 → 5.0000e+00%, 2 → 5.0000e-01%, ..., 7 → 5.0000e-05%
```

이 표는 실용적인 참고 자료이다: 반복 계산을 시작하기 전에, 필요한 유효 숫자 수를 결정하고, 해당 임계값을 찾아 중지 기준 $\varepsilon_s$로 사용한다.

### 1.2 sqrt(2) 근사 — 바빌로니아/헤론 방법

**바빌로니아 방법(Babylonian Method)**(**헤론 방법(Heron's Method)**이라고도 함)은 고대 메소포타미아로 거슬러 올라가는 가장 오래된 반복 알고리즘 중 하나이다. 다음 점화식을 사용하여 $\sqrt{S}$를 계산한다:

$$x_{n+1} = \frac{1}{2}\left(x_n + \frac{S}{x_n}\right)$$

직관은 우아하다: $x_n$이 $\sqrt{S}$의 과대추정이면 $S/x_n$은 과소추정이다(반대도 성립). 둘의 평균을 취하면 두 값보다 $\sqrt{S}$에 더 가까운 더 나은 추정값을 얻는다.

> **[선형대수]** 바빌로니아 방법은 실제로 $f(x) = x^2 - S$에 뉴턴 방법을 적용한 특수한 경우이다. 뉴턴 방법: $x_{n+1} = x_n - f(x_n)/f'(x_n) = x_n - (x_n^2 - S)/(2x_n) = (x_n + S/x_n)/2$. 이것이 이차 수렴(오차가 각 반복마다 제곱됨)을 설명한다.

다음 코드는 참 백분율 상대 오차 $\varepsilon_t$와 근사 백분율 상대 오차 $\varepsilon_a$를 모두 추적하면서 방법이 얼마나 빠르게 수렴하는지 보여준다:

```python
target = 2
x = 1.0  # initial guess
true_val = math.sqrt(target)

print(f'{"Iter":>4} {"x":>20} {"ε_t (%)":>12} {"ε_a (%)":>12}')
print('-' * 52)

for i in range(6):
    x_old = x
    x = 0.5 * (x + target / x)
    et = abs(true_val - x) / true_val * 100
    ea = abs(x - x_old) / abs(x) * 100 if i > 0 else float('inf')
    print(f'{i+1:4d} {x:20.15f} {et:12.6e} {ea:12.6e}')
```

**이차 수렴(Quadratic Convergence)**을 관찰하라: 정확한 자릿수가 각 반복마다 대략 두 배로 늘어난다. 5~6번째 반복에서 참 오차는 머신 엡실론에 도달하거나 그 근처에 있다 — 알고리즘이 $\sqrt{2}$의 최선의 배정밀도 표현으로 수렴한 것이다.

이는 핵심 원리를 보여준다: 잘 설계된 반복 방법의 경우, $\varepsilon_a$가 $\varepsilon_t$를 밀접하게 추적하므로, 참값을 모를 때에도 $\varepsilon_a$를 신뢰할 수 있는 대리 지표로 사용할 수 있다.

---

<br>

## 2. 절삭 오차와 테일러 급수

### 2.1 e^x의 매클로린 급수 — 항별 계산

$e^x$에 대한 **매클로린 급수(Maclaurin Series)**(테일러 급수를 $x = 0$ 근방에서 전개)는 다음과 같다:

$$e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots = \sum_{k=0}^{\infty} \frac{x^k}{k!}$$

이것은 모든 실수 $x$에 대해 수렴하는 급수이며, 충분한 항을 포함하면 원하는 정확도까지 $e^x$를 근사할 수 있다. 각 추가 항은 절삭 오차를 줄인다.

다음 코드는 $e^{0.5}$에 대해 항별로 근사를 구축하고, 각 추가 항에 따라 오차가 어떻게 감소하는지 보여준다:

```python
# Term-by-term approximation of e^0.5
x = 0.5
true_val = math.exp(x)

f_approx = np.zeros(6)
f_approx[0] = 1.0                           # 0th order: 1
f_approx[1] = f_approx[0] + x              # 1st: + x
f_approx[2] = f_approx[1] + x**2/2         # 2nd: + x²/2!
f_approx[3] = f_approx[2] + x**3/6         # 3rd: + x³/3!
f_approx[4] = f_approx[3] + x**4/24        # 4th: + x⁴/4!
f_approx[5] = f_approx[4] + x**5/120       # 5th: + x⁵/5!

for k in range(6):
    et = abs(true_val - f_approx[k]) / true_val * 100
    print(f'Order {k}: {f_approx[k]:.10f}, ε_t = {et:.6e}%')
```

$e^{0.5} \approx 1.6487212707$의 참값이다. 더 많은 항을 포함할수록 근사값이 빠르게 수렴한다:

- 0차 ($1$만): 오차 $\approx 39.3\%$
- 1차 ($1 + x = 1.5$): 오차 $\approx 9.02\%$
- 5차: 오차는 무시할 수 있을 정도로 작음

분모의 팩토리얼이 분자의 $x$의 거듭제곱보다 훨씬 빠르게 증가하므로, 적당한 $x$ 값에서 급수가 매우 빠르게 수렴한다.

### 2.2 중지 기준을 사용한 반복 계산

실제로는 각 항을 하드 코딩하지 않는다. 대신, 근사 오차가 임계값 아래로 떨어질 때까지 항을 반복적으로 추가하는 루프를 작성한다:

```python
def exp_series(x, es=1e-4, maxit=50):
    """Compute e^x via Taylor series with stopping criterion.

    Args:
        x: input value
        es: stopping criterion (% approximate error)
        maxit: maximum iterations
    Returns:
        (result, final_error_%, iterations)
    """
    sol = 1.0
    ea = 100.0

    for k in range(1, maxit + 1):
        sol_old = sol
        sol = sol + x**k / math.factorial(k)
        ea = abs(sol - sol_old) / abs(sol) * 100
        if ea < es:
            break

    return sol, ea, k

result, error, iters = exp_series(0.5)
print(f'e^0.5 ≈ {result:.15f} (after {iters} iterations, ε_a = {error:.2e}%)')
print(f'True:    {math.exp(0.5):.15f}')
```

이 구현의 핵심 관찰 사항:

- 중지 기준 `es=1e-4`는 $\varepsilon_s = 0.0001\%$에 해당한다. 스카보로 기준에 의하면, 이는 최소 6개의 유효 숫자를 보장한다 ($0.5 \times 10^{2-6} = 0.00005\%$이고, $0.0001\% < 0.0005\%$이므로 5 유효 숫자 제공).
- `maxit` 매개변수는 안전 장치를 제공한다: 급수가 수렴하지 않으면 (예: 매우 큰 $|x|$의 경우), 최대 반복 횟수 후 루프가 종료된다.
- 함수는 호출자에게 필요한 세 가지 정보를 모두 반환한다: 결과, 최종 오차 추정값, 사용된 반복 횟수.

### 2.3 절삭 오차 수렴의 시각화

반로그 그래프는 더 많은 항을 추가할수록 절삭 오차가 기하급수적으로 감소하는 속도를 보여준다:

```python
orders = range(1, 15)
errors = [abs(math.exp(0.5) - sum(0.5**k/math.factorial(k) for k in range(n+1))) for n in orders]

plt.semilogy(list(orders), errors, 'bo-')
plt.xlabel('Number of Terms')
plt.ylabel('|Truncation Error|')
plt.title('Taylor Series Convergence for e^0.5')
plt.grid(True)
plt.show()
```

**반로그 그래프에서의 직선**은 지수적 수렴을 확인한다: 각 추가 항이 대략 일정한 비율로 오차를 줄인다. 이것은 잘 정의된(해석적) 함수에 대한 테일러 급수의 특징이다.

$e^{0.5}$의 경우, 약 10-12개 항 이후 오차가 머신 엡실론($\approx 2.22 \times 10^{-16}$) 아래로 떨어진다. 이 지점 이후로 더 많은 항을 추가해도 더 이상 개선되지 않는데, 반올림 오차가 지배적이기 때문이다 — 이것이 절삭 오차와 반올림 오차가 만나는 지점이다.

---

<br>

## 3. 반올림 오차와 부동 소수점

### 3.1 머신 엡실론

**머신 엡실론(Machine Epsilon)**은 컴퓨터 산술에서 $1 + \varepsilon \neq 1$을 만족하는 가장 작은 부동 소수점 수 $\varepsilon$이다. 이는 부동 소수점 형식의 근본적인 정밀도 한계를 나타낸다.

```python
# Machine epsilon
eps = np.finfo(float).eps
print(f'Machine epsilon: {eps}')  # 2.220446049250313e-16

# Demonstration: 1 + eps/2 == 1 is True!
print(1 + eps/2 == 1)   # True (eps/2 is below representable precision)
print(1 + eps == 1)      # False (eps is exactly at the boundary)
```

첫 번째 테스트(`1 + eps/2 == 1`)는 `True`를 반환하는데, `eps/2`가 64비트 부동 소수점에서 1의 표현에 영향을 줄 만큼 충분히 작기 때문이다. 수 $1 + \varepsilon/2$는 $1$과 구별될 수 없다 — $1$과 $1 + \varepsilon$ 사이에 64비트 부동 소수점 수가 없다.

두 번째 테스트(`1 + eps == 1`)는 `False`를 반환하는데, 머신 엡실론은 정확히 구별 가능한 결과를 만드는 가장 작은 증분이기 때문이다. 수 $1 + \varepsilon$는 $1$ 다음에 표현 가능한 부동 소수점 수이다.

> **참고:** 머신 엡실론은 **상대적** 측정값이다. $1$ 근처에서 절대 정밀도는 $\varepsilon \approx 2.22 \times 10^{-16}$이다. $10^{10}$ 근처에서 절대 정밀도는 $10^{10} \times \varepsilon \approx 2.22 \times 10^{-6}$이다. 상대 정밀도는 항상 동일하다.

### 3.2 소수의 이진 표현

사람에게 "단순해" 보이는 많은 십진 소수(예: $0.1$ 또는 $0.3$)가 $1/3 = 0.333...$이 십진법에서 비종단적인 것처럼 **이진법에서 비종단적 표현**을 가진다. 이것은 부동 소수점 산술에서 반올림 오차의 근본적인 원인이다.

```python
def decimal_to_binary(num, bits=52):
    """Show binary representation of a decimal fraction."""
    result = ''
    for _ in range(bits):
        num *= 2
        if num >= 1:
            result += '1'
            num -= 1
        else:
            result += '0'
    return '0.' + result

print(f'0.1 in binary: {decimal_to_binary(0.1, 20)}')  # non-terminating!
print(f'0.3 in binary: {decimal_to_binary(0.3, 20)}')  # also non-terminating!
# This is why 0.1 + 0.2 != 0.3 in floating-point!
```

알고리즘은 소수를 반복적으로 2배하는 방식으로 동작한다. 결과가 1을 초과하면 다음 이진 자릿수가 1(그리고 1을 빼줌)이고, 그렇지 않으면 0이다. $0.1$의 경우:

$$0.1 \times 2 = 0.2 \to 0, \quad 0.2 \times 2 = 0.4 \to 0, \quad 0.4 \times 2 = 0.8 \to 0, \quad 0.8 \times 2 = 1.6 \to 1, \quad \ldots$$

패턴 $0.0\overline{0011}$이 무한히 반복된다. 컴퓨터가 52비트의 가수만 저장할 수 있으므로, 표현이 절삭되어 작은 오차가 도입된다. 이것이 유명한 `0.1 + 0.2 != 0.3` 놀라움이 IEEE 754 부동 소수점을 사용하는 거의 모든 프로그래밍 언어에서 발생하는 이유이다.

### 3.3 반올림 오차의 누적

많은 작은 반올림 오차가 수천 번의 연산에 걸쳐 누적되면, 총 오차가 상당해질 수 있다. 이것은 $0.0001$을 만 번 합산하여 시연된다 — 수학적 결과는 정확히 $1.0$이지만, 부동 소수점 산술은 불일치를 도입한다:

```python
# Summing 0.0001 ten thousand times
total_naive = sum(0.0001 for _ in range(10000))
print(f'Naive sum: {total_naive}')
print(f'Error: {abs(1.0 - total_naive):.2e}')  # ~9.38e-14
```

$0.0001$의 각 덧셈은 작은 반올림 오차를 도입한다($0.0001$이 이진법으로 정확히 표현될 수 없기 때문). 10,000번의 덧셈에 걸쳐 이 오차들이 약 $10^{-13}$까지 누적된다 — 절대값으로는 작지만, 높은 정확도가 요구되는 응용이나 이 합을 추가 계산에 사용하는 경우(오차가 증폭될 수 있음) 잠재적으로 문제가 된다.

### 3.4 카한 합산 알고리즘

**카한 합산 알고리즘(Kahan Summation Algorithm)**(**보상 합산(Compensated Summation)**이라고도 함)은 대규모 합산에서 반올림 오차의 누적을 극적으로 줄이는 우아한 기법이다. 핵심 아이디어는 실행 중인 반올림 오차를 추적하고 다음 덧셈에서 이를 보정하는 별도의 "보상" 변수를 유지하는 것이다:

```python
def kahan_sum(values):
    """Kahan summation for reduced round-off error."""
    total = 0.0
    compensation = 0.0
    for val in values:
        y = val - compensation          # compensate for lost low-order bits
        temp = total + y
        compensation = (temp - total) - y  # recover what was lost
        total = temp
    return total

total_kahan = kahan_sum(0.0001 for _ in range(10000))
print(f'Kahan sum: {total_kahan}')
print(f'Error: {abs(1.0 - total_kahan):.2e}')  # 0.00e+00 — perfect!
```

알고리즘이 단계별로 어떻게 동작하는지:

1. `y = val - compensation`: 누적된 보상(이전 단계의 오차)을 빼서 다음에 더할 값을 조정한다.
2. `temp = total + y`: 덧셈을 수행한다(이때 새로운 반올림 오차가 도입됨).
3. `compensation = (temp - total) - y`: 2단계에서 도입된 반올림 오차를 계산한다. 정확한 산술에서 `temp = total + y`이므로, `(temp - total) - y`는 0이어야 하지만 — 부동 소수점에서는 손실된 반올림 오차를 정확히 포착한다.
4. `total = temp`: 누적 합계를 갱신한다.

> **참고:** 카한 합산은 누적된 반올림 오차를 추적하고 다음 덧셈에서 이를 보정하는 "보상" 변수를 추가한다. 이것은 대규모 합산에서 정확도를 유지하기 위한 가장 중요한 수치 기법 중 하나이다.

결과는 놀라울 정도이다: 순진한 합산은 $\sim 10^{-13}$의 오차를 생성하는 반면, 카한 합산은 이 테스트 케이스에서 **오차가 0**을 달성한다. 보상 메커니즘은 사실상 단일 부동 소수점 변수의 두 배 정밀도를 제공한다.

### 3.5 재앙적 소거 — 이차 방정식 공식

$ax^2 + bx + c = 0$에 대한 표준 이차 방정식 공식:

$$x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$$

은 $b^2 \gg 4ac$일 때 **재앙적 소거(Catastrophic Cancellation)**를 겪는다. 이 경우, $\sqrt{b^2 - 4ac} \approx |b|$이므로, 두 근 중 하나는 거의 같은 두 수를 빼는 것을 포함한다.

구체적으로, $b > 0$일 때, 근 $x_1 = \frac{-b + \sqrt{b^2 - 4ac}}{2a}$은 $b$에서 $\sqrt{b^2 - 4ac}$를 빼는 것을 포함한다 — $b$가 클 때 거의 동일한 두 수이다.

안정적인 대안은 **비에타의 공식(Vieta's Formulas)**을 사용한다: $ax^2 + bx + c = 0$에서, 근의 곱은 $x_1 \cdot x_2 = c/a$이다. 따라서 하나의 근을 정확하게 계산할 수 있으면, 뺄셈 없이 다른 근을 $x_2 = c/(a \cdot x_1)$로 얻을 수 있다.

```python
# Standard quadratic formula for ax² + bx + c = 0:
# x = (-b ± sqrt(b²-4ac)) / (2a)
# When b >> sqrt(b²-4ac), one root suffers catastrophic cancellation

a, c = 1, 1

print(f'{"b":>12} {"x2_standard":>20} {"x2_stable":>20} {"rel_error":>15}')
print('-' * 70)

for exp in range(2, 11):
    b = 10**exp
    disc = math.sqrt(b**2 - 4*a*c)

    # Standard formula (cancellation in x2)
    x1_std = (-b + disc) / (2*a)
    x2_std = (-b - disc) / (2*a)

    # Stable alternative: x2 = c / (a * x1)
    x2_stable = c / (a * x1_std)

    rel_err = abs(x2_std - x2_stable) / abs(x2_stable) * 100
    print(f'{b:12.0e} {x2_std:20.10f} {x2_stable:20.10f} {rel_err:15.6e}%')
```

> **[미적분]** 안정적인 공식은 비에타의 공식을 사용한다: $ax^2 + bx + c = 0$에서, 근의 곱은 $x_1 \cdot x_2 = c/a$이다. 따라서 $x_1$을 정확하게 계산할 수 있으면, 뺄셈 없이 $x_2 = c/(a \cdot x_1)$로 얻을 수 있다.

$b$가 증가하면, 표준 공식은 $x_1$($-b + \sqrt{b^2 - 4ac}$ 뺄셈으로 계산되는 근)에서 점점 더 많은 유효 자릿수를 잃는 반면, 안정적인 대안은 전체 정밀도를 유지한다. $b = 10^{10}$의 경우, 표준 공식은 거의 모든 유효 자릿수를 잃었을 수 있지만, 비에타의 접근법은 정확하게 유지된다.

교훈은 일반적이다: 공식에서 거의 같은 양의 뺄셈을 만날 때마다, 뺄셈을 피하는 대수적으로 동등한 형태를 찾아라. 일반적인 기법에는 다음이 포함된다:

- **유리화(Rationalizing)**: 켤레를 곱함 (예: $\frac{(\sqrt{a} - \sqrt{b})(\sqrt{a} + \sqrt{b})}{\sqrt{a} + \sqrt{b}} = \frac{a - b}{\sqrt{a} + \sqrt{b}}$)
- **비에타의 공식**: 근의 곱 관계 사용
- **테일러 전개**: 항이 거의 소거될 때, 전개하여 소거되지 않는 선행 항을 찾기

---

<br>

## 요약

| 주제 | 핵심 개념 |
|:------|:-----------|
| 스카보로 기준(Scarborough Criterion) | $|\varepsilon_a| < 0.5 \times 10^{2-n}\%$ → n 유효 숫자 |
| 바빌로니아 방법(Babylonian Method) | $x_{n+1} = (x_n + S/x_n)/2$ — 이차 수렴 |
| 테일러 급수(Taylor Series) | $e^x = \sum x^k/k!$ — 중지 기준을 사용한 반복적 근사 |
| 머신 엡실론(Machine Epsilon) | `np.finfo(float).eps` $\approx 2.22 \times 10^{-16}$ |
| 카한 합산(Kahan Summation) | 정확한 대규모 합산을 위한 보상 기법 |
| 재앙적 소거(Catastrophic Cancellation) | 거의 같은 수의 뺄셈을 피하기 위해 대수적 대안 사용 |

---
