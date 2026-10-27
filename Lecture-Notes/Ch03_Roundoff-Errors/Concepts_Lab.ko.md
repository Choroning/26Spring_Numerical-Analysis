# 3장 실습 — 수치 오차의 원인

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **최종 수정일:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 4

> **학습 목표:**
> 1. 스카보로 기준으로 중지 허용 오차를 정하고 반복 계산에 적용할 수 있다
> 2. 테일러 급수의 절단 오차가 항의 개수에 따라 줄어드는 모습을 계산으로 확인할 수 있다
> 3. 머신 엡실론과 소수의 이진 표현으로 부동 소수점 연산의 반올림 오차를 분석할 수 있다
> 4. 카한 합산과 대수적 변형으로 반올림 오차의 누적과 재앙적 소거를 줄일 수 있다

---

## 목차

- [1. 오차 정의와 지표](#1-오차-정의와-지표)
  - [1.1 스카보로 기준 표](#11-스카보로-기준-표)
  - [1.2 sqrt(2) 근사: 바빌로니아/헤론 방법](#12-sqrt2-근사-바빌로니아헤론-방법)
- [2. 절단 오차와 테일러 급수](#2-절단-오차와-테일러-급수)
  - [2.1 e^x의 매클로린 급수: 항별 계산](#21-ex의-매클로린-급수-항별-계산)
  - [2.2 중지 기준을 사용한 반복 계산](#22-중지-기준을-사용한-반복-계산)
  - [2.3 절단 오차 수렴의 시각화](#23-절단-오차-수렴의-시각화)
- [3. 반올림 오차와 부동 소수점](#3-반올림-오차와-부동-소수점)
  - [3.1 머신 엡실론](#31-머신-엡실론)
  - [3.2 소수의 이진 표현](#32-소수의-이진-표현)
  - [3.3 반올림 오차의 누적](#33-반올림-오차의-누적)
  - [3.4 카한 합산 알고리즘](#34-카한-합산-알고리즘)
  - [3.5 재앙적 소거: 이차 방정식 공식](#35-재앙적-소거-이차-방정식-공식)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 오차 정의와 지표

### 1.1 스카보로 기준 표

**스카보로 기준(Scarborough Criterion)** 은 원하는 유효 숫자 개수에서 필요한 근사 백분율 상대 오차를 바로 구하게 해 준다. 이 기준에 따르면 $`|\varepsilon_a| < (0.5 \times 10^{2-n})\%`$일 때 결과는 최소 $`n`$개의 유효 숫자까지 정확하다.

간단한 람다 함수로 각 $`n`$ 값의 임계값을 계산할 수 있다.

```python
import numpy as np
import matplotlib.pyplot as plt
import math

# Scarborough criterion: maximum epsilon_a for n significant figures
scarborough = lambda n: 0.5 * 10**(2 - n)

# Table of criteria
for n in range(1, 8):
    print(f'{n} sig. figs → |ε_a| < {scarborough(n):.4e}%')
# 1 → 5.0000e+00%, 2 → 5.0000e-01%, ..., 7 → 5.0000e-06%
```

이 표는 실용적인 참고 자료이다. 반복 계산을 시작하기 전에 필요한 유효 숫자 개수를 정하고, 해당 임계값을 찾아 중지 기준 $`\varepsilon_s`$로 사용한다.

### 1.2 sqrt(2) 근사: 바빌로니아/헤론 방법

**바빌로니아 방법(Babylonian Method)** 은 고대 메소포타미아까지 거슬러 올라가는 가장 오래된 반복 알고리즘 중 하나이며, **헤론 방법(Heron's Method)** 이라고도 한다. 다음 점화식으로 $`\sqrt{S}`$를 계산한다.

```math
x_{n+1} = \frac{1}{2}\left(x_n + \frac{S}{x_n}\right)
```

직관은 간단하다. $`x_n`$이 $`\sqrt{S}`$의 과대추정이면 $`S/x_n`$은 과소추정이고, 그 반대도 성립한다. 둘의 평균을 취하면 두 값 각각보다 $`\sqrt{S}`$에 더 가까운 추정값을 얻는다.

> **참고:** **뉴턴-랩슨법과의 관계.** 바빌로니아 방법은 $`f(x) = x^2 - S`$에 뉴턴-랩슨법(6장)을 적용한 특수한 경우이다. 뉴턴-랩슨법의 반복식은 $`x_{n+1} = x_n - f(x_n)/f'(x_n) = x_n - (x_n^2 - S)/(2x_n) = (x_n + S/x_n)/2`$이다. 이 관계가 이차 수렴(매 반복의 오차가 대략 이전 오차의 제곱에 비례함)을 설명한다.

다음 코드는 참 백분율 상대 오차 $`\varepsilon_t`$와 근사 백분율 상대 오차 $`\varepsilon_a`$를 함께 추적하면서 방법이 얼마나 빠르게 수렴하는지 보여 준다.

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

실행 결과에서 **이차 수렴(Quadratic Convergence)** 이 나타난다. 정확한 자릿수가 매 반복 대략 두 배로 늘어난다. 5번째 반복에서 참 오차는 $`\varepsilon_t \approx 1.57 \times 10^{-14}\%`$(상대 오차로 약 $`1.6 \times 10^{-16}`$)로 머신 엡실론 수준에 이르며, 알고리즘이 $`\sqrt{2}`$의 배정밀도 최선 표현으로 수렴한다.

이는 핵심 원리를 보여 준다. 잘 설계된 반복법에서는 $`\varepsilon_a`$가 $`\varepsilon_t`$를 밀접하게 따라가므로(이 예에서 각 반복의 $`\varepsilon_a`$는 직전 반복의 $`\varepsilon_t`$와 거의 같다), 참값을 모를 때에도 $`\varepsilon_a`$를 신뢰할 수 있는 대리 지표로 사용할 수 있다.

---

<br>

## 2. 절단 오차와 테일러 급수

### 2.1 e^x의 매클로린 급수: 항별 계산

**매클로린 급수(Maclaurin Series)** 는 테일러 급수를 $`x = 0`$ 근방에서 전개한 것이며, $`e^x`$의 매클로린 급수는 다음과 같다.

```math
e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots = \sum_{k=0}^{\infty} \frac{x^k}{k!}
```

이 급수는 모든 실수 $`x`$에서 수렴하므로, 항을 충분히 포함하면 원하는 정확도까지 $`e^x`$를 근사할 수 있다. 항을 하나 더할 때마다 절단 오차가 줄어든다.

다음 코드는 $`e^{0.5}`$의 근사를 항별로 쌓아 가며 항을 더할 때마다 오차가 어떻게 줄어드는지 보여 준다.

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

참값은 $`e^{0.5} \approx 1.6487212707`$이다. 항을 더 포함할수록 근사값이 빠르게 수렴한다.

- 0차($`1`$만 사용)의 오차는 약 $`39.3\%`$이다.
- 1차($`1 + x = 1.5`$)의 오차는 약 $`9.02\%`$이다.
- 5차의 오차는 약 $`0.0014\%`$로 무시할 수 있을 만큼 작다.

분모의 팩토리얼이 분자의 $`x`$ 거듭제곱보다 훨씬 빠르게 커지므로, 적당한 $`x`$ 값에서 급수가 매우 빠르게 수렴한다.

### 2.2 중지 기준을 사용한 반복 계산

실제로는 각 항을 하드 코딩하지 않는다. 대신 근사 오차가 임계값 아래로 떨어질 때까지 항을 반복해서 더하는 루프를 작성한다.

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

이 구현에서 관찰할 핵심 사항은 다음과 같다.

- 중지 기준 `es=1e-4`는 $`\varepsilon_s = 0.0001\%`$에 해당한다. 이 값은 유효 숫자 5자리의 임계값 $`0.5 \times 10^{2-5} = 0.0005\%`$보다 작지만 6자리의 임계값 $`0.5 \times 10^{2-6} = 0.00005\%`$보다 크므로, 스카보로 기준은 최소 5자리의 유효 숫자를 보장한다. 실제로 `exp_series(0.5)`는 7번 반복 후 $`\varepsilon_a \approx 9.40 \times 10^{-5}\%`$에서 멈추며, 참 상대 오차는 약 $`6.2 \times 10^{-6}\%`$이다.
- `maxit` 매개변수는 안전 장치이다. $`|x|`$가 매우 커서 항이 많이 필요한 경우처럼 `maxit` 안에 기준을 만족하지 못하면 루프가 무한히 돌지 않고 최대 반복 횟수 후 끝난다.
- 함수는 호출자에게 필요한 세 가지 정보, 즉 결과, 최종 오차 추정값, 사용한 반복 횟수를 모두 반환한다.

### 2.3 절단 오차 수렴의 시각화

반로그 그래프는 항을 더할수록 절단 오차가 매우 빠르게 줄어드는 모습을 보여 준다.

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

**반로그 그래프에서 오차는 직선보다 가파르게 떨어지는 곡선**을 그린다. 항을 하나 더할 때 오차가 줄어드는 비율이 약 0.16에서 0.04로 점점 작아지므로, 감소 속도는 일정한 비율로 줄어드는 기하급수적 감소보다도 빠르다. 다음 항 $`x^{n+1}/(n+1)!`$의 분모가 팩토리얼로 커지기 때문이며, $`e^x`$처럼 모든 $`x`$에서 수렴하는 해석적 함수의 테일러 급수에서 나타나는 특징이다.

코드의 가로축 값은 차수 $`n`$이며 항의 개수는 $`n + 1`$이다. $`e^{0.5}`$의 경우 13차(14개 항)에서 절대 오차가 약 $`6.7 \times 10^{-16}`$으로 머신 엡실론($`\approx 2.22 \times 10^{-16}`$) 수준에 이르고, 14차에서는 계산값이 `math.exp(0.5)`와 같아져 오차가 0이 된다(반로그 그래프에는 0을 그릴 수 없으므로 이 점은 나타나지 않는다). 이 지점 이후로는 항을 더해도 더 이상 개선되지 않는데, 반올림 오차가 지배하기 때문이다. 이 지점이 절단 오차와 반올림 오차가 만나는 곳이다.

---

<br>

## 3. 반올림 오차와 부동 소수점

### 3.1 머신 엡실론

**머신 엡실론(Machine Epsilon)** 은 1과 그다음으로 큰 부동 소수점 수 사이의 간격 $`\varepsilon`$이며, 컴퓨터 산술에서 $`1 + \varepsilon \neq 1`$을 만족하는 가장 작은 2의 거듭제곱이기도 하다. 이는 부동 소수점 형식의 근본적인 정밀도 한계를 나타낸다.

```python
# Machine epsilon
eps = np.finfo(float).eps
print(f'Machine epsilon: {eps}')  # 2.220446049250313e-16

# Demonstration: 1 + eps/2 == 1 is True!
print(1 + eps/2 == 1)   # True (eps/2 is below representable precision)
print(1 + eps == 1)      # False (eps is exactly at the boundary)
```

첫 번째 검사(`1 + eps/2 == 1`)는 `True`를 반환한다. $`1 + \varepsilon/2`$는 $`1`$과 $`1 + \varepsilon`$의 정확히 중간에 있고 두 수 사이에는 64비트 부동 소수점 수가 없으므로, 최근접 짝수 반올림에 따라 1로 반올림되기 때문이다. 즉 `eps/2`는 64비트 부동 소수점에서 1의 표현에 영향을 주지 못할 만큼 작다.

두 번째 검사(`1 + eps == 1`)는 `False`를 반환한다. $`1 + \varepsilon`$는 $`1`$ 다음에 표현 가능한 부동 소수점 수이므로 반올림 없이 정확히 표현되어 1과 구별된다. 머신 엡실론은 1에 더해 구별 가능한 결과를 만드는 가장 작은 2의 거듭제곱이며, 1과 구별되는 결과를 만드는 가장 작은 수 자체는 $`\varepsilon/2`$보다 약간 큰 수이다.

> **참고:** 머신 엡실론은 **상대적** 척도이다. $`1`$ 근처에서 인접한 수의 간격(절대 정밀도)은 $`\varepsilon \approx 2.22 \times 10^{-16}`$이고, $`10^{10}`$ 근처에서는 약 $`10^{10} \times \varepsilon \approx 2.22 \times 10^{-6}`$ 수준이다(실제 간격은 약 $`1.9 \times 10^{-6}`$). 절대 정밀도는 수의 크기에 비례해 커지지만 상대 정밀도는 크기와 관계없이 거의 같다.

### 3.2 소수의 이진 표현

사람에게 단순해 보이는 많은 십진 소수(예: $`0.1`$, $`0.3`$)는 $`1/3 = 0.333...`$이 십진법에서 끝나지 않는 것처럼 **이진법에서 끝나지 않는 표현**을 가진다. 이것이 부동 소수점 산술에서 반올림 오차의 근본적인 원인이다.

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

알고리즘은 소수를 반복해서 2배 하는 방식으로 동작한다. 결과가 1 이상이면 다음 이진 자릿수는 1이고 1을 빼며, 그렇지 않으면 0이다. $`0.1`$의 경우는 다음과 같다.

```math
0.1 \times 2 = 0.2 \to 0, \quad 0.2 \times 2 = 0.4 \to 0, \quad 0.4 \times 2 = 0.8 \to 0, \quad 0.8 \times 2 = 1.6 \to 1, \quad \ldots
```

패턴 $`0.0\overline{0011}`$이 무한히 반복된다. 배정밀도는 가수 52비트(암묵적 선행 1을 포함하면 53비트)만 저장하므로 표현이 유한한 비트로 반올림되어 작은 오차가 생긴다. 예를 들어 0.1로 저장되는 값은 정확히 0.1000000000000000055511151231257827...이다. 이것이 IEEE 754 부동 소수점을 쓰는 거의 모든 프로그래밍 언어에서 `0.1 + 0.2 != 0.3`이 되는 이유이다.

> **[컴퓨터구조]** IEEE 754 배정밀도 형식은 수를 부호 1비트, 지수 11비트, 가수 52비트의 2진 비트 패턴으로 저장하므로, 유한한 비트로 정확히 나타낼 수 있는 소수는 분모가 2의 거듭제곱인 분수(0.5, 0.25, 0.375 등)로 한정된다. 0.1 = 1/10처럼 분모에 2 이외의 소인수가 있는 십진 소수는 위와 같이 2진 전개가 끝나지 않으므로 가수의 마지막 비트에서 반올림되어 저장되며, 이 저장 단계의 반올림이 이 절에서 말한 오차의 출처이다.

### 3.3 반올림 오차의 누적

많은 작은 반올림 오차가 수천 번의 연산에 걸쳐 누적되면 총 오차가 상당해질 수 있다. $`0.0001`$을 만 번 더하는 예가 이를 보여 준다. 수학적 결과는 정확히 $`1.0`$이지만 부동 소수점 산술은 불일치를 만든다.

```python
# Summing 0.0001 ten thousand times
total_naive = sum(0.0001 for _ in range(10000))
print(f'Naive sum: {total_naive}')
print(f'Error: {abs(1.0 - total_naive):.2e}')  # ~9.38e-14
```

> **참고:** **Python 버전에 따른 차이.** 위 코드의 주석 `~9.38e-14`는 Python 3.11 이하의 결과이다. Python 3.12부터 내장 함수 `sum()`은 부동 소수점 합에 보상 합산을 사용하므로, 이 코드는 `1.0`과 오차 `0.00e+00`을 출력한다. 단순 합산의 오차를 보려면 `total = 0.0`에서 시작해 `for` 루프로 `total += 0.0001`을 만 번 반복하면 되며, 이때 결과는 `0.9999999999999062`, 오차는 약 `9.38e-14`이다.

$`0.0001`$의 각 덧셈은 작은 반올림 오차를 만든다($`0.0001`$은 이진법으로 정확히 표현할 수 없기 때문이다). 10,000번의 덧셈에 걸쳐 이 오차들이 약 $`10^{-13}`$까지 누적된다. 절댓값으로는 작지만, 높은 정확도가 필요한 응용이나 이 합을 이후 계산에 사용하여 오차가 증폭될 수 있는 경우에는 문제가 될 수 있다.

### 3.4 카한 합산 알고리즘

**카한 합산 알고리즘(Kahan Summation Algorithm)** 은 대규모 합산에서 반올림 오차의 누적을 크게 줄이는 기법이며, **보상 합산(Compensated Summation)** 이라고도 한다. 핵심 아이디어는 진행 중인 반올림 오차를 추적하여 다음 덧셈에서 보정하는 별도의 보상 변수를 두는 것이다.

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

알고리즘은 단계별로 다음과 같이 동작한다.

1. `y = val - compensation`: 누적된 보상(이전 단계의 오차)을 빼서 다음에 더할 값을 조정한다.
2. `temp = total + y`: 덧셈을 수행한다(이때 새로운 반올림 오차가 생긴다).
3. `compensation = (temp - total) - y`: 2단계에서 생긴 반올림 오차를 계산한다. 정확한 산술에서는 `temp = total + y`이므로 `(temp - total) - y`가 0이어야 하지만, 부동 소수점에서는 손실된 반올림 오차를 정확히 포착한다.
4. `total = temp`: 누적 합계를 갱신한다.

> **참고:** 카한 합산은 누적된 반올림 오차를 추적하고 다음 덧셈에서 이를 보정하는 보상 변수를 추가한다. 이것은 대규모 합산에서 정확도를 유지하기 위한 가장 중요한 수치 기법 중 하나이다.

결과는 인상적이다. 이 테스트에서 단순 합산은 $`\sim 10^{-13}`$의 오차를 내지만(위 참고 상자의 `for` 루프 기준), 카한 합산은 **오차 0**을 얻는다. 보상 메커니즘은 사실상 부동 소수점 변수 하나의 두 배에 가까운 정밀도를 제공한다.

### 3.5 재앙적 소거: 이차 방정식 공식

$`ax^2 + bx + c = 0`$에 대한 표준 근의 공식은 다음과 같다.

```math
x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
```

이 공식은 $`b^2 \gg 4ac`$일 때 **재앙적 소거(Catastrophic Cancellation)**(강의 노트 6.1절의 뺄셈 소거)를 겪는다. 이 경우 $`\sqrt{b^2 - 4ac} \approx |b|`$이므로 두 근 중 하나를 계산할 때 거의 같은 두 수를 빼게 된다.

구체적으로 $`b > 0`$이면 근 $`x_1 = \frac{-b + \sqrt{b^2 - 4ac}}{2a}`$은 $`\sqrt{b^2 - 4ac}`$에서 $`b`$를 빼는 계산을 포함하며, $`b`$가 클 때 이 두 수는 거의 같다.

안정적인 대안은 **비에타의 공식(Vieta's Formulas)** 을 사용한다. $`ax^2 + bx + c = 0`$에서 근의 곱은 $`x_1 \cdot x_2 = c/a`$이다. 따라서 한 근을 정확하게 계산할 수 있으면 다른 근을 뺄셈 없이 $`x_2 = c/(a \cdot x_1)`$로 얻을 수 있다.

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

> **참고:** **위 코드의 근 선택 오류.** $`b > 0`$이면 소거가 일어나는 근은 `x1_std = (-b + disc) / (2*a)`이고, `x2_std = (-b - disc) / (2*a)`는 소거 없이 정확하다. 그런데 위 코드는 주석과 달리 부정확한 `x1_std`로 `x2_stable = c / (a * x1_std)`를 계산하므로, 표의 `x2_stable` 열이 오히려 $`b`$가 커질수록 부정확해진다($`b = 10^8`$에서 `-134217728`). 또 $`b = 10^9`$에서는 `math.sqrt(b**2 - 4*a*c)`가 $`b`$와 같아져 `x1_std`가 0이 되므로 `ZeroDivisionError`가 발생한다. 비에타의 공식은 정확한 큰 근 $`x_2`$로 작은 근을 $`x_1 = c/(a \cdot x_2)`$로 구하는 방향으로 써야 한다. 다음은 이렇게 고친 코드이다.

```python
# Corrected: compute the large root without cancellation, then the small root by Vieta
a, c = 1, 1

print(f'{"b":>12} {"x1_standard":>22} {"x1_stable":>22} {"rel_error":>15}')
print('-' * 76)

for exp in range(2, 11):
    b = 10**exp
    disc = math.sqrt(b**2 - 4*a*c)

    x1_std = (-b + disc) / (2*a)    # small root: cancellation
    x2 = (-b - disc) / (2*a)        # large root: no cancellation

    x1_stable = c / (a * x2)        # Vieta: x1 * x2 = c / a

    rel_err = abs(x1_std - x1_stable) / abs(x1_stable) * 100
    print(f'{b:12.0e} {x1_std:22.15e} {x1_stable:22.15e} {rel_err:15.6e}%')
```

> **참고:** **비에타의 공식을 쓰는 방향.** 근의 곱 $`x_1 \cdot x_2 = c/a`$를 쓰면 어느 근에서든 다른 근을 뺄셈 없이 구할 수 있지만, 기준이 되는 근은 소거 없이 정확하게 계산한 근이어야 한다. 고친 코드는 $`b > 0`$일 때 정확한 큰 근 $`x_2`$를 기준으로 삼아 작은 근을 $`x_1 = c/(a \cdot x_2)`$로 구한다.

$`b`$가 커지면 표준 공식은 $`x_1`$($`-b + \sqrt{b^2 - 4ac}`$ 뺄셈으로 계산되는 근)에서 점점 더 많은 유효 자릿수를 잃지만, 안정적인 대안은 전체 정밀도를 유지한다. 고친 코드에서 표준 공식의 $`x_1`$은 $`b = 10^8`$에서 상대 오차가 약 25%이고 $`b \geq 10^9`$에서는 0이 되어 유효 자릿수를 모두 잃지만, 비에타 방식의 $`x_1`$은 모든 $`b`$에서 약 $`-1/b`$로 정확하게 유지된다.

교훈은 일반적이다. 공식에서 거의 같은 양의 뺄셈을 만나면 그 뺄셈을 피하는 대수적으로 동등한 형태를 찾는다. 일반적인 기법은 다음과 같다.

- **유리화(Rationalizing):** 켤레를 곱한다. 예를 들어 $`\frac{(\sqrt{a} - \sqrt{b})(\sqrt{a} + \sqrt{b})}{\sqrt{a} + \sqrt{b}} = \frac{a - b}{\sqrt{a} + \sqrt{b}}`$이다.
- **비에타의 공식:** 근의 곱 관계를 사용한다.
- **테일러 전개:** 항이 거의 소거될 때 전개하여 소거되지 않는 선행 항을 찾는다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:-----|:-----|
| 스카보로 기준(Scarborough Criterion) | $`\vert\varepsilon_a\vert < (0.5 \times 10^{2-n})\%`$이면 결과가 최소 $`n`$자리 유효 숫자까지 정확하다. |
| 바빌로니아 방법(Babylonian Method) | $`x_{n+1} = (x_n + S/x_n)/2`$로 $`\sqrt{S}`$를 구하며 이차 수렴한다. |
| 테일러 급수(Taylor Series) | $`e^x = \sum x^k/k!`$의 항을 중지 기준을 만족할 때까지 더해 반복적으로 근사한다. |
| 머신 엡실론(Machine Epsilon) | `np.finfo(float).eps` $`\approx 2.22 \times 10^{-16}`$은 1과 그다음 부동 소수점 수 사이의 간격이다. |
| 카한 합산(Kahan Summation) | 보상 변수로 손실된 반올림 오차를 추적하여 대규모 합산을 정확하게 한다. |
| 재앙적 소거(Catastrophic Cancellation) | 거의 같은 수의 뺄셈을 피하는 대수적 대안(유리화, 비에타의 공식, 테일러 전개)을 사용한다. |

---

<br>

## 점검 문제

1. **스카보로 기준과 중지 기준:** `exp_series`를 `es=1e-6`, 즉 $`\varepsilon_s = 10^{-6}\%`$로 호출하면 스카보로 기준은 최소 몇 자리의 유효 숫자를 보장하는가?

   <details>
   <summary>정답</summary>

   7자리 기준은 $`0.5 \times 10^{2-7}\% = 5 \times 10^{-6}\%`$이고 8자리 기준은 $`5 \times 10^{-7}\%`$이다. $`10^{-6}\%`$는 7자리 기준보다 작고 8자리 기준보다 크므로 최소 7자리가 보장된다.

   </details>

2. **바빌로니아 방법:** $`S = 2`$, $`x_0 = 1`$에서 시작하여 $`x_1`$과 $`x_2`$를 손으로 계산하고, $`x_2`$의 참 백분율 상대 오차를 구하라.

   <details>
   <summary>정답</summary>

   $`x_1 = (1 + 2/1)/2 = 1.5`$이고 $`x_2 = (1.5 + 2/1.5)/2 = (1.5 + 1.333333)/2 \approx 1.416667`$이다. $`\sqrt{2} \approx 1.414214`$이므로 $`\varepsilon_t = |1.414214 - 1.416667| / 1.414214 \times 100\% \approx 0.173\%`$이다.

   </details>

3. **머신 엡실론:** 배정밀도에서 `1 + eps/2 == 1`은 `True`인데 `1 + 0.6*eps == 1`은 `False`인 이유를 설명하라.

   <details>
   <summary>정답</summary>

   $`1`$과 $`1 + \varepsilon`$ 사이에는 부동 소수점 수가 없으므로, 그 사이의 값은 더 가까운 쪽으로 반올림된다. $`1 + 0.6\varepsilon`$은 중간보다 $`1 + \varepsilon`$에 가까워서 $`1 + \varepsilon`$이 되므로 1과 다르다. $`1 + \varepsilon/2`$는 정확히 중간이므로 최근접 짝수 반올림에 따라 1이 된다.

   </details>

4. **재앙적 소거 회피:** $`x = 10^8`$에서 $`\sqrt{x + 1} - \sqrt{x}`$를 그대로 계산할 때 생기는 문제를 설명하고, 유리화한 식을 쓰라.

   <details>
   <summary>정답</summary>

   $`\sqrt{10^8 + 1} \approx 10000.00005`$와 $`\sqrt{10^8} = 10000`$이 거의 같아서 뺄셈에서 앞자리가 소거된다. 배정밀도로 그대로 계산하면 $`5.000000055588 \times 10^{-5}`$가 나와 참값 $`4.999999998750 \times 10^{-5}`$과 상대 오차가 약 $`1.4 \times 10^{-8}`$이며, 16자리 중 약 8자리를 잃는다. 켤레를 곱해 $`\frac{1}{\sqrt{x + 1} + \sqrt{x}}`$로 쓰면 뺄셈이 없어 거의 모든 자릿수가 정확하다.

   </details>

5. **이차 방정식의 안정적 계산:** $`x^2 + 10^8 x + 1 = 0`$의 두 근을 정확하게 구하는 절차와 그 근사값을 쓰라.

   <details>
   <summary>정답</summary>

   먼저 소거가 없는 큰 근 $`x_2 = (-b - \sqrt{b^2 - 4ac})/(2a) \approx -10^8`$을 계산한다. 그다음 비에타의 공식으로 작은 근을 $`x_1 = c/(a \cdot x_2) \approx -1.0 \times 10^{-8}`$으로 구한다. 표준 공식으로 $`x_1`$을 직접 계산하면 소거 때문에 약 $`-7.45 \times 10^{-9}`$가 나와 상대 오차가 약 25%가 된다.

   </details>
