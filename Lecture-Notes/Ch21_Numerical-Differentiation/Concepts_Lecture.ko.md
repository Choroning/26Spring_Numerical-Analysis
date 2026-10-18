# 21장 — 수치 미분(Numerical Differentiation)

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 21

> **학습 목표:**
> 1. 미분을 위해 유한 차분 공식을 적용할 수 있다
> 2. 수치 미분에서 절단 오차를 분석할 수 있다
> 3. 정확도 향상을 위해 리처드슨 외삽법을 구현할 수 있다
> 4. 비등간격 데이터와 다변수 함수에 수치 미분 공식을 적용할 수 있다

---

<br>

## 목차

- [1. 수치 미분 개요](#1-수치-미분-개요)
  - [1.1 극한으로서의 도함수](#11-극한으로서의-도함수)
  - [1.2 2계도함수](#12-2계도함수)
  - [1.3 편도함수](#13-편도함수)
  - [1.4 물리적 응용: 열유속(Heat Flux)](#14-물리적-응용-열유속heat-flux)
- [2. 1계도함수의 유한 차분 근사(Finite Difference Approximation)](#2-1계도함수의-유한-차분-근사finite-difference-approximation)
  - [2.1 전진 차분(Forward Difference) (1차)](#21-전진-차분forward-difference-1차)
  - [2.2 후진 차분(Backward Difference) (1차)](#22-후진-차분backward-difference-1차)
  - [2.3 중심 차분(Central Difference) (2차)](#23-중심-차분central-difference-2차)
- [3. 2계도함수의 유한 차분 근사](#3-2계도함수의-유한-차분-근사)
  - [3.1 테일러 급수를 이용한 중심 차분](#31-테일러-급수를-이용한-중심-차분)
  - [3.2 2계도함수의 전진 차분](#32-2계도함수의-전진-차분)
- [4. 정확도 향상: 고차 전진 차분](#4-정확도-향상-고차-전진-차분)
- [5. 미분을 위한 리처드슨 외삽법(Richardson Extrapolation)](#5-미분을-위한-리처드슨-외삽법richardson-extrapolation)
  - [5.1 기본 아이디어](#51-기본-아이디어)
  - [5.2 중심 차분의 오차 구조 재검토](#52-중심-차분의-오차-구조-재검토)
  - [5.3 첫 번째 외삽 수준](#53-첫-번째-외삽-수준)
  - [5.4 두 번째 외삽 수준](#54-두-번째-외삽-수준)
- [6. 비등간격 데이터의 도함수(Derivatives of Unequally Spaced Data)](#6-비등간격-데이터의-도함수derivatives-of-unequally-spaced-data)
  - [6.1 라그랑주 다항식 접근법](#61-라그랑주-다항식-접근법)
  - [6.2 라그랑주 기저 함수의 도함수](#62-라그랑주-기저-함수의-도함수)
- [7. 수치 편도함수(Numerical Partial Derivatives)](#7-수치-편도함수numerical-partial-derivatives)
  - [7.1 1계 편도함수 (중심 차분)](#71-1계-편도함수-중심-차분)
  - [7.2 혼합 편도함수(Mixed Partial Derivative)](#72-혼합-편도함수mixed-partial-derivative)
- [8. 미정 계수법(Method of Undetermined Coefficients)](#8-미정-계수법method-of-undetermined-coefficients)
  - [8.1 문제 설정](#81-문제-설정)
  - [8.2 테일러 전개와 매칭](#82-테일러-전개와-매칭)
  - [8.3 계수 풀기](#83-계수-풀기)
  - [8.4 절단 오차(Truncation Error)](#84-절단-오차truncation-error)
- [9. 미분 가능 프로그래밍(Differentiable Programming) (간략 개요)](#9-미분-가능-프로그래밍differentiable-programming-간략-개요)
- [10. Python 구현](#10-python-구현)
  - [10.1 전진, 후진, 중심 차분](#101-전진-후진-중심-차분)
  - [10.2 미분을 위한 리처드슨 외삽법](#102-미분을-위한-리처드슨-외삽법)
  - [10.3 비등간격 데이터의 도함수 (라그랑주)](#103-비등간격-데이터의-도함수-라그랑주)
  - [10.4 수치 편도함수](#104-수치-편도함수)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 수치 미분 개요

수치 미분(Numerical Differentiation)은 해석적으로 계산하는 대신 이산 데이터나 함수 평가로부터 도함수를 추정하는 과정이다. 이는 다음과 같은 경우에 필수적이다:

- 함수가 이산 점에서만 알려진 경우 (실험/현장 데이터)
- 해석적 도함수가 너무 복잡하거나 구할 수 없는 경우
- 더 큰 수치 알고리즘(예: ODE/PDE 풀기) 내에서 도함수가 필요한 경우

### 1.1 극한으로서의 도함수

유한 차분비(Finite Difference Ratio)는 자연스러운 출발점을 제공한다:

```math
\frac{\Delta y}{\Delta x} = \frac{f(x + \Delta x) - f(x)}{\Delta x}
```

$`\Delta x \to 0`$으로 극한을 취하면:

```math
\frac{dy}{dx} = \lim_{\Delta x \to 0} \frac{f(x + \Delta x) - f(x)}{\Delta x} =: f'(x)
```

이것이 **$`x`$에 대한 $`y`$의 1계도함수**이다. 기하학적으로 도함수는 점 $`x_i`$에서 곡선에 대한 접선의 기울기이다.

> **[미적분학]** $`\Delta x`$가 줄어들면 할선(곡선 위의 두 점을 잇는 선)은 $`x_i`$에서의 접선에 가까워지고, 할선의 기울기는 $`f'(x_i)`$에 수렴한다.

### 1.2 2계도함수

2계도함수는 1계도함수의 도함수이다:

```math
\frac{d^2 y}{dx^2} = \frac{d}{dx}\left(\frac{dy}{dx}\right)
```

이것은 **기울기가 얼마나 빠르게 변하는지**를 측정하며, 함수 그래프가 휘는 방향(오목성)과 정도를 나타낸다.

### 1.3 편도함수

두 변수의 함수 $`f(x, y)`$에 대해:

```math
\frac{\partial f}{\partial x} = \lim_{\Delta x \to 0} \frac{f(x + \Delta x, y) - f(x, y)}{\Delta x}
```

```math
\frac{\partial f}{\partial y} = \lim_{\Delta y \to 0} \frac{f(x, y + \Delta y) - f(x, y)}{\Delta y}
```

이것들은 각각 $`x`$와 $`y`$에 대한 $`f`$의 **편도함수**이다. 각 편도함수는 다른 변수를 고정한 채 $`f`$의 변화율을 측정한다.

### 1.4 물리적 응용: 열유속(Heat Flux)

대표적인 공학적 응용은 푸리에의 열전도 법칙이다. 온도 분포 $`T(x)`$가 주어지면, 열유속은:

```math
\vec{q} = -k \frac{\partial T}{\partial x}
```

여기서 $`k`$는 열전도율이다. 음의 부호는 열이 고온에서 저온으로 흐름을 나타낸다:

- **양의 기울기**($`\partial T / \partial x > 0`$)는 **음의 열유속**(열이 $`-x`$ 방향으로 흐름)을 만든다.
- **음의 기울기**($`\partial T / \partial x < 0`$)는 **양의 열유속**(열이 $`+x`$ 방향으로 흐름)을 만든다.

> **참고:** **열 흐름의 방향.** 열은 항상 온도 기울기의 반대 방향으로 흐른다. 이것은 도함수 개념의 직접적인 물리적 응용이다.

---

<br>

## 2. 1계도함수의 유한 차분 근사(Finite Difference Approximation)

이 장의 등간격 유한 차분 공식은 테일러 급수 전개(Taylor Series Expansion)로부터 유도한다. $`h`$를 등간격 점들 사이의 간격 크기로 놓으면 $`x_{i+1} = x_i + h`$이다.

### 2.1 전진 차분(Forward Difference) (1차)

$`x_i`$에서의 전진 방향 **테일러 전개**:

```math
f(x_{i+1}) = f(x_i) + f'(x_i)\,h + f''(x_i)\frac{h^2}{2} + f'''(x_i)\frac{h^3}{6} + \cdots \quad (\ast)
```

식 $`(\ast)`$에서 $`f'(x_i)`$를 풀면:

```math
f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - f''(x_i)\frac{h}{2} + O(h^2)
```

절단 오차항을 버리면 **전진 차분** 근사를 얻는다:

```math
\boxed{f'(x_i) \approx \frac{f(x_{i+1}) - f(x_i)}{h}}
```

- **절단 오차**: $`O(h)`$이며, 1차 정확도이다.
- **주도적 오차항**: $`-f''(x_i)\dfrac{h}{2}`$이다.
- 점 $`x_i`$와 $`x_{i+1}`$(해당 점과 앞의 한 점)을 사용한다.

### 2.2 후진 차분(Backward Difference) (1차)

마찬가지로, $`x_i`$에서 $`f(x_{i-1})`$을 전개하면:

```math
f(x_{i-1}) = f(x_i) - f'(x_i)\,h + f''(x_i)\frac{h^2}{2} - f'''(x_i)\frac{h^3}{6} + \cdots \quad (\ast\ast)
```

$`f'(x_i)`$를 풀면:

```math
\boxed{f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{h}}
```

- **절단 오차**: $`O(h)`$이며, 1차 정확도이다.
- 점 $`x_{i-1}`$과 $`x_i`$를 사용한다.

### 2.3 중심 차분(Central Difference) (2차)

$`(\ast)`$에서 $`(\ast\ast)`$를 빼면:

```math
f(x_{i+1}) - f(x_{i-1}) = 2f'(x_i)\,h + 2f'''(x_i)\frac{h^3}{6} + \cdots
```

$`f'(x_i)`$를 풀면:

```math
\boxed{f'(x_i) \approx \frac{f(x_{i+1}) - f(x_{i-1})}{2h}}
```

- **절단 오차**: $`O(h^2)`$이며, 2차 정확도이다.
- **주도적 오차항**: $`-f'''(x_i)\dfrac{h^2}{6}`$이다.
- 점 $`x_{i-1}`$과 $`x_{i+1}`$($`x_i`$에 대해 대칭)을 사용한다.

> **핵심:** 중심 차분은 대칭성 때문에 $`h`$의 홀수 거듭제곱 오차항이 상쇄되므로, 전진 차분과 후진 차분의 $`O(h)`$보다 높은 $`O(h^2)`$ 정확도를 달성한다. 이것은 양쪽 이웃 점을 모두 사용할 수 있을 때의 핵심적인 장점이다.

---

<br>

## 3. 2계도함수의 유한 차분 근사

### 3.1 테일러 급수를 이용한 중심 차분

$`(\ast)`$와 $`(\ast\ast)`$를 더하면:

```math
f(x_{i+1}) + f(x_{i-1}) = 2f(x_i) + f''(x_i)\,h^2 + O(h^4)
```

$`f''(x_i)`$를 풀면:

```math
\boxed{f''(x_i) = \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2} + O(h^2)}
```

이것이 **2계도함수의 중심 차분 근사**이다:

- **절단 오차**: $`O(h^2)`$이며, 2차 정확도이다.
- 등간격의 세 점 $`x_{i-1}`$, $`x_i`$, $`x_{i+1}`$을 사용한다.

> **참고:** **공식의 해석.** 이 공식의 분자는 $`2\left(\frac{f(x_{i+1}) + f(x_{i-1})}{2} - f(x_i)\right)`$이므로, $`x_i`$에서의 함수값이 양쪽 이웃의 평균으로부터 얼마나 벗어나는지를 측정한다. 양의 $`f''`$은 그래프가 아래로 볼록(concave up)함을, 음의 $`f''`$은 그래프가 위로 볼록(concave down)함을 의미한다.

### 3.2 2계도함수의 전진 차분

전진 점만 사용하여 $`f''(x_i)`$의 **전진** 공식을 유도하기 위해, $`f(x_{i+2})`$를 전개한다:

```math
f(x_{i+2}) = f(x_i) + f'(x_i)\,(2h) + f''(x_i)\frac{(2h)^2}{2} + f'''(x_i)\frac{(2h)^3}{6} + O(h^4) \quad (\ast\ast\ast)
```

$`(\ast\ast\ast) - 2(\ast)`$를 계산하면:

```math
f(x_{i+2}) - 2f(x_{i+1}) = -f(x_i) + f''(x_i)\,h^2 + O(h^3)
```

$`f''(x_i)`$를 풀면:

```math
\boxed{f''(x_i) = \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{h^2} + O(h)}
```

- **절단 오차**: $`O(h)`$이며, 1차 정확도이다(중심 공식보다 덜 정확하다).
- 세 개의 전진 점 $`x_i`$, $`x_{i+1}`$, $`x_{i+2}`$를 사용한다.

---

<br>

## 4. 정확도 향상: 고차 전진 차분

2계도함수 보정항을 포함하여 $`f'(x_i)`$의 전진 차분을 개선할 수 있다. 2.1절에서 얻은 다음 식에서 출발한다:

```math
f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \underbrace{f''(x_i)}_{\text{FD 근사 대입}}\frac{h}{2} + O(h^2)
```

$`f''(x_i)`$의 전진 차분 근사:

```math
f''(x_i) \approx \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{h^2}
```

를 오차항에 대입하면:

```math
f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{2h} + O(h^2)
```

결합하면:

```math
\boxed{f'(x_i) = \frac{-f(x_{i+2}) + 4f(x_{i+1}) - 3f(x_i)}{2h} + O(h^2)}
```

이것이 $`f'(x_i)`$에 대한 **2차 전진 차분** 공식이다.

- 추가 점 하나를 포함하여 정확도가 $`O(h)`$에서 $`O(h^2)`$로 **향상된다**.
- 확산(2계도함수) 항을 추가함으로써 정확도를 높인다.

지금까지 유도한 공식이 사용하는 점을 정리하면 그림 1과 같다.

```text
         x_{i-1}     x_i     x_{i+1}   x_{i+2}
  ──────────┼─────────┼─────────┼─────────┼────
            ●─────────●                          후진 차분: f'(x_i), O(h)
                      ●─────────●                전진 차분: f'(x_i), O(h)
            ●───────────────────●                중심 차분: f'(x_i), O(h^2)
            ●─────────●─────────●                중심 차분: f''(x_i), O(h^2)
                      ●─────────●─────────●      전진 차분: f''(x_i), O(h); 2차 전진 차분: f'(x_i), O(h^2)
```

*그림 1. 유한 차분 공식이 사용하는 점(스텐실)*

> **참고:** **보정항 대입과 리처드슨 외삽법.** 주도적 오차항에 수치 근사를 대입하는 이 아이디어는 정확도를 높이는 일반적인 기법이다. 이것은 5절의 리처드슨 외삽법과 밀접한 관련이 있다.

---

<br>

## 5. 미분을 위한 리처드슨 외삽법(Richardson Extrapolation)

### 5.1 기본 아이디어

20장의 수치 적분에서 리처드슨 외삽법은 서로 다른 간격 크기로 계산한 두 추정값을 사용하여 주도적 오차항을 제거했다. 오차가 $`O(h^2)`$인 사다리꼴 법칙의 적분값에 대해서는 다음과 같다:

```math
I = I(h_2) + \frac{1}{(h_1/h_2)^2 - 1}\bigl(I(h_2) - I(h_1)\bigr)
```

$`h_1/h_2 = 2`$일 때:

```math
I = \frac{4}{3}I(h_2) - \frac{1}{3}I(h_1)
```

마찬가지로, 중심 차분처럼 오차가 $`O(h^2)`$인 도함수 추정값에 대해서는 다음과 같다:

```math
\boxed{D = \frac{4}{3}D(h_2) - \frac{1}{3}D(h_1)}
```

여기서 $`D(h_1)`$과 $`D(h_2)`$는 간격 크기 $`h_1`$ (더 큰 것)과 $`h_2 = h_1/2`$ (더 작은 것)으로 계산한 도함수 추정값이다.

### 5.2 중심 차분의 오차 구조 재검토

테일러 전개에서 시작한다:

```math
f(x+h) = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}h^k
```

```math
f(x-h) = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}(-h)^k = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}(-1)^k h^k
```

$`f(x+h) - f(x-h)`$을 계산하면:

```math
f(x+h) - f(x-h) = 2h\,f'(x) + 2f'''(x)\frac{h^3}{3!} + 2f^{(5)}(x)\frac{h^5}{5!} + \cdots
```

$`2h`$로 나누면:

```math
f'(x) = \frac{f(x+h) - f(x-h)}{2h} - f'''(x)\frac{h^2}{3!} - f^{(5)}(x)\frac{h^4}{5!} - \cdots
```

중심 차분 연산자와 오차 계수를 정의한다:

```math
D(h) = \frac{f(x+h) - f(x-h)}{2h}
```

```math
e_2 = -\frac{f^{(3)}(x)}{3!}, \qquad e_4 = -\frac{f^{(5)}(x)}{5!}
```

그러면:

```math
f'(x) = D(h) + e_2\,h^2 + e_4\,h^4 + \cdots
```

여기서 오차 $`E(h) = e_2\,h^2 + e_4\,h^4 + \cdots`$이다.

> **핵심:** 오차 계수 $`e_2, e_4, \ldots`$는 $`h`$에 **독립**이며, 함수와 평가점 $`x`$에만 의존한다. 이것이 바로 리처드슨 외삽법을 가능하게 하는 성질이다.

### 5.3 첫 번째 외삽 수준

**목표**: $`h^2`$ 항을 제거한다.

간격 크기 $`h`$와 $`2h`$에서의 근사를 작성한다:

```math
f'(x) = D(h) + e_2\,h^2 + e_4\,h^4 + \cdots \quad \text{(4를 곱함)}
```

```math
f'(x) = D(2h) + e_2\,(2h)^2 + e_4\,(2h)^4 + \cdots
```

$`4 \times (\text{첫 번째}) - (\text{두 번째})`$를 계산하면:

```math
3f'(x) = 4D(h) - D(2h) - 12h^4 e_4 + O(h^6)
```

```math
\boxed{f'(x) = \underbrace{\frac{4}{3}D(h) - \frac{1}{3}D(2h)}_{=:\,\bar{D}(h)} - 4h^4 e_4 + O(h^6)}
```

결합된 추정값 $`\bar{D}(h)`$는 $`h^2`$ 오차가 제거되어 **$`O(h^4)`$** 정확도를 가진다.

### 5.4 두 번째 외삽 수준

같은 과정을 계속하여 $`h^4`$ 항을 제거한다. 간격 크기 $`h`$와 $`2h`$에서:

```math
f'(x) = \bar{D}(h) - 4h^4 e_4 + O(h^6)
```

```math
f'(x) = \bar{D}(2h) - 4(2h)^4 e_4 + O(h^6)
```

$`16 \times (\text{첫 번째}) - (\text{두 번째})`$를 계산하면:

```math
15f'(x) = 16\bar{D}(h) - \bar{D}(2h) + O(h^6)
```

```math
\boxed{f'(x) = \frac{16}{15}\bar{D}(h) - \frac{1}{15}\bar{D}(2h) + O(h^6)}
```

이제 **$`O(h^6)`$** 정확도를 가진다. 이 과정은 원리적으로 계속 반복할 수 있으며, 각 수준에서 $`h`$의 다음 짝수 거듭제곱을 제거한다.

> **참고:** **일반 패턴.** 리처드슨 외삽법은 수렴을 가속하는 체계적인 방법이다. 각 수준에서 서로 다른 간격 크기의 추정값을 결합하여 주도적 오차항을 상쇄한다. 오차가 $`h`$의 짝수 거듭제곱으로만 이루어진 중심 차분에서는 $`n`$번째 수준에서 비율 $`4^n`$을 사용하며, 가중치는 $`\frac{4^n}{4^n - 1}`$과 $`\frac{-1}{4^n - 1}`$이다.

간격 크기 $`h`$, $`2h`$, $`4h`$의 중심 차분으로 두 수준의 외삽을 수행하는 과정은 그림 2와 같다. 10.2절의 `richardson_diff` 함수는 이 표를 계산한다.

```text
  O(h^2)       O(h^4)                               O(h^6)

  D(4h) ──┐
          ├──> D̄(2h) = [4 D(2h) - D(4h)] / 3 ──┐
  D(2h) ──┤                                    ├──> [16 D̄(h) - D̄(2h)] / 15
          ├──> D̄(h)  = [4 D(h) - D(2h)] / 3  ──┘
  D(h)  ──┘
```

*그림 2. 리처드슨 외삽 표의 구성*

---

<br>

## 6. 비등간격 데이터의 도함수(Derivatives of Unequally Spaced Data)

실제로 실험이나 현장 조사에서 수집된 데이터는 흔히 **비등간격**이다. 위에서 유도한 유한 차분 공식은 등간격 점을 가정하므로 직접 적용할 수 없다.

### 6.1 라그랑주 다항식 접근법

비등간격 데이터를 다루는 한 가지 방법은 **라그랑주 다항식(Lagrange Polynomial)** 을 맞추고 미분하는 것이다.

세 데이터 점 $`(x_0, y_0)`$, $`(x_1, y_1)`$, $`(x_2, y_2)`$가 주어지면, 2차 라그랑주 다항식을 구성한다:

```math
f(x) = y_0\,N_0(x) + y_1\,N_1(x) + y_2\,N_2(x)
```

여기서 $`N_k(x)`$는 라그랑주 기저 함수이다. 미분하면:

```math
f'(x) = y_0\,N_0'(x) + y_1\,N_1'(x) + y_2\,N_2'(x)
```

### 6.2 라그랑주 기저 함수의 도함수

세 점에 대해, 기저 함수와 그 도함수는:

**기저 함수 $`N_0(x)`$**:

```math
N_0(x) = \frac{(x - x_1)(x - x_2)}{(x_0 - x_1)(x_0 - x_2)} = \frac{x^2 - (x_1 + x_2)x + x_1 x_2}{(x_0 - x_1)(x_0 - x_2)}
```

```math
N_0'(x) = \frac{2x - x_1 - x_2}{(x_0 - x_1)(x_0 - x_2)}
```

**기저 함수 $`N_1(x)`$**:

```math
N_1(x) = \frac{(x - x_0)(x - x_2)}{(x_1 - x_0)(x_1 - x_2)} = \frac{x^2 - (x_0 + x_2)x + x_0 x_2}{(x_1 - x_0)(x_1 - x_2)}
```

```math
N_1'(x) = \frac{2x - x_0 - x_2}{(x_1 - x_0)(x_1 - x_2)}
```

**기저 함수 $`N_2(x)`$**:

```math
N_2(x) = \frac{(x - x_0)(x - x_1)}{(x_2 - x_0)(x_2 - x_1)} = \frac{x^2 - (x_0 + x_1)x + x_0 x_1}{(x_2 - x_0)(x_2 - x_1)}
```

```math
N_2'(x) = \frac{2x - x_0 - x_1}{(x_2 - x_0)(x_2 - x_1)}
```

임의의 점 $`x`$에서의 도함수는:

```math
f'(x) = y_0 \cdot \frac{2x - x_1 - x_2}{(x_0 - x_1)(x_0 - x_2)} + y_1 \cdot \frac{2x - x_0 - x_2}{(x_1 - x_0)(x_1 - x_2)} + y_2 \cdot \frac{2x - x_0 - x_1}{(x_2 - x_0)(x_2 - x_1)}
```

> **참고:** **등간격 데이터와의 관계.** 이 접근법은 데이터 점의 간격에 관계없이 작동한다. 점들이 등간격이면(즉, $`x_1 - x_0 = x_2 - x_1 = h`$) 위 공식에 $`x = x_1`$을 대입하여 표준 중심 차분 공식 $`f'(x_1) \approx (y_2 - y_0)/(2h)`$를 얻는다.

---

<br>

## 7. 수치 편도함수(Numerical Partial Derivatives)

### 7.1 1계 편도함수 (중심 차분)

함수 $`f(x, y)`$에 대해, 편도함수의 중심 유한 차분 근사는:

```math
\boxed{\frac{\partial f}{\partial x} \approx \frac{f(x + \Delta x,\, y) - f(x - \Delta x,\, y)}{2\,\Delta x}}
```

```math
\boxed{\frac{\partial f}{\partial y} \approx \frac{f(x,\, y + \Delta y) - f(x,\, y - \Delta y)}{2\,\Delta y}}
```

각각 $`O(\Delta x^2)`$과 $`O(\Delta y^2)`$ 정확도를 가진다.

### 7.2 혼합 편도함수(Mixed Partial Derivative)

혼합 2계 편도함수 $`\dfrac{\partial^2 f}{\partial x \,\partial y}`$는 중심 차분을 두 번 적용하여 계산한다:

```math
\frac{\partial^2 f}{\partial x \,\partial y} = \frac{\partial}{\partial x}\left(\frac{\partial f}{\partial y}\right)
```

```math
\approx \frac{1}{2\,\Delta x}\left(\left.\frac{\partial f}{\partial y}\right|_{x+\Delta x,\,y} - \left.\frac{\partial f}{\partial y}\right|_{x-\Delta x,\,y}\right)
```

각 점에서 $`\partial f / \partial y`$에 중심 차분을 적용하면:

```math
\approx \frac{1}{2\,\Delta x}\left(\frac{f(x+\Delta x,\, y+\Delta y) - f(x+\Delta x,\, y-\Delta y)}{2\,\Delta y} - \frac{f(x-\Delta x,\, y+\Delta y) - f(x-\Delta x,\, y-\Delta y)}{2\,\Delta y}\right)
```

정리하면:

```math
\boxed{\frac{\partial^2 f}{\partial x \,\partial y} \approx \frac{f(x+\Delta x,\, y+\Delta y) - f(x+\Delta x,\, y-\Delta y) - f(x-\Delta x,\, y+\Delta y) + f(x-\Delta x,\, y-\Delta y)}{4\,\Delta x\,\Delta y}}
```

이 근사는 $`(x, y)`$를 둘러싼 대각선 위치의 네 점 $`(x \pm \Delta x, y \pm \Delta y)`$에서 함수를 평가해야 하며, 각 점의 가중치 부호는 그림 3과 같다.

```text
            x - Δx          x           x + Δx

  y + Δy      ● -1                        ● +1

     y                      ○ (x, y)

  y - Δy      ● +1                        ● -1
```

*그림 3. 혼합 편도함수 근사에 쓰이는 네 대각선 점과 가중치 부호*

---

<br>

## 8. 미정 계수법(Method of Undetermined Coefficients)

이것은 지정된 점 집합을 사용하여 도함수를 근사하는 유한 차분 공식을 유도하는 체계적인 기법이다.

### 8.1 문제 설정

등간격 세 점에서의 함수값 $`f(x-h)`$, $`f(x)`$, $`f(x+h)`$를 사용하여 2계도함수 $`f''(x)`$의 근사를 구한다:

```math
f''(x) \approx a\,f(x-h) + b\,f(x) + c\,f(x+h)
```

계수 $`a`$, $`b`$, $`c`$는 테일러 전개를 대입하고 각 도함수의 계수를 맞추어 결정한다.

### 8.2 테일러 전개와 매칭

$`f(x-h)`$와 $`f(x+h)`$를 테일러 급수로 전개한다(나머지항의 $`\xi_- \in [x-h, x]`$, $`\xi_+ \in [x, x+h]`$):

```math
f(x-h) = f(x) - f'(x)\,h + f''(x)\frac{h^2}{2} - f'''(x)\frac{h^3}{6} + f^{(4)}(\xi_-)\frac{h^4}{24}
```

```math
f(x+h) = f(x) + f'(x)\,h + f''(x)\frac{h^2}{2} + f'''(x)\frac{h^3}{6} + f^{(4)}(\xi_+)\frac{h^4}{24}
```

근사식에 대입하고 도함수별로 정리하면:

```math
a\,f(x-h) + b\,f(x) + c\,f(x+h) = (a + b + c)\,f(x) + (-a + c)\,h\,f'(x) + (a + c)\frac{h^2}{2}\,f''(x)
```
```math
\quad + (c - a)\frac{h^3}{6}\,f'''(x) + \frac{h^4}{24}\bigl(a\,f^{(4)}(\xi_-) + c\,f^{(4)}(\xi_+)\bigr)
```

### 8.3 계수 풀기

이 식이 $`f''(x)`$와 같으려면 다음이 필요하다:

| 조건 | 방정식 |
|:-----|:-------|
| $`f(x)`$의 계수 = 0 | $`a + b + c = 0`$ |
| $`f'(x)`$의 계수 = 0 | $`-a + c = 0`$ (즉, $`c = a`$) |
| $`f''(x)`$의 계수 = 1 | $`(a + c)\dfrac{h^2}{2} = 1`$ |

*표 1. 2계도함수 근사의 계수 조건*

$`c = a`$와 $`(a + c)\dfrac{h^2}{2} = 1`$로부터:

```math
2a \cdot \frac{h^2}{2} = 1 \implies a = \frac{1}{h^2}
```

따라서:

```math
\boxed{a = \frac{1}{h^2}, \quad c = \frac{1}{h^2}, \quad b = -\frac{2}{h^2}}
```

### 8.4 절단 오차(Truncation Error)

다시 대입하면 $`c - a = 0`$이므로 3계도함수 항이 사라진다. 절단 오차는 4계도함수를 포함하는 나머지항에서 온다:

```math
f''(x) = \frac{1}{h^2}\bigl(f(x-h) - 2f(x) + f(x+h)\bigr) - \frac{h^2}{24}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)
```

$`f^{(4)}`$가 연속이면 $`f^{(4)}(\xi_-)`$와 $`f^{(4)}(\xi_+)`$의 평균은 두 값 사이에 있으므로, **중간값 정리(Intermediate Value Theorem)** 에 의해 $`[x-h, x+h]`$에 다음을 만족하는 $`\xi`$가 존재한다:

```math
f^{(4)}(\xi) = \frac{1}{2}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)
```

따라서:

```math
\boxed{f''(x) = \frac{f(x-h) - 2f(x) + f(x+h)}{h^2} - \frac{h^2}{12}\,f^{(4)}(\xi)}
```

이는 $`f''`$의 중심 차분이 $`O(h^2)`$ 정확도이며, 명시적 절단 오차항이 $`-\dfrac{h^2}{12}f^{(4)}(\xi)`$임을 확인한다.

> **핵심:** 미정 계수법은 범용적이다. 포함할 점과 근사할 도함수를 선택하면 어떤 유한 차분 스텐실(stencil)이든 유도할 수 있다. 미지 계수에 대한 연립방정식은 항상 선형이므로 풀기가 간단하다.

---

<br>

## 9. 미분 가능 프로그래밍(Differentiable Programming) (간략 개요)

**미분 가능 프로그래밍**은 프로그램의 모든 것이 **자동으로 미분 가능**하도록 작성되는 현대적 프로그래밍 패러다임이다.

주요 특징은 다음과 같다:

- 입력의 작은 변화가 출력에 어떻게 영향을 미치는지 계산할 수 있다.
- **기울기 기반 최적화**(예: 신경망 학습)를 가능하게 한다.
- 전통적인 수치 방법과 데이터 기반 방법을 결합하는 **새로운 하이브리드 접근법**이다.
- 자동 미분(AD, Automatic Differentiation)을 사용하여 유한 차분보다 더 정확하고 기호적 미분보다 더 효율적인 **정확한 기울기**를 제공한다.

**미분 가능 프로그래밍을 위한 인기 있는 도구**:
- **PyTorch** (Python)
- **TensorFlow** (Python)
- **JAX** (Python, Google 개발)
- **Zygote** (Julia)

> **참고:** **자동 미분과 수치 미분.** 자동 미분은 수치 미분과 근본적으로 다르다. 수치 미분은 유한 차분으로 도함수를 근사하므로 절단 오차와 반올림 오차를 겪는다. 반면 AD는 프로그램의 모든 연산에 연쇄 법칙(Chain Rule)을 체계적으로 적용하여 기계 정밀도로 도함수를 계산한다.

---

<br>

## 10. Python 구현

### 10.1 전진, 후진, 중심 차분

```python
import numpy as np

def forward_diff(f, x, h):
    """First derivative using forward difference -- O(h)"""
    return (f(x + h) - f(x)) / h

def backward_diff(f, x, h):
    """First derivative using backward difference -- O(h)"""
    return (f(x) - f(x - h)) / h

def central_diff(f, x, h):
    """First derivative using central difference -- O(h^2)"""
    return (f(x + h) - f(x - h)) / (2 * h)

def central_diff_2nd(f, x, h):
    """Second derivative using central difference -- O(h^2)"""
    return (f(x + h) - 2 * f(x) + f(x - h)) / h**2

def forward_diff_2nd_order(f, x, h):
    """First derivative using 2nd-order forward difference -- O(h^2)"""
    return (-f(x + 2*h) + 4*f(x + h) - 3*f(x)) / (2 * h)
```

### 10.2 미분을 위한 리처드슨 외삽법

```python
def richardson_diff(f, x, h, levels=3):
    """
    Richardson extrapolation for central difference.
    Returns a triangular table of increasingly accurate estimates.
    """
    # Build step sizes: h, 2h, 4h, ...
    D = np.zeros((levels, levels))

    # First column: central differences with step sizes h, 2h, 4h, ...
    for i in range(levels):
        hi = h * (2 ** i)
        D[i, 0] = (f(x + hi) - f(x - hi)) / (2 * hi)

    # Richardson extrapolation (note: entries go from fine to coarse)
    # We rearrange so D[0,0] uses smallest h (most accurate base)
    # Rebuild with D[i,0] using step size h * 2^(levels-1-i)
    R = np.zeros((levels, levels))
    for i in range(levels):
        hi = h * (2 ** (levels - 1 - i))
        R[i, 0] = (f(x + hi) - f(x - hi)) / (2 * hi)

    for j in range(1, levels):
        for i in range(j, levels):
            R[i, j] = R[i, j-1] + (R[i, j-1] - R[i-1, j-1]) / (4**j - 1)

    return R

# --- Example ---
f = lambda x: np.exp(x)   # f(x) = e^x, true derivative = e^x
x0 = 1.0
h = 0.1
true_val = np.exp(x0)

print(f"True f'({x0}) = {true_val:.10f}")
print(f"Forward diff:           {forward_diff(f, x0, h):.10f}  "
      f"Error: {abs(forward_diff(f, x0, h) - true_val):.2e}")
print(f"Central diff:           {central_diff(f, x0, h):.10f}  "
      f"Error: {abs(central_diff(f, x0, h) - true_val):.2e}")
print(f"2nd-order forward diff: {forward_diff_2nd_order(f, x0, h):.10f}  "
      f"Error: {abs(forward_diff_2nd_order(f, x0, h) - true_val):.2e}")

R = richardson_diff(f, x0, h, levels=3)
print(f"Richardson (level 1):   {R[1,1]:.10f}  "
      f"Error: {abs(R[1,1] - true_val):.2e}")
print(f"Richardson (level 2):   {R[2,2]:.10f}  "
      f"Error: {abs(R[2,2] - true_val):.2e}")
```

이 예제를 실행하면 참값 $`e \approx 2.7182818285`$에 대한 오차는 전진 차분 $`1.41 \times 10^{-1}`$, 중심 차분 $`4.53 \times 10^{-3}`$, 2차 전진 차분 $`9.77 \times 10^{-3}`$, 리처드슨 1수준 `R[1,1]` $`1.46 \times 10^{-4}`$, 2수준 `R[2,2]` $`3.46 \times 10^{-8}`$이다. `R[1,1]`은 간격 $`4h`$와 $`2h`$의 추정값을 결합하므로 간격 $`2h`$와 $`h`$를 결합하는 `R[2,1]`보다 오차가 크다.

### 10.3 비등간격 데이터의 도함수 (라그랑주)

```python
def lagrange_deriv(x_pts, y_pts, x_eval):
    """
    Compute the derivative at x_eval using the derivative
    of the Lagrange interpolating polynomial.

    Parameters:
        x_pts: array of x data points
        y_pts: array of y data points (same length)
        x_eval: point at which to evaluate the derivative

    Returns:
        Approximate value of f'(x_eval)
    """
    n = len(x_pts)
    deriv = 0.0

    for i in range(n):
        # Compute N_i'(x_eval)
        Ni_prime = 0.0
        for j in range(n):
            if j == i:
                continue
            # Product of all (x_eval - x_k) / (x_i - x_k) for k != i, k != j
            prod = 1.0
            for k in range(n):
                if k == i or k == j:
                    continue
                prod *= (x_eval - x_pts[k]) / (x_pts[i] - x_pts[k])
            Ni_prime += prod / (x_pts[i] - x_pts[j])
        deriv += y_pts[i] * Ni_prime

    return deriv

# --- Example: unequally spaced data ---
x_data = np.array([1.0, 2.5, 4.0])
y_data = np.exp(x_data)  # f(x) = e^x
x_eval = 2.5

approx = lagrange_deriv(x_data, y_data, x_eval)
true = np.exp(x_eval)
print(f"\nUnequally spaced data derivative at x={x_eval}:")
print(f"  Lagrange approx: {approx:.8f}")
print(f"  True value:      {true:.8f}")
print(f"  Error:           {abs(approx - true):.2e}")
```

> **참고:** **예제 데이터의 간격.** 이 예제의 점 1.0, 2.5, 4.0은 실제로는 간격 1.5의 등간격이다. 따라서 결과 17.29329는 $`h = 1.5`$인 중심 차분 $`(e^4 - e^1)/3`$과 같고, 간격이 커서 참값 12.18249와 5.11만큼 차이가 난다. 함수 `lagrange_deriv` 자체는 비등간격 점에도 그대로 동작한다.

### 10.4 수치 편도함수

```python
def partial_x_central(f, x, y, dx):
    """Central difference for df/dx"""
    return (f(x + dx, y) - f(x - dx, y)) / (2 * dx)

def partial_y_central(f, x, y, dy):
    """Central difference for df/dy"""
    return (f(x, y + dy) - f(x, y - dy)) / (2 * dy)

def mixed_partial_xy(f, x, y, dx, dy):
    """Central difference for d^2f/(dx dy)"""
    return (f(x+dx, y+dy) - f(x+dx, y-dy)
          - f(x-dx, y+dy) + f(x-dx, y-dy)) / (4 * dx * dy)

# --- Example ---
g = lambda x, y: x**2 * y + np.sin(x * y)  # test function
x0, y0 = 1.0, 2.0
h = 1e-5

print(f"\nPartial derivatives of f(x,y) = x^2*y + sin(x*y) at ({x0},{y0}):")
print(f"  df/dx (central):      {partial_x_central(g, x0, y0, h):.8f}")
print(f"  df/dy (central):      {partial_y_central(g, x0, y0, h):.8f}")
print(f"  d^2f/(dx dy) (mixed): {mixed_partial_xy(g, x0, y0, h, h):.8f}")
```

실행 결과의 $`\partial f/\partial x`$와 $`\partial f/\partial y`$는 해석적 값 $`2xy + y\cos(xy) = 4 + 2\cos 2 \approx 3.16770633`$, $`x^2 + x\cos(xy) = 1 + \cos 2 \approx 0.58385316`$과 소수점 아래 여덟째 자리까지 일치한다. 혼합 편도함수의 결과 $`-0.23474334`$는 해석적 값 $`2x + \cos(xy) - xy\sin(xy) = 2 + \cos 2 - 2\sin 2 \approx -0.23474169`$와 약 $`1.7 \times 10^{-6}`$ 차이가 난다. 분모 $`4h^2 = 4 \times 10^{-10}`$이 매우 작아 반올림 오차가 커지기 때문이다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:-----|:-----|
| 전진 차분(Forward Difference) ($`f'`$) | $`f'(x_i) \approx \dfrac{f_{i+1} - f_i}{h}`$이며, 점 $`x_i, x_{i+1}`$을 사용하는 $`O(h)`$ 공식이다. |
| 후진 차분(Backward Difference) ($`f'`$) | $`f'(x_i) \approx \dfrac{f_i - f_{i-1}}{h}`$이며, 점 $`x_{i-1}, x_i`$를 사용하는 $`O(h)`$ 공식이다. |
| 중심 차분(Central Difference) ($`f'`$) | $`f'(x_i) \approx \dfrac{f_{i+1} - f_{i-1}}{2h}`$이며, 점 $`x_{i-1}, x_{i+1}`$을 사용하는 $`O(h^2)`$ 공식으로, 대칭성에 의한 오차 상쇄 때문에 같은 간격 크기에서 전진 차분과 후진 차분보다 정확하다. |
| 2차 전진(2nd-Order Forward) ($`f'`$) | $`f'(x_i) \approx \dfrac{-f_{i+2} + 4f_{i+1} - 3f_i}{2h}`$이며, 점 $`x_i, x_{i+1}, x_{i+2}`$를 사용하는 $`O(h^2)`$ 공식이다. |
| 중심 차분(Central Difference) ($`f''`$) | $`f''(x_i) \approx \dfrac{f_{i+1} - 2f_i + f_{i-1}}{h^2}`$이며, 점 $`x_{i-1}, x_i, x_{i+1}`$을 사용하는 $`O(h^2)`$ 공식이다. |
| 전진 차분(Forward Difference) ($`f''`$) | $`f''(x_i) \approx \dfrac{f_{i+2} - 2f_{i+1} + f_i}{h^2}`$이며, 점 $`x_i, x_{i+1}, x_{i+2}`$를 사용하는 $`O(h)`$ 공식이다. |
| 테일러 급수 전개 | 등간격 유한 차분 공식은 여러 점에서의 테일러 급수 전개를 조합하여 유도한다. |
| 리처드슨 외삽 (1수준)(Richardson Extrap. 1 level) | $`\dfrac{4}{3}D(h) - \dfrac{1}{3}D(2h)`$는 중심 차분 $`D(h)`$의 $`h^2`$ 오차를 제거하여 $`O(h^4)`$ 정확도를 얻는다. |
| 리처드슨 외삽 (2수준)(Richardson Extrap. 2 levels) | $`\dfrac{16}{15}\bar{D}(h) - \dfrac{1}{15}\bar{D}(2h)`$는 $`\bar{D}(h)`$의 $`h^4`$ 오차를 제거하여 $`O(h^6)`$ 정확도를 얻으며, 이 과정을 반복하면 정확도가 $`O(h^2) \to O(h^4) \to O(h^6) \to \cdots`$로 향상된다. |
| 비등간격 데이터 | 라그랑주 보간 다항식을 미분하여 도함수를 구하며, 점이 등간격이면 중심 차분 공식과 같아진다. |
| 혼합 편도함수(Mixed Partial) ($`f_{xy}`$) | $`\dfrac{f(x\!+\!\Delta x, y\!+\!\Delta y) - f(x\!+\!\Delta x, y\!-\!\Delta y) - f(x\!-\!\Delta x, y\!+\!\Delta y) + f(x\!-\!\Delta x, y\!-\!\Delta y)}{4\Delta x \Delta y}`$로 근사하며, 대각선 4점을 사용하는 $`O(\Delta x^2, \Delta y^2)`$ 공식이다. |
| 미정 계수법 | 테일러 전개의 계수를 맞추는 선형 연립방정식을 풀어 임의의 유한 차분 스텐실을 체계적으로 유도한다. |
| 미분 가능 프로그래밍 | PyTorch, JAX, TensorFlow 같은 도구는 자동 미분으로 정확한 기울기를 제공하며, 수치 미분의 현대적 대안이다. |

---

<br>

## 점검 문제

1. **전진 차분과 중심 차분:** $`f(x) = e^x`$, $`x = 1`$, $`h = 0.1`$일 때 전진 차분과 중심 차분으로 $`f'(1)`$을 근사하고, 참값 $`e \approx 2.718282`$에 대한 오차를 각각 구하라.

   <details>
   <summary>정답</summary>

   $`e^{1.1} \approx 3.004166`$, $`e^{0.9} \approx 2.459603`$이다. 전진 차분은 $`(3.004166 - 2.718282)/0.1 \approx 2.858842`$이고 오차는 약 0.140560이다. 중심 차분은 $`(3.004166 - 2.459603)/0.2 \approx 2.722815`$이고 오차는 약 0.004533이다. 중심 차분의 오차가 약 31배 작으며, 이는 $`O(h)`$와 $`O(h^2)`$ 공식의 차이를 보여 준다.

   </details>

2. **리처드슨 외삽:** 1번의 함수에서 중심 차분 값이 $`D(0.2) \approx 2.736440`$, $`D(0.1) \approx 2.722815`$일 때 1수준 리처드슨 외삽값과 그 오차를 구하라.

   <details>
   <summary>정답</summary>

   $`\bar{D}(0.1) = \frac{4}{3}(2.722815) - \frac{1}{3}(2.736440) \approx 3.630420 - 0.912147 = 2.718273`$이다. 참값과의 오차는 약 $`9.1 \times 10^{-6}`$이며, $`D(0.1)`$의 오차 0.004533보다 약 500배 작다.

   </details>

3. **비등간격 데이터:** 데이터 $`(0, 1)`$, $`(1, 3)`$, $`(3, 2)`$에 2차 라그랑주 다항식을 맞추어 $`x = 1`$에서의 도함수를 구하라.

   <details>
   <summary>정답</summary>

   6.2절의 식에 $`x = 1`$을 대입하면 $`N_0'(1) = \frac{2 - 1 - 3}{(0 - 1)(0 - 3)} = -\frac{2}{3}`$, $`N_1'(1) = \frac{2 - 0 - 3}{(1 - 0)(1 - 3)} = \frac{1}{2}`$, $`N_2'(1) = \frac{2 - 0 - 1}{(3 - 0)(3 - 1)} = \frac{1}{6}`$이다. 따라서 $`f'(1) = 1 \cdot \left(-\frac{2}{3}\right) + 3 \cdot \frac{1}{2} + 2 \cdot \frac{1}{6} = \frac{7}{6} \approx 1.1667`$이다.

   </details>

4. **미정 계수법:** $`f'(x) \approx a\,f(x) + b\,f(x+h) + c\,f(x+2h)`$가 2차 정확도를 갖도록 $`a`$, $`b`$, $`c`$를 구하라.

   <details>
   <summary>정답</summary>

   테일러 전개를 대입하고 계수를 맞추면 $`a + b + c = 0`$, $`(b + 2c)h = 1`$, $`(b + 4c)\frac{h^2}{2} = 0`$이다. 셋째 식에서 $`b = -4c`$이고, 둘째 식에 대입하면 $`-2ch = 1`$이므로 $`c = -\frac{1}{2h}`$, $`b = \frac{2}{h}`$, $`a = -\frac{3}{2h}`$이다. 결과 $`f'(x) \approx \frac{-3f(x) + 4f(x+h) - f(x+2h)}{2h}`$는 4절의 2차 전진 차분 공식과 같다.

   </details>

5. **간격 크기와 오차:** 간격 크기 $`h`$를 절반으로 줄이면 전진 차분, 중심 차분, 1수준 리처드슨 외삽값의 절단 오차는 대략 몇 배가 되는가? 또 $`h`$를 계속 줄일 때 주의할 점은 무엇인가?

   <details>
   <summary>정답</summary>

   절단 오차가 각각 $`O(h)`$, $`O(h^2)`$, $`O(h^4)`$이므로 대략 $`\frac{1}{2}`$, $`\frac{1}{4}`$, $`\frac{1}{16}`$배가 된다. 그러나 $`h`$가 너무 작으면 $`f(x+h) - f(x-h)`$처럼 가까운 두 값의 차를 작은 수로 나누는 과정에서 반올림 오차가 커지므로, 전체 오차가 오히려 증가할 수 있다.

   </details>
