# 제21장 강의 — 수치 미분(Numerical Differentiation)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 21

> **선수 지식**: [미적분학] 미분, 테일러 급수 (제1-20장).
>
> **학습 목표**:
> 1. 미분을 위해 유한 차분 공식을 적용할 수 있다
> 2. 수치 미분에서 절삭 오차를 분석할 수 있다
> 3. 정확도 향상을 위해 리처드슨 외삽법을 구현할 수 있다

---

<br>

## 목차

- [1. 수치 미분 개요](#1-수치-미분-개요)
  - [1.1 극한으로서의 도함수](#11-극한으로서의-도함수)
  - [1.2 이계도함수](#12-이계도함수)
  - [1.3 편도함수](#13-편도함수)
  - [1.4 물리적 응용: 열유속](#14-물리적-응용-열유속)
- [2. 1계도함수의 유한 차분 근사](#2-1계도함수의-유한-차분-근사)
  - [2.1 전진 차분 (1차)](#21-전진-차분-1차)
  - [2.2 후진 차분 (1차)](#22-후진-차분-1차)
  - [2.3 중심 차분 (2차)](#23-중심-차분-2차)
- [3. 2계도함수의 유한 차분 근사](#3-2계도함수의-유한-차분-근사)
  - [3.1 테일러 급수를 이용한 중심 차분](#31-테일러-급수를-이용한-중심-차분)
  - [3.2 2계도함수의 전진 차분](#32-2계도함수의-전진-차분)
- [4. 정확도 향상: 고차 전진 차분](#4-정확도-향상-고차-전진-차분)
- [5. 미분을 위한 리처드슨 외삽법](#5-미분을-위한-리처드슨-외삽법)
  - [5.1 기본 아이디어](#51-기본-아이디어)
  - [5.2 중심 차분의 오차 구조 재검토](#52-중심-차분의-오차-구조-재검토)
  - [5.3 첫 번째 외삽 수준](#53-첫-번째-외삽-수준)
  - [5.4 두 번째 외삽 수준](#54-두-번째-외삽-수준)
- [6. 비등간격 데이터의 도함수](#6-비등간격-데이터의-도함수)
  - [6.1 라그랑주 다항식 접근법](#61-라그랑주-다항식-접근법)
  - [6.2 라그랑주 기저 함수의 도함수](#62-라그랑주-기저-함수의-도함수)
- [7. 수치 편도함수](#7-수치-편도함수)
  - [7.1 1차 편도함수 (중심 차분)](#71-1차-편도함수-중심-차분)
  - [7.2 혼합 편도함수](#72-혼합-편도함수)
- [8. 미정 계수법](#8-미정-계수법)
  - [8.1 문제 설정](#81-문제-설정)
  - [8.2 테일러 전개와 매칭](#82-테일러-전개와-매칭)
  - [8.3 계수 풀기](#83-계수-풀기)
  - [8.4 절단 오차](#84-절단-오차)
- [9. 미분 가능 프로그래밍 (간략 개요)](#9-미분-가능-프로그래밍-간략-개요)
- [10. Python 구현](#10-python-구현)
- [요약](#요약)

---

<br>

## 1. 수치 미분 개요

수치 미분(Numerical Differentiation)은 해석적으로 계산하는 대신 이산 데이터나 함수 평가로부터 도함수를 추정하는 과정이다. 이는 다음과 같은 경우에 필수적이다:

- 함수가 이산 점에서만 알려진 경우 (실험/현장 데이터)
- 해석적 도함수가 너무 복잡하거나 구할 수 없는 경우
- 더 큰 수치 알고리즘(예: ODE/PDE 풀기) 내에서 도함수가 필요한 경우

### 1.1 극한으로서의 도함수

유한 차분비(Finite Difference Ratio)는 자연스러운 출발점을 제공한다:

$$\frac{\Delta y}{\Delta x} = \frac{f(x + \Delta x) - f(x)}{\Delta x}$$

$\Delta x \to 0$으로 극한을 취하면:

$$\frac{dy}{dx} = \lim_{\Delta x \to 0} \frac{f(x + \Delta x) - f(x)}{\Delta x} =: f'(x)$$

이것이 **$x$에 대한 $y$의 1계도함수** 이다. 기하학적으로 도함수는 점 $x_i$에서 곡선에 대한 접선의 기울기이다.

> **[미적분]** 강의의 세 그림은 $\Delta x$가 줄어들면서 할선(곡선 위의 두 점을 잇는 선)이 $x_i$에서의 접선에 어떻게 수렴하는지, 그리고 할선의 기울기가 $f'(x_i)$에 어떻게 수렴하는지를 보여준다.

### 1.2 이계도함수

이계도함수는 일계도함수의 도함수이다:

$$\frac{d^2 y}{dx^2} = \frac{d}{dx}\left(\frac{dy}{dx}\right)$$

이것은 **기울기가 얼마나 빠르게 변하는지**, 즉 함수의 곡률이나 오목성을 측정한다.

### 1.3 편도함수

두 변수의 함수 $f(x, y)$에 대해:

$$\frac{\partial f}{\partial x} = \lim_{\Delta x \to 0} \frac{f(x + \Delta x, y) - f(x, y)}{\Delta x}$$

$$\frac{\partial f}{\partial y} = \lim_{\Delta y \to 0} \frac{f(x, y + \Delta y) - f(x, y)}{\Delta y}$$

이것들은 각각 $x$와 $y$에 대한 $f$의 **편도함수** 이다. 각 편도함수는 다른 변수를 고정한 채 $f$의 변화율을 측정한다.

### 1.4 물리적 응용: 열유속(Heat Flux)

대표적인 공학적 응용은 푸리에의 열전도 법칙이다. 온도 분포 $T(x)$가 주어지면, 열유속은:

$$\vec{q} = -k \frac{\partial T}{\partial x}$$

여기서 $k$는 열전도율이다. 음의 부호는 열이 고온에서 저온으로 흐름을 나타낸다:

- **양의 기울기** ($\partial T / \partial x > 0$)는 **음의 열유동** (열이 $-x$ 방향으로 흐름)을 유발
- **음의 기울기** ($\partial T / \partial x < 0$)는 **양의 열유동** (열이 $+x$ 방향으로 흐름)을 유발

> **[미적분]** 열유동의 방향은 항상 온도 기울기의 반대이다. 이것은 도함수 개념의 직접적인 물리적 응용이다.

---

<br>

## 2. 1계도함수의 유한 차분 근사(Finite Difference Approximation)

모든 유한 차분 공식은 테일러 급수 전개(Taylor Series Expansion)로부터 유도된다. $h$를 등간격 점들 사이의 간격 크기로 놓자: $x_{i+1} = x_i + h$.

### 2.1 전진 차분(Forward Difference) (1차)

$x_i$에서의 전진 방향 **테일러 전개**:

$$f(x_{i+1}) = f(x_i) + f'(x_i)\,h + f''(x_i)\frac{h^2}{2} + f'''(x_i)\frac{h^3}{6} + \cdots \quad (\ast)$$

식 $(\ast)$에서 $f'(x_i)$를 풀면:

$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - f''(x_i)\frac{h}{2} + O(h^2)$$

절단 오차항을 버리면 **전진 차분** 근사를 얻는다:

$$\boxed{f'(x_i) \approx \frac{f(x_{i+1}) - f(x_i)}{h}}$$

- **절단 오차**: $O(h)$ -- 1차 정확도
- **주도적 오차항**: $-f''(x_i)\dfrac{h}{2}$
- 점 $x_i$와 $x_{i+1}$을 사용 (해당 점과 앞의 한 점)

### 2.2 후진 차분(Backward Difference) (1차)

마찬가지로, $x_i$에서 $f(x_{i-1})$을 전개하면:

$$f(x_{i-1}) = f(x_i) - f'(x_i)\,h + f''(x_i)\frac{h^2}{2} - f'''(x_i)\frac{h^3}{6} + \cdots \quad (\ast\ast)$$

$f'(x_i)$를 풀면:

$$\boxed{f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{h}}$$

- **절단 오차**: $O(h)$ -- 1차 정확도
- 점 $x_{i-1}$과 $x_i$를 사용

### 2.3 중심 차분(Central Difference) (2차)

$(\ast)$에서 $(\ast\ast)$를 빼면:

$$f(x_{i+1}) - f(x_{i-1}) = 2f'(x_i)\,h + 2f'''(x_i)\frac{h^3}{6} + \cdots$$

$f'(x_i)$를 풀면:

$$\boxed{f'(x_i) \approx \frac{f(x_{i+1}) - f(x_{i-1})}{2h}}$$

- **절단 오차**: $O(h^2)$ -- 2차 정확도
- **주도적 오차항**: $-f'''(x_i)\dfrac{h^2}{6}$
- 점 $x_{i-1}$과 $x_{i+1}$을 사용 ($x_i$에 대해 대칭)

> **[미적분]** 중심 차분은 대칭성으로 인해 홀수 차수의 오차항이 상쇄되기 때문에 전진/후진 차분보다 더 높은 정확도($O(h)$ 대 $O(h^2)$)를 달성한다. 이것은 양쪽 이웃 점이 모두 사용 가능할 때의 핵심적인 장점이다.

---

<br>

## 3. 2계도함수의 유한 차분 근사

### 3.1 테일러 급수를 이용한 중심 차분

$(\ast)$와 $(\ast\ast)$를 더하면:

$$f(x_{i+1}) + f(x_{i-1}) = 2f(x_i) + f''(x_i)\,h^2 + O(h^4)$$

$f''(x_i)$를 풀면:

$$\boxed{f''(x_i) = \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2} + O(h^2)}$$

이것이 **이계도함수의 중심 차분 근사** 이다:

- **절단 오차**: $O(h^2)$ -- 2차 정확도
- 등간격의 세 점 사용: $x_{i-1}$, $x_i$, $x_{i+1}$

> **[미적분]** 이 공식은 아름다운 물리적 해석을 가진다: $x_i$에서의 함수값이 양쪽 이웃의 평균으로부터 얼마나 벗어나는지를 측정한다. 큰 양의 $f''$은 함수가 위로 볼록(오목 위)함을 의미하고, 큰 음의 $f''$은 아래로 볼록(오목 아래)함을 의미한다.

### 3.2 2계도함수의 전진 차분

전진 점만 사용하여 $f''(x_i)$의 **전진** 공식을 유도하기 위해, $f(x_{i+2})$를 전개한다:

$$f(x_{i+2}) = f(x_i) + f'(x_i)\,(2h) + f''(x_i)\frac{(2h)^2}{2} + f'''(x_i)\frac{(2h)^3}{6} + O(h^4) \quad (\ast\ast\ast)$$

$(\ast\ast\ast) - 2(\ast)$를 계산하면:

$$f(x_{i+2}) - 2f(x_{i+1}) = -f(x_i) + f''(x_i)\,h^2 + O(h^3)$$

$f''(x_i)$를 풀면:

$$\boxed{f''(x_i) = \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{h^2} + O(h)}$$

- **절단 오차**: $O(h)$ -- 1차 정확도 (중심 공식보다 덜 정확)
- 세 개의 전진 점 사용: $x_i$, $x_{i+1}$, $x_{i+2}$

---

<br>

## 4. 정확도 향상: 고차 전진 차분

이계도함수 보정항을 포함하여 $f'(x_i)$의 전진 차분을 개선할 수 있다. 다음을 상기하자:

$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \underbrace{f''(x_i)}_{\text{FD 근사 대입}}\frac{h}{2} + O(h^2)$$

$f''(x_i)$의 전진 차분 근사:

$$f''(x_i) \approx \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{h^2}$$

를 오차항에 대입하면:

$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} - \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{2h} + O(h^2)$$

결합하면:

$$\boxed{f'(x_i) = \frac{-f(x_{i+2}) + 4f(x_{i+1}) - 3f(x_i)}{2h} + O(h^2)}$$

이것이 $f'(x_i)$에 대한 **2차 전진 차분** 공식이다.

- 추가 점 하나를 포함하여 **정확도가** $O(h)$에서 $O(h^2)$로 **향상**
- 확산(이계도함수) 항을 추가함으로써 정확도를 높인다

> **[미적분]** 오차의 주도항에 대해 수치 근사를 대입하는 이 아이디어는 정확도를 높이기 위한 일반적인 기법이다. 이것은 리처드슨 외삽법(다음에 논의)과 밀접한 관련이 있다.

---

<br>

## 5. 미분을 위한 리처드슨 외삽법(Richardson Extrapolation)

### 5.1 기본 아이디어

수치 적분에서 리처드슨 외삽법이 서로 다른 간격 크기로 계산한 두 추정값을 사용하여 주도적 오차항을 제거했음을 상기하자. 적분의 경우:

$$I = I(h_2) + \frac{1}{(h_1/h_2)^2 - 1}\bigl(I(h_2) - I(h_1)\bigr)$$

$h_1/h_2 = 2$일 때:

$$I = \frac{4}{3}I(h_2) - \frac{1}{3}I(h_1)$$

마찬가지로, 미분의 경우:

$$\boxed{D = \frac{4}{3}D(h_2) - \frac{1}{3}D(h_1)}$$

여기서 $D(h_1)$과 $D(h_2)$는 간격 크기 $h_1$ (더 큰 것)과 $h_2 = h_1/2$ (더 작은 것)으로 계산한 도함수 추정값이다.

### 5.2 중심 차분의 오차 구조 재검토

테일러 전개에서 시작한다:

$$f(x+h) = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}h^k$$

$$f(x-h) = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}(-h)^k = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!}(-1)^k h^k$$

$f(x+h) - f(x-h)$을 계산하면:

$$f(x+h) - f(x-h) = 2h\,f'(x) + 2f'''(x)\frac{h^3}{3!} + 2f^{(5)}(x)\frac{h^5}{5!} + \cdots$$

$2h$로 나누면:

$$f'(x) = \frac{f(x+h) - f(x-h)}{2h} - f'''(x)\frac{h^2}{3!} - f^{(5)}(x)\frac{h^4}{5!} - \cdots$$

중심 차분 연산자와 오차 계수를 정의한다:

$$D(h) = \frac{f(x+h) - f(x-h)}{2h}$$

$$e_2 = -\frac{f^{(3)}(x)}{3!}, \qquad e_4 = -\frac{f^{(5)}(x)}{5!}$$

그러면:

$$f'(x) = D(h) + e_2\,h^2 + e_4\,h^4 + \cdots$$

여기서 오차 $E(h) = e_2\,h^2 + e_4\,h^4 + \cdots$이다.

> **[미적분]** 중요한 관찰은 오차 계수 $e_2, e_4, \ldots$가 $h$에 **독립** 이라는 것이다 -- 함수와 평가점 $x$에만 의존한다. 이것이 바로 리처드슨 외삽법을 가능하게 하는 성질이다.

### 5.3 첫 번째 외삽 수준

**목표**: $h^2$ 항을 제거한다.

간격 크기 $h$와 $2h$에서의 근사를 작성한다:

$$f'(x) = D(h) + e_2\,h^2 + e_4\,h^4 + \cdots \quad \text{(4를 곱함)}$$

$$f'(x) = D(2h) + e_2\,(2h)^2 + e_4\,(2h)^4 + \cdots$$

$4 \times (\text{첫 번째}) - (\text{두 번째})$를 계산하면:

$$3f'(x) = 4D(h) - D(2h) - 12h^4 e_4 + O(h^6)$$

$$\boxed{f'(x) = \underbrace{\frac{4}{3}D(h) - \frac{1}{3}D(2h)}_{=:\,\bar{D}(h)} - 4h^4 e_4 + O(h^6)}$$

결합된 추정값 $\bar{D}(h)$는 이제 **$O(h^4)$** 정확도를 가진다 -- $h^2$ 오차가 제거되었다.

### 5.4 두 번째 외삽 수준

같은 과정을 계속하여 $h^4$ 항을 제거한다. 간격 크기 $h$와 $2h$에서:

$$f'(x) = \bar{D}(h) - 4h^4 e_4 + O(h^6)$$

$$f'(x) = \bar{D}(2h) - 4(2h)^4 e_4 + O(h^6)$$

$16 \times (\text{첫 번째}) - (\text{두 번째})$를 계산하면:

$$15f'(x) = 16\bar{D}(h) - \bar{D}(2h) + O(h^6)$$

$$\boxed{f'(x) = \frac{16}{15}\bar{D}(h) - \frac{1}{15}\bar{D}(2h) + O(h^6)}$$

이제 **$O(h^6)$** 정확도를 가진다. 이 과정은 무한히 계속할 수 있으며, 각 수준에서 $h$의 다음 짝수 거듭제곱을 제거한다.

> **[미적분]** 리처드슨 외삽법은 수렴을 가속화하는 체계적인 방법이다. 각 수준에서 서로 다른 간격 크기의 추정값을 결합하여 주도적 오차항을 상쇄한다. 중심 차분(오차가 $h$의 짝수 거듭제곱)에 대한 일반 패턴은 $n$번째 수준에서 비율 $4^n$을 사용한다: 가중치는 $\frac{4^n}{4^n - 1}$과 $\frac{-1}{4^n - 1}$이다.

---

<br>

## 6. 비등간격 데이터의 도함수(Derivatives of Unequally Spaced Data)

실제로 실험이나 현장 조사에서 수집된 데이터는 흔히 **비등간격** 이다. 위에서 유도한 유한 차분 공식은 등간격 점을 가정하므로 직접 적용할 수 없다.

### 6.1 라그랑주 다항식 접근법

비등간격 데이터를 다루는 한 가지 방법은 **라그랑주 다항식(Lagrange Polynomial)** 을 맞추고 미분하는 것이다.

세 데이터 점 $(x_0, y_0)$, $(x_1, y_1)$, $(x_2, y_2)$가 주어지면, 2차 라그랑주 다항식을 구성한다:

$$f(x) = y_0\,N_0(x) + y_1\,N_1(x) + y_2\,N_2(x)$$

여기서 $N_k(x)$는 라그랑주 기저 함수이다. 미분하면:

$$f'(x) = y_0\,N_0'(x) + y_1\,N_1'(x) + y_2\,N_2'(x)$$

### 6.2 라그랑주 기저 함수의 도함수

세 점에 대해, 기저 함수와 그 도함수는:

**기저 함수 $N_0(x)$**:

$$N_0(x) = \frac{(x - x_1)(x - x_2)}{(x_0 - x_1)(x_0 - x_2)} = \frac{x^2 - (x_1 + x_2)x + x_1 x_2}{(x_0 - x_1)(x_0 - x_2)}$$

$$N_0'(x) = \frac{2x - x_1 - x_2}{(x_0 - x_1)(x_0 - x_2)}$$

**기저 함수 $N_1(x)$**:

$$N_1(x) = \frac{(x - x_0)(x - x_2)}{(x_1 - x_0)(x_1 - x_2)} = \frac{x^2 - (x_0 + x_2)x + x_0 x_2}{(x_1 - x_0)(x_1 - x_2)}$$

$$N_1'(x) = \frac{2x - x_0 - x_2}{(x_1 - x_0)(x_1 - x_2)}$$

**기저 함수 $N_2(x)$**:

$$N_2(x) = \frac{(x - x_0)(x - x_1)}{(x_2 - x_0)(x_2 - x_1)} = \frac{x^2 - (x_0 + x_1)x + x_0 x_1}{(x_2 - x_0)(x_2 - x_1)}$$

$$N_2'(x) = \frac{2x - x_0 - x_1}{(x_2 - x_0)(x_2 - x_1)}$$

임의의 점 $x$에서의 도함수는:

$$f'(x) = y_0 \cdot \frac{2x - x_1 - x_2}{(x_0 - x_1)(x_0 - x_2)} + y_1 \cdot \frac{2x - x_0 - x_2}{(x_1 - x_0)(x_1 - x_2)} + y_2 \cdot \frac{2x - x_0 - x_1}{(x_2 - x_0)(x_2 - x_1)}$$

> **[미적분]** 이 접근법은 데이터 점의 간격에 관계없이 작동한다. 점들이 등간격인 경우(즉, $x_1 - x_0 = x_2 - x_1 = h$), 위 공식에 $x = x_1$을 대입하면 표준 중심 차분 공식 $f'(x_1) \approx (y_2 - y_0)/(2h)$를 복원한다.

---

<br>

## 7. 수치 편도함수(Numerical Partial Derivatives)

### 7.1 1차 편도함수 (중심 차분)

함수 $f(x, y)$에 대해, 편도함수의 중심 유한 차분 근사는:

$$\boxed{\frac{\partial f}{\partial x} \approx \frac{f(x + \Delta x,\, y) - f(x - \Delta x,\, y)}{2\,\Delta x}}$$

$$\boxed{\frac{\partial f}{\partial y} \approx \frac{f(x,\, y + \Delta y) - f(x,\, y - \Delta y)}{2\,\Delta y}}$$

각각 $O(\Delta x^2)$과 $O(\Delta y^2)$ 정확도를 가진다.

### 7.2 혼합 편도함수(Mixed Partial Derivative)

혼합 2차 편도함수 $\dfrac{\partial^2 f}{\partial x \,\partial y}$는 중심 차분을 두 번 적용하여 계산한다:

$$\frac{\partial^2 f}{\partial x \,\partial y} = \frac{\partial}{\partial x}\left(\frac{\partial f}{\partial y}\right)$$

$$\approx \frac{1}{2\,\Delta x}\left(\left.\frac{\partial f}{\partial y}\right|_{x+\Delta x,\,y} - \left.\frac{\partial f}{\partial y}\right|_{x-\Delta x,\,y}\right)$$

각 점에서 $\partial f / \partial y$에 중심 차분을 적용하면:

$$\approx \frac{1}{2\,\Delta x}\left(\frac{f(x+\Delta x,\, y+\Delta y) - f(x+\Delta x,\, y-\Delta y)}{2\,\Delta y} - \frac{f(x-\Delta x,\, y+\Delta y) - f(x-\Delta x,\, y-\Delta y)}{2\,\Delta y}\right)$$

정리하면:

$$\boxed{\frac{\partial^2 f}{\partial x \,\partial y} \approx \frac{f(x+\Delta x,\, y+\Delta y) - f(x+\Delta x,\, y-\Delta y) - f(x-\Delta x,\, y+\Delta y) + f(x-\Delta x,\, y-\Delta y)}{4\,\Delta x\,\Delta y}}$$

이는 $(x, y)$ 주위의 "대각선" 점에서 네 번의 함수 평가를 필요로 한다.

---

<br>

## 8. 미정 계수법(Method of Undetermined Coefficients)

이것은 지정된 점 집합을 사용하여 도함수를 근사하기 위한 유한 차분 공식을 유도하는 체계적인 기법이다.

### 8.1 문제 설정

등간격 세 점 $f(x-h)$, $f(x)$, $f(x+h)$를 기반으로 이계도함수 $f''(x)$의 근사를 구한다:

$$f''(x) \approx a\,f(x-h) + b\,f(x) + c\,f(x+h)$$

**질문**: $a$, $b$, $c$를 어떻게 결정하는가?

### 8.2 테일러 전개와 매칭

$f(x-h)$와 $f(x+h)$를 테일러 급수로 전개한다(나머지항은 $[x-h, x+h]$에서의 $\xi_-$와 $\xi_+$를 포함):

$$f(x-h) = f(x) - f'(x)\,h + f''(x)\frac{h^2}{2} - f'''(x)\frac{h^3}{6} + f^{(4)}(\xi_-)\frac{h^4}{24}$$

$$f(x+h) = f(x) + f'(x)\,h + f''(x)\frac{h^2}{2} + f'''(x)\frac{h^3}{6} + f^{(4)}(\xi_+)\frac{h^4}{24}$$

근사식에 대입하고 도함수별로 정리하면:

$$a\,f(x-h) + b\,f(x) + c\,f(x+h) = (a + b + c)\,f(x) + (-a + c)\,h\,f'(x) + (a + c)\frac{h^2}{2}\,f''(x)$$
$$\quad + (c - a)\frac{h^3}{6}\,f'''(x) + \frac{h^4}{24}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)$$

### 8.3 계수 풀기

이 식이 $f''(x)$와 같으려면 다음이 필요하다:

| 조건 | 방정식 |
|:-----|:-------|
| $f(x)$의 계수 = 0 | $a + b + c = 0$ |
| $f'(x)$의 계수 = 0 | $-a + c = 0$ (즉, $c = a$) |
| $f''(x)$의 계수 = 1 | $(a + c)\dfrac{h^2}{2} = 1$ |

$c = a$와 $(a + c)\dfrac{h^2}{2} = 1$로부터:

$$2a \cdot \frac{h^2}{2} = 1 \implies a = \frac{1}{h^2}$$

따라서:

$$\boxed{a = \frac{1}{h^2}, \quad c = \frac{1}{h^2}, \quad b = -\frac{2}{h^2}}$$

### 8.4 절단 오차(Truncation Error)

다시 대입하면, $c - a = 0$이므로 삼계항이 소멸한다. 절단 오차는 사계 나머지항에서 온다:

$$f''(x) = \frac{1}{h^2}\bigl(f(x-h) - 2f(x) + f(x+h)\bigr) - \frac{h^4}{24}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)$$

**중간값 정리(Intermediate Value Theorem)** 에 의해, $f^{(4)}$가 연속이고 $f^{(4)}(\xi_-)$와 $f^{(4)}(\xi_+)$가 $f^{(4)}$의 두 값이므로, $[x-h, x+h]$에 어떤 $\xi$가 존재하여:

$$f^{(4)}(\xi) = \frac{1}{2}\bigl(f^{(4)}(\xi_-) + f^{(4)}(\xi_+)\bigr)$$

따라서:

$$\boxed{f''(x) = \frac{f(x-h) - 2f(x) + f(x+h)}{h^2} - \frac{h^2}{12}\,f^{(4)}(\xi)}$$

이는 $f''$의 중심 차분이 $O(h^2)$ 정확도이며, 명시적 절단 오차항이 $-\dfrac{h^2}{12}f^{(4)}(\xi)$임을 확인한다.

> **[미적분]** 미정 계수법은 매우 범용적이다. 포함할 점과 근사할 도함수를 선택함으로써 어떤 유한 차분 스텐실(stencil)이든 유도할 수 있다. 미지 계수에 대한 연립방정식은 항상 선형이므로 풀기가 간단하다.

---

<br>

## 9. 미분 가능 프로그래밍(Differentiable Programming) (간략 개요)

**미분 가능 프로그래밍** 은 프로그램의 모든 것이 **자동으로 미분 가능** 하도록 작성되는 현대적 프로그래밍 패러다임이다.

핵심 사항:

- 입력의 작은 변화가 출력에 어떻게 영향을 미치는지 계산할 수 있다
- **기울기 기반 최적화**(예: 신경망 학습)를 가능하게 한다
- 전통적인 수치 방법과 데이터 기반 방법을 결합하는 **새로운 하이브리드 접근법** 을 나타낸다
- 자동 미분(AD, Automatic Differentiation)을 사용하여 유한 차분보다 더 정확하고 기호적 미분보다 더 효율적인 **정확한 기울기** 를 제공한다

**미분 가능 프로그래밍을 위한 인기 있는 도구**:
- **PyTorch** (Python)
- **TensorFlow** (Python)
- **JAX** (Python, Google 개발)
- **Zygote** (Julia)

> **[미적분]** 자동 미분은 수치 미분과 근본적으로 다르다. 수치 미분이 유한 차분을 사용하여 도함수를 근사하고 (절단 오차와 반올림 오차를 겪는 반면), AD는 프로그램의 모든 연산에 연쇄 법칙(Chain Rule)을 체계적으로 적용하여 기계 정밀도로 도함수를 계산한다.

---

<br>

## 10. Python 구현

### 전진, 후진, 중심 차분

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

### 미분을 위한 리처드슨 외삽법

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

### 비등간격 데이터의 도함수 (라그랑주)

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

### 수치 편도함수

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

---

<br>

## 요약

| 방법 | 공식 | 정확도 | 사용 점 |
|:-----|:-----|:-------|:--------|
| **전진 차분(Forward Difference)** ($f'$) | $\dfrac{f_{i+1} - f_i}{h}$ | $O(h)$ | $x_i, x_{i+1}$ |
| **후진 차분(Backward Difference)** ($f'$) | $\dfrac{f_i - f_{i-1}}{h}$ | $O(h)$ | $x_{i-1}, x_i$ |
| **중심 차분(Central Difference)** ($f'$) | $\dfrac{f_{i+1} - f_{i-1}}{2h}$ | $O(h^2)$ | $x_{i-1}, x_{i+1}$ |
| **2차 전진(2nd-Order Forward)** ($f'$) | $\dfrac{-f_{i+2} + 4f_{i+1} - 3f_i}{2h}$ | $O(h^2)$ | $x_i, x_{i+1}, x_{i+2}$ |
| **중심 차분(Central Difference)** ($f''$) | $\dfrac{f_{i+1} - 2f_i + f_{i-1}}{h^2}$ | $O(h^2)$ | $x_{i-1}, x_i, x_{i+1}$ |
| **전진 차분(Forward Difference)** ($f''$) | $\dfrac{f_{i+2} - 2f_{i+1} + f_i}{h^2}$ | $O(h)$ | $x_i, x_{i+1}, x_{i+2}$ |
| **리처드슨 외삽 (1수준)(Richardson Extrap. 1 level)** | $\dfrac{4}{3}D(h) - \dfrac{1}{3}D(2h)$ | $O(h^4)$ | $D(h)$로부터 |
| **리처드슨 외삽 (2수준)(Richardson Extrap. 2 levels)** | $\dfrac{16}{15}\bar{D}(h) - \dfrac{1}{15}\bar{D}(2h)$ | $O(h^6)$ | $\bar{D}(h)$로부터 |
| **혼합 편도함수(Mixed Partial)** ($f_{xy}$) | $\dfrac{f(x\!+\!\Delta x, y\!+\!\Delta y) - f(x\!+\!\Delta x, y\!-\!\Delta y) - f(x\!-\!\Delta x, y\!+\!\Delta y) + f(x\!-\!\Delta x, y\!-\!\Delta y)}{4\Delta x \Delta y}$ | $O(\Delta x^2, \Delta y^2)$ | 대각선 4점 |

**핵심 내용**:
- 모든 유한 차분 공식은 **테일러 급수 전개** 로부터 유도된다
- **중심 차분** 은 대칭성에 의한 오차 상쇄로 인해 같은 간격 크기에서 전진/후진 차분보다 더 정확하다
- **리처드슨 외삽법** 은 서로 다른 간격 크기의 추정값을 결합하여 주도적 오차항을 체계적으로 제거하며, 정확도를 $O(h^2) \to O(h^4) \to O(h^6) \to \cdots$로 향상시킨다
- **비등간격 데이터** 에는 라그랑주 보간 다항식을 미분한다
- **미정 계수법** 은 임의의 유한 차분 스텐실을 유도하는 체계적인 방법을 제공한다
- **미분 가능 프로그래밍**(PyTorch, JAX, TensorFlow)은 자동 미분을 통해 정확한 기울기를 제공하며, 수치 미분의 현대적 대안을 나타낸다
