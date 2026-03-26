# 제17장 강의 — 다항식 보간(Polynomial Interpolation)

> **최종 수정일:** 2026-03-26

---

<br>

## 목차

- [1. 다항식 보간(Polynomial Interpolation)](#1-다항식-보간polynomial-interpolation)
  - [1.1 보간이란 무엇인가?](#11-보간이란-무엇인가)
  - [1.2 다항식 형태와 미지수](#12-다항식-형태와-미지수)
  - [1.3 계수 결정](#13-계수-결정)
  - [1.4 반데르몬드 행렬(Vandermonde Matrix)](#14-반데르몬드-행렬vandermonde-matrix)
  - [1.5 표준화(Standardization)](#15-표준화standardization)
- [2. 뉴턴 보간 다항식(Newton Interpolating Polynomial)](#2-뉴턴-보간-다항식newton-interpolating-polynomial)
  - [2.1 선형 보간(1차)](#21-선형-보간1차)
  - [2.2 이차 보간(2차)](#22-이차-보간2차)
  - [2.3 유한 분할 차분(Finite Divided Differences)](#23-유한-분할-차분finite-divided-differences)
  - [2.4 유한 차분의 재귀적 성질](#24-유한-차분의-재귀적-성질)
  - [2.5 뉴턴 보간 다항식의 일반 형태](#25-뉴턴-보간-다항식의-일반-형태)
  - [2.6 예제 17.4 — ln 2 추정](#26-예제-174--ln-2-추정)
- [3. 라그랑주 보간 다항식(Lagrange Interpolating Polynomial)](#3-라그랑주-보간-다항식lagrange-interpolating-polynomial)
  - [3.1 선형 라그랑주 보간(1차)](#31-선형-라그랑주-보간1차)
  - [3.2 기저 함수 N1(x)과 N2(x)](#32-기저-함수-n1x과-n2x)
  - [3.3 2차 라그랑주 보간 다항식](#33-2차-라그랑주-보간-다항식)
  - [3.4 일반 (n-1)차 라그랑주 다항식](#34-일반-n-1차-라그랑주-다항식)
- [4. 역보간(Inverse Interpolation)](#4-역보간inverse-interpolation)
  - [4.1 개념](#41-개념)
  - [4.2 예제](#42-예제)
- [5. 외삽과 진동(Extrapolation and Oscillations)](#5-외삽과-진동extrapolation-and-oscillations)
  - [5.1 외삽(Extrapolation)](#51-외삽extrapolation)
  - [5.2 진동(룽게 현상, Runge Phenomenon)](#52-진동룽게-현상-runge-phenomenon)
- [요약 표](#요약-표)

---

<br>

## 1. 다항식 보간(Polynomial Interpolation)

### 1.1 보간이란 무엇인가?

보간(Interpolation)은 이산 데이터 점들이 주어졌을 때 중간 지점에서 함수의 값을 추정하는 기법이다. 핵심 질문은: **알려진 두 데이터 점 사이의 값을 어떻게 채우는가?**

답: **다항식 보간**을 사용한다. 핵심 아이디어는 알려진 데이터 점을 통과하는 다항식을 구성한 다음, 원하는 중간 위치에서 그 다항식을 평가하는 것이다.

사용하는 데이터 점의 수에 따라 세 가지 일반적인 유형의 다항식 보간이 있다:

| 유형 | 필요한 데이터 점 | 다항식 차수 |
|------|---------------------|-----------------|
| 선형 보간 | 2 | 1차 |
| 이차 보간 | 3 | 2차 |
| 삼차 보간 | 4 | 3차 |

> **동기 부여 예제:** 온도 $T$ [C]와 공기의 밀도 $\rho$ [kg/m$^3$] 표가 주어졌다고 하자:
>
> | $T$ [C] | $\rho$ [kg/m$^3$] |
> |---------|-------------------|
> | $-40$   | $1.52$            |
> | $0$     | $1.29$            |
> | $20$    | $1.20$            |
> | $50$    | $1.09$            |
> | $\vdots$ | $\vdots$         |
> | $500$   | $0.46$            |
>
> $T = 30$ C에서의 밀도는 얼마인가? 보간이 이 질문에 답한다.

### 1.2 다항식 형태와 미지수

$n$차 다항식의 일반 형태:

$$f(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_n x^n$$

이 다항식은 **$n + 1$개의 미지수**(계수 $a_0, a_1, \ldots, a_n$)를 가지므로, 이를 유일하게 결정하려면 **$n + 1$개의 데이터 점**이 필요하다.

> **핵심 원리:** $n$차 다항식을 적합하려면 정확히 $n + 1$개의 데이터 점이 필요하다.

### 1.3 계수 결정

**예제 17.1:** 세 개의 데이터 점을 통과하는 2차 다항식을 구하고자 한다:

$$f(x) = p_0 x^2 + p_1 x + p_2 \qquad (\ast)$$

주어진 데이터:
- $f(300) = 0.616$
- $f(400) = 0.525$
- $f(500) = 0.457$

각 데이터 점을 $(\ast)$에 대입하면:

$$0.616 = p_0 (300)^2 + p_1 (300) + p_2$$

$$0.525 = p_0 (400)^2 + p_1 (400) + p_2$$

$$0.457 = p_0 (500)^2 + p_1 (500) + p_2$$

이것은 3개의 미지수($p_0, p_1, p_2$)에 대한 3개의 방정식으로 이루어진 연립방정식이다.

### 1.4 반데르몬드 행렬(Vandermonde Matrix)

1.3절의 연립방정식은 **벡터/행렬 형태**로 다시 쓸 수 있다:

$$\begin{pmatrix} 0.616 \\ 0.525 \\ 0.457 \end{pmatrix} = \begin{pmatrix} 300^2 & 300 & 1 \\ 400^2 & 400 & 1 \\ 500^2 & 500 & 1 \end{pmatrix} \begin{pmatrix} p_0 \\ p_1 \\ p_2 \end{pmatrix}$$

일반적으로, $x_1, x_2, x_3$에서의 데이터로 2차 다항식을 구하면:

$$\begin{pmatrix} f(x_1) \\ f(x_2) \\ f(x_3) \end{pmatrix} = \underbrace{\begin{pmatrix} x_1^2 & x_1 & 1 \\ x_2^2 & x_2 & 1 \\ x_3^2 & x_3 & 1 \end{pmatrix}}_{V \text{ (반데르몬드 행렬)}} \begin{pmatrix} p_0 \\ p_1 \\ p_2 \end{pmatrix}$$

계수 행렬 $V$는 **반데르몬드 행렬(Vandermonde Matrix)**로 알려져 있다. $V \mathbf{p} = \mathbf{f}$를 풀어 다항식 계수를 구한다.

> **문제:** 반데르몬드 행렬은 $x$ 값이 클 때 종종 **조건 수가 나쁘다(Ill-Conditioned)**. 위 예제의 경우 $\text{cond}(V) \sim 10^6$으로, 수치적 불안정성을 유발한다.

### 1.5 표준화(Standardization)

조건 수를 개선하기 위해, **새로운 표준화 변수**를 정의한다:

$$z := \frac{x - 400}{100}$$

이것은 원래의 $x$ 값 $\{300, 400, 500\}$을 $\{-1, 0, 1\}$로 매핑하여, 새로운 반데르몬드 행렬을 생성한다:

$$V' = \begin{pmatrix} 1 & -1 & 1 \\ 0 & 0 & 1 \\ 1 & 1 & 1 \end{pmatrix}$$

이제 $\text{cond}(V') = 3.26$으로, $10^6$보다 극적으로 개선된다.

> **표준화**는 반데르몬드 행렬의 조건 수를 개선하고 계산을 안정화하기 위해 필요하다.

---

<br>

## 2. 뉴턴 보간 다항식(Newton Interpolating Polynomial)

뉴턴의 접근법은 보간 다항식을 **점진적으로** 구축하며, 한 번에 하나의 항을 추가한다. 각 새 항은 더 높은 차수의 보정(곡률)을 도입한다.

### 2.1 선형 보간(1차)

두 데이터 점 $(x_1, f(x_1))$과 $(x_2, f(x_2))$가 주어졌을 때, **뉴턴 선형 보간 공식**은:

$$\hat{f}(x) = f(x_1) + \frac{f(x_2) - f(x_1)}{x_2 - x_1}(x - x_1)$$

이것은 단순히 두 점을 통과하는 직선(할선의 기울기를 사용한 기울기-절편 형태)이다.

> **정확도**는 구간 $|x_2 - x_1|$이 감소하고 $x$가 데이터 점에 가까울수록 증가한다.

### 2.2 이차 보간(2차)

세 데이터 점 $f(x_1), f(x_2), f(x_3)$이 주어졌을 때, 이차 보간 공식을 구성한다:

$$\hat{f}(x) = b_1 + b_2(x - x_1) + b_3(x - x_1)(x - x_2)$$

**계수 결정:**

**단계 1:** $\hat{f}(x_1) = f(x_1)$로 놓으면:

$$f(x_1) = b_1 \quad \Longrightarrow \quad b_1 = f(x_1)$$

**단계 2:** $\hat{f}(x_2) = f(x_2)$로 놓으면:

$$f(x_2) = f(x_1) + b_2(x_2 - x_1)$$

$$\therefore \quad b_2 = \frac{f(x_2) - f(x_1)}{x_2 - x_1} =: f[x_2, x_1]$$

이것이 **1차 유한 분할 차분**(기울기)이다.

**단계 3:** $\hat{f}(x_3) = f(x_3)$으로 놓고 $b_3$에 대해 풀면:

$$f(x_3) = f(x_1) + \frac{f(x_2) - f(x_1)}{x_2 - x_1}(x_3 - x_1) + b_3(x_3 - x_1)(x_3 - x_2)$$

대수적 조작 후:

$$b_3 = \frac{\dfrac{f(x_3) - f(x_2)}{x_3 - x_2} - \dfrac{f(x_2) - f(x_1)}{x_2 - x_1}}{x_3 - x_1} = \frac{f[x_3, x_2] - f[x_2, x_1]}{x_3 - x_1} =: f[x_3, x_2, x_1]$$

이것이 **2차 유한 분할 차분**이다.

### 2.3 유한 분할 차분(Finite Divided Differences)

유한 분할 차분은 재귀적으로 정의된다:

**0차** (함수값):

$$f[x_i] = f(x_i)$$

**1차** (기울기):

$$f[x_i, x_j] = \frac{f(x_i) - f(x_j)}{x_i - x_j}$$

**2차** (곡률):

$$f[x_i, x_j, x_k] = \frac{f[x_i, x_j] - f[x_j, x_k]}{x_i - x_k}$$

**$n$차** (일반):

$$f[x_n, x_{n-1}, \ldots, x_1] = \frac{f[x_n, x_{n-1}, \ldots, x_2] - f[x_{n-1}, x_{n-2}, \ldots, x_1]}{x_n - x_1}$$

### 2.4 유한 차분의 재귀적 성질

분할 차분은 **삼각 표**로 정리할 수 있다:

| $x_i$ | $f(x_i)$ | 1차 분할 차분 | 2차 분할 차분 | 3차 분할 차분 |
|--------|-----------|-------------------|--------------------|---------------------|
| $x_1$ | $f(x_1) \; {}^{\prime\prime}b_1{}^{\prime\prime}$ | | | |
| | | $f[x_2, x_1] \; {}^{\prime\prime}b_2{}^{\prime\prime}$ | | |
| $x_2$ | $f(x_2)$ | | $f[x_3, x_2, x_1] \; {}^{\prime\prime}b_3{}^{\prime\prime}$ | |
| | | $f[x_3, x_2]$ | | $f[x_4, x_3, x_2, x_1] \; {}^{\prime\prime}b_4{}^{\prime\prime}$ |
| $x_3$ | $f(x_3)$ | | $f[x_4, x_3, x_2]$ | |
| | | $f[x_4, x_3]$ | | |
| $x_4$ | $f(x_4)$ | | | |

표의 **위쪽 대각선**이 뉴턴 다항식의 계수 $b_1, b_2, b_3, b_4, \ldots$를 제공한다.

### 2.5 뉴턴 보간 다항식의 일반 형태

$n$개의 데이터 점에 대해, **(n-1)차** 뉴턴 보간 다항식은:

$$\hat{f}_{n-1}(x) = f(x_1) + f[x_2, x_1](x - x_1) + f[x_3, x_2, x_1](x - x_1)(x - x_2) + \cdots$$

$$+ f[x_n, x_{n-1}, \ldots, x_1](x - x_1)(x - x_2) \cdots (x - x_{n-1})$$

간결한 표기법:

$$\hat{f}_{n-1}(x) = \sum_{k=1}^{n} b_k \prod_{j=1}^{k-1}(x - x_j)$$

여기서

- $b_1 = f(x_1)$
- $b_2 = f[x_2, x_1]$
- $b_3 = f[x_3, x_2, x_1]$
- $b_k = f[x_k, x_{k-1}, \ldots, x_1]$

> **장점:** 뉴턴 형태는 **점진적(Incremental)**이다 — 새 데이터 점을 추가할 때 하나의 추가 분할 차분만 계산하고 항 하나를 추가하면 된다. 이전에 계산된 계수는 변경되지 않는다.

### 2.6 예제 17.4 — ln 2 추정

**문제:** 3차 뉴턴 보간 다항식을 사용하여 $\ln 2$를 추정하라.

**데이터 점** ($f(x) = \ln(x)$ 사용):

| $i$ | $x_i$ | $f(x_i) = \ln(x_i)$ |
|-----|--------|----------------------|
| 1   | 1      | 0                    |
| 2   | 4      | 1.386294             |
| 3   | 6      | 1.791759             |
| 4   | 5      | 1.609438             |

**단계 1:** 분할 차분 표 작성:

| | $f(x_i)$ | 1차 | 2차 | 3차 |
|---|----------|-----|-----|-----|
| $x_1 = 1$ | $0 = b_1$ | | | |
| | | $f[x_2, x_1]$ | | |
| $x_2 = 4$ | $1.386294$ | | $f[x_3, x_2, x_1]$ | |
| | | $f[x_3, x_2]$ | | $f[x_4, x_3, x_2, x_1] = b_4$ |
| $x_3 = 6$ | $1.791759$ | | $f[x_4, x_3, x_2]$ | |
| | | $f[x_4, x_3]$ | | |
| $x_4 = 5$ | $1.609438$ | | | |

**단계 2:** 각 분할 차분 계산:

$$f[x_2, x_1] = \frac{1.386294 - 0}{4 - 1} = 0.462098 = b_2$$

$$f[x_3, x_2] = \frac{1.791759 - 1.386294}{6 - 4} = 0.202733$$

$$f[x_4, x_3] = \frac{1.609438 - 1.791759}{5 - 6} = 0.182321$$

$$f[x_3, x_2, x_1] = \frac{0.202733 - 0.462098}{6 - 1} = -0.051873 = b_3$$

$$f[x_4, x_3, x_2] = \frac{0.182321 - 0.202733}{5 - 4} = -0.020412$$

$$f[x_4, x_3, x_2, x_1] = \frac{-0.020412 - (-0.051873)}{5 - 1} = 0.007865 = b_4$$

**단계 3:** 3차 다항식 구성:

$$\hat{f}_3(x) = 0 + 0.462098(x - 1) + (-0.051873)(x - 1)(x - 4) + 0.007865(x - 1)(x - 4)(x - 6)$$

**단계 4:** $x = 2$에서 평가:

$$\hat{f}_3(2) = 0.462098(1) + (-0.051873)(1)(-2) + 0.007865(1)(-2)(-4)$$

$$= 0.462098 + 0.103746 + 0.062920 = 0.628764$$

정확한 값은 $\ln 2 = 0.693147$이므로, 데이터 점의 간격으로 인해 추정에 다소 오차가 있다.

```python
import numpy as np

# Data points
x = np.array([1, 4, 6, 5])
f = np.array([0, 1.386294, 1.791759, 1.609438])

# Divided differences
def divided_diff(x, f):
    n = len(x)
    table = np.zeros((n, n))
    table[:, 0] = f
    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (x[i+j] - x[i])
    return table[0]  # top row = coefficients b1, b2, b3, ...

# Newton evaluation
def newton_eval(coeffs, x_data, x_eval):
    n = len(coeffs)
    result = coeffs[-1]
    for k in range(n - 2, -1, -1):
        result = result * (x_eval - x_data[k]) + coeffs[k]
    return result

coeffs = divided_diff(x, f)
print(f"Coefficients: {coeffs}")
print(f"f(2) = {newton_eval(coeffs, x, 2):.6f}")
print(f"Exact ln(2) = {np.log(2):.6f}")
```

---

<br>

## 3. 라그랑주 보간 다항식(Lagrange Interpolating Polynomial)

라그랑주 방법은 연립방정식을 풀거나 분할 차분 표를 만들 필요 없이 보간 다항식의 **대안적** 형태를 제공한다. 대신, **기저 함수(Basis Functions)**(또는 **형상 함수, Shape Functions**)를 사용하여 다항식을 직접 구성한다.

### 3.1 선형 라그랑주 보간(1차)

두 데이터 점 $(x_1, f(x_1))$과 $(x_2, f(x_2))$가 주어졌을 때, 이를 통과하는 직선은:

$$y - f(x_1) = \frac{f(x_2) - f(x_1)}{x_2 - x_1}(x - x_1)$$

정리하면:

$$y = \frac{x_2 - x}{x_2 - x_1} f(x_1) + \frac{x - x_1}{x_2 - x_1} f(x_2)$$

**기저 함수(형상 함수)**를 정의한다:

$$N_1(x) = \frac{x_2 - x}{x_2 - x_1}, \qquad N_2(x) = \frac{x - x_1}{x_2 - x_1}$$

그러면 **선형 라그랑주 보간 다항식**은:

$$\hat{f}_1(x) = N_1(x) \, f(x_1) + N_2(x) \, f(x_2)$$

### 3.2 기저 함수 N1(x)과 N2(x)

기저 함수는 중요한 **단위 분할(Partition of Unity)** 성질을 갖는다:

- $N_1(x_1) = 1, \quad N_1(x_2) = 0$
- $N_2(x_1) = 0, \quad N_2(x_2) = 1$

$x = x_1$일 때: 보간은 $f(x_1)$을 정확히 준다($N_1$에 의해 100% 가중).
$x = x_2$일 때: 보간은 $f(x_2)$를 정확히 준다($N_2$에 의해 100% 가중).
$x_1$과 $x_2$ 사이: 결과는 $f(x_1)$과 $f(x_2)$의 **가중 평균**이다.

> **기하학적 해석:** $N_1(x)$과 $N_2(x)$는 중간 지점에서 교차하는 선형 함수이다. 보간값 $y = N_1(x) f(x_1) + N_2(x) f(x_2)$는 데이터 점에서의 함수값의 가중 결합이다.

### 3.3 2차 라그랑주 보간 다항식

세 데이터 점 $(x_1, f(x_1)), (x_2, f(x_2)), (x_3, f(x_3))$에 대해, 세 개의 라그랑주 기저 함수를 정의한다:

$$L_1(x) = \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)}$$

$$L_2(x) = \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)}$$

$$L_3(x) = \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)}$$

2차 라그랑주 보간 다항식은:

$$\hat{f}_2(x) = L_1(x) \, f(x_1) + L_2(x) \, f(x_2) + L_3(x) \, f(x_3) = \sum_{i=1}^{3} L_i(x) \, f(x_i)$$

각 $L_i(x)$는 다음을 만족한다:

$$L_i(x_j) = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{if } i \neq j \end{cases}$$

> **참고:** 2차 다항식에는 3개의 데이터 점이 필요하다.

### 3.4 일반 (n-1)차 라그랑주 다항식

$n$개의 데이터 점 $x_1, x_2, \ldots, x_n$에 대해, **(n-1)차 라그랑주 보간 다항식**은:

$$\hat{f}_{n-1}(x) = \sum_{i=1}^{n} L_i(x) \, f(x_i)$$

여기서 라그랑주 기저 함수는:

$$L_i(x) = \prod_{\substack{j=1 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}$$

여기서:
- $n$은 **데이터 점의 수**
- $\prod$는 $j = i$를 제외한 모든 $j$에 대한 **곱**

```python
import numpy as np

def lagrange_interp(x_data, f_data, x_eval):
    """Lagrange interpolation at a single point x_eval."""
    n = len(x_data)
    result = 0.0
    for i in range(n):
        # Compute L_i(x_eval)
        Li = 1.0
        for j in range(n):
            if j != i:
                Li *= (x_eval - x_data[j]) / (x_data[i] - x_data[j])
        result += Li * f_data[i]
    return result

# Example: Estimate ln(2) using the same 4 data points
x_data = np.array([1, 4, 6, 5], dtype=float)
f_data = np.log(x_data)

estimate = lagrange_interp(x_data, f_data, 2.0)
print(f"Lagrange f(2) = {estimate:.6f}")
print(f"Exact ln(2) = {np.log(2):.6f}")
```

> **뉴턴 vs 라그랑주:** 두 방법 모두 **동일한** 보간 다항식을 생성한다 — 대수적으로 동치이다. 뉴턴 형태는 새 데이터 점을 **점진적으로 추가**할 때 더 효율적이다. 라그랑주 형태는 개념적으로 더 단순하고 구현이 더 직관적이다.

---

<br>

## 4. 역보간(Inverse Interpolation)

### 4.1 개념

표준 보간에서는 $x$가 주어지고 $f(x)$를 추정한다. **역보간(Inverse Interpolation)**에서는 목표값 $y$가 주어지고 $f(x) = y$인 $x$를 찾고자 한다.

데이터 점 $(x_0, f(x_0)), (x_1, f(x_1)), (x_2, f(x_2)), \ldots$이 주어졌을 때, $f(x) = y$인 $x$를 구하고자 한다.

**핵심 아이디어:** 명시적인 역함수 $f^{-1}(x)$가 필요하지 않다. 대신, $x$와 $f(x)$의 **역할을 교체**하여 — $f(x_i)$를 독립 변수로, $x_i$를 종속 변수로 취급한 다음, 표준 보간을 적용한다.

$$x(y) = \sum_{i=0}^{n} x_i \, L_i(y)$$

여기서 $L_i(y)$는 $f(x_i)$ 값을 "x-데이터"로 사용하여 구성한 라그랑주 기저 함수이다.

### 4.2 예제

세 데이터 점이 주어졌다:

$$
(x_0, f(x_0)) = (2, 1), \quad
(x_1, f(x_1)) = (4, 2), \quad
(x_2, f(x_2)) = (8, 3)
$$

**목표:** $f(x) = 5$인 $x$를 구하라.

**단계 1:** 쌍을 반전:

$$(f(x_0), x_0) = (1, 2), \quad (f(x_1), x_1) = (2, 4), \quad (f(x_2), x_2) = (3, 8)$$

**단계 2:** $y$ 값을 독립 변수로 하여 라그랑주 보간을 적용:

$$x(y) = \sum_{i=0}^{2} x_i \, L_i(y)$$

**단계 3:** $y = 5$에서 평가:

$$x(5) = \sum_{i=0}^{2} x_i \, L_i(5)$$

각 기저 함수를 $y = 5$에서 계산:

$$L_0(5) = \frac{(5-2)(5-3)}{(1-2)(1-3)} = \frac{(3)(2)}{(-1)(-2)} = 3$$

$$L_1(5) = \frac{(5-1)(5-3)}{(2-1)(2-3)} = \frac{(4)(2)}{(1)(-1)} = -8$$

$$L_2(5) = \frac{(5-1)(5-2)}{(3-1)(3-2)} = \frac{(4)(3)}{(2)(1)} = 6$$

$$x(5) = 2 \cdot 3 + 4 \cdot (-8) + 8 \cdot 6 = 6 - 32 + 48 = 22$$

```python
import numpy as np

# Inverse interpolation example
y_data = np.array([1, 2, 3], dtype=float)   # f(x) values as independent variable
x_data = np.array([2, 4, 8], dtype=float)   # x values as dependent variable

def lagrange_interp(ind, dep, target):
    n = len(ind)
    result = 0.0
    for i in range(n):
        Li = 1.0
        for j in range(n):
            if j != i:
                Li *= (target - ind[j]) / (ind[i] - ind[j])
        result += Li * dep[i]
    return result

x_at_y5 = lagrange_interp(y_data, x_data, 5.0)
print(f"x such that f(x)=5: {x_at_y5:.2f}")
```

---

<br>

## 5. 외삽과 진동(Extrapolation and Oscillations)

### 5.1 외삽(Extrapolation)

**외삽(Extrapolation)**은 알려진 데이터 점 범위 $[x_1, x_n]$ **밖**에 있는 $f(x)$의 값을 추정하는 것이다.

- **보간(Interpolation):** $x \in [x_1, x_n]$ (알려진 점들 사이)
- **외삽(Extrapolation):** $x < x_1$ 또는 $x > x_n$ (알려진 점들 밖)

> **경고:** 다항식 보간을 이용한 외삽은 일반적으로 **신뢰할 수 없다**. 데이터 범위 밖에서는 다항식의 동작을 제약할 데이터가 없으므로, 다항식이 급격히 발산할 수 있다. 오차가 극적으로 증가할 수 있다.

### 5.2 진동(룽게 현상, Runge Phenomenon)

보간에 **고차 다항식**을 사용하면 특히 데이터 구간의 가장자리 근처에서 **가짜 진동**이 발생할 수 있다. 이를 **룽게 현상(Runge's Phenomenon)**이라 한다.

**예제 17.7 — 룽게 함수:**

$$f(x) = \frac{1}{1 + 25x^2}, \qquad x \in [-1, 1]$$

이 함수를 등간격 점을 사용하여 고차 다항식(예: 4차 이상)으로 보간하면, 모든 데이터 점을 정확히 통과함에도 불구하고 경계 $x = \pm 1$ 근처에서 다항식이 큰 진동을 보인다.

> **룽게 현상:** 등간격 노드에서 다항식 차수가 증가하면, 보간 오차가 끝점 근처에서 감소하지 않고 오히려 **증가**할 수 있다. 이것은 고차 다항식 보간의 근본적인 한계이다.

**진동에 대한 해결책:**
1. 단일 고차 다항식 대신 **저차** 구간별 다항식(스플라인)을 사용한다
2. 경계 근처에 밀집되는 **비등간격** 노드(예: 체비셰프 노드, Chebyshev Nodes)를 사용한다
3. 다항식 차수를 합리적으로 낮게 유지하고 더 많은 구간을 사용한다

```python
import numpy as np
import matplotlib.pyplot as plt

# Runge function demonstration
f_runge = lambda x: 1 / (1 + 25 * x**2)

x_fine = np.linspace(-1, 1, 500)
y_true = f_runge(x_fine)

for n in [5, 9, 13]:
    x_nodes = np.linspace(-1, 1, n)
    y_nodes = f_runge(x_nodes)
    # Lagrange interpolation
    y_interp = np.zeros_like(x_fine)
    for i in range(n):
        Li = np.ones_like(x_fine)
        for j in range(n):
            if j != i:
                Li *= (x_fine - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        y_interp += Li * y_nodes[i]
    plt.plot(x_fine, y_interp, label=f"Order {n-1}")

plt.plot(x_fine, y_true, 'k--', linewidth=2, label="True f(x)")
plt.ylim(-0.5, 1.5)
plt.legend()
plt.title("Runge Phenomenon")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.show()
```

---

<br>

## 요약 표

| 주제 | 핵심 수식 / 개념 | 필요한 점 수 | 차수 |
|-------|----------------------|---------------|-------|
| **다항식 보간** | $f(x) = a_0 + a_1 x + \cdots + a_n x^n$ | $n+1$ | $n$ |
| **반데르몬드 행렬** | $V \mathbf{p} = \mathbf{f}$ | 계수를 풀어줌 | 조건 수가 나쁜 경우가 많음 |
| **표준화** | $z = \frac{x - \bar{x}}{\Delta x}$ | $\text{cond}(V)$를 개선 | 동일한 차수 |
| **뉴턴 선형** | $\hat{f}(x) = f(x_1) + f[x_2,x_1](x - x_1)$ | 2 | 1 |
| **뉴턴 이차** | $+ \, f[x_3,x_2,x_1](x-x_1)(x-x_2)$ | 3 | 2 |
| **뉴턴 일반** | $\hat{f}_{n-1}(x) = \sum b_k \prod_{j=1}^{k-1}(x - x_j)$ | $n$ | $n-1$ |
| **분할 차분** | $f[x_i,x_j] = \frac{f(x_i)-f(x_j)}{x_i - x_j}$ | 재귀적 | 표를 구축 |
| **라그랑주 선형** | $\hat{f}_1 = N_1(x)f(x_1) + N_2(x)f(x_2)$ | 2 | 1 |
| **라그랑주 일반** | $\hat{f}_{n-1}(x) = \sum L_i(x) f(x_i)$ | $n$ | $n-1$ |
| **라그랑주 기저** | $L_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}$ | 크로네커 델타 성질 | -- |
| **역보간** | $x$와 $f(x)$의 역할을 교체한 후 보간 | 순방향과 동일 | 동일 |
| **외삽** | $[x_1, x_n]$ 밖에서의 평가 | -- | 신뢰 불가 |
| **진동(룽게)** | $f(x) = \frac{1}{1+25x^2}$ | 고차 + 등간격 | 가장자리에서 발산 |
