# 15장 — 비선형 회귀

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 15

> **학습 목표:**
> 1. 지수, 거듭제곱, 포화 성장률 모델을 선형화하여 선형 회귀로 적합할 수 있다
> 2. 다항 회귀와 다중 선형 회귀의 정규 방정식을 유도하고 적용할 수 있다
> 3. 일반 선형 최소제곱 모델의 정규 방정식을 행렬로 유도하고 적합도 통계량을 계산할 수 있다
> 4. 매개변수에 대해 본질적으로 비선형인 모델을 반복 최적화로 적합할 수 있다

---

<br>

## 목차

- [1. 복습: 표준 오차와 결정 계수](#1-복습-표준-오차와-결정-계수)
  - [1.1 표본 표준 편차(SD)](#11-표본-표준-편차sd)
  - [1.2 평균의 표준 오차(SE)](#12-평균의-표준-오차se)
  - [1.3 회귀의 표준 오차](#13-회귀의-표준-오차)
  - [1.4 결정 계수(R-제곱)](#14-결정-계수r-제곱)
- [2. 비선형 관계의 선형화(Linearization)](#2-비선형-관계의-선형화linearization)
  - [2.1 지수 모델(Exponential Model)](#21-지수-모델exponential-model)
  - [2.2 거듭제곱 모델(Power Model)](#22-거듭제곱-모델power-model)
  - [2.3 포화 성장률 모델(Saturation-Growth-Rate Model)](#23-포화-성장률-모델saturation-growth-rate-model)
- [3. 다항 회귀(Polynomial Regression)](#3-다항-회귀polynomial-regression)
  - [3.1 문제 정의](#31-문제-정의)
  - [3.2 정규 방정식 유도](#32-정규-방정식-유도)
  - [3.3 행렬 형태의 정규 방정식](#33-행렬-형태의-정규-방정식)
  - [3.4 추정의 표준 오차](#34-추정의-표준-오차)
- [4. 다중 선형 회귀(Multiple Linear Regression)](#4-다중-선형-회귀multiple-linear-regression)
  - [4.1 문제 정의](#41-문제-정의)
  - [4.2 정규 방정식](#42-정규-방정식)
  - [4.3 추정의 표준 오차](#43-추정의-표준-오차)
  - [4.4 다변수 거듭제곱 모델로의 확장](#44-다변수-거듭제곱-모델로의-확장)
- [5. 일반 선형 최소제곱 회귀(General Linear Least Square Regression)](#5-일반-선형-최소제곱-회귀general-linear-least-square-regression)
  - [5.1 일반 모델 정의](#51-일반-모델-정의)
  - [5.2 행렬 정의](#52-행렬-정의)
  - [5.3 정규 방정식 유도(행렬 미적분)](#53-정규-방정식-유도행렬-미적분)
  - [5.4 해 계산](#54-해-계산)
  - [5.5 적합도 통계량](#55-적합도-통계량)
  - [5.6 매개변수 추정치의 공분산 행렬](#56-매개변수-추정치의-공분산-행렬)
- [6. 비선형 회귀(가우스-뉴턴, Gauss-Newton)](#6-비선형-회귀가우스-뉴턴-gauss-newton)
  - [6.1 동기: 앙투안 방정식(Antoine Equation)](#61-동기-앙투안-방정식antoine-equation)
  - [6.2 scipy.optimize를 이용한 Python 구현](#62-scipyoptimize를-이용한-python-구현)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 복습: 표준 오차와 결정 계수

### 1.1 표본 표준 편차(SD)

**표본 표준 편차**(SD)는 개별 데이터 점이 데이터셋의 평균으로부터 얼마나 변동(퍼짐)하는지를 측정한다.

```math
s_y = \sqrt{\frac{1}{n-1} \sum_{j=1}^{n} (y_j - \bar{y})^2}
```

여기서:
- $`n`$ = 데이터 점의 수
- $`y_j`$ = 개별 관측값
- $`\bar{y}`$ = 표본 평균

**예제:** 키 = [170, 172, 168, 174, 176] cm

- 표본 평균 $`\bar{y} = 172`$ cm
- $`s_y = 3.16`$ cm

> **참고:** **해석.** 개별 키 측정값은 평균과 대략 3.16 cm 차이가 난다.

### 1.2 평균의 표준 오차(SE)

**표준 오차**(SE)는 표본 평균이 참 모집단 평균을 얼마나 정밀하게 추정하는지를 측정한다.

```math
SE = \frac{s_y}{\sqrt{n}}
```

**예제 (계속):**
- $`SE = \frac{3.16}{\sqrt{5}} = 1.41`$ cm

> **참고:** **해석.** 평균 키는 172 cm로 추정되며, 표준 오차는 $`\pm 1.41`$ cm이다.

### 1.3 회귀의 표준 오차

**회귀의 표준 오차**($`s_e`$)는 관측값이 선형 회귀선으로부터 얼마나 벗어나는지를 측정한다.

```math
s_e = \sqrt{\frac{1}{n-2} \sum_{j=1}^{n} (y_j - \hat{y}_j)^2} = \sqrt{\frac{SSE}{n-2}}
```

여기서:
- $`\hat{y}_j`$ = 회귀선에서의 예측값
- $`n - 2`$ = 자유도(두 매개변수 $`\beta_0, \beta_1`$이 추정됨)

**예제:** 키 = [170, 172, 168, 174, 176] cm, $`x = [1, 2, 3, 4, 5]`$

- $`s_e = \sqrt{\frac{SSE}{n-2}} = 2.61`$ cm

> **참고:** **해석.** 관측된 $`y`$ 값은 회귀선으로부터 약 2.61 cm 벗어난다.

### 1.4 결정 계수(R-제곱)

**결정 계수** $`R^2`$는 회귀 모델이 종속 변수의 변동성을 얼마나 잘 설명하는지를 측정한다.

총 제곱합은 다음과 같이 분해된다:

```math
SST = SSR + SSE
```

양변을 $`SST`$로 나누면:

```math
1 = \frac{SSR}{SST} + \frac{SSE}{SST}
```

따라서:

```math
R^2 = \frac{SSR}{SST} = 1 - \frac{SSE}{SST}
```

여기서:
- $`SST = \sum (y_i - \bar{y})^2`$: 총 제곱합(Total Sum of Squares)
- $`SSR = \sum (\hat{y}_i - \bar{y})^2`$: 회귀 제곱합(설명된 변동)
- $`SSE = \sum (y_i - \hat{y}_i)^2`$: 오차 제곱합(설명되지 않은 변동)

| $`R^2`$ 값 | 의미 |
|---|---|
| $`R^2 = 1`$ | 완벽한 적합 |
| $`R^2 = 0`$ | 적합 없음 |
| $`0 < R^2 < 1`$ | 부분적 적합 |

*표 1. 결정 계수 값에 따른 적합의 의미*

**예제 (교재 예제 14.6):** $`R^2 = 88\%`$

> **참고:** **해석.** 두 변수 간의 관계는 데이터 변동의 88%(약 90%)를 설명한다.

---

<br>

## 2. 비선형 관계의 선형화(Linearization)

핵심 아이디어는 **비선형 함수를 직선 형태로 변환** 하여 선형 회귀 기법을 적용할 수 있게 하는 것이다.

### 2.1 지수 모델(Exponential Model)

**원래 모델:**

```math
y = \alpha_1 e^{\beta_1 x}
```

**선형화:** 양변에 자연로그를 취한다:

```math
\ln y = \ln(\alpha_1 e^{\beta_1 x}) = \ln \alpha_1 + \ln e^{\beta_1 x} = \ln \alpha_1 + \beta_1 x
```

이로써 다음 선형 형태를 얻는다: $`\ln y = \underbrace{\ln \alpha_1}_{\text{절편}} + \underbrace{\beta_1}_{\text{기울기}} \cdot x`$

> **참고:** **그래프 효과.** $`x`$에 대한 $`y`$의 그래프는 지수 곡선이다. $`x`$에 대한 $`\ln y`$의 그래프는 기울기 $`\beta_1`$, 절편 $`\ln \alpha_1`$인 직선이 된다.

### 2.2 거듭제곱 모델(Power Model)

**원래 모델:**

```math
y = \alpha_2 x^{\beta_2}
```

**선형화:** 양변에 자연로그를 취한다:

```math
\ln y = \ln(\alpha_2 x^{\beta_2}) = \ln \alpha_2 + \ln x^{\beta_2} = \ln \alpha_2 + \beta_2 \ln x
```

이로써 다음 선형 형태를 얻는다: $`\ln y = \underbrace{\ln \alpha_2}_{\text{절편}} + \underbrace{\beta_2}_{\text{기울기}} \cdot \ln x`$

> **참고:** **그래프 효과.** $`x`$에 대한 $`y`$의 그래프는 거듭제곱 곡선이다. $`\ln x`$에 대한 $`\ln y`$의 그래프는 기울기 $`\beta_2`$, 절편 $`\ln \alpha_2`$인 직선이 된다.

### 2.3 포화 성장률 모델(Saturation-Growth-Rate Model)

**원래 모델:**

```math
y = \frac{\alpha_3 x}{\beta_3 + x}
```

**선형화:** 양변의 역수를 취한다:

```math
\frac{1}{y} = \frac{\beta_3 + x}{\alpha_3 x} = \frac{\beta_3}{\alpha_3} \cdot \frac{1}{x} + \frac{1}{\alpha_3}
```

이로써 다음 선형 형태를 얻는다: $`\frac{1}{y} = \underbrace{\frac{\beta_3}{\alpha_3}}_{\text{기울기}} \cdot \frac{1}{x} + \underbrace{\frac{1}{\alpha_3}}_{\text{절편}}`$

> **참고:** **그래프 효과.** $`x`$에 대한 $`y`$의 그래프는 포화 곡선이다. $`\frac{1}{x}`$에 대한 $`\frac{1}{y}`$의 그래프는 기울기 $`\frac{\beta_3}{\alpha_3}`$, 절편 $`\frac{1}{\alpha_3}`$인 직선이 된다.

| 모델 | 원래 형태 | 변환된 변수 | 기울기 | 절편 |
|---|---|---|---|---|
| 지수 모델 | $`y = \alpha_1 e^{\beta_1 x}`$ | $`x`$에 대한 $`\ln y`$ | $`\beta_1`$ | $`\ln \alpha_1`$ |
| 거듭제곱 모델 | $`y = \alpha_2 x^{\beta_2}`$ | $`\ln x`$에 대한 $`\ln y`$ | $`\beta_2`$ | $`\ln \alpha_2`$ |
| 포화 성장 모델 | $`y = \frac{\alpha_3 x}{\beta_3 + x}`$ | $`\frac{1}{x}`$에 대한 $`\frac{1}{y}`$ | $`\frac{\beta_3}{\alpha_3}`$ | $`\frac{1}{\alpha_3}`$ |

*표 2. 선형화 모델의 변환 변수, 기울기, 절편*

---

<br>

## 3. 다항 회귀(Polynomial Regression)

### 3.1 문제 정의

데이터가 곡선 형태(예: 포물선)를 보일 때, 직선만으로는 충분하지 않다. 선형 회귀의 아이디어를 **다항 회귀**로 확장한다.

다항 회귀로 추정한 매개변수의 신뢰 구간은 5.6절의 공분산 행렬로 구한다.

데이터셋 $`\{x_i, y_i\}_{i=1}^{n}`$이 주어졌을 때, 다음을 만족하는 $`\beta_0, \beta_1, \beta_2`$를 구한다:

```math
\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2
```

목적은 오차 제곱합을 최소화하는 것이다:

```math
\min_{\beta_0, \beta_1, \beta_2} SSE
```

여기서:

```math
SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)^2
```

### 3.2 정규 방정식 유도

$`\beta_0, \beta_1, \beta_2`$를 결정하기 위해 편도함수를 취하고 0으로 놓는다:

```math
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-1) = 0
```

```math
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i) = 0
```

```math
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i^2) = 0
```

정리하면, **정규 방정식**은 다음과 같다:

```math
\sum y_i = \beta_0 \cdot n + \beta_1 \sum x_i + \beta_2 \sum x_i^2
```

```math
\sum y_i x_i = \beta_0 \sum x_i + \beta_1 \sum x_i^2 + \beta_2 \sum x_i^3
```

```math
\sum y_i x_i^2 = \beta_0 \sum x_i^2 + \beta_1 \sum x_i^3 + \beta_2 \sum x_i^4
```

### 3.3 행렬 형태의 정규 방정식

정규 방정식은 다음과 같이 간결하게 표현할 수 있다:

```math
\begin{bmatrix} \sum y_i \\ \sum y_i x_i \\ \sum y_i x_i^2 \end{bmatrix}
=
\begin{bmatrix} n & \sum x_i & \sum x_i^2 \\ \sum x_i & \sum x_i^2 & \sum x_i^3 \\ \sum x_i^2 & \sum x_i^3 & \sum x_i^4 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
```

> **핵심:** 이것은 9장의 가우스 소거법이나 10장의 LU 분해 같은 연립 일차 방정식의 풀이법으로 풀 수 있는 $`3 \times 3`$ 선형 시스템이다.

### 3.4 추정의 표준 오차

3개의 매개변수($`\beta_0, \beta_1, \beta_2`$)를 갖는 다항 회귀의 경우:

```math
e_i = y_i - \hat{y}_i
```

```math
SSE = \sum e_i^2
```

```math
s_e = \sqrt{\frac{SSE}{n - 3}}
```

> **참고:** **자유도.** 분모가 $`n - 3`$인 이유는 세 개의 매개변수($`\beta_0, \beta_1, \beta_2`$)가 $`\hat{y}`$를 구성하는 데 사용되었기 때문이다. 일반적으로 $`m`$차 다항식의 자유도는 $`n - (m+1)`$이다.

---

<br>

## 4. 다중 선형 회귀(Multiple Linear Regression)

다중 선형 회귀는 1차원에서 2차원(및 그 이상)으로 개념을 확장한다: **두 개 이상의 독립 변수**를 사용하여 데이터를 적합한다.

### 4.1 문제 정의

모델은 데이터에 평면(또는 초평면)을 적합한다:

```math
\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2
```

목적은:

```math
SSE = \sum e_i^2 = \sum (y_i - \hat{y}_i)^2 = \sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})^2
```

### 4.2 정규 방정식

$`\min_{\beta_0, \beta_1, \beta_2} SSE`$를 만족하는 $`\beta_0, \beta_1, \beta_2`$를 구한다:

```math
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-1) = 0
```

```math
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{1i}) = 0
```

```math
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{2i}) = 0
```

행렬 형태:

```math
\begin{bmatrix} \sum y_i \\ \sum y_i x_{1i} \\ \sum y_i x_{2i} \end{bmatrix}
=
\begin{bmatrix} n & \sum x_{1i} & \sum x_{2i} \\ \sum x_{1i} & \sum x_{1i}^2 & \sum x_{1i} x_{2i} \\ \sum x_{2i} & \sum x_{1i} x_{2i} & \sum x_{2i}^2 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
```

### 4.3 추정의 표준 오차

두 개의 독립 변수(3개의 매개변수)를 갖는 다중 선형 회귀의 경우:

```math
s_e = \sqrt{\frac{SSE}{n - 3}}
```

### 4.4 다변수 거듭제곱 모델로의 확장

다중 선형 회귀는 **다변수 거듭제곱 모델**로 확장될 수 있다:

```math
y = \beta_0 x_1^{\beta_1} x_2^{\beta_2} \cdots x_m^{\beta_m} + \varepsilon
```

오차 항 $`\varepsilon`$을 무시하고 자연로그를 취하면:

```math
\ln y = \ln \beta_0 + \beta_1 \ln x_1 + \beta_2 \ln x_2 + \cdots + \beta_m \ln x_m
```

> **핵심:** 로그 변환 후, 이것은 변환된 변수 $`\ln y`$, $`\ln x_1`$, $`\ln x_2`$, ..., $`\ln x_m`$에 대한 표준 다중 선형 회귀 문제가 된다.

---

<br>

## 5. 일반 선형 최소제곱 회귀(General Linear Least Square Regression)

### 5.1 일반 모델 정의

일반 선형 최소제곱 모델은 이전의 모든 경우를 통합한다. 모델은 매개변수 $`\beta_j`$에 대해 선형이지만, 기저 함수 $`f_j`$는 임의적(예: 삼각 함수)일 수 있다:

**예제(푸리에 유형 모델):**

```math
y = \beta_0 + \beta_1 \cos(\omega x) + \beta_2 \sin(\omega x)
```

**일반 형태:**

```math
y = \beta_0 f_0(x_1, x_2, \ldots, x_p) + \beta_1 f_1(x_1, x_2, \ldots, x_p) + \cdots + \beta_m f_m(x_1, x_2, \ldots, x_p)
```

> **핵심:** 여기서 "선형"은 $`x`$가 아니라 **매개변수** $`\beta_j`$에 대해 선형이라는 의미이다. 함수 $`f_j`$는 $`x`$에 대해 비선형일 수 있다.

### 5.2 행렬 정의

$`n`$개의 데이터 점 $`\{y_i, x_{1i}, x_{2i}, \ldots, x_{pi}\}_{i=1}^{n}`$을 수집한다.

각 데이터 점 $`i = 1, 2, \ldots, n`$에 대해:

```math
y_i = \beta_0 f_0(x_{1i}, x_{2i}, \ldots, x_{pi}) + \beta_1 f_1(x_{1i}, x_{2i}, \ldots, x_{pi}) + \cdots + \beta_m f_m(x_{1i}, x_{2i}, \ldots, x_{pi})
```

행렬 형태:

```math
\underset{(n \times 1)}{\mathbf{y}}
=
\underset{(n \times (m+1))}{\mathbf{X}}
\underset{((m+1) \times 1)}{\boldsymbol{\beta}}
+ \underset{(n \times 1)}{\mathbf{e}}
```

여기서:

```math
\mathbf{y} = \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{bmatrix}, \quad
\mathbf{X} = \begin{bmatrix} f_{01} & f_{11} & f_{21} & \cdots & f_{m1} \\ f_{02} & f_{12} & f_{22} & \cdots & f_{m2} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ f_{0n} & f_{1n} & f_{2n} & \cdots & f_{mn} \end{bmatrix}, \quad
\boldsymbol{\beta} = \begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \\ \vdots \\ \beta_m \end{bmatrix}, \quad
\mathbf{e} = \begin{bmatrix} e_1 \\ e_2 \\ \vdots \\ e_n \end{bmatrix}
```

여기서 $`f_{ji}`$는 데이터 점 $`i`$에서 평가된 기저 함수 $`f_j`$의 값이고, $`\mathbf{e}`$는 잔차(오차) 벡터이다.

> **참고:** **표기법.** $`\mathbf{X}`$는 $`n`$개의 행(데이터 점)과 $`m+1`$개의 열(매개변수 수, $`\beta_0`$ 포함)을 갖는다.

### 5.3 정규 방정식 유도(행렬 미적분)

행렬 형태의 오차 제곱합:

```math
SSE = \sum e_i^2 = \mathbf{e}^T \mathbf{e} = (\mathbf{y} - \hat{\mathbf{y}})^T (\mathbf{y} - \hat{\mathbf{y}})
```

최소화하기 위해 $`\boldsymbol{\beta}`$에 대한 도함수를 취하고 0으로 놓는다:

```math
\frac{\partial SSE}{\partial \boldsymbol{\beta}} = \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{e}^T \mathbf{e})
= \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})
= (-2\mathbf{X})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta}) = \mathbf{0}
```

> **[미적분학]** 다변수 함수의 기울기(gradient)는 각 변수에 대한 편도함수를 모은 벡터이며, 미분 가능한 함수의 극솟값에서는 모든 편도함수가 0이다. SSE는 $`\boldsymbol{\beta}`$의 이차 함수이므로 3.2절에서 성분마다 편미분한 계산을 벡터로 한꺼번에 쓰면 위의 기울기가 되고, 이를 0으로 놓은 식이 정규 방정식이다.

이로부터 **정규 방정식**을 얻는다:

```math
\mathbf{X}^T \mathbf{X} \boldsymbol{\beta} = \mathbf{X}^T \mathbf{y}
```

$`\boldsymbol{\beta}`$에 대해 풀면:

```math
\boxed{\boldsymbol{\beta} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}}
```

> **핵심:** 이것은 이 장의 가장 중요한 결과이다. 임의의 일반 선형 모델에 적용 가능한, 오차 제곱합을 최소화하는 매개변수 벡터 $`\boldsymbol{\beta}`$에 대한 닫힌 형태의 해를 제공한다.

### 5.4 해 계산

$`\boldsymbol{\beta}`$가 계산되면:

| 양 | 수식 |
|---|---|
| 잔차 벡터 | $`\mathbf{e} = \mathbf{y} - \mathbf{X}\boldsymbol{\beta}`$ |
| 오차 제곱합 | $`SSE = \mathbf{e}^T \mathbf{e}`$ |
| 예측값 | $`\hat{\mathbf{y}} = \mathbf{X}\boldsymbol{\beta}`$ |

*표 3. 매개변수 벡터를 구한 뒤 계산하는 양*

### 5.5 적합도 통계량

```math
SST = \sum (y_i - \bar{y})^2
```

```math
SSR = \sum (\hat{y}_i - \bar{y})^2
```

```math
R^2 = \frac{SSR}{SST}
```

```math
SSE = SST - SSR
```

**추정의 표준 오차** ($`m+1`$개의 매개변수를 갖는 일반적인 경우):

```math
s_e = \sqrt{\frac{SSE}{n - (m+1)}}
```

> **참고:** **자유도.** 분모 $`n - (m+1)`$은 모델 적합에 사용된 $`m+1`$개의 매개변수를 고려한다.

### 5.6 매개변수 추정치의 공분산 행렬

추정된 매개변수의 공분산 행렬:

```math
\text{cov}(\hat{\boldsymbol{\beta}}) = \sigma_e^2 (\mathbf{X}^T \mathbf{X})^{-1}
```

여기서 $`\sigma_e^2`$는 $`s_e^2 = \frac{SSE}{n-(m+1)}`$로 추정된다.

> **참고:** **실용적 활용.** $`\text{cov}(\hat{\boldsymbol{\beta}})`$의 대각 원소는 각 개별 매개변수 추정치의 분산을 제공하며, 그 제곱근은 매개변수의 표준 오차를 제공한다. 이를 사용하여 신뢰 구간을 구성할 수 있다.

---

<br>

## 6. 비선형 회귀(가우스-뉴턴, Gauss-Newton)

### 6.1 동기: 앙투안 방정식(Antoine Equation)

일부 모델은 매개변수에 대해 **본질적으로 비선형**이며 변환으로 선형화할 수 없다. 예를 들어, 증기압에 대한 **앙투안 방정식(Antoine Equation)**:

```math
\log_{10} P_v = A - \frac{B}{C + T}
```

여기서:
- $`P_v`$ = 증기압
- $`T`$ = 온도
- $`A, B, C`$ = 조정 가능한 매개변수

$`C`$가 비선형 표현(분모의 $`C + T`$) 안에 나타나므로, 이 방정식은 일반 선형 모델로 표현할 수 **없다**. 이를 적합하려면 비선형 회귀(최적화)가 필요하다.

> **참고:** **가우스-뉴턴 방법.** 비선형 회귀 전용 기법인 가우스-뉴턴 방법(Gauss-Newton method)은 현재 매개변수 추정값 주위에서 모델 $`f(x_i; a_0, \ldots, a_m)`$을 1차 테일러 급수로 선형화한다. 그러면 잔차 벡터 $`\mathbf{d}`$($`d_i = y_i - f(x_i)`$)와 매개변수 보정량 $`\Delta\mathbf{a}`$ 사이에 편도함수 $`\partial f / \partial a_k`$를 원소로 하는 행렬 $`\mathbf{Z}`$에 대한 선형 관계가 생기므로, 일반 선형 최소제곱의 정규 방정식 $`\mathbf{Z}^T \mathbf{Z} \, \Delta\mathbf{a} = \mathbf{Z}^T \mathbf{d}`$를 풀어 보정량을 구한다. 매개변수를 보정한 뒤 선형화를 반복하고, 보정량이 충분히 작아지면 멈춘다. 6.2절의 코드는 이 방법 대신 범용 최적화 함수로 SSE를 직접 최소화한다.

### 6.2 scipy.optimize를 이용한 Python 구현

```python
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Sample vapor pressure data (e.g., for water) in mmHg at different temperatures (deg C)
T_data = np.array([40, 50, 60, 70, 80])    # Temperature in Celsius
P_data = np.array([55.3, 92.5, 149.4, 233.7, 355.1])  # Vapor pressure in mmHg

# Convert pressure to log base 10
logP_data = np.log10(P_data)

# Antoine equation model: log10(P) = A - B / (C + T)
def antoine_model(params, T):
    A, B, C = params
    return A - B / (C + T)

# Objective function: sum of squared errors between predicted and actual logP
def objective(params):
    predictions = antoine_model(params, T_data)
    return np.sum((logP_data - predictions) ** 2)

# Initial guess for A, B, C
initial_guess = [8, 1500, 200]

# Minimize the objective function
result = minimize(objective, initial_guess, method='L-BFGS-B')

# Extract optimized parameters
A_opt, B_opt, C_opt = result.x
```

> **참고:** **방법.** `scipy.optimize.minimize` 함수의 `method='L-BFGS-B'`는 유사 뉴턴(Quasi-Newton) 최적화 알고리즘을 사용한다. 초기 추정값이 필요하며, 알고리즘이 $`SSE`$를 최소화하도록 매개변수를 반복적으로 개선한다.

> **주의:** **초기 추정값 민감도.** 비선형 회귀 방법은 초기 추정값에 의존한다. 부적절한 초기 추정값은 전역 최솟값이 아닌 지역 최솟값으로의 수렴을 유발할 수 있다. 위 코드에서도 매개변수의 크기 차이가 커서 L-BFGS-B는 $`B`$가 초기값 1500 근처에 머문 지점에서 수렴 판정을 내린다. 같은 초기값으로 가우스-뉴턴 방법을 개선한 레벤버그-마쿼트 방법(`scipy.optimize.curve_fit`의 기본값)을 쓰면 SSE가 더 작은 $`A \approx 8.01`$, $`B \approx 1697`$, $`C \approx 230.6`$을 얻는다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|---|---|
| 선형 회귀 | 모델 $`\hat{y} = \beta_0 + \beta_1 x`$는 매개변수가 2개이고 14장의 닫힌 형태로 풀며, 추정의 표준 오차는 $`s_e = \sqrt{\frac{SSE}{n-2}}`$이다. |
| 다항 회귀 | 모델 $`\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2`$는 매개변수가 3개이고 $`3 \times 3`$ 정규 방정식으로 풀며, $`s_e = \sqrt{\frac{SSE}{n-3}}`$이다. |
| 다중 선형 회귀 | 모델 $`\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2`$는 매개변수가 3개이고 $`3 \times 3`$ 정규 방정식으로 풀며, $`s_e = \sqrt{\frac{SSE}{n-3}}`$이다. |
| 일반 선형 최소제곱 | 모델 $`\hat{y} = \sum \beta_j f_j(\mathbf{x})`$는 매개변수가 $`m+1`$개이고 $`\boldsymbol{\beta} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\mathbf{y}`$로 풀며, $`s_e = \sqrt{\frac{SSE}{n-(m+1)}}`$이다. |
| 선형화 모델 | 지수, 거듭제곱, 포화 성장 모델은 매개변수가 2개이고 변수를 변환한 뒤 선형 회귀로 풀며, $`s_e`$는 변환에 따라 다르다. |
| 비선형 회귀 | 앙투안 방정식 같은 모델은 매개변수 수가 모델마다 다르고 반복 최적화(예: L-BFGS-B)로 풀며, 매개변수가 $`p`$개이면 $`s_e = \sqrt{\frac{SSE}{n-p}}`$이다. |

---

<br>

## 점검 문제

1. **회귀의 표준 오차와 결정 계수:** 1.3절의 예제(키 = [170, 172, 168, 174, 176] cm, $`x = [1, 2, 3, 4, 5]`$)에서 최소제곱 직선, $`SSE`$, $`s_e`$, $`R^2`$를 구하라.

   <details>
   <summary>정답</summary>

   $`\bar{x} = 3`$, $`\bar{y} = 172`$이고 $`\sum (x_i - \bar{x})(y_i - \bar{y}) = 14`$, $`\sum (x_i - \bar{x})^2 = 10`$이므로 기울기는 $`\beta_1 = 1.4`$, 절편은 $`\beta_0 = 172 - 1.4 \cdot 3 = 167.8`$이다. 예측값은 169.2, 170.6, 172.0, 173.4, 174.8이고 잔차는 0.8, 1.4, -4.0, 0.6, 1.2이므로 $`SSE = 20.4`$이다. 따라서 $`s_e = \sqrt{20.4/3} \approx 2.61`$ cm이다. $`SST = 40`$이므로 $`R^2 = 1 - 20.4/40 = 0.49`$이며, 직선은 키 변동의 49%만 설명한다.

   </details>

2. **포화 성장률 모델의 선형화:** $`1/x`$에 대한 $`1/y`$의 선형 회귀 결과 기울기가 0.5, 절편이 0.2였다. 포화 성장률 모델 $`y = \frac{\alpha_3 x}{\beta_3 + x}`$의 $`\alpha_3`$와 $`\beta_3`$를 구하라.

   <details>
   <summary>정답</summary>

   절편이 $`1/\alpha_3 = 0.2`$이므로 $`\alpha_3 = 5`$이다. 기울기가 $`\beta_3/\alpha_3 = 0.5`$이므로 $`\beta_3 = 0.5 \cdot 5 = 2.5`$이다. 따라서 모델은 $`y = \frac{5x}{2.5 + x}`$이며, $`x`$가 커지면 $`y`$는 포화값 5에 가까워진다.

   </details>

3. **자유도:** 데이터 점 10개에 2차 다항식을 적합했더니 $`SSE = 3.5`$였다. 추정의 표준 오차를 구하고, 같은 $`SSE`$를 직선 적합에서 얻었다면 값이 어떻게 달라지는지 설명하라.

   <details>
   <summary>정답</summary>

   2차 다항식은 매개변수가 3개이므로 자유도는 $`n - 3 = 7`$이고 $`s_e = \sqrt{3.5/7} = \sqrt{0.5} \approx 0.707`$이다. 직선은 매개변수가 2개이므로 자유도가 $`n - 2 = 8`$이고 $`s_e = \sqrt{3.5/8} \approx 0.661`$이다. 매개변수가 많을수록 자유도가 줄어 같은 $`SSE`$에서도 $`s_e`$가 커진다.

   </details>

4. **매개변수에 대한 선형성:** 다음 모델을 일반 선형 최소제곱으로 직접 적합할 수 있는지, 변환 후 선형 회귀로 적합할 수 있는지, 비선형 회귀가 필요한지 판단하라. (a) $`\omega`$가 알려진 $`y = \beta_0 + \beta_1 \cos(\omega x) + \beta_2 \sin(\omega x)`$ (b) $`y = \alpha_1 e^{\beta_1 x}`$ (c) $`\log_{10} P_v = A - \frac{B}{C + T}`$

   <details>
   <summary>정답</summary>

   (a)는 기저 함수 $`1`$, $`\cos(\omega x)`$, $`\sin(\omega x)`$의 계수가 매개변수이고 매개변수에 대해 선형이므로 일반 선형 최소제곱으로 직접 적합한다. (b)는 매개변수에 대해 비선형이지만 $`\ln y = \ln \alpha_1 + \beta_1 x`$로 변환하면 선형 회귀로 적합할 수 있다. (c)는 $`C`$가 분모의 $`C + T`$ 안에 있어 $`\log_{10} P_v`$를 매개변수의 선형결합으로 쓸 수 없으므로 비선형 회귀가 필요하다.

   </details>
